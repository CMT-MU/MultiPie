"""
Utility functions for reading Wannier90 files used by ModelAnalyzer.
"""

from contextlib import contextmanager
from io import TextIOWrapper
from pathlib import Path
import gzip
import re
import tarfile

import numpy as np

BOHR2ANG = 0.529177249


# ==================================================
@contextmanager
def _open_text(path):
    """Open a plain, ``.gz``, or ``.tar.gz`` text file."""
    path = Path(path)

    if path.exists():
        with path.open("r", encoding="utf-8") as f:
            yield f
        return

    gz_path = Path(str(path) + ".gz")
    if gz_path.exists():
        with gzip.open(gz_path, "rt", encoding="utf-8") as f:
            yield f
        return

    tar_path = Path(str(path) + ".tar.gz")
    if tar_path.exists():
        with tarfile.open(tar_path, "r:gz") as tf:
            files = [member for member in tf.getmembers() if member.isfile()]
            exact = [member for member in files if Path(member.name).name == path.name]

            if exact:
                member = exact[0]
            elif len(files) == 1:
                member = files[0]
            else:
                raise FileNotFoundError(f"{path.name} is not found in {tar_path}")

            raw = tf.extractfile(member)
            if raw is None:
                raise FileNotFoundError(f"failed to read {member.name} in {tar_path}")

            with TextIOWrapper(raw, encoding="utf-8") as f:
                yield f
        return

    raise FileNotFoundError(f"failed to read file: {path}")


# ==================================================
def _read_lines(path):
    """Read a plain, ``.gz``, or ``.tar.gz`` text file as a list of lines."""
    with _open_text(path) as f:
        return f.read().splitlines()


# ==================================================
def _blocks(lines):
    """Return ``{block_name: lines_between_begin_and_end}``."""
    blocks = {}
    i = 0

    while i < len(lines):
        line = lines[i].strip().lower()
        if not line.startswith("begin "):
            i += 1
            continue

        name = line[6:].strip()
        end = f"end {name}"

        if name in blocks:
            raise ValueError(f"block '{name}' is defined more than once")

        for j in range(i + 1, len(lines)):
            if lines[j].strip().lower() == end:
                blocks[name] = lines[i + 1 : j]
                i = j + 1
                break
        else:
            raise ValueError(f"missing '{end}'")

    return blocks


# ==================================================
def _array(lines, dtype=float):
    """Read numerical rows and always return a two-dimensional array."""
    return np.atleast_2d(np.loadtxt(lines, dtype=dtype))


# ==================================================
def _convert_w90_orbital(l, m, r, s):
    """
    Convert a Wannier90 projection label to a MultiPie orbital label.

    Args:
        l (int): Angular-momentum index.
        m (int): Wannier90 angular-function index.
        r (int): Wannier90 radial-function index.
        s (int): Spin component: ``1`` for up, ``-1`` for down, and ``0``
                 for a non-spinor projection.

    Returns:
        str: MultiPie orbital label, such as ``s``, ``px``, ``du``,
             ``(px,u)``, or ``(px,d)``.

    Raises:
        ValueError: If ``(l,m)`` or ``s`` is not supported.

    Notes:
        The radial index ``r`` is currently not encoded in the returned
        MultiPie orbital label, but it is retained in the argument list for
        consistency with the Wannier90 projection format.
    """
    orbitals = {
        0: {1: (0, "s")},
        1: {1: (2, "pz"), 2: (0, "px"), 3: (1, "py")},
        2: {1: (4, "du"), 2: (2, "dxz"), 3: (3, "dyz"), 4: (0, "dv"), 5: (1, "dxy")},
        3: {
            1: (6, "faz"),
            2: (4, "f3x"),
            3: (5, "f3y"),
            4: (2, "fbz"),
            5: (3, "f3"),
            6: (0, "f2"),
            7: (1, "f1"),
        },
    }

    try:
        orbital = orbitals[l][m]
    except KeyError as exc:
        raise ValueError(f"invalid orbital projection: l={l}, m={m}, r={r}, s={s}") from exc

    if s == 1:
        return orbital[0], f"({orbital[1]},u)"
    if s == -1:
        return orbital[0], f"({orbital[1]},d)"
    if s == 0:
        return orbital

    raise ValueError(f"invalid spin component: s={s}")


# ==================================================
def _projection_data(rows, spinors):
    """Read Wannier projection metadata from an ``nnkp`` projection block."""
    num_wann = int(rows[0])
    step = 3 if spinors else 2

    if len(rows) < 1 + step * num_wann:
        raise ValueError("incomplete projection block")

    projections = [rows[1 + step * iw] for iw in range(num_wann)]
    orbital = np.zeros((num_wann, 3), dtype=int)
    spin = np.zeros(num_wann, dtype=int)
    spin_axis = np.zeros((num_wann, 3), dtype=float)

    for iw, row in enumerate(projections):
        values = row.split()
        if len(values) < 6:
            raise ValueError(f"invalid projection entry: {row}")

        # Columns 4--6 contain l, m, and the radial-function index.
        orbital[iw] = values[3:6]

        if spinors:
            values = rows[1 + step * iw + 2].split()
            if len(values) < 4:
                raise ValueError("invalid spinor projection entry")
            spin[iw] = int(values[0])
            spin_axis[iw] = values[1:4]

    # Wannier90 uses fixed-width projection strings.  These slices retain the
    # original convention for identifying orbital specifications and centers.
    orb_key = [row[:40] for row in projections]
    pos_key = [row[:35] for row in projections]
    orb_unique = list(dict.fromkeys(orb_key))
    pos_unique = list(dict.fromkeys(pos_key))

    atom_orb = [[iw for iw, value in enumerate(orb_key) if value == key] for key in orb_unique]
    atom_pos = [[iorb for iorb, value in enumerate(orb_unique) if key in value] for key in pos_unique]
    atom_pos_r = [[float(x) for x in key.split()[:3]] for key in pos_unique]

    # Map each Wannier function to its zero-based projection-center index.
    nw2n = np.zeros(num_wann, dtype=int)
    for ia, positions in enumerate(atom_pos):
        for ipos in positions:
            nw2n[atom_orb[ipos]] = ia

    return {
        "num_wann": num_wann,
        "num_atom": len(atom_pos_r),
        "nw2n": nw2n.tolist(),
        "nw2l": orbital[:, 0].tolist(),
        "nw2m": orbital[:, 1].tolist(),
        "nw2r": orbital[:, 2].tolist(),
        "nw2s": spin.tolist(),
        "nw2saxis": spin_axis.tolist(),
        "atom_orb": atom_orb,
        "atom_pos": atom_pos,
        "atom_pos_r": atom_pos_r,
    }


# ==================================================
def read_win(topdir, seedname):
    """
    Read ``seedname.win``.

    Args:
        topdir (str): Top directory containing ``seedname.win``.
        seedname (str): Wannier90 seedname.

    Returns:
        dict:
            - num_k       : Number of k points inferred from ``mp_grid`` (int).
            - num_bands   : Number of Bloch bands passed to Wannier90 (int).
            - num_wann    : Number of Wannier functions (int).
            - spinors     : Whether the Wannier functions are spinors (bool).
            - fermi_energy: Fermi energy in eV (float).
            - A           : Real lattice vectors ``[a1,a2,a3]`` in angstrom,
                            shape ``[3][3]`` (list).
            - atoms_frac  : Atomic positions in fractional coordinates,
                            ``{(symbol, occurrence): [r1,r2,r3]}`` (dict).
            - atoms_cart  : Atomic positions in Cartesian coordinates in
                            angstrom, ``{(symbol, occurrence): [x,y,z]}`` (dict).
            - mp_grid     : Monkhorst-Pack grid ``[nk1,nk2,nk3]`` (list).
            - kpoints     : Explicit k points in crystal coordinates, reduced
                            to ``0 <= k_i < 1``, shape ``[num_k][3]`` (list).
            - kpoint      : Labelled representative k points,
                            ``{label: [k1,k2,k3]}`` (dict or None).
            - kpoint_path : Labelled high-symmetry path such as
                            ``G-X-M-G|R-X`` (str or None).

    Raises:
        FileNotFoundError: If the input file cannot be found.
        ValueError: If the file contents are inconsistent with the expected
                    Wannier90 format.

    Notes:
        If only one of ``atoms_frac`` and ``atoms_cart`` is present, the other
        is calculated using ``A``.  Cartesian lengths are returned in angstrom.
    """
    lines = _read_lines(Path(topdir) / f"{seedname}.win")
    blocks = _blocks(lines)

    def param(name, default, cast):
        """Read a scalar ``name = value`` or ``name : value`` parameter."""
        pattern = re.compile(rf"^\s*{re.escape(name)}\s*[=:]\s*([^!]*)", re.IGNORECASE)
        values = [match.group(1).strip() for line in lines if (match := pattern.match(line))]

        if len(values) > 1:
            raise ValueError(f"{name} is defined more than once")
        if not values:
            return default

        value = values[0]
        if cast is bool:
            value = value.lower()
            if value in ("true", ".true."):
                return True
            if value in ("false", ".false."):
                return False
            raise ValueError(f"invalid boolean value for {name}: {values[0]}")

        return cast(value)

    def vectors(name):
        """Read a three-vector block and convert bohr to angstrom."""
        rows = blocks.get(name)
        if not rows:
            return None

        unit = rows[0].strip().lower()
        start = int(unit in ("ang", "bohr"))
        scale = BOHR2ANG if unit == "bohr" else 1.0
        data = _array([row for row in rows[start:] if row.strip()]) * scale

        if data.shape != (3, 3):
            raise ValueError(f"{name} must contain three three-dimensional vectors")

        return data

    def atoms(name):
        """Read an atoms block as ``{(symbol, occurrence): position}``."""
        rows = blocks.get(name)
        if not rows:
            return None

        unit = rows[0].strip().lower()
        start = int(unit in ("ang", "bohr"))
        scale = BOHR2ANG if unit == "bohr" else 1.0
        count, out = {}, {}

        for row in rows[start:]:
            if not row.strip():
                continue

            values = row.split()
            if len(values) != 4:
                raise ValueError(f"invalid entry in {name}: {row}")

            symbol, *xyz = values
            count[symbol] = count.get(symbol, 0) + 1
            out[symbol, count[symbol]] = (np.asarray(xyz, dtype=float) * scale).tolist()

        return out

    # --- scalar parameters ---
    mp_grid = [int(x) for x in param("mp_grid", "1 1 1", str).split()]
    if len(mp_grid) != 3:
        raise ValueError("mp_grid must contain three integers")

    # --- lattice vectors and atomic positions ---
    A = vectors("unit_cell_cart")
    atoms_frac = atoms("atoms_frac")
    atoms_cart = atoms("atoms_cart")

    if A is not None:
        if atoms_cart is not None and atoms_frac is None:
            atoms_frac = {key: (np.asarray(pos) @ np.linalg.inv(A)).tolist() for key, pos in atoms_cart.items()}
        elif atoms_frac is not None and atoms_cart is None:
            atoms_cart = {key: (np.asarray(pos) @ A).tolist() for key, pos in atoms_frac.items()}

    d = {
        "num_k": int(np.prod(mp_grid)),
        "num_bands": param("num_bands", 0, int),
        "num_wann": param("num_wann", 0, int),
        "spinors": param("spinors", False, bool),
        "fermi_energy": param("fermi_energy", 0.0, float),
        "A": None if A is None else A.tolist(),
        "atoms_frac": atoms_frac,
        "atoms_cart": atoms_cart,
        "mp_grid": mp_grid,
        "kpoints": [[0.0, 0.0, 0.0]],
        "kpoint": None,
        "kpoint_path": None,
    }

    # --- explicit k points ---
    if rows := blocks.get("kpoints"):
        kpoints = np.mod(_array([row for row in rows if row.strip()])[:, :3], 1)
        if len(kpoints) != d["num_k"]:
            raise ValueError(
                f"inconsistent number of k points in {seedname}.win: "
                f"mp_grid gives {d['num_k']}, but the kpoints block contains {len(kpoints)}"
            )
        d["kpoints"] = kpoints.tolist()

    # --- labelled high-symmetry k-point path ---
    if rows := blocks.get("kpoint_path"):
        points, paths = {}, []

        for row in rows:
            if not row.strip():
                continue

            values = row.split()
            if len(values) != 8:
                raise ValueError(f"invalid kpoint_path entry: {row}")

            X, x1, x2, x3, Y, y1, y2, y3 = values
            points.setdefault(X, list(map(float, (x1, x2, x3))))
            points.setdefault(Y, list(map(float, (y1, y2, y3))))

            if paths and paths[-1][-1] == X:
                paths[-1].append(Y)
            else:
                paths.append([X, Y])

        d["kpoint"] = points
        d["kpoint_path"] = "|".join("-".join(path) for path in paths)

    return d


# ==================================================
def read_nnkp(topdir, seedname):
    """
    Read ``seedname.nnkp`` generated by Wannier90.

    Args:
        topdir (str): Top directory containing ``seedname.nnkp``.
        seedname (str): Wannier90 seedname.

    Returns:
        dict:
            - A                : Real lattice vectors ``[a1,a2,a3]`` in
                                 Cartesian coordinates, shape ``[3][3]`` (list).
            - B                : Reciprocal lattice vectors ``[b1,b2,b3]`` in
                                 Cartesian coordinates, shape ``[3][3]`` (list).
            - num_k            : Number of k points (int).
            - num_wann         : Number of Wannier projections (int).
            - num_atom         : Number of distinct projection centers (int).
            - num_b            : Number of b-vectors per k point (int).
            - kpoints          : k points in crystal coordinates, reduced to
                                 ``0 <= k_i < 1``, shape ``[num_k][3]`` (list).
            - kpoints_wo_shift : Original k points before modulo reduction,
                                 shape ``[num_k][3]`` (list).
            - nnkpts           : Nearest-neighbour table
                                 ``[ik,ikb,G1,G2,G3]``, shape
                                 ``[num_k][num_b][5]`` (list).
            - nw2n             : Zero-based projection-center index of each
                                 Wannier function, shape ``[num_wann]`` (list).
            - nw2l, nw2m       : Angular indices ``l`` and ``m`` of each
                                 projection, shape ``[num_wann]`` (list).
            - nw2r             : Radial-function index of each projection,
                                 shape ``[num_wann]`` (list).
            - nw2s             : Spin component of each projection; zero for
                                 non-spinor projections (list).
            - nw2saxis         : Spin-quantization axis of each projection,
                                 shape ``[num_wann][3]`` (list).
            - atom_orb         : Wannier-function indices grouped by orbital
                                 projection (list).
            - atom_pos         : Orbital-projection indices grouped by
                                 projection center (list).
            - atom_pos_r       : Projection-center positions in fractional
                                 coordinates, shape ``[num_atom][3]`` (list).
            - bvec_cart        : b-vectors at Gamma in Cartesian reciprocal
                                 coordinates, shape ``[num_b][3]`` (list).
            - bvec_crys        : The same b-vectors in crystal coordinates,
                                 shape ``[num_b][3]`` (list).
            - wb               : Finite-difference weights of ``bvec_cart``,
                                 shape ``[num_b]`` (list).
            - bveck            : Cartesian b-vector for each ``(k,b)`` pair,
                                 shape ``[num_k][num_b][3]`` (list).
            - wk               : Finite-difference weight for each ``(k,b)``
                                 pair, shape ``[num_k][num_b]`` (list).
            - kb2k             : Zero-based index of ``k+b`` for each pair,
                                 shape ``[num_k][num_b]`` (list).

    Raises:
        FileNotFoundError: If the input file cannot be found.
        ValueError: If required blocks are missing or the file contents are
                    inconsistent with the expected Wannier90 format.

    Notes:
        ``num_atom`` and ``atom_pos_r`` describe Wannier projection centers,
        which need not coincide with the physical atoms in ``seedname.win``.
        For ``auto_projections``, only ``num_wann`` is available; the detailed
        projection entries are returned as ``None``.
    """
    blocks = _blocks(_read_lines(Path(topdir) / f"{seedname}.nnkp"))

    required = ("real_lattice", "recip_lattice", "kpoints", "nnkpts")
    missing = [name for name in required if name not in blocks]
    if missing:
        raise ValueError(f"missing blocks in {seedname}.nnkp: {', '.join(missing)}")

    # --- lattice vectors ---
    A = _array(blocks["real_lattice"])
    B = _array(blocks["recip_lattice"])
    if A.shape != (3, 3) or B.shape != (3, 3):
        raise ValueError("real_lattice and recip_lattice must have shape [3][3]")

    # --- k points ---
    rows = blocks["kpoints"]
    num_k = int(rows[0])
    kpoints_wo_shift = _array(rows[1 : 1 + num_k])
    if len(kpoints_wo_shift) != num_k:
        raise ValueError("invalid number of k points")
    kpoints = np.mod(kpoints_wo_shift, 1)

    # --- nearest-neighbour k points ---
    rows = blocks["nnkpts"]
    num_b = int(rows[0])
    nnkpts = _array(rows[1 : 1 + num_k * num_b], dtype=int)
    if len(nnkpts) != num_k * num_b:
        raise ValueError("invalid number of nearest-neighbour k-point entries")
    nnkpts = nnkpts.reshape(num_k, num_b, 5)

    d = {
        "A": A.tolist(),
        "B": B.tolist(),
        "num_k": num_k,
        "num_wann": 0,
        "num_atom": 0,
        "num_b": num_b,
        "kpoints": kpoints.tolist(),
        "kpoints_wo_shift": kpoints_wo_shift.tolist(),
        "nnkpts": nnkpts.tolist(),
        "nw2n": None,
        "nw2l": None,
        "nw2m": None,
        "nw2r": None,
        "nw2s": None,
        "nw2saxis": None,
        "atom_orb": None,
        "atom_pos": None,
        "atom_pos_r": None,
        "bvec_cart": None,
        "bvec_crys": None,
        "wb": None,
        "bveck": None,
        "wk": None,
        "kb2k": None,
    }

    # --- projection information ---
    projection_block = "spinor_projections" if "spinor_projections" in blocks else "projections"
    projection_rows = blocks.get(projection_block)

    if projection_rows and int(projection_rows[0]) > 0:
        d.update(
            _projection_data(
                projection_rows,
                spinors=projection_block.startswith("spinor"),
            )
        )
    elif "auto_projections" in blocks:
        # With automatic projections, nnkp stores num_wann but does not store
        # the per-Wannier projection metadata listed above.
        d["num_wann"] = int(blocks["auto_projections"][0])
    else:
        raise ValueError(f"projection information is missing in {seedname}.nnkp")

    # b = k_ikb + G - k_ik in crystal coordinates, followed by conversion to
    # Cartesian reciprocal coordinates using the row-vector convention.
    b_crys = kpoints_wo_shift[nnkpts[..., 1] - 1] + nnkpts[..., 2:5] - kpoints_wo_shift[nnkpts[..., 0] - 1]
    bveck = b_crys @ B

    # Wannier90 finite-difference shells are identified using the Gamma point.
    gamma = np.flatnonzero(np.all(np.isclose(kpoints_wo_shift, 0), axis=1))
    if not len(gamma):
        raise ValueError("Gamma point must be included")

    bvec_cart = bveck[gamma[0]]
    bvec_crys = bvec_cart @ A.T / (2 * np.pi)

    # Solve sum_b w_b b_a b_c = delta_ac using the pseudoinverse.
    bbmat = np.einsum("bi,bj->bij", bvec_cart, bvec_cart).reshape(num_b, 9)
    wb = np.eye(3).ravel() @ np.linalg.pinv(bbmat)

    # Assign the corresponding Gamma-shell weight to every (k,b) pair.
    match = np.isclose(
        bveck[..., None, :],
        bvec_cart[None, None, ...],
        rtol=1e-5,
        atol=1e-5,
    ).all(axis=-1)
    if not match.any(axis=-1).all():
        raise ValueError("unknown b-vector")

    d.update(
        bvec_cart=bvec_cart.tolist(),
        bvec_crys=bvec_crys.tolist(),
        wb=wb.tolist(),
        bveck=bveck.tolist(),
        wk=wb[match.argmax(axis=-1)].tolist(),
        kb2k=(nnkpts[..., 1] - 1).tolist(),
    )

    return d


# ==================================================
def merge_wannier_info(win, nnkp, seedname):
    """
    Check duplicated entries in ``win`` and ``nnkp``, then merge them.

    Args:
        win (dict): Information returned by :func:`read_win`.
        nnkp (dict): Information returned by :func:`read_nnkp`.
        seedname (str): Wannier90 seedname used in error messages.

    Returns:
        dict:
            Merged Wannier90 information.  Every key present in both input
            dictionaries is checked for consistency, and the value from
            ``win`` is retained for the duplicated key.  With the current
            readers, the duplicated entries are ``A``, ``num_k``,
            ``num_wann``, and ``kpoints``.

    Raises:
        ValueError: If any duplicated entry is inconsistent between
                    ``seedname.win`` and ``seedname.nnkp``.
    """
    common_keys = win.keys() & nnkp.keys()

    for key in sorted(common_keys):
        value_win = win[key]
        value_nnkp = nnkp[key]

        if isinstance(value_win, (list, tuple, np.ndarray)) or isinstance(value_nnkp, (list, tuple, np.ndarray)):
            a = np.asarray(value_win)
            b = np.asarray(value_nnkp)

            if a.shape != b.shape:
                raise ValueError(
                    f"inconsistent {key} between {seedname}.win and " f"{seedname}.nnkp: shape {a.shape} != {b.shape}"
                )

            if np.issubdtype(a.dtype, np.number) and np.issubdtype(b.dtype, np.number):
                same = np.allclose(a, b, rtol=1e-7, atol=1e-8)
                detail = f"maximum difference = {np.max(np.abs(a - b)):.3e}"
            else:
                same = np.array_equal(a, b)
                detail = "values differ"
        else:
            same = value_win == value_nnkp
            detail = f"win={value_win!r}, nnkp={value_nnkp!r}"

        if not same:
            raise ValueError(f"inconsistent {key} between {seedname}.win and " f"{seedname}.nnkp: {detail}")

    # Duplicated entries have already been checked; retain the win values.
    return win | {key: value for key, value in nnkp.items() if key not in common_keys}


# ==================================================
def matrix_dict_r(Or, rpoints, diagonal=False):
    """
    dictionary form of an arbitrary operator matrix in real-space representation.

    Args:
        Or (ndarray): real-space representation of the given operator, O_{ab}(R) = <φ_{a}(0)|O|φ_{b}(R)>.
        rpoints (ndarray): lattice points (crystal coordinate, [[n1,n2,n3]], nj: integer).
        diagonal (bool, optional): diagonal matrix ?

    Returns:
        dict: real-space representation of the given operator, {(n2,n2,n3,a,b) = O_{ab}(R)}.
    """
    # number of pseudo atomic orbitals
    dim_r = len(Or[0])
    if not diagonal:
        dim_c = len(Or[0][0])

    Or_dict = {}

    r_list = [[r, round(n1), round(n2), round(n3)] for r, (n1, n2, n3) in enumerate(rpoints)]

    if diagonal:
        Or_dict = {(n1, n2, n3, a, a): Or[r][a] for r, n1, n2, n3 in r_list for a in range(dim_r)}
    else:
        Or_dict = {(n1, n2, n3, a, b): Or[r][a][b] for r, n1, n2, n3 in r_list for a in range(dim_r) for b in range(dim_c)}

    return Or_dict


# ==================================================
def read_hr(topdir, hr_file):
    """
    Read ``seedname_hr.dat`` generated by Wannier90.

    Args:
        topdir (str): Top directory containing ``seedname_hr.dat``.
        hr_file (str): Wannier90 Hamiltonian file.

    Returns:
        dict:
            - HH_R   : Real-space Hamiltonian matrices ``H_mn(R)`` in eV,
                       shape ``[num_R][num_wann][num_wann]`` (ndarray).
            - irvec  : Lattice vectors ``R = n1*a1 + n2*a2 + n3*a3`` in
                       crystal coordinates, shape ``[num_R][3]`` (ndarray).
            - ndegen : Wigner-Seitz degeneracy of each ``R`` vector, shape
                       ``[num_R]`` (ndarray).

    Raises:
        FileNotFoundError: If the input file cannot be found.
        ValueError: If the file is incomplete, contains invalid indices, has
                    inconsistent ``R`` blocks, or violates
                    ``H(R) = H(-R)^dagger``.

    Notes:
        Wannier indices in the file are one-based and are converted to
        zero-based array indices.  ``ndegen`` is the denominator used in the
        Wannier90 Fourier interpolation over ``R`` vectors.  The file is read
        sequentially to avoid storing the full text representation in memory.
    """

    def next_nonempty(f):
        """Return the next non-empty line or raise a format error."""
        for line in f:
            if line.strip():
                return line
        raise ValueError(f"unexpected end of {hr_file}")

    with _open_text(Path(topdir) / f"{hr_file}") as f:
        # The first line is a comment written by Wannier90.
        try:
            next(f)
        except StopIteration as exc:
            raise ValueError(f"invalid hr file: {hr_file}") from exc

        num_wann = int(next_nonempty(f))
        num_R = int(next_nonempty(f))
        if num_wann <= 0 or num_R <= 0:
            raise ValueError("num_wann and num_R must be positive")

        # Degeneracies are normally written with at most 15 integers per line.
        ndegen = []
        while len(ndegen) < num_R:
            ndegen.extend(map(int, next_nonempty(f).split()))

        if len(ndegen) != num_R or any(value <= 0 for value in ndegen):
            raise ValueError(f"invalid degeneracies in {hr_file}")

        ndegen = np.asarray(ndegen, dtype=int)
        irvec = np.empty((num_R, 3), dtype=int)
        HH_R = np.zeros((num_R, num_wann, num_wann), dtype=complex)
        r_index = {}

        # Each R block contains num_wann**2 Hamiltonian matrix elements.
        for iR in range(num_R):
            R_block = None
            assigned = np.zeros((num_wann, num_wann), dtype=bool)

            for _ in range(num_wann * num_wann):
                values = next_nonempty(f).split()
                if len(values) != 7:
                    raise ValueError(f"invalid Hamiltonian entry: {' '.join(values)}")

                n1, n2, n3, m, n = map(int, values[:5])
                R = (n1, n2, n3)

                if R_block is None:
                    R_block = R
                elif R != R_block:
                    raise ValueError(f"inconsistent R vectors within block {iR}")

                if not (1 <= m <= num_wann and 1 <= n <= num_wann):
                    raise ValueError(f"invalid Wannier indices m={m}, n={n}")
                if assigned[m - 1, n - 1]:
                    raise ValueError(f"duplicated Hamiltonian entry for R={R}, m={m}, n={n}")

                assigned[m - 1, n - 1] = True
                HH_R[iR, m - 1, n - 1] = float(values[5]) + 1j * float(values[6])

            if not assigned.all():
                raise ValueError(f"incomplete Hamiltonian block for R={R_block}")
            if R_block in r_index:
                raise ValueError(f"duplicated R vector: {R_block}")

            irvec[iR] = R_block
            r_index[R_block] = iR

        # No additional non-empty Hamiltonian rows should remain.
        if any(line.strip() for line in f):
            raise ValueError(f"too many Hamiltonian entries in {hr_file}")

    # Check H(R) = H(-R)^dagger.
    for iR, R in enumerate(map(tuple, irvec)):
        minus_R = tuple(-np.asarray(R))
        if minus_R not in r_index:
            raise ValueError(f"missing -R vector for R={R}")

        jR = r_index[minus_R]
        if not np.allclose(HH_R[iR], HH_R[jR].conj().T, rtol=1e-8, atol=1e-8):
            diff = np.max(np.abs(HH_R[iR] - HH_R[jR].conj().T))
            raise ValueError(f"Hermiticity is violated for R={R}: maximum difference = {diff:.3e}")

    hr_dict = matrix_dict_r(HH_R, irvec)

    return hr_dict, irvec, ndegen


# ==================================================
def build_ket_wannier(nnkp, site_dict, rtol=1e-4, atol=1e-4):
    """
    Construct the MultiPie ket basis from Wannier90 projection information.

    Args:
        nnkp (dict): Information returned by :func:`read_nnkp`.
        site_dict (dict): MultiPie sites in the form
                          ``{(name, sublattice): fractional_position}``.
        rtol (float, optional): Relative tolerance for comparing fractional
                               positions.
        atol (float, optional): Absolute tolerance for comparing fractional
                               positions.

    Returns:
        list:
            Wannier ket basis with one entry per Wannier function.  Each entry
            has the form ``[name, sublattice, l, orbital]``.

    Raises:
        ValueError: If detailed projection information is unavailable, a
                    projection center cannot be matched to exactly one
                    MultiPie site, or the projection metadata are inconsistent.

    Notes:
        Fractional positions are compared modulo lattice translations, so
        coordinates differing by an integer lattice vector are treated as the
        same site.
    """
    keys = ("nw2n", "nw2l", "nw2m", "nw2r", "nw2s", "atom_pos_r")
    missing = [key for key in keys if nnkp.get(key) is None]
    if missing:
        raise ValueError(
            "ket_wannier cannot be constructed because detailed projection " f"information is unavailable: {', '.join(missing)}"
        )

    arrays = [nnkp[key] for key in ("nw2n", "nw2l", "nw2m", "nw2r", "nw2s")]
    lengths = {len(values) for values in arrays}
    if lengths != {nnkp["num_wann"]}:
        raise ValueError("inconsistent number of Wannier projection entries")

    site_items = [((name, sublattice), np.asarray(position, dtype=float)) for (name, sublattice), position in site_dict.items()]

    ket_wannier = []
    for iw, (pos_idx, l, m, r, s) in enumerate(zip(*arrays)):
        if not 0 <= pos_idx < len(nnkp["atom_pos_r"]):
            raise ValueError(f"invalid projection-center index for Wannier function {iw}: {pos_idx}")

        position = np.asarray(nnkp["atom_pos_r"][pos_idx], dtype=float)
        matches = []

        for site, site_position in site_items:
            delta = position - site_position
            delta -= np.rint(delta)
            if np.allclose(delta, 0.0, rtol=rtol, atol=atol):
                matches.append(site)

        if len(matches) == 0:
            raise ValueError(f"no MultiPie site matches Wannier projection center {position.tolist()}")
        if len(matches) > 1:
            raise ValueError(f"multiple MultiPie sites match Wannier projection center " f"{position.tolist()}: {matches}")

        name, sublattice = matches[0]
        comp, orbital = _convert_w90_orbital(l, m, r, s)
        ket_wannier.append([name, sublattice, l, comp, orbital])

    return ket_wannier


# ==================================================
def sort_ket_matrix(Ok, ket1, ket2):
    """
    sort ket to align with the MultiPie definition (ket2).

    Args:
        Ok (ndarray): arbitrary operator in k-space representation, O_{ab}(k) = <φ_{a}(k)|H|φ_{b}(k)>.
        ket1 (list): ket basis list, [[atom name, sublattice, rank, orbital]]
        ket2 (list): ket basis list, orbital@site.

    Returns:
        Ok (ndarray): operator.
    """
    idx_list = [ket1.index(o) for o in ket2]
    Ok = Ok[:, idx_list, :]
    Ok = Ok[:, :, idx_list]

    return Ok


# ==================================================
def sort_ket_matrix_dict(Or_dict, ket1, ket2):
    """
    sort ket to align with the MultiPie definition (ket2).

    Args:
        Or_dict (dict): dictionary form of an arbitrary operator matrix in reak-space/k-space representation.
        ket1 (list): ket basis list, [[atom name, sublattice, rank, orbital]]
        ket2 (list): ket basis list, orbital@site.

    Returns:
        Or_dict (dict):  dictionary form of an arbitrary operator matrix.
    """
    idx_list = [ket2.index(o) for o in ket1]
    Or_dict = {(n1, n2, n3, idx_list[a], idx_list[b]): v for (n1, n2, n3, a, b), v in Or_dict.items()}

    return Or_dict


# ==================================================
def sort_ket_list(lst, ket1, ket2):
    """
    sort ket to align with the MultiPie definition (ket2).

    Args:
        lst (list): arbitrary list that has the same dimensions as the ket.
        ket1 (list): ket basis list, [[atom name, sublattice, rank, orbital]]
                ket2 (list): ket basis list, orbital@site.

    Returns:
        lst (list): arbitrary list that has the same dimensions as the ket.
    """
    idx_list = [ket1.index(o) for o in ket2]

    lst = np.array(lst)

    if lst.ndim == 1:
        lst = list(np.array(lst)[idx_list])
    elif lst.ndim == 2:
        lst = list(np.array(lst)[idx_list, :])
    elif lst.ndim == 3:
        lst = list(np.array(lst)[idx_list, :, :])
    else:
        raise Exception(f"invalid dimension of lst = {lst.ndim} was given.")

    return list(lst)


# ==================================================
def decompose_operator_by_SAMB(Or_dict, Zr_dict):
    """
    decompose arbitrary operator into linear combination of SAMBs.

    Args:
        Or_dict (dict): dictionary form of an arbitrary operator matrix in real-space/k-space representation.
        Zr_dict (dict): dictionary form of SAMBs.

    Returns:
        z (dict): parameter set, {zj: coeff}.
    """
    z = {
        zj: np.real(np.sum([v * Or_dict.get((-k[0], -k[1], -k[2], k[4], k[3]), 0) for k, v in d.items()]))
        for zj, d in Zr_dict.items()
    }

    return z


# ==================================================
def read_umat(topdir, seedname):
    """Read ``seedname_u.mat``, ``seedname_u_dis.mat``.  Not implemented yet."""
    pass


# ==================================================
def read_eig(topdir, seedname):
    """Read ``seedname.eig``.  Not implemented yet."""
    pass


# ==================================================
def read_mmn(topdir, seedname):
    """Read ``seedname.mmn``.  Not implemented yet."""
    pass


# ==================================================
def read_spn(topdir, seedname):
    """Read ``seedname.spn``.  Not implemented yet."""
    pass


# ==================================================
def read_uHu(topdir, seedname):
    """Read ``seedname.uHu``.  Not implemented yet."""
    pass


# ==================================================
def read_uIu(topdir, seedname):
    """Read ``seedname.uIu``.  Not implemented yet."""
    pass
