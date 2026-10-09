"""
Regression tests for reading Wannier90 data in ModelAnalyzer (wannier/symcw mode).

The symcw tests write synthetic Wannier90 files (win, nnkp, hr.dat) from the H(R) of a MultiPie model,
expressed in another primitive cell, origin and orbital order, and check that the SAMB parameters are recovered.
"""

import gzip
import copy
import os
import shutil
import tarfile
import warnings

import numpy as np
import pytest
import seekpath

from multipie import Group, MaterialModel
from multipie.core.cmd import create_model
from multipie.core.model_analyzer import ModelAnalyzer, _join_kpath
from multipie.util.util_model_analyzer import fourier_r_to_k
from multipie.util.util_wannier import (
    convert_hr_to_model,
    convert_to_conventional_frac,
    create_ket_wannier_multipie,
    find_lattice_transformation,
    get_or_add_vector,
    map_wannier_to_model,
    model_primitive_vector,
    apply_ws_degeneracy,
    read_hr,
    read_nnkp,
    read_win,
    read_wsvec,
)

# Mn2Au-type, I4/mmm (Mn 4e, Au 2a).
A_LAT, C_LAT, Z_MN = 3.328, 8.539, 1 / 3
# primitive cell of Quantum ESPRESSO, ibrav=7.
A_QE = [[A_LAT / 2, -A_LAT / 2, C_LAT / 2], [A_LAT / 2, A_LAT / 2, C_LAT / 2], [-A_LAT / 2, -A_LAT / 2, C_LAT / 2]]
# Wannier90 (l, m) of orbitals.
W90_LM = {"s": (0, 1), "pz": (1, 1), "px": (1, 2), "py": (1, 3)}


# ==================================================
@pytest.fixture(autouse=True)
def unchanged_cwd():
    """
    Check that the working directory is not changed by MultiPie (restored for the following tests).
    """
    cwd = os.getcwd()
    yield
    changed = os.getcwd()
    os.chdir(cwd)
    assert changed == cwd


# ==================================================
def wannier_info(A, cart):
    """
    Wannier info. for create_ket_wannier_multipie (one s projection per atom).
    """
    A = np.asarray(A, dtype=float)
    frac = {k: list(np.linalg.solve(A.T, np.asarray(v, dtype=float))) for k, v in cart.items()}
    n = list(range(len(cart)))
    return dict(
        A=A,
        atoms_frac=frac,
        atoms_cart=cart,
        fermi_energy=0.0,
        nw2n=n,
        nw2l=[0] * len(n),
        nw2m=[1] * len(n),
        nw2r=[1] * len(n),
        nw2s=[0] * len(n),
        atom_pos_r=list(frac.values()),
    )


# ==================================================
def test_get_or_add_vector_different_length():
    lst = [np.zeros(6)]
    assert get_or_add_vector(lst, np.zeros(12)) == 1
    assert get_or_add_vector(lst, np.zeros(6)) == 0


# ==================================================
def test_ket_different_multiplicity():
    # rutile TiO2 (P4_2/mnm), Ti 2a and O 4f.
    a, c, u = 4.594, 2.959, 0.305
    cart = {
        ("Ti", 1): [0, 0, 0],
        ("Ti", 2): [a / 2, a / 2, c / 2],
        ("O", 1): [u * a, u * a, 0],
        ("O", 2): [-u * a, -u * a, 0],
        ("O", 3): [(0.5 + u) * a, (0.5 - u) * a, c / 2],
        ("O", 4): [(0.5 - u) * a, (0.5 + u) * a, c / 2],
    }
    _, _, ket, _, _ = create_ket_wannier_multipie(wannier_info(np.diag([a, a, c]), cart))
    assert len(ket) == 6


# ==================================================
def test_wyckoff_centred_lattice_in_primitive_cell():
    cart = {("Mn", 1): [0, 0, Z_MN * C_LAT], ("Mn", 2): [0, 0, -Z_MN * C_LAT], ("Au", 1): [0, 0, 0]}
    info = wannier_info(A_QE, cart)
    no, conv = convert_to_conventional_frac(info["A"], info["atoms_frac"])
    assert no == 139
    assert Group(139).find_wyckoff_site(conv[("Mn", 1)])[0] == "4e"
    assert Group(139).find_wyckoff_site(conv[("Au", 1)])[0] == "2a"
    _, _, ket, _, _ = create_ket_wannier_multipie(info)
    assert len(ket) == 3


# ==================================================
def test_ket_projection_centres_not_in_atom_order():
    # rutile TiO2: projections on O(4), Ti(1) and O(2) shifted by a lattice vector, not in the order of atoms.
    a, c, u = 4.594, 2.959, 0.305
    cart = {
        ("Ti", 1): [0, 0, 0],
        ("Ti", 2): [a / 2, a / 2, c / 2],
        ("O", 1): [u * a, u * a, 0],
        ("O", 2): [-u * a, -u * a, 0],
        ("O", 3): [(0.5 + u) * a, (0.5 - u) * a, c / 2],
        ("O", 4): [(0.5 - u) * a, (0.5 + u) * a, c / 2],
    }
    info = wannier_info(np.diag([a, a, c]), cart)
    frac = info["atoms_frac"]
    centres = [frac[("O", 4)], frac[("Ti", 1)], (np.asarray(frac[("O", 2)]) + [1, 0, -1]).tolist()]
    info.update(nw2n=[0, 1, 2], nw2l=[0] * 3, nw2m=[1] * 3, nw2r=[1] * 3, nw2s=[0] * 3, atom_pos_r=centres)
    w2m, m2w, ket, pos, cart_pos = create_ket_wannier_multipie(info)
    assert [ket[w2m[w]].split("@")[1][:2] for w in range(3)] == ["O(", "Ti", "O("]
    assert np.allclose([pos[w2m[w]] for w in range(3)], centres)
    assert np.allclose(cart_pos, np.asarray(pos) @ info["A"])

    # a projection centre not on an atom (e.g., bond centre) is named X1.
    info["atom_pos_r"] = [[0.1, 0.2, 0.3]] + centres[1:]
    w2m, _, ket, pos, _ = create_ket_wannier_multipie(info)
    assert ket[w2m[0]] == "s@X1(1)"
    assert np.allclose(pos[w2m[0]], [0.1, 0.2, 0.3])


# ==================================================
def test_ket_without_projection_centres():
    # without atom_pos_r, the n-th projection centre is the n-th atom in seedname.win.
    cart = {("Mn", 1): [0, 0, Z_MN * C_LAT], ("Mn", 2): [0, 0, -Z_MN * C_LAT], ("Au", 1): [0, 0, 0]}
    info = wannier_info(A_QE, cart)
    legacy = {k: v for k, v in info.items() if k != "atom_pos_r"}
    result, result_legacy = create_ket_wannier_multipie(info), create_ket_wannier_multipie(legacy)
    assert result[:3] == result_legacy[:3]
    assert np.allclose(result[3], result_legacy[3])
    assert np.allclose(result[4], result_legacy[4])


# ==================================================
def write_wannier(path, seed, A, centers, wann, HR, ndegen=None, wsvec=None, species=None, atoms=None):
    """
    Write seed.win, seed.nnkp and seed_hr.dat (and seed_wsvec.dat if wsvec is given).

    Args:
        wann (list): (centre index, l, m) of each Wannier function.
        HR (dict): H(R), dict[((n1,n2,n3), a, b), value], written as it is.
        ndegen (dict, optional): dict[(n1,n2,n3), ndegen], 1 if not given.
        wsvec (dict, optional): dict[((n1,n2,n3), a, b), [T]], [(0,0,0)] if not given.
        species (list, optional): element of each centre, "X" if not given.
        atoms (list, optional): atoms in seedname.win, [(element, position)], (species, centers) if not given.
    """
    species = species or ["X"] * len(centers)
    atoms = atoms or list(zip(species, centers))
    ndegen = ndegen or {}
    os.makedirs(path, exist_ok=True)
    nw = len(wann)
    A = np.asarray(A)
    B = 2 * np.pi * np.linalg.inv(A).T
    with open(os.path.join(path, f"{seed}.win"), "w") as f:
        f.write(f"num_wann = {nw}\nnum_bands = {nw}\nmp_grid = 1 1 1\n\nbegin unit_cell_cart\nAng\n")
        f.writelines("  %.10f %.10f %.10f\n" % tuple(v) for v in A)
        f.write("end unit_cell_cart\n\nbegin atoms_frac\n")
        f.writelines("%s %.10f %.10f %.10f\n" % (e, *c) for e, c in atoms)
        f.write("end atoms_frac\n")
    with open(os.path.join(path, f"{seed}.nnkp"), "w") as f:
        f.write("synthetic\n\nbegin real_lattice\n")
        f.writelines("  %.7f %.7f %.7f\n" % tuple(v) for v in A)
        f.write("end real_lattice\n\nbegin recip_lattice\n")
        f.writelines("  %.7f %.7f %.7f\n" % tuple(v) for v in B)
        f.write("end recip_lattice\n\nbegin kpoints\n  1\n  0.0 0.0 0.0\nend kpoints\n\n")
        f.write(f"begin projections\n {nw:5d}\n")
        for ic, l, m in wann:
            c = centers[ic]
            f.write(f"{c[0]:10.5f}{c[1]:11.5f}{c[2]:11.5f}{l:6d}{m:3d}{1:3d}\n")
            f.write("    0.0000000  0.0000000  1.0000000   1.0000000  0.0000000  0.0000000    1.00\n")
        f.write("end projections\n\nbegin nnkpts\n  6\n")
        f.writelines("  1  1  %d %d %d\n" % G for G in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)])
        f.write("end nnkpts\n")
    Rs = sorted({k[0] for k in HR})
    with open(os.path.join(path, f"{seed}_hr.dat"), "w") as f:
        f.write(f"synthetic\n{nw}\n{len(Rs)}\n")
        for i in range(0, len(Rs), 15):
            f.write(" ".join(str(ndegen.get(R, 1)) for R in Rs[i : i + 15]) + "\n")
        for R in Rs:
            for b in range(nw):
                for a in range(nw):
                    v = HR.get((R, a, b), 0.0)
                    f.write(f"{R[0]:5d}{R[1]:5d}{R[2]:5d}{a+1:5d}{b+1:5d}{v.real:20.12f}{v.imag:20.12f}\n")
    if wsvec is not None:
        with open(os.path.join(path, f"{seed}_wsvec.dat"), "w") as f:
            f.write("## synthetic with use_ws_distance=T\n")
            for R in Rs:
                for a in range(nw):
                    for b in range(nw):
                        T = wsvec.get((R, a, b), [(0, 0, 0)])
                        f.write(f"{R[0]:5d}{R[1]:5d}{R[2]:5d}{a+1:5d}{b+1:5d}\n{len(T):5d}\n")
                        f.writelines(f"{t[0]:5d}{t[1]:5d}{t[2]:5d}\n" for t in T)


# ==================================================
def model_to_wannier(model, HR, A, t, site_order, translation=None):
    """
    Express H(R) of model in the cell A with origin shift t (x_model = x_wannier U + t).
    If translation is given, the projection centres are wrapped into the home cell of A, and then translated by translation[i].

    Returns:
        - (list) -- projection centres.
        - (list) -- (centre index, l, m) of each Wannier function.
        - (dict) -- H(R) in Wannier index.
        - (list) -- Wannier index of each MultiPie ket.
    """
    U = np.asarray(A) @ np.linalg.inv(model_primitive_vector(model))
    assert np.allclose(U, np.rint(U))
    Ui = np.linalg.inv(np.rint(U))
    ket = [tuple(k) for k in model["full_matrix"]["ket"]]
    pos = np.asarray(list(model.get_ket_site().values()), dtype=float)
    site_pos = {(k[0], k[1]): p for k, p in zip(ket, pos)}
    centers = [(site_pos[s] - t) @ Ui for s in site_order]
    if translation is not None:
        centers = [c - np.floor(c) + np.asarray(n) for c, n in zip(centers, translation)]

    wann, m2w = [], [None] * len(ket)
    for ic, s in enumerate(site_order):
        for lm, i in sorted((W90_LM[k[4]], i) for i, k in enumerate(ket) if (k[0], k[1]) == s):
            m2w[i] = len(wann)
            wann.append((ic, *lm))

    HR_w = {}
    for (n1, n2, n3, m, n), v in HR.items():
        a, b = m2w[m], m2w[n]
        R = (pos[n] + np.array([n1, n2, n3], dtype=float) - pos[m]) @ Ui - (centers[wann[b][0]] - centers[wann[a][0]])
        assert np.allclose(R, np.rint(R))
        key = (tuple(int(i) for i in np.rint(R)), a, b)
        HR_w[key] = HR_w.get(key, 0) + complex(v)

    return centers, wann, HR_w, m2w


# ==================================================
@pytest.fixture(scope="module")
def centred_model(tmp_path_factory):
    """
    Create I4/mmm model with H(R) from random SAMB parameters (once per module).
    """
    topdir = str(tmp_path_factory.mktemp("wannier"))
    name = "Mn2Au"
    model = {
        "model": name,
        "group": 139,
        "cell": {"a": A_LAT, "c": C_LAT},
        "site": {"Mn": (f"[0,0,{Z_MN}]", ["px", "py", "pz"]), "Au": ("[0,0,0]", "s")},
        "bond": [("Mn", "Mn", [1, 2]), ("Mn", "Au", [1])],
        "spinful": False,
        "pdf": {"create": False},
        "qtdraw": {"create": False},
    }
    cwd = os.getcwd()
    try:
        create_model(model, topdir=topdir)
        mm = MaterialModel(topdir=topdir)
        mm.load(name)
    finally:
        changed = os.getcwd()
        os.chdir(cwd)
    assert changed == cwd  # the working directory is not changed (checked here, before the autouse fixture).
    Zr = mm.get_samb_matrix({})["matrix"]
    rng = np.random.default_rng(1)
    parameter = {z: float(rng.normal()) for z in Zr}
    HR = {k: complex(v) for k, v in mm.get_hr(parameter, Zr).items()}

    return topdir, name, mm, parameter, HR


# ==================================================
def run_symcw(topdir, name, ket_wannier=None):
    control = {
        "mode": "symcw",
        "samb": {"model": name},
        "wannier": {"seedname": name, "ket_wannier": ket_wannier or []},
        "output": {"dispersion": {"k_path": None}},
    }
    ma = ModelAnalyzer(topdir)
    ma.analyze(control)
    return ma


# ==================================================
def test_find_lattice_transformation(centred_model):
    _, _, mm, _, _ = centred_model
    U, _ = find_lattice_transformation(A_QE, model_primitive_vector(mm))
    assert U is not None
    assert np.allclose(np.asarray(A_QE), U @ model_primitive_vector(mm))


# ==================================================
@pytest.mark.parametrize(
    "lattice, group, site",
    [("I", 139, "[0,0,0.3]"), ("F", 225, "[0.3,0,0]"), ("C", 65, "[0.3,0,0]"), ("A", 38, "[0,0.3,0.1]"), ("R", 166, "[0,0,0.3]")],
)
def test_model_primitive_vector(lattice, group, site):
    # Cartesian positions from conventional and primitive coordinates agree modulo primitive lattice vectors.
    mm = MaterialModel()
    mm.analyze({"model": "test", "group": group, "cell": {"a": 3.0, "b": 4.0, "c": 5.0}, "site": {"X": (site, "s")}})
    assert mm.group.info.lattice == lattice
    A = np.asarray(mm["unit_vector"], dtype=float)
    Ap = model_primitive_vector(mm)
    n_lattice = {"I": 2, "F": 4, "C": 2, "A": 2, "R": 3}[lattice]
    assert np.isclose(abs(np.linalg.det(Ap)) * n_lattice, abs(np.linalg.det(A)))
    for i in mm["site"]["cell"]["X"]:
        d = np.asarray(i.position, dtype=float) @ A - np.asarray(i.position_primitive, dtype=float) @ Ap
        n = np.linalg.solve(Ap.T, d)
        assert np.allclose(n, np.rint(n), atol=1e-8)


# ==================================================
SHEAR = np.array([[1, 0, 0], [3, 1, 0], [0, -2, 1]])
# large shear amplifies the rounding of projection centres in nnkp (5 decimals).
LARGE_SHEAR = np.array([[1, 0, 0], [40, 1, 0], [0, 40, 1]])


@pytest.mark.parametrize(
    "cell, shift, site_order, translation",
    [
        ("qe", [0.0, 0.0, 0.0], [("Mn", 1), ("Au", 1), ("Mn", 2)], None),
        ("qe", [0.13, -0.21, 0.07], [("Au", 1), ("Mn", 2), ("Mn", 1)], [(0, 0, 0), (1, -1, 0), (0, 0, 2)]),
        ("multipie", [0.0, 0.5, 0.25], [("Mn", 2), ("Au", 1), ("Mn", 1)], [(0, 0, 0), (0, 0, 0), (0, 0, 0)]),
        ("sheared", [0.31, 0.0, -0.17], [("Mn", 1), ("Mn", 2), ("Au", 1)], [(0, 1, 0), (-1, 0, 0), (0, 0, 0)]),
        ("large_shear", [0.0, 0.0, 0.0], [("Mn", 1), ("Mn", 2), ("Au", 1)], [(0, 0, 0), (0, 0, 0), (0, 0, 0)]),
    ],
)
def test_symcw_centred_lattice(centred_model, cell, shift, site_order, translation):
    topdir, name, mm, parameter, HR = centred_model
    Ap = model_primitive_vector(mm)
    A = {"qe": np.asarray(A_QE), "multipie": Ap, "sheared": SHEAR @ Ap, "large_shear": LARGE_SHEAR @ Ap}[cell]
    centers, wann, HR_w, m2w = model_to_wannier(mm, HR, A, np.asarray(shift), site_order, translation)
    path = os.path.join(topdir, name, "wannier")
    shutil.rmtree(path, ignore_errors=True)
    write_wannier(path, name, A, centers, wann, HR_w)

    # H(R) is recovered in the primitive cell of the model.
    win, nnkp = read_win(name, path), read_nnkp(name, path)
    hr_dict, _, _ = read_hr(f"{name}_hr.dat", path)
    HR_m = {
        k: v for k, (v, _) in convert_hr_to_model(hr_dict, map_wannier_to_model(nnkp, win["A"], mm)).items() if abs(v) > 1e-10
    }
    assert HR_m.keys() == {k for k, v in HR.items() if abs(v) > 1e-10}
    assert max(abs(HR_m[k] - HR[k]) for k in HR_m) < 1e-10

    # eigenvalues agree at the same (Cartesian) k points.
    k_w = np.random.default_rng(0).random((5, 3))
    k_m = k_w @ np.linalg.inv(np.rint(A @ np.linalg.inv(Ap))).T
    zero = np.zeros((len(m2w), 3))
    Hk_w = fourier_r_to_k({(R, a, b): v for (R, a, b), v in HR_w.items()}, zero, k_w, s=False)
    Hk_m = fourier_r_to_k({((int(n1), int(n2), int(n3)), m, n): v for (n1, n2, n3, m, n), v in HR.items()}, zero, k_m, s=False)
    assert np.allclose(np.linalg.eigvalsh(Hk_w), np.linalg.eigvalsh(Hk_m), atol=1e-8)

    ma = run_symcw(topdir, name)
    assert ma["wannier"]["multipie_to_wannier"] == m2w
    # primitive cell of analyzer is the same as that used for the mapping.
    assert np.allclose(ma["info"]["A"], Ap)
    assert np.allclose(np.asarray(ma["info"]["A"]) @ np.asarray(ma["info"]["B"]).T, 2 * np.pi * np.eye(3))
    assert max(abs(ma.parameter[z] - v) for z, v in parameter.items()) < 1e-8

    # user-given correspondence gives the same result.
    ket = [tuple(k) for k in mm["full_matrix"]["ket"]]
    w2m = [m2w.index(w) for w in range(len(m2w))]
    ket_wannier = [f"{ket[i][4]}@{ket[i][0]}({ket[i][1]})" for i in w2m]
    ma = run_symcw(topdir, name, ket_wannier)
    assert max(abs(ma.parameter[z] - v) for z, v in parameter.items()) < 1e-8

    # wrong correspondence is rejected.
    ket_wannier[0], ket_wannier[-1] = ket_wannier[-1], ket_wannier[0]
    with pytest.raises(ValueError):
        run_symcw(topdir, name, ket_wannier)


# ==================================================
def test_symcw_centre_not_on_site(centred_model):
    topdir, name, mm, _, HR = centred_model
    centers, wann, HR_w, _ = model_to_wannier(mm, HR, A_QE, np.zeros(3), [("Mn", 1), ("Au", 1), ("Mn", 2)])
    centers[0] = centers[0] + np.array([0.0, 0.0, 0.01])
    path = os.path.join(topdir, name, "wannier")
    shutil.rmtree(path, ignore_errors=True)
    write_wannier(path, name, A_QE, centers, wann, HR_w)

    with pytest.raises(ValueError, match="cannot be mapped"):
        run_symcw(topdir, name)


# ==================================================
def test_symcw_ambiguous_mapping(tmp_path):
    # P1 with two s sites related by a half translation: the mapping depends on the origin shift.
    topdir = str(tmp_path)
    name = "ambiguous"
    model = {
        "model": name,
        "group": 1,
        "cell": {"a": 4.0, "b": 3.0, "c": 5.0},
        "site": {"A": ("[0,0,0]", "s"), "B": ("[1/2,0,0]", "s")},
        "pdf": {"create": False},
        "qtdraw": {"create": False},
    }
    create_model(model, topdir=topdir)
    mm = MaterialModel(topdir=topdir)
    mm.load(name)
    Ap = model_primitive_vector(mm)
    path = os.path.join(topdir, name, "wannier")
    write_wannier(
        path,
        name,
        Ap,
        [np.zeros(3), np.array([0.5, 0, 0])],
        [(0, 0, 1), (1, 0, 1)],
        {((0, 0, 0), 0, 0): 1.0, ((0, 0, 0), 1, 1): 2.0},
    )
    win, nnkp = read_win(name, path), read_nnkp(name, path)

    with pytest.warns(UserWarning, match="several ways"):
        mapping = map_wannier_to_model(nnkp, win["A"], mm)
    assert np.allclose(mapping["t"], 0.0)

    # explicit correspondence: no warning.
    mapping = map_wannier_to_model(nnkp, win["A"], mm, ["s@B(1)", "s@A(1)"])
    assert [mm["full_matrix"]["ket"][i][0] for i in mapping["w2m"]] == ["B", "A"]


# ==================================================
def test_symcw_rounding_of_centres(tmp_path):
    # with a large shear, five-decimal rounding of two centres in nnkp gives opposite errors (about 2e-4 each).
    topdir = str(tmp_path)
    name = "rounding"
    U = np.array([[1, 0, 0], [40, 1, 0], [0, 0, 1]])
    c = [np.array([0.1, 0.1000049, 0.0]), np.array([0.2, 0.1999951, 0.0])]
    model = {
        "model": name,
        "group": 1,
        "cell": {"a": 4.0, "b": 3.0, "c": 5.0},
        "site": {"A": ("[0.100196,0.1000049,0]", "s"), "B": ("[0.199804,0.1999951,0]", "s")},  # c U modulo 1.
        "pdf": {"create": False},
        "qtdraw": {"create": False},
    }
    create_model(model, topdir=topdir)
    mm = MaterialModel(topdir=topdir)
    mm.load(name)
    A = U @ model_primitive_vector(mm)
    path = os.path.join(topdir, name, "wannier")
    HR_w = {((0, 0, 0), 0, 0): 1.0, ((0, 0, 0), 1, 1): 2.0, ((0, 0, 0), 0, 1): 0.5, ((0, 0, 0), 1, 0): 0.5}
    write_wannier(path, name, A, c, [(0, 0, 1), (1, 0, 1)], HR_w)
    win, nnkp = read_win(name, path), read_nnkp(name, path)
    hr_dict, _, _ = read_hr(f"{name}_hr.dat", path)

    # bond (c_B - c_A) U = r_B - r_A + (4, 0, 0).
    for ket_wannier in [None, ["s@A(1)", "s@B(1)"]]:
        mapping = map_wannier_to_model(nnkp, win["A"], mm, ket_wannier)
        assert mapping["U"].tolist() == U.tolist()
        m = mapping["w2m"]
        HR = convert_hr_to_model(hr_dict, mapping)
        assert set(HR) == {(0, 0, 0, m[0], m[0]), (0, 0, 0, m[1], m[1]), (4, 0, 0, m[0], m[1]), (-4, 0, 0, m[1], m[0])}


# ==================================================
def test_symcw_periodic_images_of_site(tmp_path):
    # s and px projections on the same site, but at centres differing by a lattice vector.
    topdir = str(tmp_path)
    name = "image"
    model = {
        "model": name,
        "group": 1,
        "cell": {"a": 4.0, "b": 3.0, "c": 5.0},
        "site": {"A": ("[0.1,0.2,0.3]", ["s", "px"])},
        "pdf": {"create": False},
        "qtdraw": {"create": False},
    }
    create_model(model, topdir=topdir)
    mm = MaterialModel(topdir=topdir)
    mm.load(name)
    path = os.path.join(topdir, name, "wannier")
    c = [np.array([0.1, 0.2, 0.3]), np.array([1.1, 0.2, 0.3])]
    HR_w = {((0, 0, 0), 0, 0): 1.0, ((0, 0, 0), 1, 1): 2.0, ((0, 0, 0), 0, 1): 0.5, ((0, 0, 0), 1, 0): 0.5}
    write_wannier(path, name, model_primitive_vector(mm), c, [(0, 0, 1), (1, 1, 2)], HR_w)
    win, nnkp = read_win(name, path), read_nnkp(name, path)
    hr_dict, _, _ = read_hr(f"{name}_hr.dat", path)

    for ket_wannier in [None, ["s@A(1)", "px@A(1)"], [["A", 1, "s"], ["A", 1, "px"]]]:
        mapping = map_wannier_to_model(nnkp, win["A"], mm, ket_wannier)
        m = mapping["w2m"]
        assert [mm["full_matrix"]["ket"][i][4] for i in m] == ["s", "px"]
        assert set(convert_hr_to_model(hr_dict, mapping)) == {
            (0, 0, 0, m[0], m[0]),
            (0, 0, 0, m[1], m[1]),
            (1, 0, 0, m[0], m[1]),
            (-1, 0, 0, m[1], m[0]),
        }


# ==================================================
def test_symcw_read_ks_not_implemented(centred_model):
    topdir, name, mm, _, HR = centred_model
    centers, wann, HR_w, _ = model_to_wannier(mm, HR, A_QE, np.zeros(3), [("Mn", 1), ("Au", 1), ("Mn", 2)])
    path = os.path.join(topdir, name, "wannier")
    shutil.rmtree(path, ignore_errors=True)
    write_wannier(path, name, A_QE, centers, wann, HR_w)
    control = {"mode": "symcw", "samb": {"model": name}, "wannier": {"seedname": name, "read_KS": True}}
    with pytest.raises(NotImplementedError):
        ModelAnalyzer(topdir).analyze(control)


# ==================================================
def test_point_group_model(tmp_path):
    # the model of a point group has no centring.
    topdir = str(tmp_path)
    model = {"model": "mol", "group": "D3h", "site": {"A": ("[1,0,0]", "s")}, "bond": [("A", "A", [1])]}
    create_model(model | {"pdf": {"create": False}, "qtdraw": {"create": False}}, topdir=topdir)
    mm = MaterialModel(topdir=topdir)
    mm.load("mol")
    assert np.allclose(model_primitive_vector(mm), mm["unit_vector"])
    ma = ModelAnalyzer(topdir)
    ma.analyze({"samb": {"model": "mol", "parameter": {"z1": 1.0}}})
    assert np.allclose(ma["info"]["A"], mm["unit_vector"])


def test_apply_ws_degeneracy(tmp_path):
    # one orbital: the hopping at R=+-2 is on the boundary of the Wigner-Seitz supercell (ndegen=2).
    hr = {(0, 0, 0, 0, 0): 1.0, (1, 0, 0, 0, 0): 0.5, (-1, 0, 0, 0, 0): 0.5, (2, 0, 0, 0, 0): 0.2, (-2, 0, 0, 0, 0): 0.2}
    irvec = [(-2, 0, 0), (-1, 0, 0), (0, 0, 0), (1, 0, 0), (2, 0, 0)]
    ndegen = [2, 1, 1, 1, 2]
    expected = {(0, 0, 0, 0, 0): 1.0, (1, 0, 0, 0, 0): 0.5, (-1, 0, 0, 0, 0): 0.5, (2, 0, 0, 0, 0): 0.1, (-2, 0, 0, 0, 0): 0.1}
    assert apply_ws_degeneracy(hr, irvec, ndegen) == pytest.approx(expected)

    # use_ws_distance: H(R=1) is distributed over R+T, T=(0,0,0) and (0,1,0).
    (tmp_path / "x_wsvec.dat").write_text(
        "## header\n"
        + "".join(f"{R[0]:5d}{R[1]:5d}{R[2]:5d}    1    1\n    1\n    0    0    0\n" for R in irvec if R[0] not in (1, -1))
        + "    1    0    0    1    1\n    2\n    0    0    0\n    0    1    0\n"
        + "   -1    0    0    1    1\n    2\n    0    0    0\n    0   -1    0\n"
    )
    wsvec = read_wsvec("x", str(tmp_path))
    assert wsvec[(1, 0, 0, 0, 0)] == [(0, 0, 0), (0, 1, 0)]
    expected = {
        (0, 0, 0, 0, 0): 1.0,
        (1, 0, 0, 0, 0): 0.25,
        (1, 1, 0, 0, 0): 0.25,
        (-1, 0, 0, 0, 0): 0.25,
        (-1, -1, 0, 0, 0): 0.25,
        (2, 0, 0, 0, 0): 0.1,
        (-2, 0, 0, 0, 0): 0.1,
    }
    assert apply_ws_degeneracy(hr, irvec, ndegen, wsvec) == pytest.approx(expected)

    assert read_wsvec("missing", str(tmp_path)) is None

    # compressed file is read in the same way, and an archive without the file is an error.
    with gzip.open(tmp_path / "y_wsvec.dat.gz", "wt") as f:
        f.write((tmp_path / "x_wsvec.dat").read_text())
    assert read_wsvec("y", str(tmp_path)) == wsvec
    with tarfile.open(tmp_path / "z_wsvec.dat.tar.gz", "w:gz") as tf:
        tf.add(tmp_path / "x_wsvec.dat", arcname="a.dat")
        tf.add(tmp_path / "x_wsvec.dat", arcname="b.dat")
    with pytest.raises(FileNotFoundError):
        read_wsvec("z", str(tmp_path))
    (tmp_path / "bad_wsvec.dat").write_text("    0    0    0    1    1\n    2\n    0    0    0\n")
    with pytest.raises(ValueError, match="invalid format"):
        read_wsvec("bad", str(tmp_path))
    (tmp_path / "dup_wsvec.dat").write_text("    0    0    0    1    1\n    1\n    0    0    0\n" * 2)
    with pytest.raises(ValueError, match="invalid format"):
        read_wsvec("dup", str(tmp_path))

    # inconsistent degeneracies or shifts of R and -R break the Hermiticity.
    with pytest.raises(ValueError, match="Hermiticity"):
        apply_ws_degeneracy(hr, irvec, [1, 1, 1, 1, 2])
    bad = {k: [(0, 0, 0)] for k in hr}
    bad[(1, 0, 0, 0, 0)] = [(0, 1, 0)]
    with pytest.raises(ValueError, match="Hermiticity"):
        apply_ws_degeneracy(hr, irvec, ndegen, bad)


# ==================================================
def test_apply_ws_degeneracy_fourier():
    # two orbitals, complex H(R), orbital-dependent shifts with nT > 1, and shifted elements landing on existing R.
    rng = np.random.default_rng(3)
    irvec = [(i, j, 0) for i in range(-2, 3) for j in range(-1, 2)]
    nw = 2
    H = {}
    for R in irvec:
        mR = tuple(-i for i in R)
        for a in range(nw):
            for b in range(nw):
                if (mR, b, a) in H:
                    H[(R, a, b)] = np.conj(H[(mR, b, a)])
                else:
                    H[(R, a, b)] = complex(rng.normal(), rng.normal()) if R != (0, 0, 0) or a != b else rng.normal()
    ndegen = [1 + (abs(R[0]) == 2) + (abs(R[1]) == 1) for R in irvec]
    shifts = [[(0, 0, 0)], [(0, 0, 0), (-1, 0, 0)], [(0, 1, 0), (0, 0, 0), (1, -1, 0)]]
    wsvec = {}
    for R in irvec:
        mR = tuple(-i for i in R)
        for a in range(nw):
            for b in range(nw):
                if (mR + (b, a)) in wsvec:  # T(-R,b,a) = -T(R,a,b).
                    wsvec[R + (a, b)] = [tuple(-t for t in T) for T in wsvec[mR + (b, a)]]
                else:
                    wsvec[R + (a, b)] = shifts[(sum(map(abs, R)) + 2 * a + b) % 3]
    hr = {R + (a, b): v for (R, a, b), v in H.items()}

    HR = apply_ws_degeneracy(hr, irvec, ndegen, wsvec)
    assert len(HR) < sum(len(wsvec[k]) for k in hr)  # some elements are accumulated.

    for k in rng.random((4, 3)):
        Hk_w90 = np.zeros((nw, nw), dtype=complex)
        for (R, a, b), v in H.items():
            nd = ndegen[irvec.index(R)]
            for T in wsvec[R + (a, b)]:
                Hk_w90[a, b] += v * np.exp(2j * np.pi * np.dot(k, np.add(R, T))) / (nd * len(wsvec[R + (a, b)]))
        Hk = np.zeros((nw, nw), dtype=complex)
        for (n1, n2, n3, a, b), v in HR.items():
            Hk[a, b] += v * np.exp(2j * np.pi * np.dot(k, (n1, n2, n3)))
        assert np.allclose(Hk, Hk_w90)
        assert np.allclose(Hk, Hk.conj().T)


# ==================================================
@pytest.mark.parametrize("use_wsvec", [False, True])
def test_symcw_ws_degeneracy(centred_model, use_wsvec):
    # hr.dat with ndegen > 1, and (with wsvec) the elements of R0 and -R0 written at another image R0+S.
    topdir, name, mm, parameter, HR = centred_model
    centers, wann, HR_w, _ = model_to_wannier(mm, HR, A_QE, np.zeros(3), [("Mn", 1), ("Au", 1), ("Mn", 2)])
    Rs = sorted({k[0] for k, v in HR_w.items() if abs(v) > 1e-10})
    ndegen = {R: 1 + sum(map(abs, R)) % 3 for R in Rs}
    ndegen = {R: max(n, ndegen.get(tuple(-i for i in R), 1)) for R, n in ndegen.items()}  # n(R) = n(-R).
    assert max(ndegen.values()) > 1

    R0 = next(R for R in Rs if R != (0, 0, 0))
    S = (5, 0, 0)
    HR_file, wsvec = {}, {}
    for (R, a, b), v in HR_w.items():
        n = ndegen.get(R, 1)
        if use_wsvec and R in (R0, tuple(-i for i in R0)):
            sign = 1 if R == R0 else -1
            Rf = tuple(r + sign * s for r, s in zip(R, S))
            ndegen[Rf] = n
            wsvec[(Rf, a, b)] = [tuple(-sign * s for s in S)]
            R = Rf
        HR_file[(R, a, b)] = v * n

    path = os.path.join(topdir, name, "wannier")
    shutil.rmtree(path, ignore_errors=True)
    write_wannier(path, name, A_QE, centers, wann, HR_file, ndegen, wsvec if use_wsvec else None)

    ma = run_symcw(topdir, name)
    assert max(abs(ma.parameter[z] - v) for z, v in parameter.items()) < 1e-8


# ==================================================
def test_wannier_mode(centred_model):
    # wannier mode with and without model, for seedname.win in the primitive cell of QE (ibrav=7),
    # with the atoms of seedname.win not in the order of projection centres, and seedname different from the model name.
    topdir, name, mm, parameter, HR = centred_model
    seed = name + "_w"
    site_order = [("Mn", 1), ("Au", 1), ("Mn", 2)]
    species = ["Mn", "Au", "Mn"]
    centers, wann, HR_w, _ = model_to_wannier(mm, HR, A_QE, np.zeros(3), site_order)
    atoms = [("Au", centers[1]), ("Mn", centers[2]), ("Mn", centers[0])]
    path = os.path.join(topdir, seed, "wannier")
    shutil.rmtree(path, ignore_errors=True)
    write_wannier(path, seed, A_QE, centers, wann, HR_w, atoms=atoms)

    # the same Cartesian k path, in the reciprocal basis of the model and of seedname.win.
    U = np.rint(np.asarray(A_QE) @ np.linalg.inv(model_primitive_vector(mm)))
    k_m = {"Γ": [0, 0, 0], "X": [0.5, 0, 0], "Z": [0.5, 0.5, -0.5], "P": [0.25, 0.25, 0.25]}
    k_w = {k: (np.asarray(v) @ U.T).tolist() for k, v in k_m.items()}

    ma = ModelAnalyzer(topdir)  # the same analyzer is used for all runs.

    def run(mode, model, k_point, tb_gauge=True, samb_parameter=None):
        control = {
            "mode": mode,
            "wannier": {"seedname": seed},
            "output": {
                "fourier": {"tb_gauge": tb_gauge},
                "dispersion": {"k_path": "Γ-X-P-Z|X-Γ", "k_point": {k: str(v) for k, v in k_point.items()}},
            },
        }
        if model:
            control["samb"] = {"model": name}
            if samb_parameter:
                control["samb"]["parameter"] = samb_parameter
        ma.analyze(control)
        out = name if model else seed  # output directory.
        disp = np.loadtxt(os.path.join(topdir, out, "output", f"{out}_dispersion.txt"))
        assert ma["info"]["name"] == out
        HR = {k: complex(v[0] if isinstance(v, tuple) else v) for k, v in ma.HR.items()}
        return copy.deepcopy(dict(ma)), HR, disp

    # model-free -> model -> samb with the same analyzer: no state is carried over.
    without_model, HR_without, disp_without = run("wannier", False, k_w)
    with_model, _, disp_with = run("wannier", True, k_m)
    _, _, disp_samb = run("samb", True, k_m, samb_parameter=parameter)
    _, _, disp_gauge = run("wannier", False, k_w, tb_gauge=False)
    symcw, _, disp_symcw = run("symcw", True, k_m)
    # the same k path (Cartesian distance) and bands.
    for disp in [disp_with, disp_without, disp_gauge, disp_samb]:
        assert np.allclose(disp, disp_symcw, atol=1e-8)
    B_m = 2 * np.pi * np.linalg.inv(model_primitive_vector(mm)).T
    segments = [("Γ", "X"), ("X", "P"), ("P", "Z"), ("X", "Γ")]
    length = sum(np.linalg.norm((np.asarray(k_m[e]) - np.asarray(k_m[s])) @ B_m) for s, e in segments)
    assert disp_symcw[:, 0].max() == pytest.approx(length)

    # with model: lattice of the model. without model: those of seedname.win.
    assert np.allclose(with_model["info"]["A"], model_primitive_vector(mm))
    assert np.allclose(without_model["info"]["A"], A_QE)
    w2m = without_model["wannier"]["wannier_to_multipie"]
    HR_ma = {k: v for k, v in HR_without.items() if abs(v) > 1e-10}
    assert HR_ma == pytest.approx({(*R, w2m[a], w2m[b]): v for (R, a, b), v in HR_w.items() if abs(v) > 1e-10})
    # ket names and positions (projection centres) of each Wannier function; the names are those of the model.
    ket = without_model["wannier"]["ket"]
    assert sorted(ket) == sorted(symcw["wannier"]["ket"])
    for w, (ic, _, _) in enumerate(wann):
        assert ket[w2m[w]].split("@")[1].startswith(species[ic] + "(")
        assert np.allclose(without_model["wannier"]["atoms_frac"][w2m[w]], centers[ic])

    # projection centres in other periodic images, with the atoms of seedname.win unchanged:
    # the positions are the projection centres, and the dispersion is the same in both gauges.
    centers_s, wann, HR_w, _ = model_to_wannier(mm, HR, A_QE, np.zeros(3), site_order, [(0, 1, 0), (-1, 0, 2), (0, 0, 0)])
    assert not np.allclose(centers_s, centers)
    shutil.rmtree(path, ignore_errors=True)
    write_wannier(path, seed, A_QE, centers_s, wann, HR_w, atoms=atoms)
    for tb_gauge in [True, False]:
        d, _, disp = run("wannier", False, k_w, tb_gauge)
        assert np.allclose(disp, disp_symcw, atol=1e-8)
        w2m = d["wannier"]["wannier_to_multipie"]
        for w, (ic, _, _) in enumerate(wann):
            assert np.allclose(d["wannier"]["atoms_frac"][w2m[w]], centers_s[ic])

    # explicit ket_wannier (not the automatic correspondence) without model: H(R), positions and H(k) in the
    # tight-binding gauge follow it.
    d_auto = without_model
    w2m_auto = d_auto["wannier"]["wannier_to_multipie"]
    ket_auto = d_auto["wannier"]["ket"]
    perm = list(range(len(ket_auto)))
    perm[0], perm[-1] = perm[-1], perm[0]  # exchange kets on different centres.
    ket_wannier = [ket_auto[perm[w2m_auto[w]]] for w in range(len(ket_auto))]
    ma.analyze(
        {"mode": "wannier", "wannier": {"seedname": seed, "ket_wannier": ket_wannier}, "output": {"dispersion": {"k_path": None}}}
    )
    w2m = ma["wannier"]["wannier_to_multipie"]
    m2w = ma["wannier"]["multipie_to_wannier"]
    assert [ket_wannier[w] for w in m2w] == ma["wannier"]["ket"]
    atom = np.asarray(ma["wannier"]["atoms_frac"])
    for m, w in enumerate(m2w):
        assert np.allclose(atom[m], centers_s[wann[w][0]], atol=1e-5)
    k = np.random.default_rng(2).random((4, 3))
    HR_m = {((n1, n2, n3), m, n): complex(v[0]) for (n1, n2, n3, m, n), v in ma.HR.items()}
    Hk = fourier_r_to_k(HR_m, atom, k, s=True)
    Hk_ref = np.zeros((len(k), len(wann), len(wann)), dtype=complex)
    for (R, a, b), v in HR_w.items():
        r = np.asarray(R) + centers_s[wann[b][0]] - centers_s[wann[a][0]]
        Hk_ref[:, a, b] += v * np.exp(2j * np.pi * k @ r)
    assert np.allclose(Hk, Hk_ref[:, m2w][:, :, m2w], atol=1e-3)  # centres in seedname.nnkp have 5 decimals.
    with pytest.raises(ValueError, match="permutation"):
        ma.analyze({"mode": "wannier", "wannier": {"seedname": seed, "ket_wannier": ket_wannier[:-1]}})

    # without projections in seedname.nnkp (auto_projections).
    nnkp = os.path.join(path, f"{seed}.nnkp")
    text = open(nnkp).read()
    i, j = text.index("begin projections"), text.index("end projections") + len("end projections")
    with open(nnkp, "w") as f:
        f.write(text[:i] + f"begin auto_projections\n {len(wann)}\n 0\nend auto_projections" + text[j:])
    with pytest.raises(ValueError, match="projection information is missing"):
        ma.analyze({"mode": "wannier", "wannier": {"seedname": seed}})
    with open(nnkp, "w") as f:
        f.write(text)

    # default k path for seedname.win (body-centred tetragonal, in the reciprocal basis of A_QE).
    with warnings.catch_warnings():
        warnings.simplefilter("error", seekpath.SupercellWarning)
        ma.analyze({"mode": "wannier", "wannier": {"seedname": seed}, "grid": (10, 10, 10)})
    disp = ma["output"]["dispersion"]
    # the same path and Cartesian k points as in the standardised cell of seekpath (no rotation for this cell).
    info = seekpath.get_path((A_QE, [p for _, p in atoms], [2, 1, 1]))
    assert np.allclose(info["rotation_matrix"], np.eye(3))
    assert disp["k_path"] == _join_kpath(info["path"])
    B = np.asarray(ma["info"]["B"])
    B_std = np.asarray(info["reciprocal_primitive_lattice"])
    assert len(disp["k_point"]) == len(info["point_coords"])
    for label, k in disp["k_point"].items():
        k_std = np.asarray(info["point_coords"]["GAMMA" if label == "Γ" else label])
        k = np.asarray([float(i) for i in k.strip("[]").split(",")])
        assert np.allclose(k @ B, k_std @ B_std)


# ==================================================
@pytest.mark.parametrize(
    "A, cart, ket",
    [
        # primitive cell (QE ibrav=7) of I4/mmm: sublattices as in a model.
        (
            A_QE,
            {("Mn", 1): [0, 0, Z_MN * C_LAT], ("Mn", 2): [0, 0, -Z_MN * C_LAT], ("Au", 1): [0, 0, 0]},
            ["s@Au(1)", "s@Mn(1)", "s@Mn(2)"],
        ),
        # conventional cell of bcc: atoms related by the centring translation.
        (np.diag([3.0, 3.0, 3.0]), {("Fe", 1): [0, 0, 0], ("Fe", 2): [1.5, 1.5, 1.5]}, ["s@Fe(1)", "s@Fe(2)"]),
        # cubic cell of a rhombohedral structure, species interleaved.
        (
            np.diag([3.0, 3.0, 3.0]),
            {("Fe", 1): [0, 0, 0], ("Co", 1): [0.75, 0.75, 0.75], ("Fe", 2): [1.5, 1.5, 1.5], ("Co", 2): [2.25, 2.25, 2.25]},
            ["s@Co(1)", "s@Co(2)", "s@Fe(1)", "s@Fe(2)"],
        ),
        # the same element in two Wyckoff orbits.
        (
            np.diag([3.0, 3.0, 3.0]),
            {("Fe", 1): [0, 0, 0], ("Fe", 2): [1.5, 1.5, 1.5], ("O", 1): [1.5, 0, 0]},
            ["s@Fe1(1)", "s@Fe2(1)", "s@O(2)"],
        ),
    ],
)
def test_ket_names_unique(A, cart, ket):
    assert create_ket_wannier_multipie(wannier_info(A, cart))[2] == ket


# ==================================================
def test_ket_off_atom_label_not_used_by_atoms():
    # atoms labelled X in two Wyckoff orbits (X1, X2) and a projection centre not on an atom: the latter is X3.
    cart = {("X", 1): [0, 0, 0], ("X", 2): [1.5, 1.5, 1.5], ("O", 1): [1.5, 0, 0]}
    info = wannier_info(np.diag([3.0, 3.0, 3.0]), cart)
    info["atom_pos_r"] = info["atom_pos_r"] + [[0.25, 0, 0]]
    info.update(nw2n=[0, 1, 2, 3], nw2l=[0] * 4, nw2m=[1] * 4, nw2r=[1] * 4, nw2s=[0] * 4)
    w2m, _, ket, _, _ = create_ket_wannier_multipie(info)
    assert [ket[w2m[w]] for w in range(4)] == ["s@X1(1)", "s@X2(1)", "s@O(2)", "s@X3(1)"]


# ==================================================
def test_ket_wannier_not_unique(tmp_path):
    # two s projections on one atom have the same name: ket_wannier cannot be used, but the automatic correspondence works.
    topdir, seed = str(tmp_path), "two_s"
    A = np.diag([3.0, 3.0, 3.0])
    HR = {((0, 0, 0), 0, 0): 1.0, ((0, 0, 0), 1, 1): -1.0, ((1, 0, 0), 0, 0): 0.1, ((-1, 0, 0), 0, 0): 0.1}
    write_wannier(os.path.join(topdir, seed, "wannier"), seed, A, [[0, 0, 0]], [(0, 0, 1), (0, 0, 1)], HR, species=["Fe"])
    ma = ModelAnalyzer(topdir)
    ma.analyze({"mode": "wannier", "wannier": {"seedname": seed}, "output": {"dispersion": {"k_path": None}}})
    assert ma["wannier"]["ket"] == ["s@Fe(1)", "s@Fe(1)"]
    with pytest.raises(ValueError, match="not unique"):
        ma.analyze({"mode": "wannier", "wannier": {"seedname": seed, "ket_wannier": ["s@Fe(1)", "s@Fe(1)"]}})


# ==================================================
def test_wannier_mode_bond_centre(tmp_path):
    # model-free wannier mode with a projection centre on a bond (not on an atom).
    topdir, seed = str(tmp_path), "bond"
    A = np.diag([3.0, 6.0, 6.0])
    centers = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]
    HR = {((0, 0, 0), 0, 0): 1.0, ((0, 0, 0), 1, 1): -1.0, ((0, 0, 0), 0, 1): 0.3, ((0, 0, 0), 1, 0): 0.3}
    HR |= {((-1, 0, 0), 0, 1): 0.3, ((1, 0, 0), 1, 0): 0.3}
    write_wannier(os.path.join(topdir, seed, "wannier"), seed, A, centers, [(0, 0, 1), (1, 0, 1)], HR, atoms=[("Fe", centers[0])])
    ma = ModelAnalyzer(topdir)
    ma.analyze({"mode": "wannier", "wannier": {"seedname": seed}, "grid": (10, 1, 1)})
    w2m = ma["wannier"]["wannier_to_multipie"]
    assert [ma["wannier"]["ket"][w2m[w]] for w in range(2)] == ["s@Fe(1)", "s@X1(1)"]
    assert np.allclose([ma["wannier"]["atoms_frac"][w2m[w]] for w in range(2)], centers)
    assert os.path.isfile(os.path.join(topdir, seed, "output", f"{seed}_dispersion.txt"))
