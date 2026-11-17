# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Conversion of a model (.pkl) of MultiPie Ver.2 into Ver.3.

The model is created again with Ver.3 from the input rebuilt from the model of Ver.2, and the weights z_j of the SAMB
of Ver.2 are converted into those of Ver.3 by projection.
"""

import re
import copy
import gzip
import pickle
from fractions import Fraction
import numpy as np
from scipy import sparse

from multipie.core.group import Group
from multipie.core.material_model import MaterialModel
from multipie.core.default_model import default_model
from multipie.util.util_crystal import get_cell_info, create_igrid
from multipie.util.util_material_model import (
    get_basis_type,
    parse_representative_site,
    create_site_grid,
    create_site_so,
    parse_representative_bond,
    parse_samb_select,
    SEARCH_CELL_RANGE,
    RENAMED_ORBITAL,
    REMOVED_MODEL_KEY,
)

# tolerance of the residual of the projection (the SAMB are normalized to 1).
CONVERT_TOL = 1e-8
# tolerance to identify the distance of a bond (Angstrom, or the unit of the cell).
DISTANCE_TOL = 1e-5
# denominators tried to write a coordinate as a fraction.
MAX_DENOMINATOR = 24


# ==================================================
def model_version(filename):
    """
    Version of MultiPie that created a model file.

    Args:
        filename (str): model file (.pkl).

    Returns:
        - (str) -- version, e.g., "2.3.4".
        - (dict) -- data of the model.
    """
    with gzip.open(filename, "rb") as f:
        dic = pickle.load(f)
    data = dic["data"] if isinstance(dic, dict) and "data" in dic else dic
    return str(data.get("version", "0")), data


# ==================================================
def major_version(version):
    """
    Major version.

    Args:
        version (str): version, e.g., "2.3.4".

    Returns:
        - (int) -- major version (0 if unknown).
    """
    m = re.match(r"(\d+)", str(version))
    return 0 if m is None else int(m[1])


# ==================================================
def _coordinate(x):
    """
    Coordinate as a str, a fraction if it is close to one with a small denominator.

    Args:
        x (float): coordinate.

    Returns:
        - (str) -- coordinate.
    """
    f = Fraction(float(x)).limit_denominator(MAX_DENOMINATOR)
    if abs(float(f) - x) < 1e-6:
        return str(f)
    return repr(round(float(x), 10))


# ==================================================
def _orbital(orbital, basis_type):
    """
    Orbitals of a site in the input format of Ver.3.

    Args:
        orbital ([[str]]): orbitals of the representative site of Ver.2, in the order of rank.
        basis_type (str): basis type, "jml/lgs/lg".

    Returns:
        - ([str]) -- orbitals.
        - ([(str, str)]) -- renamed orbitals, [(old, new)].
    """
    lst = []
    renamed = []
    for rank, block in enumerate(orbital):
        for o in block:
            if basis_type == "lgs":  # (px,u) -> px, the spin is given by "spinful".
                o = o.strip("()").split(",")[0]
            elif basis_type == "jml":  # (j,m) of rank 1 -> (j,p), a full multiplet of j.
                o = f"({o.strip('()').split(',')[0]},{'spdf'[rank]})"
            if o in RENAMED_ORBITAL:
                renamed.append((o, RENAMED_ORBITAL[o]))
                o = RENAMED_ORBITAL[o]
            if o not in lst:
                lst.append(o)
    return lst, list(dict.fromkeys(renamed))


# ==================================================
def renamed_orbital(orbital):
    """
    Orbital with the f orbitals renamed in Ver.3, e.g., "(f3x,u)" -> "(fx,u)".

    Args:
        orbital (str): orbital of Ver.2.

    Returns:
        - (str) -- orbital of Ver.3.
    """
    return re.sub(r"\b(" + "|".join(RENAMED_ORBITAL) + r")\b", lambda m: RENAMED_ORBITAL[m[1]], orbital)


# ==================================================
def rebuild_model_input(data):
    """
    Rebuild the model input from the data of a model of Ver.2.

    Args:
        data (dict): data of the model of Ver.2.

    Returns:
        - (dict) -- model input (the neighbors of "bond" are those of Ver.2).
        - (dict) -- notes, {"renamed_orbital": {site: [(old, new)]}, "removed_key": {key: value}}.
    """
    basis_type = data["basis_type"]
    site = {}
    renamed = {}
    for name, r in data["site"]["representative"].items():
        orbital, ren = _orbital(r.orbital, basis_type)
        site[name] = ("[" + ",".join(_coordinate(x) for x in r.position) + "]", orbital)
        if ren:
            renamed[name] = ren

    bond = {}
    for b in data["bond"]["representative"].values():
        bond.setdefault((b.tail, b.head), set()).add(b.neighbor)
    bond = [(t, h, sorted(n)) for (t, h), n in bond.items()]

    def option(key, src):  # the keys of Ver.3 with a value other than the default.
        return {k: v for k, v in data.get(src, {}).items() if k in default_model[key] and v != default_model[key][k]}

    irrep = list(Group(data["group"]).character["table"].keys())

    def select(key):  # the conditions other than the default.
        default = parse_samb_select(copy.deepcopy(default_model[key]), irrep)
        return {k: v for k, v in data[key].items() if k in ("X", "l", "Gamma", "s") and v != default.get(k)}

    model = {
        "model": data["model"],
        "group": data["group"],
        "cell": {k: float(v) for k, v in data["cell_info"]["cell"].items()},
        "spinful": basis_type in ("lgs", "jml"),
        "site": site,
        "bond": bond,
        "SAMB_select": select("SAMB_select"),
        "atomic_select": select("atomic_select"),
        "site_select": select("site_select"),
        "bond_select": select("bond_select"),
        "toroidal_priority": data["toroidal_priority"],
        "qtdraw": option("qtdraw", "qtdraw_prop"),
        "pdf": option("pdf", "pdf_ctrl"),
    }
    model = {k: v for k, v in model.items() if v != {} or k in ("site",)}
    removed = {k: data[k] for k in REMOVED_MODEL_KEY if k in data}

    return model, {"renamed_orbital": renamed, "removed_key": removed}


# ==================================================
def bond_shell_distance(model, tail, head, max_neighbor):
    """
    Distances of the neighbor shells of a bond in Ver.3.

    Args:
        model (dict): model input.
        tail (str): tail site.
        head (str): head site.
        max_neighbor (int): max. neighbor (all shells for a point group).

    Returns:
        - ([float]) -- distance of each neighbor, 1, 2, ...
    """
    group = Group(model["group"])
    cell_info = get_cell_info(group.info.crystal, model["cell"])
    G = cell_info["G"][0:3, 0:3]
    basis_type = get_basis_type(model["site"], model["spinful"])
    basis_info = {k: group.atomic_basis(k) for k in ["jml", "lgs", "lg"]}
    igrid = None if group.is_point_group else create_igrid(SEARCH_CELL_RANGE)
    site = {k: model["site"][k] for k in dict.fromkeys([tail, head])}
    site_dict = parse_representative_site(group, site, basis_type, basis_info)
    site_grid = create_site_grid(site_dict, igrid)
    site_so = create_site_so(group, site_dict)
    bond = [(tail, head, list(range(1, max_neighbor + 1)))]
    rep = parse_representative_bond(group, G, site_grid, site_so, site_dict, bond, False)["representative"]
    distance = {}
    for b in rep.values():
        distance.setdefault(b.neighbor, b.distance)
    return [distance[n] for n in sorted(distance)]


# ==================================================
def remap_neighbor(model, data, max_neighbor=2000):
    """
    Choose the neighbors of Ver.3 by the distances of the bonds of Ver.2.

    Args:
        model (dict): model input rebuilt by rebuild_model_input (modified in place).
        data (dict): data of the model of Ver.2.
        max_neighbor (int, optional): max. neighbor to search.

    Returns:
        - (dict) -- changed neighbors, {(tail, head): (neighbors of Ver.2, neighbors of Ver.3)}.

    Note:
        - Ver.2 searched the bonds in a finite range of cells (search_cell_range), which may miss a part of the neighbor
          shells and shift the numbering of the farther neighbors. Ver.3 finds all the shells, so the neighbors are chosen
          as the shells at the distances of the bonds of Ver.2.
    """
    dist2 = {}
    for b in data["bond"]["representative"].values():
        dist2.setdefault((b.tail, b.head), set()).add(b.distance)

    changed = {}
    bond = []
    for tail, head, neighbor in model["bond"]:
        d2 = sorted(dist2[(tail, head)])
        n = max(neighbor)
        while True:
            shell = bond_shell_distance(model, tail, head, n)
            if Group(model["group"]).is_point_group or shell[-1] >= d2[-1] - DISTANCE_TOL or len(shell) < n:
                break
            if n >= max_neighbor:
                raise ValueError(f"bond {tail};{head}: the bonds of Ver.2 are not found within {max_neighbor} neighbors.")
            n = min(2 * n, max_neighbor)
        new = [no + 1 for no, d in enumerate(shell) if any(abs(d - i) < DISTANCE_TOL for i in d2)]
        missing = [i for i in d2 if all(abs(d - i) >= DISTANCE_TOL for d in shell)]
        if missing:
            raise ValueError(f"bond {tail};{head}: no neighbor shell of Ver.3 at the distance(s) {missing} of Ver.2.")
        if new != neighbor:
            changed[(tail, head)] = (neighbor, new)
        bond.append((tail, head, new))
    model["bond"] = bond

    return changed


# ==================================================
def samb_matrix_rows(matrices, column):
    """
    SAMB as the rows of a sparse matrix.

    Args:
        matrices (dict): SAMB matrices, dict[zj, dict[(n1, n2, n3, m, n), value]].
        column (dict): column index of (n1, n2, n3, m, n), updated in place.

    Returns:
        - ([str]) -- names of the SAMB (rows).
        - (tuple) -- (data, row, col) of the sparse matrix.
    """
    names = list(matrices.keys())
    data, row, col = [], [], []
    for i, z in enumerate(names):
        for Rmn, v in matrices[z].items():
            j = column.setdefault(Rmn, len(column))
            data.append(complex(v))
            row.append(i)
            col.append(j)
    return names, (data, row, col)


# ==================================================
class ModelConverter:
    """
    Projection of the SAMB of a model of Ver.2 onto those of the model of Ver.3.
    """

    # ==================================================
    def __init__(self, data2, model3):
        """
        Projection of the SAMB of Ver.2 onto those of Ver.3.

        Args:
            data2 (dict): data of the model of Ver.2.
            model3 (MaterialModel): model of Ver.3.

        Raises:
            ValueError: if the kets differ, or the SAMB of Ver.2 are not in the space of those of Ver.3.
        """
        mm2 = MaterialModel()
        mm2.clear()
        mm2.update(data2)

        ket2, ket3 = mm2["full_matrix"]["ket"], model3["full_matrix"]["ket"]
        if [list(k[:4]) + [renamed_orbital(k[4])] for k in ket2] != [list(k) for k in ket3]:
            raise ValueError(f"the kets of Ver.2 and Ver.3 differ:\n  Ver.2: {ket2}\n  Ver.3: {ket3}")

        Z2 = mm2._combined_samb_matrix(mm2["combined_id"], fmt="value", bond=False)
        Z3 = model3._combined_samb_matrix(model3["combined_id"], fmt="value", bond=False)
        column = {}
        self.name2, coo2 = samb_matrix_rows(Z2, column)
        self.name3, coo3 = samb_matrix_rows(Z3, column)
        A2 = sparse.csr_matrix((coo2[0], (coo2[1], coo2[2])), shape=(len(self.name2), len(column)))
        A3 = sparse.csr_matrix((coo3[0], (coo3[1], coo3[2])), shape=(len(self.name3), len(column)))
        self._A2, self._A3 = A2, A3

        # the SAMB are orthonormal under Tr[Z_i^dagger Z_j], unless the orbitals of a site are not a complete set of
        # irreps of its site symmetry, e.g., only a part of the f orbitals: then some SAMB vanish or are linearly dependent.
        G3 = A3.conj() @ A3.T - sparse.identity(len(self.name3))
        self.orthonormal = bool(G3.nnz == 0 or abs(G3).max() < CONVERT_TOL)

        norm2 = np.asarray(abs(A2).power(2).sum(axis=1)).ravel()  # 1, or 0 for a vanishing SAMB.
        if self.orthonormal:
            # P[i, j] = Tr[Z3_i^dagger Z2_j]: z_j(Ver.3) = sum_k P[j, k] z_k(Ver.2).
            self._P = (A3.conj() @ A2.T).tocsc()
            proj2 = np.asarray(abs(self._P).power(2).sum(axis=0)).ravel()
            residual = np.abs(norm2 - proj2)
        else:
            # least squares (min. norm): A3^T z3 = A2^T z2.
            A3d = A3.T.toarray()
            X, *_ = np.linalg.lstsq(A3d, A2.T.toarray(), rcond=None)
            self._P = sparse.csc_matrix(X)
            residual = np.abs(A3d @ X - A2.T.toarray()).max(axis=0) if len(self.name2) > 0 else np.zeros(0)

        # zero if Z2_j is in the space spanned by the SAMB of Ver.3.
        self.span_residual = float(residual.max(initial=0.0))
        if self.span_residual > CONVERT_TOL:
            j = int(np.argmax(residual))
            raise ValueError(
                f"the SAMB of Ver.2 are not in the space of those of Ver.3 (residual {self.span_residual:.3g} for "
                f"{self.name2[j]} of Ver.2)."
            )

    # ==================================================
    def convert(self, parameter):
        """
        Convert the weights z_j of Ver.2 into those of Ver.3.

        Args:
            parameter (dict): weights of Ver.2, dict[zj, float].

        Returns:
            - (dict) -- weights of Ver.3, dict[zj, float] (only the finite ones).
            - (float) -- max. residual of the matrix elements of H = sum_j z_j Z_j.
        """
        unknown = [z for z in parameter if z not in self.name2]
        if unknown:
            raise ValueError(f"{unknown} are not SAMB of the model of Ver.2.")
        idx = {z: i for i, z in enumerate(self.name2)}
        z2 = np.zeros(len(self.name2))
        for z, v in parameter.items():
            z2[idx[z]] = float(v)

        z3 = self._P @ z2
        if np.abs(z3.imag).max(initial=0.0) > CONVERT_TOL:
            raise ValueError("the converted weights are not real.")
        z3 = z3.real

        H2 = self._A2.T @ z2
        H3 = self._A3.T @ z3
        residual = float(np.abs(H2 - H3).max(initial=0.0))
        if residual > CONVERT_TOL:
            raise ValueError(f"the weights are not converted: residual of H is {residual:.3g}.")

        scale = max(np.abs(z3).max(initial=0.0), 1.0)
        # 12 significant digits, enough within the round-off of the projection.
        return {z: float(f"{v:.12g}") for z, v in zip(self.name3, z3) if abs(v) > CONVERT_TOL * scale}, residual
