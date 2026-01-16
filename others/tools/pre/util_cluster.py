import numpy as np
import sympy as sp

from multipie.util.util import str_to_sympy, replace
from multipie.util.util_crystal import shift_bond


# ==================================================
def remove_equivalent_vectors(vectors):
    """
    Remove equivalnet nondirectional vectors.

    Args:
        vectors (ndarray): vectors.

    Returns:
        - (ndarray) -- removed vectors.
    """

    def is_equivalent(v1, v2):
        return all((a - b).equals(0) for a, b in zip(v1, v2)) or all((a + b).equals(0) for a, b in zip(v1, v2))

    result = []
    for v in vectors:
        if any(is_equivalent(v, w) for w in result):
            continue
        result.append(v)
    result = np.asarray(result, dtype=object)

    return result


# ==================================================
def create_unique_bonds(bonds, plus_set):
    """
    Create unique cluster bonds.

    Args:
        bonds (ndarrays): bonds.
        plus_set (ndarray): plus set (None for point group).

    Returns:
        - (tuple) -- unique bonds.
    """
    sub = {"X": 7, "Y": 3, "Z": 1}  # dummy values in order to evaluate sign of each component.

    def is_regular_vector(v):
        vv = replace(v, sub)
        return next((i > 0 for i in vv if i != 0), True)

    def regular_bond(bond):
        if not is_regular_vector(bond[0:3]):
            bond = np.concatenate([-bond[0:3], bond[3:6]])
        if plus_set is not None:
            bond = shift_bond(bond)
        return bond

    if plus_set is not None:
        vs, cs = bonds[:, 0:3], bonds[:, 3:6]
        cs = np.concatenate([cs + i for i in plus_set])
        # cs = shift_site(cs)
        vs = np.tile(vs, (len(plus_set), 1))
        bonds = np.hstack((vs[:, None], cs[:, None])).reshape(-1, 6)

    bonds = tuple(sorted(list(map(tuple, [regular_bond(b) for b in bonds])), key=lambda x: str(x)))

    return bonds


# ==================================================
def remove_sub_vector(vectors, sub_vector):
    """
    Remove sub-level vectors.

    Args:
        vectors (ndarray): vectors.
        sub_vector (dict): sub-level dict.

    Returns:
        - (ndarray) -- removed vectors.
    """
    result = []
    seen = set()

    for i, v in enumerate(vectors):
        vs = str(v.tolist()).replace(" ", "").replace("*", "")
        if vs in seen:
            continue
        result.append(v)
        seen.add(vs)
        d = set(sub_vector.get(vs, []))
        for u in vectors[i + 1 :]:
            us = str(u.tolist()).replace(" ", "").replace("*", "")
            if us in d:
                seen.add(us)
    result = np.asarray(result)

    return result


# ==================================================
def create_transformed_vector_first(pg_so, mapping):
    """
    Create transformed vectors at 1st Wyckoff site.

    Args:
        pg_so (ndarray): point-group symmetry operations.
        mapping (list): SO mapping.

    Returns:
        - (ndarray) -- transformed vectors.
    """
    v = str_to_sympy("[X,Y,Z]")
    vectors = pg_so[mapping[0]] @ v
    return vectors


# ==================================================
def create_unique_vector_source(vectors, vector_source, sub_vector):
    """
    Create unique vector source.

    Args:
        vectors (ndarray): transformed vectors.
        vector_source (list): vector source candidates.
        sub_vector (dict): sub-level vector dict.

    Returns:
        - (dict) -- unique vector source.
    """
    unique_vector_source = {}
    for source in vector_source:
        sub = dict(zip(["X", "Y", "Z"], source))
        subvectors = replace(vectors, sub)
        subvectors = remove_equivalent_vectors(subvectors)
        n = len(subvectors)
        unique_vector_source[n] = unique_vector_source.get(n, []) + [source]

    unique_vector_source = dict(
        sorted({n: remove_sub_vector(np.asarray(v), sub_vector) for n, v in unique_vector_source.items()}.items())
    )

    return unique_vector_source


# ==================================================
def create_unique_vector_set(unique_vector_source, pg_so, mapping):
    """
    Create unique vector set at 1st Wyckoff site.

    Args:
        unique_vector_source (dict): unique vector source.
        pg_so (ndarray): point-group symmetry operations.
        mapping (list): SO mapping.

    Returns:
        - (dict) -- unique vector set.
    """
    unique_vector_set = {
        n: np.asarray([remove_equivalent_vectors(pg_so[mapping[0]] @ vi) for vi in v]) for n, v in unique_vector_source.items()
    }

    return unique_vector_set


# ==================================================
def remove_equivalent_bond_set(bonds, plus_set):
    """
    Remove equivalent bond set.

    Args:
        bonds (ndarray): bonds.
        plus_set (ndarray): plus set (None for point group).

    Returns:
        - (ndarray) -- unique bond set.
    """
    X, Y, Z = str_to_sympy("[X,Y,Z]")
    x, y, z = str_to_sympy("[x,y,z]")
    sub_pg = [
        {"x": -x},
        {"y": -y},
        {"z": -z},
        {"x": x, "y": -y},
        {"x": -x, "y": y},
        {"x": -x, "y": -y},
        {"x": y, "y": x},
        {"x": y, "y": -x},
        {"x": -y, "y": x},
        {"x": -y, "y": -x},
        {"X": -X},
        {"Y": -Y},
        {"Z": -Z},
        {"X": X, "Y": -Y},
        {"X": -X, "Y": Y},
        {"X": -X, "Y": -Y},
        {"X": Y, "Y": X},
        {"X": Y, "Y": -X},
        {"X": -Y, "Y": X},
        {"X": -Y, "Y": -X},
    ]
    sub_sg = [
        {"x": +x + sp.S(1) / 2},
        {"x": -x + sp.S(1) / 2},
        {"x": -x + sp.S(1) / 4},
        {"x": -x + sp.S(3) / 4},
        {"z": +z + sp.S(1) / 2},
        {"z": -z + sp.S(1) / 2},
        {"z": -z + sp.S(3) / 4},
        {"x": +x + sp.S(1) / 2, "X": Y, "Y": X},
    ]

    unique_bond_set = []
    seen = []
    for bs in bonds:
        sub = sub_pg if plus_set is None else sub_sg + sub_pg
        bss = [create_unique_bonds(replace(bs, d), plus_set) for d in sub]
        unique = create_unique_bonds(bs, plus_set)
        if unique not in seen:
            seen.append(unique)
            for i in bss:
                seen.append(i)
            unique_bond_set.append(bs)
    unique_bond_set = np.asarray(unique_bond_set)
    return unique_bond_set


# ==================================================
def create_unique_bond_set(unique_vector_set, centers, pg_so, mapping, plus_set):
    """
    Create unique bond set.

    Args:
        unique_vector_set (ndarray): unique vector set.
        centers (ndarray): centers.
        pg_so (ndarray): point-group symmetry operations.
        mapping (list): SO mapping.
        plus_set (bool): plus set (None for point group).

    Returns:
        - (dict) -- unique bond set.
    """
    unique_bond_set = {}
    for n, vector_set in unique_vector_set.items():
        all_centers = np.repeat(centers, n, axis=0)
        n_bonds = []
        for v in vector_set:
            all_vectors = np.asarray([v @ pg_so[si[0]].T for si in mapping]).reshape(-1, 3)
            all_bonds = np.hstack((all_vectors[:, None], all_centers[:, None])).reshape(-1, 6)
            n_bonds.append(all_bonds)
        n_bonds = remove_equivalent_bond_set(n_bonds, plus_set)
        unique_bond_set[n] = n_bonds

    return unique_bond_set


# ==================================================
def create_bond_cluster_mapping(pg_so, wyckoff_site_wp, cluster_bonds):
    """
    Create bond cluster mapping.

    Args:
        pg_so (ndarray): point-group symmetry operations.
        wyckoff_site_wp (dict): wyckoff site for s_wp dict.
        cluster_bonds (ndarray): cluster bonds.

    Returns:
        - (str) -- vectors.
        - (str) -- centers.
        - (list) -- mapping.
    """
    vectors = list(map(tuple, pg_so @ cluster_bonds[0, 0:3]))
    site_dict = dict(zip(list(map(tuple, wyckoff_site_wp["conventional"])), wyckoff_site_wp["mapping"]))

    bond_mapping = []
    for bond in cluster_bonds:
        vector = bond[0:3]
        center = tuple(bond[3:6])
        site_mapping = site_dict[center]
        mapping = []
        for no, x in enumerate(vectors):
            if x == tuple(vector) and no + 1 in site_mapping:
                mapping.append(no + 1)
            elif x == tuple(-vector) and no + 1 in site_mapping:
                mapping.append(-(no + 1))
        mapping = sorted(mapping, key=lambda x: abs(x))
        bond_mapping.append(mapping)

    dic = {}
    for bond, mapping in zip(cluster_bonds, bond_mapping):
        if mapping[0] < 0:
            mapping = (-np.asarray(mapping, dtype=int)).tolist()
            bond[0:3] = -bond[0:3]
        dic[tuple(mapping)] = bond.tolist()
    bond_mapping = dict(sorted(dic.items()))

    mapping = list(map(list, bond_mapping.keys()))
    bonds = np.asarray(list(bond_mapping.values()))
    vectors, centers = bonds[:, 0:3], bonds[:, 3:6]

    return str(vectors.tolist()).replace(" ", ""), str(centers.tolist()).replace(" ", ""), mapping


# ==================================================
def create_wyckoff_bond_wp(crystal, pg_so, wyckoff_site, s_wp, plus_set, debug=False):
    """
    Create Wyckoff bond for given s_wp.

    Args:
        crystal (str): crystal.
        pg_so (ndarray): point-group symmetry operations.
        wyckoff_site (dict): Wyckoff site dict.
        s_wp (str): Wyckoff site tag.
        plus_set (ndarray): plus set (None for point group).
        debug (bool, optional): debug print ?

    Returns:
        - (dict) -- Wyckoff bond dict for given s_wp.
    """
    vector_source = vector_source_dict[crystal]
    sub_vector = sub_vector_dict[crystal]

    mapping = np.asarray(wyckoff_site[s_wp]["mapping"], dtype=int) - 1
    centers = wyckoff_site[s_wp]["conventional"]
    ns_wp = int(s_wp[:-1])

    vectors = create_transformed_vector_first(pg_so, mapping)
    unique_vector_source = create_unique_vector_source(vectors, vector_source, sub_vector)
    unique_vector_set = create_unique_vector_set(unique_vector_source, pg_so, mapping)
    unique_bond_set = create_unique_bond_set(unique_vector_set, centers, pg_so, mapping, plus_set)

    wyckoff_bond_wp = {}
    no = 0
    for n, unique_bond in unique_bond_set.items():
        for bond in unique_bond:
            key = f"{n*ns_wp}{chr(97+no)}@{s_wp}"
            wyckoff_bond_wp[key] = create_bond_cluster_mapping(pg_so, wyckoff_site[s_wp], bond)
            no += 1

    if debug:
        print("vectors at 1st WS =", vectors.tolist())
        print("--- unique vector source ----")
        for n, v in unique_vector_source.items():
            print(n, v.tolist())
        print("--- unique vector set ----")
        for n, vs in unique_vector_set.items():
            print(n, vs.tolist())
        print("--- unique bond set ----")
        for n, n_bonds in unique_bond_set.items():
            print("---", n, "---")
            for i in n_bonds:
                print(i.tolist())

    return wyckoff_bond_wp


# ==================================================
# vector direction source.
vector_source_dict = {
    # Ci : remove opposite direction.
    "triclinic": ["[X,Y,Z]"],
    # C2h : remove opposite direction.
    "monoclinic": ["[X,Y,Z]", "[X,0,Z]", "[0,Y,0]"],  # Y=0, X=Z=0.
    # D2h : remove opposite direction, x:-x, y:-y, z:-z.
    "orthorhombic": [
        "[X,Y,Z]",
        "[X,Y,0]",  # Z=0.
        "[X,0,Z]",  # Y=0.
        "[0,Y,Z]",  # X=0.
        "[0,0,Z]",  # X=Y=0.
        "[0,Y,0]",  # X=Z=0.
        "[X,0,0]",  # Y=Z=0.
    ],
    # D4h : remove opposite direction, x:-x, y:-y, z:-z, x:y.
    "tetragonal": [
        "[X,Y,Z]",
        #
        "[X,0,Z]",  # Y=0.
        "[0,X,Z]",
        #
        "[X,X,Z]",  # Y=X.
        "[X,-X,Z]",
        #
        "[X,Y,0]",  # Z=0.
        #
        "[X,0,0]",  # Y=Z=0.
        "[0,X,0]",
        #
        "[X,X,0]",  # Z=0.
        "[X,-X,0]",
        #
        "[0,0,Z]",  # X=Y=0.
    ],
    # D3d, D3d-1 : remove opposite direction, x:-x, y:-y, z:-z.
    "trigonal": [
        "[X,Y,Z]",
        #
        "[X,0,Z]",  # Y=0, D3d.
        "[0,X,Z]",
        "[X,X,-Z]",
        #
        "[X,-X,Z]",  # Y=-X, D3d-1.
        "[X,2X,Z]",
        "[2X,X,-Z]",
        #
        "[X,-X,0]",  # Y=-X, Z=0, D3d.
        "[X,2X,0]",
        "[2X,X,0]",
        #
        "[X,0,0]",  # Y=Z=0, D3d-1.
        "[0,X,0]",
        "[X,X,0]",
        #
        "[0,0,Z]",  # X=Y=0.
    ],
    # D6h : remove opposite direction, x:-x, y:-y, z:-z, x:y.
    "hexagonal": [
        "[X,Y,Z]",
        #
        "[X,Y,0]",  # Z=0.
        "[X-Y,-Y,0]",
        "[X,X-Y,0]",
        #
        "[X,2X,Z]",  # Y=2X.
        "[X,-X,Z]",
        "[2X,X,Z]",
        #
        "[X,0,Z]",  # Y=0.
        "[0,X,Z]",
        "[X,X,Z]",
        #
        "[X,2X,0]",  # Y=2X, Z=0.
        "[X,-X,0]",
        "[2X,X,0]",
        #
        "[X,0,0]",  # Y=Z=0.
        "[0,X,0]",
        "[X,X,0]",
        #
        "[0,0,Z]",
    ],
    # Oh : remove opposite direction, x:-x, y:-y, z:-z, y:z.
    "cubic": [
        "[X,Y,Z]",
        #
        "[X,X,Y]",  # Y=X, Z=Y.
        "[X,Y,X]",
        "[Y,X,X]",
        "[X,-X,Y]",
        "[-X,Y,X]",
        "[Y,X,-X]",
        #
        "[X,Y,0]",  # Z=0.
        "[Y,0,X]",
        "[0,X,Y]",
        #
        "[X,X,X]",  # Y=Z=X.
        "[X,X,-X]",
        "[X,-X,X]",
        "[-X,X,X]",
        #
        "[X,X,0]",  # Y=X, Z=0.
        "[X,0,X]",
        "[0,X,X]",
        "[X,-X,0]",
        "[-X,0,X]",
        "[0,X,-X]",
        #
        "[X,0,0]",  # Y=Z=0.
        "[0,0,X]",
        "[0,X,0]",
    ],
}
vector_source_dict = {crystal: np.asarray([str_to_sympy(i) for i in data]) for crystal, data in vector_source_dict.items()}

# ==================================================
# sub vector.
sub_vector_dict = {
    "triclinic": {"[X,Y,Z]": []},
    "monoclinic": {"[X,Y,Z]": ["[X,0,Z]", "[0,Y,0]"], "[X,0,Z]": [], "[0,Y,0]": []},
    "orthorhombic": {
        "[X,Y,Z]": ["[X,Y,0]", "[X,0,Z]", "[0,Y,Z]", "[0,0,Z]", "[0,Y,0]", "[X,0,0]"],
        "[X,Y,0]": ["[X,0,0]", "[0,Y,0]"],
        "[X,0,Z]": ["[X,0,0]", "[0,0,Z]"],
        "[0,Y,Z]": ["[0,Y,0]", "[0,0,Z]"],
        "[0,0,Z]": [],
        "[0,Y,0]": [],
        "[X,0,0]": [],
    },
    "tetragonal": {
        "[X,Y,Z]": [
            "[X,0,Z]",
            "[0,X,Z]",
            "[X,X,Z]",
            "[X,-X,Z]",
            "[X,Y,0]",
            "[X,0,0]",
            "[0,X,0]",
            "[X,X,0]",
            "[X,-X,0]",
            "[0,0,Z]",
        ],
        #
        "[X,0,Z]": ["[X,0,0]", "[0,0,Z]"],
        "[0,X,Z]": ["[0,X,0]", "[0,0,Z]"],
        #
        "[X,X,Z]": ["[X,X,0]", "[0,0,Z]"],
        "[X,-X,Z]": ["[X,-X,0]", "[0,0,Z]"],
        #
        "[X,Y,0]": ["[X,0,0]", "[0,X,0]", "[X,X,0]", "[X,-X,0]"],
        #
        "[X,0,0]": [],
        "[0,X,0]": [],
        #
        "[X,X,0]": [],
        "[X,-X,0]": [],
        #
        "[0,0,Z]": [],
    },
    "trigonal": {
        "[X,Y,Z]": [
            "[X,0,Z]",
            "[0,X,Z]",
            "[X,X,-Z]",
            "[X,-X,Z]",
            "[X,2X,Z]",
            "[2X,X,-Z]",
            "[X,-X,0]",
            "[X,2X,0]",
            "[2X,X,0]",
            "[X,0,0]",
            "[0,X,0]",
            "[X,X,0]",
            "[0,0,Z]",
        ],
        #
        "[X,0,Z]": ["[X,0,0]", "[0,0,Z]"],
        "[0,X,Z]": ["[0,X,0]", "[0,0,Z]"],
        "[X,X,-Z]": ["[X,X,0]", "[0,0,Z]"],
        #
        "[X,-X,Z]": ["[X,-X,0]", "[0,0,Z]"],
        "[X,2X,Z]": ["[X,2X,0]", "[0,0,Z]"],
        "[2X,X,-Z]": ["[2X,X,0]", "[0,0,Z]"],
        #
        "[X,-X,0]": [],
        "[X,2X,0]": [],
        "[2X,X,0]": [],
        #
        "[X,0,0]": [],
        "[0,X,0]": [],
        "[X,X,0]": [],
        #
        "[0,0,Z]": [],
    },
    "hexagonal": {
        "[X,Y,Z]": [
            "[X,Y,0]",
            "[X-Y,-Y,0]",
            "[X,X-Y,0]",
            "[X,2X,Z]",
            "[X,-X,Z]",
            "[2X,X,Z]",
            "[X,0,Z]",
            "[0,X,Z]",
            "[X,X,Z]",
            "[X,2X,0]",
            "[X,-X,0]",
            "[2X,X,0]",
            "[X,0,0]",
            "[0,X,0]",
            "[X,X,0]",
            "[0,0,Z]",
        ],
        #
        "[X,Y,0]": ["[X-Y,-Y,0]", "[X,X-Y,0]", "[X,2X,0]", "[X,-X,0]", "[2X,X,0]", "[X,0,0]", "[0,X,0]", "[X,X,0]"],
        "[X-Y,-Y,0]": ["[X,Y,0]", "[X,X-Y,0]", "[X,2X,0]", "[X,-X,0]", "[2X,X,0]", "[X,0,0]", "[0,X,0]", "[X,X,0]"],
        "[X,X-Y,0]": ["[X,Y,0]", "[X-Y,-Y,0]", "[X,2X,0]", "[X,-X,0]", "[2X,X,0]", "[X,0,0]", "[0,X,0]", "[X,X,0]"],
        #
        "[X,2X,Z]": ["[X,2X,0]", "[0,0,Z]"],
        "[X,-X,Z]": ["[X,-X,0]", "[0,0,Z]"],
        "[2X,X,Z]": ["[2X,X,0]", "[0,0,Z]"],
        #
        "[X,0,Z]": ["[X,0,0]", "[0,0,Z]"],
        "[0,X,Z]": ["[0,X,0]", "[0,0,Z]"],
        "[X,X,Z]": ["[X,X,0]", "[0,0,Z]"],
        #
        "[X,2X,0]": [],
        "[X,-X,0]": [],
        "[2X,X,0]": [],
        #
        "[X,0,0]": [],
        "[0,X,0]": [],
        "[X,X,0]": [],
        #
        "[0,0,Z]": [],
    },
    "cubic": {
        "[X,Y,Z]": [
            "[X,X,Y]",
            "[X,Y,X]",
            "[Y,X,X]",
            "[X,-X,Y]",
            "[-X,Y,X]",
            "[Y,X,-X]",
            "[X,Y,0]",
            "[Y,0,X]",
            "[0,X,Y]",
            "[X,X,X]",
            "[X,X,-X]",
            "[X,-X,X]",
            "[-X,X,X]",
            "[X,X,0]",
            "[X,0,X]",
            "[0,X,X]",
            "[X,-X,0]",
            "[-X,0,X]",
            "[0,X,-X]",
            "[X,0,0]",
            "[0,0,X]",
            "[0,X,0]",
        ],
        #
        "[X,X,Y]": ["[X,X,X]", "[X,X,-X]", "[X,X,0]", "[0,0,X]"],
        "[X,Y,X]": ["[X,X,X]", "[X,-X,X]", "[X,0,X]", "[0,X,0]"],
        "[Y,X,X]": ["[X,X,X]", "[-X,X,X]", "[0,X,X]", "[X,0,0]"],
        "[X,-X,Y]": ["[X,-X,X]", "[-X,X,X]", "[X,-X,0]", "[0,0,X]"],
        "[-X,Y,X]": ["[X,X,-X]", "[-X,X,X]", "[-X,0,X]", "[0,X,0]"],
        "[Y,X,-X]": ["[X,X,-X]", "[X,-X,X]", "[0,X,-X]", "[X,0,0]"],
        #
        "[X,Y,0]": ["[X,X,0]", "[X,-X,0]", "[X,0,0]", "[0,X,0]"],
        "[Y,0,X]": ["[X,0,X]", "[-X,0,X]", "[X,0,0]", "[0,0,X]"],
        "[0,X,Y]": ["[0,X,X]", "[0,X,-X]", "[0,X,0]", "[0,0,X]"],
        #
        "[X,X,X]": [],
        "[X,X,-X]": [],
        "[X,-X,X]": [],
        "[-X,X,X]": [],
        #
        "[X,X,0]": [],
        "[X,0,X]": [],
        "[0,X,X]": [],
        "[X,-X,0]": [],
        "[-X,0,X]": [],
        "[0,X,-X]": [],
        #
        "[X,0,0]": [],
        "[0,0,X]": [],
        "[0,X,0]": [],
    },
}
