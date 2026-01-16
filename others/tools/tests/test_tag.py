import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)

from multipie.util.util_tag import TagSymmetryOperation, TagIrrep, TagMultipole, TagBasis

from others.tools.data.data_definition import atomic_basis


# ==================================================
def test_tag_symmetry_operation():
    print("=== tag_symmetry_operation ===")

    tag = [
        "1",
        "2[001]",
        "2[-101]",
        "3+[-1-11]",
        "3-[-11-1]",
        "3-[1-1-1]",
        "-1",
        "m[010]",
        "m[01-1]",
        "-3+[-1-11]",
        "2[010]:[0,0,0]",
        "2[010]:[0,1/2,1/2]",
        "m[010]:[0,1/2,1/2]",
        "m'[010]",
        "6'+[001]",
        "-4'-[001]:[1/4,1/4,3/4]",
        "m'[100]:[0,0,1/2]",
        "m'[120]:[2/3,1/2,1/2]",
    ]

    for t in tag:
        d = TagSymmetryOperation.parse(t)
        t1 = TagSymmetryOperation.str(d)
        print(f"'{t}' = '{t1}'")
        print(TagSymmetryOperation.latex(t))
        assert f"'{t}'" == f"'{t1}'"


# ==================================================
def test_tag_irrep():
    print("=== tag_irrep ===")

    tag = ["Tg", "B1", "A'", "A2", "T2u", "Ega", "Eb"]

    for t in tag:
        d = TagIrrep.parse(t)
        t1 = TagIrrep.str(d)
        print(f"'{t}' = '{t1}'")
        print(TagIrrep.latex(t))


# ==================================================
def test_tag_multipole():
    print("=== tag_multipole ===")

    # (idx, mp_type, component).
    lst = [
        (("Q", 1), "spherical", None),
        (("Q", 1, 1, 1), "spherical", None),
        (("Q", 2, "E", 1, -1), "point_group", 1),
        (("Q", 2, "E", -1, 1), "point_group", 1),
        (("Q", 2, "E", -1, -1), "point_group", 1),
        (("Q", 2, "E", 2, 2), "point_group", 1),
    ]
    tag = {"spherical": "a", "point_group": "c"}

    for idx, mp_type, component in lst:
        d = TagMultipole.parse(idx, mp_type, component)
        t1 = TagMultipole.str(d)
        print(f"'({idx}, {mp_type}, {component})' => '{t1}'")
        print(TagMultipole.latex(idx, mp_type, component, tag=tag[mp_type]))


# ==================================================
def test_tag_basis():
    print("=== tag_basis ===")
    basis = atomic_basis["spinless"] | atomic_basis["spinful"]
    for grp, bs in basis.items():
        print("===", grp, "===")
        for L, basis in bs.items():
            for t in basis:
                d = TagBasis.parse(t)
                t1 = TagBasis.str(d, L)
                print(f"'{t}' = '{t1}'")
                print(TagBasis.latex(t, L))


# ================================================== main
# test_tag_symmetry_operation()
# test_tag_irrep()
test_tag_multipole()
# test_tag_basis()
