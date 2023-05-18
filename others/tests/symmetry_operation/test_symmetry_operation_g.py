from gcoreutils.nsarray import NSArray
from multipie.tag.tag_group import TagGroup
from multipie.symmetry_operation.symmetry_operation_g import SymmetryOperationG


# ==================================================
def test_symmetry_operation_g():
    print("=== symmetry_operation_g ===")

    pg = TagGroup.create()
    sg = TagGroup.create(space_group=True)

    for i in pg[:2] + sg[-2:]:
        print("---", i, "---")
        sog = SymmetryOperationG(i)
        print("polar", *sog.mat())
        print("axial", *sog.mat(axial=True))

    print("--- product ---")
    sog = SymmetryOperationG("C3v")
    for i in sog.keys():
        for j in sog.keys():
            print(i, "*", j, "=", sog.product(i, j))

    print("--- equivalent vector (point group) ---")
    sog = SymmetryOperationG("C3v")
    v = NSArray.vector3d(head="Q")
    gv = NSArray.vector3d(head="G")
    print(*sog._equivalent_vector(v=v, cc_only=True, axial=False))
    print(*sog._equivalent_vector(v=gv, cc_only=True, axial=True))

    print("--- equivalent bond (point group) ---")
    bond = NSArray("[0,1/2,0];[1/2,0,0]")
    print(*sog._equivalent_bond(bond, cc_only=True))

    print("--- equivalent vector (space group) ---")
    sog = SymmetryOperationG("C3v^5")
    print(*sog._equivalent_vector(v=v, cc_only=True, axial=False, plus_set=True))
    print(*sog._equivalent_vector(v=gv, cc_only=True, axial=True, plus_set=True))

    print("--- equivalent bond (space group) ---")
    print(*sog._equivalent_bond(bond, cc_only=True))


# ================================================== main
test_symmetry_operation_g()
