from gcoreutils.nsarray import NSArray
from multipie.symmetry_operation.util.symmetry_operation_util import to_reduced, to_cartesian, to_primitive, to_conventional


# ==================================================
def test_symmetry_operation_util():
    print("=== symmetry_operation_util ===")

    print("--- point group hexagonal ---")
    v = NSArray("[1,2,1]")
    m = NSArray("[[0,1,0],[1,0,0],[0,0,1]]")
    crystal = "hexagonal"

    print(v)
    print(to_cartesian(crystal, v))
    print(to_reduced(crystal, to_cartesian(crystal, v)))
    print(m)
    print(to_cartesian(crystal, m))
    print(to_reduced(crystal, to_cartesian(crystal, m)))

    print("--- space group cubic (I) ---")
    v = NSArray("[1/2,1/2,1/2]")
    m = NSArray("[[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]]")
    lattice = "I"

    print(v)
    print(to_primitive(lattice, v))
    print(to_conventional(lattice, to_primitive(lattice, v)))
    print(*to_conventional(lattice, to_primitive(lattice, v), plus_set=True))
    print(m)
    print(to_primitive(lattice, m))
    print(to_conventional(lattice, to_primitive(lattice, m)))
    print(*to_conventional(lattice, to_primitive(lattice, m), plus_set=True))


# ================================================== main
test_symmetry_operation_util()
