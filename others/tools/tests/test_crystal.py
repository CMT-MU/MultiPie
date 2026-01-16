import numpy as np
import sympy as sp

from multipie.util.util import str_to_sympy
from multipie.util.util_crystal import (
    convert_to_cartesian_hexagonal,
    convert_to_cartesian_hexagonal_matrix,
    convert_to_fractional_hexagonal,
    convert_to_fractional_hexagonal_matrix,
    convert_to_primitive,
    convert_to_primitive_matrix,
    convert_to_conventional,
    convert_to_conventional_matrix,
    shift_site,
)


# ==================================================
def test_crystal_convert():
    # vector(s) : to cartesian.
    v = str_to_sympy("[1,2,1]")
    vc = convert_to_cartesian_hexagonal(v)
    print("v =", v.tolist())
    print("v_c =", vc.tolist())
    print("v =", convert_to_fractional_hexagonal(vc).tolist())
    print()
    v2 = str_to_sympy("[[1,2,1],[2/3,1/3,1/2]]")
    v2c = convert_to_cartesian_hexagonal(v2)
    print("v2 =", v2.tolist())
    print("v2_c =", v2c.tolist())
    print("v2 =", convert_to_fractional_hexagonal(v2c).tolist())
    print()

    # matrix(s) : to cartesian.
    m = str_to_sympy("[[0,-1,0],[1,-1,0],[0,0,1]]")
    mc = convert_to_cartesian_hexagonal_matrix(m)
    print("m =", m.tolist())
    print("m_c =", mc.tolist())
    print("m =", convert_to_fractional_hexagonal_matrix(mc).tolist())
    print()
    m2 = np.asarray([str_to_sympy(i) for i in ["[[0,-1,0],[1,-1,0],[0,0,1]]", "[[-1,0,0],[0,-1,0],[0,0,-1]]"]])
    m2c = convert_to_cartesian_hexagonal_matrix(m2)
    print("m2 =", m2.tolist())
    print("m2_c =", m2c.tolist())
    print("m2 =", convert_to_fractional_hexagonal_matrix(m2c).tolist())
    print()

    # vector(s) : to primitive (I).
    v = str_to_sympy("[1/10,-3/10,7/5]")
    vc = convert_to_primitive("I", v)
    print("v =", v.tolist())
    print("v_p(I) =", vc.tolist())
    print("v =", convert_to_conventional("I", vc).tolist(), "- [1/2,1/2,-1/2]")
    print()
    v2 = str_to_sympy("[[1,0,0],[0,1,0],[0,0,1]]")
    v2c = convert_to_primitive("I", v2, shift=False)
    print("v2 =", v2.tolist())
    print("v2_p(I) =", v2c.tolist())
    print("v2 =", convert_to_conventional("I", v2c, plus_set=True, shift=False).tolist())
    print()

    # matrix(s) : to primitive (I).
    m = str_to_sympy("[[0,-1,0],[1,-1,0],[0,0,1]]")
    mc = convert_to_primitive_matrix("I", m)
    print("m =", m.tolist())
    print("m_p(I) =", mc.tolist())
    print("m =", convert_to_conventional_matrix("I", mc).tolist())
    print()
    m2 = np.asarray([str_to_sympy(i) for i in ["[[0,-1,0],[1,-1,0],[0,0,1]]", "[[-1,0,0],[0,-1,0],[0,0,-1]]"]])
    m2c = convert_to_primitive_matrix("I", m2)
    print("m2 =", m2.tolist())
    print("m2_p =", m2c.tolist())
    print("m2 =", convert_to_conventional_matrix("I", m2c).tolist())


# ================================================== test.
def test_shift():
    # in the 3rd example, all elements are treated as int!
    test_site = [
        np.array([0.0, 0.5, 0.999999999, 1.000000001, 1.0]),
        str_to_sympy("[0,1/2,1]"),
        np.array([0, sp.S(1) / 2, 1]),
        str_to_sympy("[x+3/2,1/2-y,1-z]"),
        str_to_sympy("[[x+3/2,1/2-y,1-z],[1/3,3/2,-1/4]]"),
    ]

    print("=== test shift site ===")
    for site in test_site:
        sm = np.mod(site, 1)
        s = shift_site(site)
        print(site.tolist(), "=>", s.tolist(), ": simple mod =", sm.tolist())


# ==================================================
test_crystal_convert()
test_shift()
