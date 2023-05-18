from gcoreutils.nsarray import NSArray
from multipie.harmonics.harmonics import Harmonics
from multipie.data.data_harmonics import _data_harmonics_polar


# ==================================================
def test_harmonics():
    print("=== harmonics ===")

    pg_tag = "C3"

    hs = {tag: Harmonics(tag, *h) for tag, h in _data_harmonics_polar[pg_tag].items()}

    h = hs["Qh(2,Ea,1,)"]
    v1 = NSArray.vector3d("Q")
    v2 = NSArray.vector3d("G")

    print(h.definition())
    print(h.expression(v=v1))
    print(h.u_matrix())
    print(h.equivalent_operator(j="2"))

    h = hs["Qh(2,Ea,1,)"]

    print(h.definition())
    print(h.expression(v=v2))
    print(h.u_matrix())
    print(h.equivalent_operator(j="2"))

    print("--- equivalent operator ---")
    o40 = "3 * x**4 + 3 * y**4 - 24 * y**2 * z**2 + 8 * z**4 + 6 * x**2 * (y**2 - 4 * z**2)"
    o44 = "x**4 - 6 * x**2 * y**2 + y**4"
    j = "2"
    m40 = Harmonics.create_equivalent_operator(o40, j)
    m44 = Harmonics.create_equivalent_operator(o44, j)
    print(m40 + 5 * m44)


# ================================================== main
test_harmonics()
