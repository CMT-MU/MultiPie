from multipie.harmonics.util.equivalent_operator import equivalent_operator_from_poly


# ==================================================
def test_equivalent_operator():
    o40 = "3 * x**4 + 3 * y**4 - 24 * y**2 * z**2 + 8 * z**4 + 6 * x**2 * (y**2 - 4 * z**2)"
    o44 = "x**4 - 6 * x**2 * y**2 + y**4"
    j = "2"
    m40 = equivalent_operator_from_poly(o40, j)
    m44 = equivalent_operator_from_poly(o44, j)
    print(m40 + 5 * m44)


# ================================================== main
test_equivalent_operator()
