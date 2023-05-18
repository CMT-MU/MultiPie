from multipie.symmetry_operation.symmetry_operation_g_set import SymmetryOperationGSet


# ==================================================
def test_symmetry_operation_g_set():
    print("=== symmetry_operation_g_set ===")

    soset = SymmetryOperationGSet()

    lst = ["C3", "D3^2"]

    for i in lst:
        print("---", i, "---")
        sog = soset[i]
        print("polar", sog.mat())
        print("axial", sog.mat(axial=True))


# ================================================== main
test_symmetry_operation_g_set()
