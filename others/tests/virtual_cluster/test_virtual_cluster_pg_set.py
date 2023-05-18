from multipie.virtual_cluster.virtual_cluster_pg_set import VirtualClusterPGSet


# ==================================================
def test_virtual_cluster_pg_set():
    print("=== virtual_cluster_pg_set ===")

    vcset = VirtualClusterPGSet()

    print("--- C3 ---")
    vc = vcset["C3"]
    for tag, basis in vc.items():
        print(tag, list(basis))
    print("site =", vc.site)


# ================================================== main
test_virtual_cluster_pg_set()
