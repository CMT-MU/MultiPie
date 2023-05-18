from multipie.virtual_cluster.virtual_cluster_pg import VirtualClusterPG
from multipie.tag.tag_group import TagGroup


# ==================================================
def test_virtual_cluster_pg():
    print("=== virtual_cluster_pg ===")

    pg = TagGroup.create()

    for i in pg[:2]:
        print("---", i, "---")
        vc = VirtualClusterPG(str(i))
        for tag, basis in vc.items():
            print(tag, list(basis))
        print("site =", vc.site)


# ================================================== main
test_virtual_cluster_pg()
