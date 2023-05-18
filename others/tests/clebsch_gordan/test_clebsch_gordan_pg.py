from multipie.clebsch_gordan.clebsch_gordan_pg import ClebschGordanPG
from multipie.group.point_group import PointGroup


# ==================================================
def test_clebsch_gordan_pg():
    print("=== clebsch_gordan_pg ===")

    cg = ClebschGordanPG("D6h")
    pg = PointGroup("D6h")

    tag_list = [h.tag for h in pg.harmonics.select(head="Q")]
    tag1_list = [h.tag for h in pg.harmonics.select(rank=0, head="Q")]
    tag2_list = [h.tag for h in pg.harmonics.select(rank=3, head="Q")]

    print("--- CG coefficients ---")
    d = {}
    for tag in tag_list:
        lst = []
        for tag1 in tag1_list:
            for tag2 in tag2_list:
                c = cg.cg(tag1, tag2, tag)
                if c != 0:
                    lst.append((c, tag1, tag2))
        if lst != []:
            d[tag] = lst

    for i, (tag, lst) in enumerate(d.items()):
        print(i, ":", tag)
        for c, tag1, tag2 in lst:
            print(f" <{str(tag1)};{str(tag2)}|{str(tag)}> = {c}")

    print(len(d), len(tag1_list), len(tag2_list))


# ================================================== main
test_clebsch_gordan_pg()
