from multipie.clebsch_gordan.clebsch_gordan_pg_set import ClebschGordanPGSet
from multipie.group.point_group import PointGroup


# ==================================================
def test_clebsch_gordan_pg_set():
    print("=== clebsch_gordan_pg_set ===")

    cgs32 = ClebschGordanPGSet()

    cg = cgs32["D3"]
    pg = PointGroup("D3")

    tag_list = [h.tag for h in pg.harmonics.select()]
    tag1_list = [h.tag for h in pg.harmonics.select(rank=1)]
    tag2_list = [h.tag for h in pg.harmonics.select(rank=1, head="Q")]

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
test_clebsch_gordan_pg_set()
