from multipie.group.point_group import PointGroup


# ==================================================
def test_point_group():
    print("=== point_group ===")

    pg = PointGroup("C3v", verbose=False)

    print("--- harmonics ---")
    for h in pg.harmonics.select(rank=1, head="G"):
        print(h, h.definition(), h.expression(), h.u_matrix())

    print("--- character ---")
    print(*pg.character.irrep_list)
    print(pg.character.symmetric_product_decomposition(("E", "E"), ret_ex=True))

    print("--- symmetry_operation ---")
    print("symmetry operations:", *pg.symmetry_operation.keys())
    print("polar:", pg.symmetry_operation.mat(axial=False))
    print("axial:", pg.symmetry_operation.mat(axial=True))

    print("--- wyckoff ---")
    for wp in pg.wyckoff.keys():
        print(wp, pg.wyckoff.position(wp))

    print("--- virtual_cluster ---")
    for tag, basis in pg.virtual_cluster.items():
        print(tag, basis)
    print("site =", pg.virtual_cluster.site)

    print("--- response_tensor ---")
    pg.response._dump()

    print("--- transform_matrix (site) ---")
    site = "[1/sqrt(3),2/sqrt(3),0]"
    m, b = pg.transform_matrix_site(site)
    print("matrix:", m)
    print("basis:", b)

    print("--- transform_matrix (bond) ---")
    bond = "[0,0,0];[1/sqrt(3),2/sqrt(3),0]"
    m, b = pg.transform_matrix_bond(bond)
    print("matrix:", m)
    print("basis:", b)

    print("--- transform_matrix (vector) ---")
    m, b = pg.transform_matrix_vector(axial=True)
    print("matrix:", m)
    print("basis:", b)

    print("--- transform_matrix (orbital) ---")
    m, b = pg.transform_matrix_orbital(rank=1)
    print("matrix:", m)
    print("basis:", b)

    print("--- transform (site) ---")
    print(pg.transform_site(site))

    print("--- transform (bond) ---")
    print(pg.transform_bond(bond))

    print("--- transform (vector) ---")
    print(pg.transform_vector(site, axial=True))

    print("--- transform (orbital) ---")
    print(pg.transform_orbital("x"))

    print("--- site mapping ---")
    d = pg.site_mapping(site)
    for k, v in d.items():
        print(k, v)

    print("--- bond mapping ---")
    d, nd = pg.bond_mapping(bond)
    for k, v in d.items():
        print(k, v)
    print(nd)

    print("--- find wyckoff position ---")
    print(pg.find_wyckoff_position(site))

    print("--- virtual cluster basis ---")
    d, s = pg.virtual_cluster_basis(ortho=True)
    print("sites:", s)
    for i, b in d.items():
        print(i, b)


# ================================================== main
test_point_group()
