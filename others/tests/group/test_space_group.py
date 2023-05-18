from multipie.group.space_group import SpaceGroup


# ==================================================
def test_space_group():
    print("=== space_group ===")

    sg = SpaceGroup("C3v^5")

    print("--- symmetry_operation ---")
    print("symmetry operations:", *sg.symmetry_operation.keys())
    print("polar:", sg.symmetry_operation.mat(axial=False))
    print("axial:", sg.symmetry_operation.mat(axial=True))

    print("--- wyckoff ---")
    print("9b", sg.wyckoff.position("9b"))

    print("--- transform_matrix (site) ---")
    site = "[1/2,1/3,0]"
    m, b = sg.transform_matrix_site(site)
    print("matrix:", m)
    print("basis:", b)

    print("--- transform_matrix (bond) ---")
    bond = "[1/2,1/3,0]@[1/2,1/2,0]"
    m, b = sg.transform_matrix_bond(bond)
    print("matrix:", m)
    print("basis:", b)

    print("--- transform (site) ---")
    print(sg.transform_site(site))

    print("--- transform (bond) ---")
    print(sg.transform_bond(bond))

    print("--- site mapping ---")
    d = sg.site_mapping(site)
    for k, v in d.items():
        print(k, v)

    print("--- bond mapping ---")
    d, nd = sg.bond_mapping(bond)
    for k, v in d.items():
        print(k, v)
    print(nd)

    print("--- find wyckoff position ---")
    print(sg.find_wyckoff_position(site))


# ================================================== main
test_space_group()
