from multipie.character.character_pg import CharacterPG


# ==================================================
def test_character_pg():
    print("=== character_pg ===")

    c_pg = CharacterPG("D6h")

    print(c_pg)
    print(*c_pg.irrep_list)

    print(*c_pg.symmetry_operation(ret_num=False))
    for irrep in c_pg.keys():
        print(irrep)
        print(c_pg.character(irrep, all_so=True), c_pg[irrep])

    cr = c_pg.compatibility_relation("C3")
    print("D6h => C3\n", *cr.keys(), "\n", *cr.values())
    pc = c_pg.parity_conversion()
    print("parity conv.\n", *pc.keys(), "\n", *pc.values())

    print(c_pg.symmetric_product_decomposition(("E1g", "E2u"), ret_ex=True))
    print(c_pg.anti_symmetric_product_decomposition("E2g", ret_ex=True))


# ================================================== main
test_character_pg()
