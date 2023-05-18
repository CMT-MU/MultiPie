from multipie.character.character_pg_set import CharacterPGSet


# ==================================================
def test_character_pg_set():
    print("=== character_pg_set ===")

    char_set = CharacterPGSet()
    c_pg = char_set["D6h"]

    print(c_pg)
    print(*c_pg.irrep_list)

    print(*c_pg.symmetry_operation(ret_num=False))
    for irrep, ct in c_pg.items():
        print(irrep)
        print(c_pg.character(irrep), c_pg[irrep])

    cr = c_pg.compatibility_relation("C3")
    print("D6h => C3\n", *cr.keys(), "\n", *cr.values())
    pc = c_pg.parity_conversion()
    print("parity conv.\n", *pc.keys(), "\n", *pc.values())

    print(c_pg.symmetric_product_decomposition(("E1g", "E2u"), ret_ex=True))
    print(c_pg.anti_symmetric_product_decomposition("E2g", ret_ex=True))


# ================================================== main
test_character_pg_set()
