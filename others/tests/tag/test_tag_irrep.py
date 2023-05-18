from multipie.tag.tag_irrep import TagIrrep


# ==================================================
def test_tag_irrep():
    print("=== tag_irrep ===")

    tags = TagIrrep.create()

    print("--- repr ---")
    print(tags[:3])

    print("--- latex ---")
    print(tags[:3].latex())

    print("--- symbol ---")
    print(tags[:3].symbol())

    print("--- str ---")
    print(tags[:3].str_list())

    print("--- first and last 3 ---")
    print(*(tags[:3] + tags[-3:]))

    print("--- parity g ---")
    print(*tags.select(parity="g"))

    print("--- dim 2 ---")
    print(*tags.select(dim=2))

    print("--- index ---")
    print(tags.index(tags[2]))


# ================================================== main
test_tag_irrep()
