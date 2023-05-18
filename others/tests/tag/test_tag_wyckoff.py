from multipie.tag.tag_wyckoff import TagWyckoff


# ==================================================
def test_tag_wyckoff():
    print("=== tag_wyckoff ===")

    tags = TagWyckoff.create("Oh")

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

    print("--- n=24 ---")
    print(*tags.select(n=24))


# ================================================== main
test_tag_wyckoff()
