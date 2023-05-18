from multipie.tag.tag_group import TagGroup
from multipie.const import __def_dict__


# ==================================================
def test_tag_group():
    print("=== tag_group ===")

    tags = TagGroup.create(space_group=False)

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

    print("--- crystal hexagonal ---")
    print(tags.select(crystal="hexagonal"))

    tags = TagGroup.create(space_group=True)

    print("--- crystal cubic ---")
    print(tags.select(crystal="cubic"))

    print("--- subgroup ---")
    print(*TagGroup.create(crystal="cubic_series"))
    print(*TagGroup.create(crystal="hexagonal_series"))

    print("--- crystal ---")
    for crystal in __def_dict__["crystal"]:
        print(crystal, ":", *TagGroup.create(space_group=True, crystal=crystal))


# ================================================== main
test_tag_group()
