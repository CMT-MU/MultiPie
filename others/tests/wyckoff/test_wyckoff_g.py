from multipie.wyckoff.wyckoff_g import WyckoffG
from multipie.tag.tag_group import TagGroup


# ==================================================
def test_wyckoff_g():
    print("=== wyckoff_g ===")

    pg = TagGroup.create()
    sg = TagGroup.create(space_group=True)

    for i in pg[:2]:
        print("---", i, "---")
        sog = WyckoffG(str(i))
        for wp in sog.keys():
            print(wp, *sog.position(wp), *sog.position(wp, default=True))

    for i in sg[:2]:
        print("---", i, "---")
        sog = WyckoffG(str(i))
        for wp in sog.keys():
            print(wp, *sog.position(wp))


# ================================================== main
test_wyckoff_g()
