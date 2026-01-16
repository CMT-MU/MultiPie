from multipie.core.group import Group


# ==================================================
def test_tag_group():
    print("=== tag_group ===")

    tag = ["C3v", "D3-1", "D6h^3", "Oh^4"]

    for t in tag:
        g = Group(t)
        t1 = str(g)
        print(f"'{t}' = '{t1}'")
        print(g.latex())
        print(g.latex(detail=True))


# ================================================== main
test_tag_group()
