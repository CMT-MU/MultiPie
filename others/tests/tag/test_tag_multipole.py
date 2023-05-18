from multipie.tag.tag_multipole import TagMultipole


# ==================================================
def test_tag_multipole():
    print("=== tag_multipole ===")

    tags = TagMultipole.create(axial=True) + TagMultipole.create(axial=False)

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

    print("--- rank 1 ---")
    print(*tags.select(rank=1))

    print("--- rank 1 and G")
    print(*(tags.select(rank=1, head="G")))

    print("--- A1g ---")
    print(*tags.select(irrep="A1g"))

    print("--- atomic multipole ---")
    t = TagMultipole("Ga(3,,,-2|1,-1)")
    print(t, ": ", t.latex())

    print("--- create spherical ---")
    t = TagMultipole.create_spherical(rank=1, head="Q")
    print(t)
    print(t.replace(head="T"))

    print("--- convert to harmonics ---")
    lst = [
        "Qa(1,E,,1)",
        "Qs(1,E,,0)",
        "Qa(1,E,,0)",
        "Qs(1,E,,1)",
        "Ma(1,E,,0)",
        "Qs(1,E,,0)",
        "Ma(1,E,,1)",
        "Qs(1,E,,1)",
        "Ma(1,A2,,)",
        "Ts(1,A2,,)",
    ]
    for i in lst:
        ctag = TagMultipole(i).to_harmonics()
        print(i, "=>", ctag)


# ================================================== main
test_tag_multipole()
