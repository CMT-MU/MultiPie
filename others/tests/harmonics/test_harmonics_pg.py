from multipie.harmonics.harmonics_pg import HarmonicsPG


# ==================================================
def test_harmonics_pg():
    print("=== harmonics_pg ===")

    hs = HarmonicsPG("C3v")
    print(hs)

    print(*hs.group(rank=1, axial=True))

    for h in hs.select(rank=1, head="G"):
        print(h, h.definition(), h.expression(), h.u_matrix())


# ================================================== main
test_harmonics_pg()
