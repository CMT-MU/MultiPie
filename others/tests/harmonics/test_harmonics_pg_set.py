from multipie.harmonics.harmonics_pg_set import HarmonicsPGSet


# ==================================================
def test_harmonics_pg_set():
    print("=== harmonics_pg_set ===")

    hs32 = HarmonicsPGSet()

    hs = hs32["Oh"]
    print(hs)

    for h in hs.select(rank=1, head="G"):
        print(h, h.definition(), h.expression(), h.u_matrix())


# ================================================== main
test_harmonics_pg_set()
