from multipie.wyckoff.wyckoff_g_set import WyckoffGSet


# ==================================================
def test_wyckoff_g_set():
    print("=== wyckoff_g_set ===")

    wset = WyckoffGSet()

    print("--- C3 ---")
    sog = wset["C3"]
    for wp in sog.keys():
        print(wp, *sog.position(wp))

    print("--- D3^2 ---")
    sog = wset["D3^2"]
    for wp in sog.keys():
        print(wp, *sog.position(wp))


# ================================================== main
test_wyckoff_g_set()
