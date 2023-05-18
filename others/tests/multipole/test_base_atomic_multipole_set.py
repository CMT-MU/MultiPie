from multipie.multipole.base.base_atomic_multipole_set import BaseAtomicMultipoleSet


# ==================================================
def test_base_atomic_multipole_set():
    print("=== base_atomic_multipole_set ===")

    bam = BaseAtomicMultipoleSet("lm")
    print(*bam.key_list()[:4])
    print(*bam.key_list()[-4:])


# ================================================== main
test_base_atomic_multipole_set()
