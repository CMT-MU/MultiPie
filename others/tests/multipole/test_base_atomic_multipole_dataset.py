from multipie.multipole.base.base_atomic_multipole_dataset import BaseAtomicMultipoleDataset


# ==================================================
def test_base_atomic_multipole_dataset():
    print("=== base_atomic_multipole_dataset ===")

    bam = BaseAtomicMultipoleDataset()["jm"]
    print(*bam.key_list()[:4])
    print(*bam.key_list()[-4:])


# ================================================== main
test_base_atomic_multipole_dataset()
