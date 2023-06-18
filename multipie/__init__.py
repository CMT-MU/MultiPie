import os
import sys

__version__ = "1.1.10"
__top_dir__ = os.path.normpath(os.path.dirname(__file__) + "/../") + "/"
__bin_dir__ = __top_dir__ + "multipie/binary_data/"


# ==================================================
from gcoreutils.string_util import class_name
from gcoreutils.binary_manager import BinaryManager

from multipie.character.character_pg_set import CharacterPGSet
from multipie.harmonics.harmonics_pg_set import HarmonicsPGSet
from multipie.wyckoff.wyckoff_g_set import WyckoffGSet
from multipie.symmetry_operation.symmetry_operation_g_set import SymmetryOperationGSet
from multipie.virtual_cluster.virtual_cluster_pg_set import VirtualClusterPGSet
from multipie.clebsch_gordan.clebsch_gordan_pg_set import ClebschGordanPGSet
from multipie.response_tensor.response_tensor_pg_set import ResponseTensorPGSet
from multipie.multipole.base.base_atomic_multipole_dataset import BaseAtomicMultipoleDataset


# ==================================================
def create_binary(clean=False, verbose=True):
    """
    create binary data.

    Args:
        clean (bool, optional): remove all binaries at beginning.
        verbose (bool, optional): verbose progress ?
    """
    dataset = [
        CharacterPGSet,
        HarmonicsPGSet,
        WyckoffGSet,
        SymmetryOperationGSet,
        VirtualClusterPGSet,
        ClebschGordanPGSet,
        ResponseTensorPGSet,
        BaseAtomicMultipoleDataset,
    ]

    def dprint(s, verbose):
        if verbose:
            print(s, file=sys.stderr)

    if clean:
        for ds in dataset:
            name = class_name(ds) + ".pkl"
            f = __bin_dir__ + name
            if os.path.exists(f):
                os.remove(f)
                dprint("removed " + f, verbose)

    dprint("creating binaries ...", verbose)
    core = BinaryManager(__bin_dir__, verbose=True)
    for ds in dataset:
        core[ds]


# ==================================================
def get_binary(load=False, verbose=False):
    """
    get binary data.

    Args:
        load (bool, optional): load all binaries ?
        verbose (bool, optional): verbose progress ?

    Returns:
        BinaryManager: binary manager.
    """
    core = BinaryManager(__bin_dir__, verbose=verbose)
    if load:
        core.load()
    return core
