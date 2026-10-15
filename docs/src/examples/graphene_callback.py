"""
Example of callbacks for graphene_post_ctrl.py.

The functions receive the eigenvalues Ek[k,n], the eigenvectors Uk[k,:,n], and the info dict,
and write their results in the output directory of the model.
"""

import os
import numpy as np


# ==================================================
def output_dir(info):
    # "output" is the default of output/dir; change it if output/dir is changed.
    return os.path.join(info["_topdir"], info["name"], "output")


# ==================================================
def callback_k_path(Ek, Uk, info):
    """
    Band energies measured from the Fermi energy along the k path, with the position along the path.
    """
    data = np.column_stack([info["k_linear"], Ek - info["fermi_energy"]])
    header = "k " + " ".join(f"E{n + 1}" for n in range(Ek.shape[1]))
    np.savetxt(os.path.join(output_dir(info), "my_band.txt"), data, header=header)


# ==================================================
def callback_k_full(Ek, Uk, info):
    """
    Occupation of each ket at T = 0, n_m = (1/Nk) sum_{k,n} theta(E_F - E_kn) |U_k[m,n]|^2.
    """
    occupied = Ek < info["fermi_energy"]
    n_m = np.einsum("kn,kmn->m", occupied, np.abs(Uk) ** 2) / len(Ek)
    with open(os.path.join(output_dir(info), "my_occupation.txt"), "w", encoding="utf-8") as f:
        for ket, n in zip(info["ket"], n_m):
            print(f"{ket} {n:.6f}", file=f)
