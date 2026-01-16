"""
Create binary data concerning all of crystallographic point and space groups.
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import os
import logging

from create_binary.create_data_info import create_info
from create_binary.create_data_so_matrix import create_so_matrix
from create_binary.create_data_harmonics_spherical import create_harmonics_spherical
from create_binary.create_data_group import create_group
from create_binary.create_data_harmonics import create_harmonics
from create_binary.create_data_harmonics_multipole import create_harmonics_multipole
from create_binary.create_data_rep_matrix import create_rep_matrix
from create_binary.create_data_root_cluster import create_root_cluster
from create_binary.create_data_atomic_multipole import create_atomic_multipole
from create_binary.create_data_atomic_multipole_group import create_atomic_multipole_group
from create_binary.create_data_harmonics_root_cluster import create_harmonics_root_cluster
from create_binary.create_data_cluster_samb import create_cluster_samb


# ==================================================
def create_data():
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    os.makedirs(BIN_DIR, exist_ok=True)

    create_info()  # info: 120 KB, 0.24 s.
    create_so_matrix()  # so_matrix: 93 KB, 730 s.
    create_harmonics_spherical()  # harmonics_spherical: 1.3 MB, 130 s.
    create_group()  # group: 1.1 MB, 1200 s.
    create_harmonics()  # harmonics: 970 KB, 93 s.
    create_harmonics_multipole()  # harmonics_multipole: 59 MB, 20000 s.
    create_rep_matrix()  # rep_matrix: 3.5 KB, 0.91 s.
    create_root_cluster()  # root_cluster: 0.98 KB, 0.19 s.
    create_atomic_multipole()  # atomic_multipole: 79 KB, 55 s.
    create_atomic_multipole_group()  # atomic_multipole_group: 6.3 MBytes, 6100 s.
    create_harmonics_root_cluster()  # harmonics_root_cluster: 8.6 MBytes, 1800 s.
    create_cluster_samb()  # cluster_samb: 3.7 MBytes, 17000 s.


# ================================================== main
if __name__ == "__main__":
    create_data()
