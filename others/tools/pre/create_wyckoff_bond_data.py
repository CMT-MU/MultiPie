"""
Create Wyckoff bond data for all groups.
Created "output/data_wyckoff_bond.py" should be moved to "multipie/data".
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import os
import subprocess

from multipie.util.util_binary import BinaryManager

from util_cluster import create_wyckoff_bond_wp

header_wb = """
Wyckoff bond for all space and point groups (fractional coordinate, conventional cell, no plus set).
- "wyckoff_bond_data" (dict): bond data.
    {group: { Wyckoff: (vector, center, mapping) } }
"""


# ==================================================
def create_wyckoff_bond_data():
    """
    Create wyckoff bond.
    """
    info = BinaryManager("info", topdir=BIN_DIR)
    group = BinaryManager("group", topdir=BIN_DIR)

    cur_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(cur_dir, "output/data_wyckoff_bond.py")

    wyckoff_bond = {}
    for no in info["id_set"]["PG"]["all"]:
        tag = info["tag"][no]
        gp = group[no]
        print("creating wyckoff bond for", tag)

        crystal = gp["info"].crystal
        pg_so = gp["symmetry_operation"]["fractional"][:, 0:3, 0:3]
        wyckoff_site = gp["wyckoff"]["site"]
        plus_set = None

        dic = {}
        for s_wp in wyckoff_site.keys():
            dic.update(create_wyckoff_bond_wp(crystal, pg_so, wyckoff_site, s_wp, plus_set))
        wyckoff_bond[tag] = dic

    for no in info["id_set"]["SG"]["all"]:
        tag = info["tag"][no]
        gp = group[no]
        print("creating wyckoff bond for", tag)

        crystal = gp["info"].crystal
        pg_so = gp["symmetry_operation"]["fractional"][:, 0:3, 0:3]
        wyckoff_site = gp["wyckoff"]["site"]
        plus_set = gp["symmetry_operation"]["plus_set"]

        dic = {}
        for s_wp in wyckoff_site.keys():
            dic.update(create_wyckoff_bond_wp(crystal, pg_so, wyckoff_site, s_wp, plus_set))
        wyckoff_bond[tag] = dic

    # output file.
    data = '"""' + header_wb + '"""\n\n' + "wyckoff_bond_data = " + str(wyckoff_bond)
    with open(out_file, "w", encoding="utf-8") as f:
        print(data, file=f, end="\n")


# ==================================================
if __name__ == "__main__":
    create_wyckoff_bond_data()

    cmd = "black --line-length=130 data_wyckoff_bond.py"
    cur_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    try:
        subprocess.run(cmd, capture_output=True, check=True, cwd=cur_dir, shell=True)
    except subprocess.CalledProcessError:
        raise Exception("Formatting by black is failed.")
