"""
Create binary data (information of bond for each groups).
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import logging
import numpy as np

from multipie.util.util import str_to_sympy
from multipie.util.util_binary import BinaryManager
from multipie.util.util_crystal import convert_to_primitive

from others.tools.data.data_wyckoff_bond import wyckoff_bond_data


# ==================================================
def bond_mapping_g(g_info, group_tag, s1c, s1h):
    """
    Bond mapping for give group.

    Args:
        g_info (dict): group data.
        group_tag (str): group tag.
        s1c (ndarray): rep_site for cubic (fractional).
        s1h (ndarray): rep_site for hexagonal (fractional).

    Returns:
        - (dict) -- wyckoff bond dict, {wp: {"conventional/fractional/mapping": data}}.
        - (dict) -- wyckoff site bond dict, {s_wp: [b_wp]}.
    """
    gp = g_info
    wyckoff = gp["wyckoff"]["site"]
    hexagonal = gp["info"].hexagonal_g
    s = s1h if hexagonal else s1c
    point_group = group_tag.count("^") == 0
    lattice = gp["info"].lattice
    if not point_group:
        plus_set = gp["symmetry_operation"]["plus_set"]

    wyckoff_bond_dict = {}
    wyckoff_site_bond = {}
    for wp, (vector, center, mapping) in wyckoff_bond_data[group_tag].items():
        s_wp = wp.split("@")[1]
        wyckoff_site_bond[s_wp] = wyckoff_site_bond.get(s_wp, []) + [wp]
        wyckoff_wp = wyckoff[s_wp]
        site = [str(i) for i in wyckoff_wp["conventional"].tolist()]
        frac = wyckoff_wp["reference"]
        site_mp = dict(zip(site, frac))

        center = str_to_sympy(center)
        bond_v = str_to_sympy(vector.replace("X", f"({s[0]})").replace("Y", f"({s[1]})").replace("Z", f"({s[2]})"))
        bond_c = np.asarray([site_mp[str(i)] for i in center.tolist()])
        vector = str_to_sympy(vector)

        if not point_group:
            bond_v = np.tile(bond_v, (len(plus_set), 1))
            bond_c = np.mod(np.concatenate([bond_c + t for t in plus_set]), 1)  # shift by mod works fine due to symbolic.

        position = np.hstack((vector[:, None], center[:, None])).reshape(-1, 6)
        fractional = np.hstack((bond_v[:, None], bond_c[:, None])).reshape(-1, 6)

        if point_group:
            prim = position
            wyckoff_bond_dict[wp] = {
                "conventional": position,
                "reference": fractional,
                "mapping": mapping,
            }
        else:
            center = convert_to_primitive(lattice, center, shift=True)
            prim = np.hstack((vector[:, None], center[:, None])).reshape(-1, 6)
            wyckoff_bond_dict[wp] = {
                "conventional": position,
                "primitive": prim,
                "reference": fractional,
                "mapping": mapping,
            }

    return wyckoff_bond_dict, wyckoff_site_bond


# ==================================================
def create_group_bond():
    group = BinaryManager("group", topdir=BIN_DIR)
    info = BinaryManager("info", topdir=BIN_DIR)

    s1c = info["root_cluster"]["rep_vector_c"]
    s1h = info["root_cluster"]["rep_vector_h"]

    wp_bond = BinaryManager("group", topdir=BIN_DIR)

    for no in info["id_set"]["PG"]["all"]:
        group_tag = info["tag"][no]
        print("creating", group_tag, flush=True)
        dic = group[no]
        bond, site_bond = bond_mapping_g(dic, group_tag, s1c, s1h)
        dic["wyckoff"]["bond"] = bond
        for wp, lst in site_bond.items():
            dic["wyckoff"]["site"][wp]["bond"] = lst
        wp_bond[no] = dic

    for no in info["id_set"]["SG"]["all"]:
        group_tag = info["tag"][no]
        print("creating", group_tag, flush=True)
        dic = group[no]
        bond, site_bond = bond_mapping_g(dic, group_tag, s1c, s1h)
        dic["wyckoff"]["bond"] = bond
        for wp, lst in site_bond.items():
            dic["wyckoff"]["site"][wp]["bond"] = lst
        wp_bond[no] = dic

    wp_bond.save_binary("group")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_group_bond()
