"""
Create binary data (information of SO mapping for each groups).
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import logging
import numpy as np

from multipie.util.util_binary import BinaryManager


# ==================================================
def create_so_mapping(m_tag, m_prod, c_tag):
    """
    Create SO mapping from mother group.

    Args:
        m_tag (list): SO tag of mother group.
        m_prod (dict): product table of mother group.
        c_tag (list): SO tag of target group.

    Returns:
        - (dict) -- SO mapping.
    """
    diff = [so for so in m_tag if so not in c_tag]
    lst = []
    for i in diff:
        s = []
        for j in c_tag:
            mc_tag = m_prod[(j, i)]
            if mc_tag not in sum(lst, []):
                s.append(mc_tag)
        if s:
            lst.append(s)

    lst = [c_tag] + lst
    lst = np.asarray(lst).T.tolist()
    dic = {t: mp for t, mp in zip(c_tag, lst)}

    return dic


# ==================================================
def create_group_mapping():
    info = BinaryManager("info", topdir=BIN_DIR)
    group = BinaryManager("group", topdir=BIN_DIR)

    Oh = info["id"]["Oh"]
    D6h = info["id"]["D6h"]

    oh_tag = group[Oh]["symmetry_operation"]["tag"]
    oh_prod = group[Oh]["symmetry_operation"]["product"]
    d6h_tag = group[D6h]["symmetry_operation"]["tag"]
    d6h_prod = group[D6h]["symmetry_operation"]["product"]

    g_map = BinaryManager("group", topdir=BIN_DIR)

    for no in info["id_set"]["PG"]["all"]:
        dic = group[no]
        hexagonal = dic["info"].hexagonal_g
        m_tag = d6h_tag if hexagonal else oh_tag
        m_prod = d6h_prod if hexagonal else oh_prod
        c_tag = dic["symmetry_operation"]["tag"]
        dic["symmetry_operation"]["mapping"] = create_so_mapping(m_tag, m_prod, c_tag)
        g_map[no] = dic

    g_map.save_binary("group")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_group_mapping()
