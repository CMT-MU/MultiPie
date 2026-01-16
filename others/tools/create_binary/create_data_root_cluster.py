"""
Create binary data (information for root cluster).
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import logging
import numpy as np
import sympy as sp

from multipie.util.util import timer, str_to_sympy
from multipie.util.util_binary import BinaryManager

# ==================================================
h_root_cluster = """
* Information of root cluster.
- "position" (str): (ndarray) list of root-cluster sites (cartesian).
- "op_mapping" (str): (dict) operation-site mapping.
  - "Oh/D6h" (str): (dict) data.
    - SO_tag (str): (int) rc_site_index.
"""


# ==================================================
@timer
def create_root_cluster():
    group_data = BinaryManager("group", topdir=BIN_DIR)
    info = BinaryManager("info", topdir=BIN_DIR)
    rep_point = info["root_cluster"]["rep_site_c"]
    pg_set = info["root_cluster"]["point_group"]

    # create cluster site.
    cluster_info = {}
    cluster_tag = {}
    for group in pg_set:
        no = info["id"][group]
        g_info = group_data[no]
        ops = g_info["symmetry_operation"]["cartesian"]
        ops_tag = g_info["symmetry_operation"]["tag"]
        pos = (np.vectorize(sp.expand)(ops @ rep_point)).tolist()
        cluster_info[group] = pos
        cluster_tag[group] = dict(zip(ops_tag, pos))

    rs = [str(i) for i in sum(cluster_info.values(), [])]
    rs = np.asarray([str_to_sympy(i) for i in list(dict.fromkeys(rs))])

    # create SO and site mapping.
    op_site_mapping = {}
    for pg, dic in cluster_tag.items():
        mapping_pg = {}
        for name, p in dic.items():
            idx = int(np.where(np.all(rs == p, axis=1))[0][0])
            mapping_pg[name] = idx
        op_site_mapping[pg] = mapping_pg

    root_cluster = BinaryManager(verbose=True, topdir=BIN_DIR)
    root_cluster["position"] = rs
    root_cluster["op_mapping"] = op_site_mapping
    root_cluster.add_comment(h_root_cluster)

    root_cluster.save_binary("root_cluster")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    create_root_cluster()
