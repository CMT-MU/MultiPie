"""
create virtual-cluster basis.
"""
import os
from multipie.group.point_group import PointGroup
from multipie.tag.tag_group import TagGroup
from tools.data.util.output_util import write_data

from multipie import __top_dir__, get_binary

create_dir = __top_dir__ + "tools/data"
core = get_binary()

# ==================================================
header_basis = """
basis of virtual cluster for other wyckoff positions (orthonormalized).
    { point-group tag : { wp: {multipole-tag: basis} } }
real version
"""

header_site = """
site positions of virtual cluster for other wyckoff positions (reduced coordinate).
    { point-group tag : { wp: sites} }
real version
"""


# ==================================================
def create_vc_basis(output_dir=None):
    if output_dir is not None:
        os.chdir(output_dir)

    ofile = "data_virtual_cluster_real.py"

    print("=== create data_virtual_cluster_real ===")
    basis = {}
    site = {}
    for t in TagGroup.create():
        print("creating ...", t)
        pg = PointGroup(t, core=core)
        wp = list(map(str, pg.wyckoff.keys()))[::-1]
        basis[str(pg)] = {}
        site[str(pg)] = {}
        for wpi in wp:
            b = pg.virtual_cluster_basis(wyckoff=wpi, ortho=True, create=True)
            basis[str(pg)][wpi] = {str(k).replace("h", "s"): str(v).replace(" ", "") for k, v in b[0].items()}
            site[str(pg)][wpi] = str(b[1]).replace(" ", "")

    write_data(ofile, basis, header_basis, "_data_virtual_cluster_basis_real", mode="w")
    write_data(ofile, site, header_site, "_data_virtual_cluster_site_real", mode="a")


# ================================================== main
create_vc_basis(output_dir=create_dir)
