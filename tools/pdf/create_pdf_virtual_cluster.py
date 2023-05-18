"""
create virtual cluster basis set.
"""
from gcoreutils.pdf_via_latex import PDFviaLaTeX
from gcoreutils.nsarray import NSArray
from multipie.tag.tag_group import TagGroup
from multipie.virtual_cluster.virtual_cluster_pg_set import VirtualClusterPGSet
from multipie import __top_dir__, get_binary

point_group_dir = __top_dir__ + "docs/pdf/point_group"


# ==================================================
def create_virtual_cluster(core):
    vcset = core[VirtualClusterPGSet]

    pdf = PDFviaLaTeX("vc_basis", pt=9, landscape=True, dir=point_group_dir)

    ppg = TagGroup.create()

    for tag in ppg:
        pg = vcset[tag]
        sites = [pg.site.latex()]

        pdf.title(tag.info(latex=True) + " (orthonormalized)")

        pdf.text(r"\vspace{5mm}")
        pdf.text("* site positions (reduced coordinate)")
        pdf.text(r"\begin{quote}")
        pdf.simple_table(sites, cols=5, tmath=True)
        pdf.text(r"\end{quote}")

        rc = "type"
        col = [str(i) for i in range(1, 11)]
        cpos = "l" + "c" * len(col)

        row = []
        tbl = []
        for i, (info, basis) in enumerate(pg.items()):
            row.append("\#" + str(i + 1) + r"\quad " + info.latex())
            basis = NSArray(basis.tolist(), "scalar")
            tbl.append(basis.latex())

        pdf.text("* site basis")
        pdf.table(tbl, row, col, rc, hl=True, cpos=cpos, rmath=True, tmath=True, stretch=1.3, long=True)

        pdf.text("\\newpage")


# ================================================== main
core = get_binary()
create_virtual_cluster(core)
