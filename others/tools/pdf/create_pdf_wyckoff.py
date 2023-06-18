"""
create Wyckoff positions.
"""
from gcoreutils.pdf_via_latex import PDFviaLaTeX
from multipie.tag.tag_group import TagGroup
from multipie.wyckoff.wyckoff_g_set import WyckoffGSet
from multipie.symmetry_operation.symmetry_operation_g_set import SymmetryOperationGSet
from multipie import get_binary

point_group_dir = __file__[: __file__.rfind("/")] + "/../../../docs/pdf/point_group"
space_group_dir = __file__[: __file__.rfind("/")] + "/../../../docs/pdf/space_group"


# ==================================================
def create_wyckoff(wset, soset, t_pg):
    outdir = point_group_dir if t_pg else space_group_dir
    outstr = "point" if t_pg else "space"

    pdf = PDFviaLaTeX("wyckoff", landscape=True, pt=9, dir=outdir)
    pdf.title("Wyckoff position (" + outstr + " group)")

    if t_pg:
        pgset = TagGroup.create()
    else:
        pgset = TagGroup.create(space_group=True)

    for tag in pgset:
        p = wset[tag]
        name = tag.info(latex=True)
        if not t_pg:
            so = soset[tag]
            name += r"\quad" + ",\quad ".join(["$+" + i + "$" for i in so.plus_set.latex()])
        row = []
        tbl = []
        for w in p.keys():
            sym = p.site_symmetry(w)
            row.append(r"{\tt " + str(w) + "}" + r" ({\tt " + sym + "})")
            tbl.append(p.position(w).latex())

        pdf.text(name)
        col = [str(i + 1) for i in range(6)]
        pdf.table(tbl, row, col, "WL (SS)", hl=True, stretch=1.2, tmath=True, long=True)

        if t_pg:
            tbl = []
            for w in p.keys():
                wp = p.position(w, default=True)[0].latex()
                tbl += [r"{\tt " + str(w) + r"}", wp]

            pdf.text("* default ($x$,$y$,$z$)")
            pdf.text(r"\begin{quote}")
            pdf.simple_table(tbl, cols=10, tmath=True)
            pdf.text(r"\end{quote}")

        pdf.text(r"\newpage")


# ================================================== main
core = get_binary()
wset = core[WyckoffGSet]
soset = core[SymmetryOperationGSet]

create_wyckoff(wset, soset, t_pg=True)
create_wyckoff(wset, soset, t_pg=False)
