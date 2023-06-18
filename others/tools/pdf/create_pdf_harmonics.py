"""
create harmonics.
"""
from gcoreutils.pdf_via_latex import PDFviaLaTeX
from gcoreutils.nsarray import NSArray
from multipie.harmonics.harmonics_pg_set import HarmonicsPGSet
from multipie import get_binary

harmonics_dir = __file__[: __file__.rfind("/")] + "/../../../docs/pdf/harmonics"


# ==================================================
def create_harmonics(hs32, axial):
    i_type = "axial" if axial else "polar"
    h_dir = harmonics_dir + "/" + i_type

    if axial:
        head = "G"
    else:
        head = "Q"

    for tag, pg in hs32.items():
        name = tag.info(latex=True)

        pdf = PDFviaLaTeX("harmonics_" + str(tag), dir=h_dir, pt=8, style="narrowest")
        pdf.title(name + " (" + i_type + ")")
        pdf2 = PDFviaLaTeX("harmonics_" + str(tag) + "_ex", dir=h_dir, pt=8, landscape=True, style="narrowest")
        pdf2.title(name + " (" + i_type + ")")

        no = 1
        for rank in range(12):
            tbl = []
            tbl2 = []
            row = []
            hs = pg.select(rank=rank, head=head)
            for h in hs:
                mul = h.tag.mul
                comp = h.tag.comp
                irrep = h.tag.tag_irrep()
                dim = irrep.dim
                if mul < 1:
                    mul = "-"
                if dim == 1:
                    comp = "-"
                tbl.append(
                    [
                        irrep.latex(),
                        r"{\tt " + str(irrep) + "}",
                        mul,
                        comp,
                        h.latex(),
                        r"{\tt " + str(h) + "}",
                        h.definition().latex(),
                    ]
                )
                row.append(str(no))
                tbl2.append([h.expression(v=NSArray.vector3d(head)).latex()])
                no += 1

            cap = "rank " + str(rank)
            rc = "No."
            col = ["irrep.", "(tag)", "mul.", "comp.", "harmonics", "(tag)", "definition"]
            pdf.table(tbl, row, col, rc, caption=cap, stretch=1.3, rmath=True, tmath=True)
            col2 = ["expression"]
            pdf2.table(tbl2, row, col2, rc, cpos="cl", caption=cap, stretch=1.3, rmath=True, tmath=True)


# ================================================== main
core = get_binary()
hs32 = core[HarmonicsPGSet]

create_harmonics(hs32, axial=False)
create_harmonics(hs32, axial=True)
