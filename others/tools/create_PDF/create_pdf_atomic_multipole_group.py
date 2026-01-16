"""
Create PDF for atomic multipole for group.
"""

import os

from multipie.core.multipie_info import __top_dir__
from multipie.core.group import Group
from multipie.util.util_binary import BinaryManager
from multipie.util.util import timer, to_latex
from multipie.util.util_pdf_latex import PDFviaLaTeX
from multipie.util.util_tag import TagMultipole

pdfdir = os.path.join(__top_dir__, "others/pdf/")


# ==================================================
def create_atomic_multipole_each(tag, b_type, h_dir):
    rank = {0: "s", 1: "p", 2: "d", 3: "f"}
    group = Group(tag)
    file = f"atomic_multipole_" + b_type
    title = group.latex(True) + " (" + b_type + " basis)"

    pdf = PDFviaLaTeX(file, dir=h_dir, pt=8, style="narrowest", landscape=True, english=True)
    pdf.title(title)

    pg = group._group()["atomic_samb"][b_type]
    basis = group.atomic_basis(b_type)

    cnt = 1
    for (bra, ket), samb in pg.items():
        samb1 = {}
        for (head, l, Gamma, n, p, s, k, _), (mat, ex) in samb.items():
            if len(mat) == 1:
                samb1[(head, l, Gamma, n, p, s, k, "q", -1)] = (ex[0], mat[0])
            else:
                for m, (e, mm) in enumerate(zip(ex, mat)):
                    samb1[(head, l, Gamma, n, p, s, k, "q", m)] = (e, mm)
        tbl1 = [
            [TagMultipole.latex(m_tag[:-1], "point_group", m_tag[-1], tag="a", vector=False), to_latex(mat, "matrix")]
            for m_tag, (ex, mat) in samb1.items()
        ]
        tbl2 = [[r"\text{symmetry}", to_latex(ex, "scalar")] for ex, mat in samb1.values()]
        d = len(samb1)
        tbl = []
        for i, j in zip(tbl2, tbl1):
            tbl.append(i)
            tbl.append(j)
        row1 = list(range(cnt, cnt + d))
        row = []
        for i in row1:
            row.append(i)
            row.append("")
        cnt += d
        hl = list(range(1, len(row), 2))

        bras = ", ".join([group.tag_atomic_basis(i, bra, latex=True, ket=False) for i in basis[bra]])
        kets = ", ".join([group.tag_atomic_basis(i, ket, latex=True, ket=True) for i in basis[ket]])
        pdf.text("bra: $" + "=" + bras + "$\n\n\\noindent")
        pdf.text("ket: $" + "=" + kets + "$")
        pdf.table(
            tbl,
            row,
            ["multipole", "matrix"],
            "No.",
            rmath=True,
            tmath=True,
            long=True,
            stretch=1.6,
            hl=hl,
            caption=f"({rank[bra]},{rank[ket]}) block.",
        )

    pdf.build()


# ==================================================
@timer
def create_atomic_multipole_group():
    info = BinaryManager("info")
    for id_s in info["id_set"]["PG"]["all"]:
        no = int(id_s.split(":")[1])
        tag = info["tag"][id_s]
        print(f"creating for {tag}")
        h_dir = os.path.join(__top_dir__, f"{pdfdir}PG/{no:03d}-{tag}")
        create_atomic_multipole_each(tag, "lg", h_dir)
        create_atomic_multipole_each(tag, "lgs", h_dir)
        create_atomic_multipole_each(tag, "jml", h_dir)


# ================================================== main
if __name__ == "__main__":
    create_atomic_multipole_group()
