"""
Create PDF for harmonics.
"""

import os
import sympy as sp

from multipie.core.multipie_info import __top_dir__
from multipie.core.group import Group
from multipie.util.util import timer
from multipie.util.util_binary import BinaryManager
from multipie.util.util_pdf_latex import PDFviaLaTeX

pdfdir = os.path.join(__top_dir__, "others/pdf/")


# ==================================================
def create_harmonics_each(tag, h_dir, axial, definition=False):
    if tag.count("^") > 0:
        return

    group = Group(tag)
    p_type = "axial" if axial else "polar"
    X = "G" if axial else "Q"
    file = f"harmonics_{p_type}"
    if definition:
        file = f"harmonics_{tag}_{p_type}"
    else:
        file = f"harmonics_{p_type}"
    title = group.latex(True)

    pdf = PDFviaLaTeX(file, dir=h_dir, pt=8, style="narrowest", landscape=not definition, english=True)
    pdf.title(title)

    rc = "No."
    if definition:
        col = ["multipole", "definition"]
    else:
        col = ["multipole", "expression"]

    no = 1
    for l in range(12):
        cap = f"Harmonics for rank {l}."
        harm = group.harmonics.select(X=X, l=l)
        tbl = []
        row = []
        for idx, (cset, _, tset) in harm.items():
            m_tags = group.tag_multipole(idx, latex=True)
            for m_tag, t, c in zip(m_tags, tset, cset):
                if definition:
                    tbl.append([m_tag, sp.latex(t)])
                else:
                    tbl.append([m_tag, sp.latex(c)])
                row.append(str(no))
                no += 1
        pdf.table(tbl, row, col, rc, cpos="ccl", caption=cap, stretch=1.3, rmath=True, tmath=True)

    pdf.build()


# ==================================================
@timer
def create_harmonics():
    # definition of Oh and D6h.
    h_dir = f"{pdfdir}misc"
    create_harmonics_each("Oh", h_dir, axial=False, definition=True)
    create_harmonics_each("Oh", h_dir, axial=True, definition=True)
    create_harmonics_each("D6h", h_dir, axial=False, definition=True)
    create_harmonics_each("D6h", h_dir, axial=True, definition=True)

    # expression.
    info = BinaryManager("info")

    for id_s in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        no = int(id_s.split(":")[1])
        tag = info["tag"][id_s]
        h_dir = os.path.join(__top_dir__, f"{pdfdir}PG/{no:03d}-{tag}")
        create_harmonics_each(tag, h_dir, axial=False)
        create_harmonics_each(tag, h_dir, axial=True)


# ================================================== main
if __name__ == "__main__":
    create_harmonics()
