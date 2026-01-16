"""
Create PDF for representation matrix.
"""

import os

from multipie.core.multipie_info import __top_dir__
from multipie.core.group import Group
from multipie.util.util import timer, to_latex
from multipie.util.util_binary import BinaryManager
from multipie.util.util_pdf_latex import PDFviaLaTeX

pdfdir = os.path.join(__top_dir__, "others/pdf/")


# ==================================================
def create_representation_matrix_each(tag, h_dir):
    group = Group(tag, with_opt=True)
    file = f"representation_matrix"
    title = group.latex(True)
    rep_matrix = group.opt["representation_matrix"]
    rep = rep_matrix["matrix"]
    so_tag = rep_matrix["tag"]

    row = []
    tbl = []
    for irrep, mat in rep.items():
        row.append(group.tag_irrep(irrep, True))
        tbl.append([group.tag_symmetry_operation(t, True) + ":" + to_latex(i, "matrix") for t, i in zip(so_tag, mat)])

    pdf = PDFviaLaTeX(file, dir=h_dir, pt=8, landscape=True, english=True)
    pdf.title(title)
    col = ["" for i in range(6)]

    pdf.table(tbl, row, col, "Irrep.", hl=True, stretch=1.2, rmath=True, tmath=True, long=True, caption="Representation matrices")

    pdf.build()


# ==================================================
@timer
def create_representation_matrix():
    info = BinaryManager("info")

    for id_s in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        no = int(id_s.split(":")[1])
        tag = info["tag"][id_s]
        h_dir = os.path.join(__top_dir__, f"{pdfdir}PG/{no:03d}-{tag}")
        create_representation_matrix_each(tag, h_dir)


# ================================================== main
if __name__ == "__main__":
    create_representation_matrix()
