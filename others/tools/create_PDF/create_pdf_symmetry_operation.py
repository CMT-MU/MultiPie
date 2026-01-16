"""
Create PDF for symmetry operation.
"""

import os

from multipie.core.multipie_info import __top_dir__
from multipie.core.group import Group
from multipie.util.util import timer, to_latex
from multipie.util.util_binary import BinaryManager
from multipie.util.util_pdf_latex import PDFviaLaTeX

pdfdir = os.path.join(__top_dir__, "others/pdf/")


# ==================================================
def create_symmetry_operation_each(tag, h_dir):
    group = Group(tag)
    file = f"symmetry_operation"
    title = group.latex(True)
    SO = group.symmetry_operation
    g_type = group.group_type

    pdf = PDFviaLaTeX(file, dir=h_dir, pt=8, style="narrowest", english=True)
    pdf.title(title)

    if g_type in ["PG", "SG"]:
        gen = [group.tag_symmetry_operation(i, True) for i in SO["generator"]]
        generator = "$" + r",\,\,".join(gen) + "$"
        pdf.text("\n* generator : " + generator + "\n")

    if g_type in ["PG"]:
        character = group.character["conjugacy"]
        cc = [[group.tag_symmetry_operation(i, True) for i in j] for j in character]

        irop = [["[" + j[0] + r"]:", r",\,\,\,".join(j)] for j in cc]
        pdf.text("* conjugacy class")
        pdf.text(r"\begin{quote}")
        pdf.simple_table(irop, tmath=True, cpos="rl")
        pdf.text(r"\end{quote}" + "\n")

    ss = r"* symmetry operation"
    if g_type in ["SG"]:
        ss += r"\quad" + r",\quad ".join(["$+" + to_latex(i, "vector") + "$" for i in SO["plus_set"]])
    pdf.text(ss)

    ops = [group.tag_symmetry_operation(i, True) for i in SO["tag"]]
    if g_type in ["PG", "MPG"]:
        mat = [to_latex(i, "matrix") for i in SO["fractional"]]
    else:
        mat = [to_latex(i[0:3, :], "matrix") for i in SO["fractional"]]
    det = [str(i) for i in SO["det"]]
    if g_type in ["MPG", "MSG"]:
        tr = [str(i) for i in SO["tr_sign"]]

    row = [str(i + 1) for i in range(len(ops))]
    if g_type in ["PG", "SG"]:
        tbl = [list(i) for i in zip(ops, mat, det)]
        col = ["tag", "matrix (polar)", "det"]
        cpos = "cccc"
    else:
        tbl = [list(i) for i in zip(ops, mat, det, tr)]
        col = ["tag", "matrix (polar)", "det", "TR"]
        cpos = "ccccc"

    pdf.table(
        tbl,
        row,
        col,
        "No.",
        rmath=True,
        tmath=True,
        stretch=1.3,
        caption="Symmetry operations for 3d polar vector.",
        cpos=cpos,
        long=True,
        hl=True,
    )

    pdf.build()


# ==================================================
@timer
def create_symmetry_operation():
    info = BinaryManager("info")

    for id_s in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        no = int(id_s.split(":")[1])
        tag = info["tag"][id_s]
        h_dir = os.path.join(__top_dir__, f"{pdfdir}PG/{no:03d}-{tag}")
        create_symmetry_operation_each(id_s, h_dir)

    for id_s in info["id_set"]["SG"]["all"]:
        no = int(id_s.split(":")[1])
        tag = info["tag"][id_s]
        h_dir = os.path.join(__top_dir__, f"{pdfdir}SG/{no:03d}-{tag.replace("^", "_")}")
        create_symmetry_operation_each(id_s, h_dir)

    for id_s in info["id_set"]["MPG"]["all"]:
        no = id_s.split(":")[1]
        h_dir = os.path.join(__top_dir__, f"{pdfdir}MPG/{no}")
        create_symmetry_operation_each(id_s, h_dir)

    for id_s in info["id_set"]["MSG"]["all"]:
        no = id_s.split(":")[1]
        h_dir = os.path.join(__top_dir__, f"{pdfdir}MSG/{no}")
        create_symmetry_operation_each(id_s, h_dir)


# ================================================== main
if __name__ == "__main__":
    create_symmetry_operation()
