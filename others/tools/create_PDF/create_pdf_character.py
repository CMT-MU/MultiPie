"""
Create PDF for character.
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
def create_character_each(tag, h_dir):
    group = Group(tag)
    file = f"character_table"
    title = group.latex(True)
    character = group.character

    rc = group.latex()
    col = [group.tag_symmetry_operation(i[0], True) + f"({len(i)})" for i in character["conjugacy"]]
    row = [group.tag_irrep(i, True) for i in character["table"].keys()]
    wp, wm = sp.symbols(r"\omega \omega^*")
    sub = {"wp": wp, "wm": wm}
    wpc = 0
    for i in character["table"].values():
        for c in i:
            wpc += str(c).count("wp")
    tbl = [[sp.latex(c.subs(sub)) for c in i] for i in character["table"].values()]
    align = "c" + "r" * len(col)

    pdf = PDFviaLaTeX(file, dir=h_dir, pt=8, style="narrowest", landscape=True, english=True)
    pdf.title(title)

    tt = "\n* character table"
    if wpc > 0:
        tt += r" ($\omega=e^{2\pi i/3}$)"
    pdf.text(tt)
    pdf.text(r"\begin{quote}")
    pdf.table(
        tbl,
        row,
        col,
        rc,
        rcmath=True,
        rmath=True,
        cmath=True,
        tmath=True,
        cpos=align,
    )
    pdf.text(r"\end{quote}")

    pdf.text("\n" + r"* polar $\leftrightarrow$ axial conversion")
    pdf.text(r"\begin{quote}")
    tbl = [
        [
            group.tag_irrep(i, True) + "\\,\\,(" + group.tag_irrep(j, True) + ")"
            for i, j in character["polar_axial_conversion"].items()
        ]
    ]
    pdf.simple_table(tbl, tmath=True)
    pdf.text(r"\end{quote}")

    tbl = []
    for i, ir1 in enumerate(character["table"].keys()):
        t1 = []
        for j, ir2 in enumerate(character["table"].keys()):
            v = (
                sp.latex(sum([c * sp.symbols(group.tag_irrep(e, True)) for c, e in character["symmetric_product"][(ir1, ir2)]]))
                if j >= i
                else ""
            )
            t1.append(v)
        tbl.append(t1)

    cp = "c|" + "c" * len(row)
    pdf.text("\n* symmetric product")
    pdf.text(r"\begin{quote}")
    pdf.table(tbl, row, row, rmath=True, cmath=True, tmath=True, cpos=cp)
    pdf.text(r"\end{quote}")

    tbl = []
    for ir in character["table"].keys():
        v = sp.latex(sum([c * sp.symbols(group.tag_irrep(e, True)) for c, e in character["anti_symmetric_product"][ir]]))
        if v == "0":
            v = "-"
        tbl.append(v)
    pdf.text("\n* anti-symmetric product")
    pdf.text(r"\begin{quote}")
    pdf.table([tbl], [""], row, rmath=True, cmath=True, tmath=True)
    pdf.text(r"\end{quote}")

    pdf.build()


# ==================================================
@timer
def create_character():
    info = BinaryManager("info")

    for id_s in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        no = int(id_s.split(":")[1])
        tag = info["tag"][id_s]
        h_dir = os.path.join(__top_dir__, f"{pdfdir}PG/{no:03d}-{tag}")
        create_character_each(tag, h_dir)


# ================================================== main
if __name__ == "__main__":
    create_character()
