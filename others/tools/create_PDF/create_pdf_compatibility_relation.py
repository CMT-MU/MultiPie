"""
Create PDF for compatibility relation.
"""

import os

from multipie.core.multipie_info import __top_dir__
from multipie.util.util import timer
from multipie.util.util_binary import BinaryManager
from multipie.util.util_pdf_latex import PDFviaLaTeX
from multipie.util.util_tag import TagIrrep

h_dir = os.path.join(__top_dir__, "others/pdf/info")


# ==================================================
@timer
def create_compatibility_relation():
    def to_latex(s):
        if s.count("("):
            ir, p = s.split("(")
            if ir.count("-"):
                s = "-" + TagIrrep.latex(ir[1:]) + f"({p}"
            else:
                s = TagIrrep.latex(ir) + f"({p}"
        else:
            s = TagIrrep.latex(s)
        return "$" + s + "$"

    def tag_latex(s):
        if s.count("-"):
            t, p = s.split("-")
            s = "$" + t[0] + r"_{\rm " + t[1:] + "}$(" + p + ")"
        else:
            s = "$" + s[0] + r"_{\rm " + s[1:] + "}$"
        return s

    info = BinaryManager("info")
    cr = info["character"]["compatibility_relation"]
    data = [
        (
            "Oh",
            [
                "s",
                "i",
                "du",
                "dv",
                "gax",
                "gay",
                "gaz",
                "dyz",
                "dxz",
                "dxy",
                "l",
                "f3",
                "hu",
                "hv",
                "px",
                "py",
                "pz",
                "fbx",
                "fby",
                "fbz",
            ],
            [
                "l",
                "f3",
                "hu",
                "hv",
                "px",
                "py",
                "pz",
                "fbx",
                "fby",
                "fbz",
                "s",
                "i",
                "du",
                "dv",
                "gax",
                "gay",
                "gaz",
                "dyz",
                "dxz",
                "dxy",
            ],
            [1, 8, 19, 22, 25],
            r"$O_{\rm h}$ subgroup.",
            "c|cccccccccc|cccccccccc",
        ),
        (
            "D6h",
            [
                "s/du",
                "i0",
                "ga",
                "gb",
                "dyz",
                "dxz'",
                "dv",
                "dxy'",
                "j",
                "pz",
                "f1",
                "f2",
                "px",
                "py",
                "f3",
                "fbz",
            ],
            [
                "j",
                "pz",
                "f1",
                "f2",
                "px",
                "py",
                "f3",
                "fbz",
                "s/du",
                "i0",
                "ga",
                "gb",
                "dyz",
                "dxz'",
                "dv",
                "dxy'",
            ],
            [1, 12],
            r"$D_{\rm 6h}$ subgroup.",
            "c|cccccccc|cccccccc",
        ),
    ]

    pdf = PDFviaLaTeX("compatibility_relation", dir=h_dir, pt=9, landscape=True, style="narrowest", english=True)
    pdf.title("Compatibility Relation for Point Groups")

    for rc, ex_q, ex_g, hl, cap, cpos in data:
        col = [to_latex(i[0]) for i in cr[(rc, rc)]]
        row = ["Q", "G"]
        tbl = [[r"\texttt{" + i + "}" for i in ex_q], [r"\texttt{" + i + "}" for i in ex_g]]
        for k, v in cr.items():
            if k[0] == rc:
                row.append(tag_latex(k[1]))
                tbl.append([to_latex(i[1]) for i in v])

        pdf.table(tbl, row, col, tag_latex(rc), caption=cap, stretch=1.0, cmath=False, tmath=False, hl=hl, cpos=cpos)
        pdf.text(r"\newpage")

    pdf.build()


# ================================================== main
if __name__ == "__main__":
    create_compatibility_relation()
