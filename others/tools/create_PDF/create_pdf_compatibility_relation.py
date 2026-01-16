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
            [0, 7, 18, 21, 24],
            r"$O_{\rm h}$ subgroup.",
        ),
        (
            "D6h",
            [
                "s",
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
            [0, 11],
            r"$D_{\rm 6h}$ subgroup.",
        ),
    ]

    pdf = PDFviaLaTeX("compatibility_relation", dir=h_dir, pt=9, landscape=True, style="narrowest", english=True)
    pdf.title("Compatibility Relation for Point Groups")

    for rc, ex, hl, cap in data:
        col = [to_latex(i[0]) for i in cr[(rc, rc)]]
        row = [""]
        tbl = [[r"\texttt{" + i + "}" for i in ex]]
        for k, v in cr.items():
            if k[0] == rc:
                row.append(tag_latex(k[1]))
                tbl.append([to_latex(i[1]) for i in v])

        pdf.text(r"\begin{center}")
        pdf.table(
            tbl,
            row,
            col,
            tag_latex(rc),
            caption=cap,
            stretch=1.0,
            long=True,
            cmath=False,
            tmath=False,
            hl=hl,
        )
        pdf.text(r"\end{center}")
        pdf.text(r"\newpage")

    pdf.build()


# ================================================== main
if __name__ == "__main__":
    create_compatibility_relation()
