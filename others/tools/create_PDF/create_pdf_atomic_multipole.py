"""
Create PDF for atomic multipole.
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import os

from multipie.util.util_binary import BinaryManager
from multipie.util.util import to_latex, timer
from multipie.util.util_pdf_latex import PDFviaLaTeX
from multipie.util.util_tag import TagMultipole, TagBasis

pdfdir = os.path.join(__top_dir__, "others/pdf/misc")


# ==================================================
def create_atomic_multipole_each(am, info, b_type):
    rank = {0: "s", 1: "p", 2: "d", 3: "f"}
    if b_type == "lm":
        title = "Atomic Multipoles (spinless LM basis)"
        ofile = "atomic_multipole_lm"
        basis = info["harmonics"]["atomic_basis"]["spinless"][b_type]
    else:
        title = "Atomic Multipoles (spinful JML basis)"
        ofile = "atomic_multipole_jml"
        basis = info["harmonics"]["atomic_basis"]["spinful"][b_type]

    pdf = PDFviaLaTeX(ofile, pt=8, dir=pdfdir, landscape=True, english=True)
    pdf.title(title)

    cnt = 1
    for (bra, ket), samb in am[b_type].items():
        samb1 = {}
        for (head, l, s, k, _), mat in samb.items():
            for m, mm in zip(range(l, -l - 1, -1), mat):
                samb1[(head, l, s, k, "q", m)] = mm
        tbl = [
            [TagMultipole.latex(tag[:-1], "spherical", tag[-1], tag="a", vector=False), to_latex(mat, "matrix")]
            for tag, mat in samb1.items()
        ]
        row = list(range(cnt, cnt + len(samb1)))
        cnt += len(samb1)

        bras = ", ".join([TagBasis.latex(i, bra) for i in basis[bra]])
        kets = ", ".join([TagBasis.latex(i, ket) for i in basis[ket]])
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
            hl=True,
            caption=f"({rank[bra]},{rank[ket]}) block.",
        )

    pdf.build()


# ==================================================
@timer
def create_atomic_multipole():
    info = BinaryManager("info", topdir=BIN_DIR)
    am = BinaryManager("atomic_multipole", topdir=BIN_DIR)
    create_atomic_multipole_each(am, info, "lm")
    create_atomic_multipole_each(am, info, "jml")


# ================================================== main
if __name__ == "__main__":
    create_atomic_multipole()
