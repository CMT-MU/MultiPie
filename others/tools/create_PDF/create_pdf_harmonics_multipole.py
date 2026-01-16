"""
Create PDF for multipolar harmonics.
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
def create_harmonics_each(tag, h_dir, irank, axial, iaxial, max_rank):
    group = Group(tag, with_opt=True)

    i_type = "axial" if axial else "polar"
    ii_type = "axial" if iaxial else "polar"
    head = "G" if axial else "Q"
    ihead = "g" if iaxial else "q"
    mrank = {1: "dipole", 2: "quadrupole", 3: "octupole"}
    file = f"harmonics_s{irank}_" + i_type + "_" + ihead
    title = group.latex(True) + " (" + i_type + ", internal " + ii_type + " " + mrank[irank] + ")"
    harmonics = group.opt["harmonics_multipole"]

    pdf = PDFviaLaTeX(file, dir=h_dir, pt=8, style="narrowest", english=True)
    pdf.title(title)

    for rank in range(max_rank + 1):
        pdf.text("\n" + r"\noindent * Harmonics for rank " + str(rank) + "\n")
        harm = harmonics.select(s=irank, x=ihead, X=head, l=rank)
        for idx, (cset, _, fset) in harm.items():
            if len(fset) > 0:
                m_tags = group.tag_multipole(idx, latex=True, vector=True, internal=True)
                pdf.text(r"\noindent\quad " + "$" + "$, $".join(m_tags) + "$\n")
                pdf.text(r"\noindent\quad " + f"** symmetry")
                for t in cset:
                    pdf.equation(sp.latex(t), long=True)
                pdf.text(r"\noindent\quad " + f"** expression")
                for t in fset:
                    pdf.equation(sp.latex(t), long=True)

    pdf.build()


# ==================================================
@timer
def create_harmonics_multipole():
    max_s_rank = 3
    max_rank = 4

    info = BinaryManager("info")

    for id_s in info["id_set"]["PG"]["all"]:
        no = int(id_s.split(":")[1])
        tag = info["tag"][id_s]
        h_dir = os.path.join(__top_dir__, f"{pdfdir}PG/{no:03d}-{tag}")

        for s in range(1, max_s_rank + 1):
            create_harmonics_each(tag, h_dir, irank=s, axial=False, iaxial=False, max_rank=max_rank)
            create_harmonics_each(tag, h_dir, irank=s, axial=False, iaxial=True, max_rank=max_rank)
            create_harmonics_each(tag, h_dir, irank=s, axial=True, iaxial=False, max_rank=max_rank)
            create_harmonics_each(tag, h_dir, irank=s, axial=True, iaxial=True, max_rank=max_rank)


# ================================================== main
if __name__ == "__main__":
    create_harmonics_multipole()
