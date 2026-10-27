# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create PDF for atomic multipole for group.
"""

import os
import numpy as np
import sympy as sp

from multipie.core.multipie_info import __top_dir__
from multipie.core.group import Group
from multipie.util.util import timer, to_latex
from multipie.util.util_binary import BinaryManager
from multipie.util.util_pdf_latex import PDFviaLaTeX

pdfdir = os.path.join(__top_dir__, "others/pdf/")


# ==================================================
def create_response_tensor_each(tag, h_dir):
    group = Group(tag)

    if group.info.type == "II":
        xl = ["Q", "G"]
    else:
        xl = ["Q", "G", "T", "M"]

    for X in xl:
        file = f"response_tensor_{X}"
        title = group.latex(True) + f" [{X} tensor]"
        pdf = PDFviaLaTeX(file, dir=h_dir, pt=8, style="narrowest", landscape=True, english=True)
        pdf.title(title)
        d = group.response_tensor_all(X)
        empty = True
        for t, (m, ex) in d.items():
            X, rank, opt = t
            if opt == "":
                pdf.text(f"* Rank {rank} tensor.")
            else:
                pdf.text(f"* Rank {rank} tensor ({opt}).")
            if not np.all(m == sp.S(0)):
                empty = False
                ml = to_latex(m, "matrix")
                pdf.equation(ml)
                eqs = [sp.latex(sp.Eq(i, j)) for i, j in ex.items()]
                pdf.equation(eqs)
        if empty:
            pdf.text()

        pdf.build()


# ==================================================
@timer
def create_response_tensor():
    info = BinaryManager("info")

    for id_s in info["id_set"]["MPG"]["all"]:
        no = id_s.split(":")[1]
        h_dir = os.path.join(__top_dir__, f"{pdfdir}MPG/{no}")
        create_response_tensor_each(id_s, h_dir)


# ================================================== main
if __name__ == "__main__":
    create_response_tensor()
