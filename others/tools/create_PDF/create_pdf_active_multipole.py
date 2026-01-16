"""
Create PDF for magnetic point group info.
"""

import os

from multipie.core.multipie_info import __top_dir__
from multipie.core.group import Group
from multipie.util.util import timer, to_latex
from multipie.util.util_binary import BinaryManager
from multipie.util.util_pdf_latex import PDFviaLaTeX

h_dir = os.path.join(__top_dir__, "others/pdf/info")


# ==================================================
def active_list(lst, c, basis):
    order = basis["hexagonal"] if c in ["trigonal", "hexagonal"] else basis["cubic"]

    check = ["" for _ in range(len(order))]
    for ri in lst:
        i = order.index(ri)
        check[i] = r"$\bullet$"
    return check


# ==================================================
@timer
def create_active_multipole():
    info = BinaryManager("info")
    basis = info["response_tensor"]["multipole"]
    basis_orbital = Group.global_info()["harmonics"]["basis_function"]

    pdf = PDFviaLaTeX("active_multipole", dir=h_dir, landscape=True, style="narrow", pt=9, english=True)

    pdf.title("Active Multipoles in Magnetic Point Groups")

    row = []
    tbl = []
    hl = [0, 3, 8, 19]
    for key, b in basis_orbital.items():
        row.append(r"\texttt{" + key + "}")
        tbl.append(["$" + to_latex(b[0]) + "$"])

    pdf.table(tbl, row, ["symmetry"], "tag", cpos="ll", long=True, stretch=1.6, caption="Multipoles in response tensor.", hl=hl)

    cap = "Active Multipole (cubic subgroups)"
    lst = [r"\texttt{" + i + "}" for i in basis["cubic"]]
    lbl = ["MPG", "Type", "Crystal", "X"] + lst
    cpos = "lllll|l|lll|ll|lll|l|lll|lll|l|ll|lll|lll"

    row = []
    tbl = []
    hl = []
    no = 0
    for id_s in info["id_set"]["MPG"]["all"]:
        group = Group(id_s)
        name = "$" + group.info.international + "$"
        active = group.active_multipole
        pg = group.info.PG.split(":")[1]
        mpg = id_s.split(":")[1]
        SS = "$" + Group(group.info.PG).info.schoenflies + "$"
        tp = group.info.type
        c = group.info.crystal
        if c in ["trigonal", "hexagonal"]:
            continue
        row.append(mpg)
        row.append("(" + pg + ")")
        row.append("")
        row.append("")
        tbl0 = [r"\texttt{" + name + "}", tp, c]
        d = {}
        for i in active:
            X = i[0]
            d[X] = d.get(X, []) + [i[1:]]
        tbl.append(tbl0 + ["Q"] + active_list(d.get("Q", []), c, basis))
        tbl.append([SS, "", "", "G"] + active_list(d.get("G", []), c, basis))
        tbl.append(["", "", "", "T"] + active_list(d.get("T", []), c, basis))
        tbl.append(["", "", "", "M"] + active_list(d.get("M", []), c, basis))
        no += 4
        hl.append(no - 1)

    pdf.table(tbl, row, lbl, r"\#", cpos=cpos, long=True, stretch=1.6, hl=hl[:-1], caption=cap)

    cap = "Active Multipole (hexagonal subgroups)"
    lst = [r"\texttt{" + i + "}" for i in basis["hexagonal"]]
    lbl = ["MPG", "Type", "Crystal", "X"] + lst
    cpos = "lllll|l|l|l|l|l|l|l|l|l|ll|ll|ll|ll|ll|ll|ll|ll"

    row = []
    tbl = []
    hl = []
    no = 0
    for id_s in info["id_set"]["MPG"]["all"]:
        group = Group(id_s)
        name = "$" + group.info.international + "$"
        active = group.active_multipole
        pg = group.info.PG.split(":")[1]
        mpg = id_s.split(":")[1]
        SS = "$" + group.info.schoenflies + "$"
        tp = group.info.type
        c = group.info.crystal
        if c not in ["trigonal", "hexagonal"]:
            continue
        row.append(mpg)
        row.append("(" + pg + ")")
        row.append("")
        row.append("")
        tbl0 = [r"\texttt{" + name + "}", tp, c]
        d = {}
        for i in active:
            X = i[0]
            d[X] = d.get(X, []) + [i[1:]]
        tbl.append(tbl0 + ["Q"] + active_list(d.get("Q", []), c, basis))
        tbl.append([SS, "", "", "G"] + active_list(d.get("G", []), c, basis))
        tbl.append(["", "", "", "T"] + active_list(d.get("T", []), c, basis))
        tbl.append(["", "", "", "M"] + active_list(d.get("M", []), c, basis))
        no += 4
        hl.append(no - 1)

    pdf.table(tbl, row, lbl, r"\#", cpos=cpos, long=True, stretch=1.6, hl=hl[:-1], caption=cap)

    pdf.build()


# ================================================== main
if __name__ == "__main__":
    create_active_multipole()
