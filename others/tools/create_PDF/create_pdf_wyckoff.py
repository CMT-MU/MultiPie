"""
Create PDF for wyckoff.
"""

import os

from multipie.core.multipie_info import __top_dir__
from multipie.core.group import Group
from multipie.util.util import timer, to_latex
from multipie.util.util_binary import BinaryManager
from multipie.util.util_pdf_latex import PDFviaLaTeX

pdfdir = os.path.join(__top_dir__, "others/pdf/")


# ==================================================
def chunk(lst, n=24):
    if len(lst) > n:
        return [lst[i : i + n] for i in range(0, len(lst), n)]
    else:
        return [lst]


# ==================================================
def create_wyckoff_site_each(tag, h_dir):
    group = Group(tag)
    file = f"wyckoff_site"
    title = group.latex(True)
    pset = group.symmetry_operation.get("plus_set", None)

    pdf = PDFviaLaTeX(file, dir=h_dir, pt=8, style="narrowest", english=True)
    pdf.title(title)

    wp = group.wyckoff["site"]
    if not group.is_point_group:
        ps = r"* plus set:\,\,\," + r",\quad ".join(["$+" + to_latex(i, "vector") + "$" for i in pset])
        pdf.text(ps)

    for w, val in wp.items():
        row = []
        tbl = []
        sym = val["symmetry"]
        pos = val["conventional"]
        mp = val["mapping"]
        cap = r"Wyckoff site: {\tt " + str(w) + "}" + r", site symmetry: {\tt " + sym + "}"
        for no, (i, m) in enumerate(zip(pos, mp)):
            mlst = chunk(m)
            if len(mlst) == 1:
                row.append(no + 1)
                ms = str(mlst[0]).replace(" ", "")
                tbl.append(["$" + to_latex(i, "vector") + "$", r"{\tt " + ms + "}"])
            else:
                row.append(no + 1)
                ms = str(mlst[0])[:-1].replace(" ", "") + ","
                tbl.append(["$" + to_latex(i, "vector") + "$", r"{\tt " + ms + "}"])
                for mi in mlst[1:-1]:
                    row.append("")
                    ms = str(mi)[1:-1].replace(" ", "") + ","
                    tbl.append(["", r"{\tt " + ms + "}"])
                row.append("")
                ms = str(mlst[-1])[1:].replace(" ", "")
                tbl.append(["", r"{\tt " + ms + "}"])
        pdf.table(tbl, row, ["position", "mapping"], "No.", hl=True, stretch=1.2, long=True, caption=cap)

    pdf.build()


# ==================================================
def create_wyckoff_bond_each(tag, h_dir):
    group = Group(tag)
    file = f"wyckoff_bond"
    title = group.latex(True)
    pset = group.symmetry_operation.get("plus_set", None)
    n_so = len(group.symmetry_operation["tag"])

    pdf = PDFviaLaTeX(file, dir=h_dir, pt=8, style="narrowest", english=True)
    pdf.title(title)

    wp = group.wyckoff
    if not group.is_point_group:
        ps = r"* plus set:\,\,\," + r",\quad ".join(["$+" + to_latex(i, "vector") + "$" for i in pset])
        pdf.text("\n" + ps)
    pdf.text("\n")

    for s_wp, val in wp["site"].items():
        if s_wp == "1o":
            continue
        sym = val["symmetry"]
        txt = r"* Wyckoff site: {\tt " + str(s_wp) + "}" + r", site symmetry: {\tt " + sym + "}"
        pdf.text(txt)
        for b_wp in val["bond"]:
            bond = wp["bond"][b_wp]["conventional"]
            mp = wp["bond"][b_wp]["mapping"]
            vector, center = bond[:, 0:3], bond[:, 3:6]
            row = []
            tbl = []
            cap = r"Wyckoff bond: {\tt " + str(b_wp) + "}"
            for no, (v, c, m) in enumerate(zip(vector, center, mp)):
                mlst = chunk(m)
                if len(mlst) == 1:
                    row.append(no + 1)
                    ms = str(mlst[0]).replace(" ", "")
                    tbl.append(["$" + to_latex(v, "vector") + "$", "$" + to_latex(c, "vector") + "$", r"{\tt " + ms + "}"])
                else:
                    row.append(no + 1)
                    ms = str(mlst[0])[:-1].replace(" ", "") + ","
                    tbl.append(["$" + to_latex(v, "vector") + "$", "$" + to_latex(c, "vector") + "$", r"{\tt " + ms + "}"])
                    for mi in mlst[1:-1]:
                        row.append("")
                        ms = str(mi)[1:-1].replace(" ", "") + ","
                        tbl.append(["", "", r"{\tt " + ms + "}"])
                    row.append("")
                    ms = str(mlst[-1])[1:].replace(" ", "")
                    tbl.append(["", "", r"{\tt " + ms + "}"])
                    row.append(no + 1)
            pdf.table(tbl, row, ["vector", "center", "mapping"], "No.", hl=True, stretch=1.2, long=True, caption=cap)

    pdf.build()


# ==================================================
def create_wyckoff_m_site_each(tag, h_dir):
    group = Group(tag)
    file = f"wyckoff_site"
    title = group.latex(True)

    pdf = PDFviaLaTeX(file, dir=h_dir, pt=8, style="narrowest", english=True)
    pdf.title(title)

    wp = group.wyckoff["site"]
    for w, val in wp.items():
        row = []
        tbl = []
        sym = val["symmetry"]
        pos = val["conventional"]
        mp = val["mapping"]
        cap = r"Wyckoff site: {\tt " + str(w) + "}" + r", site symmetry: {\tt " + sym + "}"
        for no, (i, m) in enumerate(zip(pos, mp)):
            mlst = chunk(m)
            if len(mlst) == 1:
                row.append(no + 1)
                ms = str(mlst[0]).replace(" ", "")
                tbl.append(["$" + to_latex(i, "vector") + "$", r"{\tt " + ms + "}"])
            else:
                row.append(no + 1)
                ms = str(mlst[0])[:-1].replace(" ", "") + ","
                tbl.append(["$" + to_latex(i, "vector") + "$", r"{\tt " + ms + "}"])
                for mi in mlst[1:-1]:
                    row.append("")
                    ms = str(mi)[1:-1].replace(" ", "") + ","
                    tbl.append(["", r"{\tt " + ms + "}"])
                row.append("")
                ms = str(mlst[-1])[1:].replace(" ", "")
                tbl.append(["", r"{\tt " + ms + "}"])
        pdf.table(tbl, row, ["position", "mapping"], "No.", hl=True, stretch=1.2, long=True, caption=cap)

    pdf.build()


# ==================================================
@timer
def create_wyckoff():
    info = BinaryManager("info")

    for id_s in info["id_set"]["PG"]["all"]:
        no = int(id_s.split(":")[1])
        tag = info["tag"][id_s]
        h_dir = os.path.join(__top_dir__, f"{pdfdir}PG/{no:03d}-{tag}")
        create_wyckoff_site_each(tag, h_dir)
        create_wyckoff_bond_each(tag, h_dir)

    for id_s in info["id_set"]["SG"]["all"]:
        no = int(id_s.split(":")[1])
        tag = info["tag"][id_s]
        h_dir = os.path.join(__top_dir__, f"{pdfdir}SG/{no:03d}-{tag.replace("^", "_")}")
        create_wyckoff_site_each(tag, h_dir)
        create_wyckoff_bond_each(tag, h_dir)

    for id_s in info["id_set"]["MPG"]["all"]:
        no = id_s.split(":")[1]
        tag = info["tag"][id_s]
        h_dir = os.path.join(__top_dir__, f"{pdfdir}MPG/{no}")
        create_wyckoff_m_site_each(tag, h_dir)

    for id_s in info["id_set"]["MSG"]["all"]:
        no = id_s.split(":")[1]
        tag = info["tag"][id_s]
        h_dir = os.path.join(__top_dir__, f"{pdfdir}MSG/{no}")
        create_wyckoff_m_site_each(tag, h_dir)


# ================================================== main
if __name__ == "__main__":
    create_wyckoff()
