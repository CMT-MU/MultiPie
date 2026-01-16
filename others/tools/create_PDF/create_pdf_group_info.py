"""
Create PDF for group info.
"""

import os

from multipie.core.multipie_info import __top_dir__
from multipie.core.group import Group
from multipie.util.util import timer
from multipie.util.util_binary import BinaryManager
from multipie.util.util_pdf_latex import PDFviaLaTeX

h_dir = os.path.join(__top_dir__, "others/pdf/info")


# ==================================================
@timer
def create_group_info_pg(info):
    pdf = PDFviaLaTeX("PG", dir=h_dir, pt=10, style="narrowest", english=True)

    rc = "ID"
    col = ["tag", "", r"Sch\"onflies", "international", "setting", "crystal"]
    hl = [1, 4, 7, 14, 19, 26, 31, 32, 35, 36, 39, 41, 44]

    row = []
    tbl = []
    for i in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        g = Group(i)
        g_info = g.info
        no = r"\texttt{" + g.ID + "}"
        tag = r"\texttt{" + g_info.tag + "}"
        tag1 = r"\texttt{PG:" + g.ID + "}"
        SN = "$" + g_info.schoenflies + "$"
        IN = "$" + g_info.international + "$"
        crystal = g_info.crystal
        setting = g_info.setting

        row.append(no)
        tbl.append([tag, tag1, SN, IN, setting, crystal])

    pdf.text(r"\begin{center}")
    pdf.table(tbl, row, col, rc, caption="Point group information.", stretch=1.0, long=True, hl=hl, cpos="lllllll")
    pdf.text(r"\end{center}")

    pdf.build()


# ==================================================
@timer
def create_group_info_sg(info):
    pdf = PDFviaLaTeX("SG", dir=h_dir, pt=10, style="narrowest", english=True)

    rc = "ID"
    col = ["tag", "", "", r"Sch\"onflies", "international", "setting", "crystal", "PG"]
    hl = [
        1,
        4,
        8,
        14,
        23,
        45,
        73,
        79,
        81,
        87,
        97,
        109,
        121,
        141,
        145,
        147,
        154,
        160,
        166,
        172,
        173,
        175,
        181,
        185,
        189,
        193,
        198,
        205,
        213,
        219,
    ]

    row = []
    tbl = []
    for i in info["id_set"]["SG"]["all"]:
        g = Group(i)
        g_info = g.info
        no = r"\texttt{" + g.ID + "}"
        tag = r"\texttt{" + g_info.tag + "}"
        tag1 = r"\texttt{SG:" + g.ID + "}"
        tag2 = r"\texttt{" + g.ID + "}"
        SN = "$" + g_info.schoenflies + "$"
        IN = "$" + g_info.international + "$"
        crystal = g_info.crystal
        setting = g_info.setting
        pg = "$" + Group(g_info.PG).info.schoenflies + "$"

        row.append(no)
        tbl.append([tag, tag1, tag2, SN, IN, setting, crystal, pg])

    pdf.text(r"\begin{center}")
    pdf.table(tbl, row, col, rc, caption="Space group information.", stretch=1.0, long=True, hl=hl, cpos="lllllllll")
    pdf.text(r"\end{center}")

    pdf.build()


# ==================================================
@timer
def create_group_info_mpg(info):
    pdf = PDFviaLaTeX("MPG", dir=h_dir, pt=10, style="narrowest", english=True)

    rc = "ID"
    col = ["tag", "international", "setting", "crystal", "PG", "type"]

    row = []
    tbl = []
    hl = []
    pg_no = 1
    n = -1
    for i in info["id_set"]["MPG"]["all"]:
        g = Group(i)
        g_info = g.info
        no = r"\texttt{" + g.ID + "}"
        tag = r"\texttt{" + "MPG:" + g_info.tag + "}"
        IN = "$" + g_info.international + "$"
        crystal = g_info.crystal
        setting = "$" + g_info.setting + "$"
        tp = g_info.type
        pg = "$" + Group(g_info.PG).info.schoenflies + "$"
        if int(g.ID.split(".")[0]) > pg_no:
            pg_no += 1
            hl.append(n)
        n += 1

        row.append(no)
        tbl.append([tag, IN, setting, crystal, pg, tp])

    pdf.text(r"\begin{center}")
    pdf.table(tbl, row, col, rc, caption="Magnetic point group information.", stretch=1.0, long=True, hl=hl, cpos="lllllll")
    pdf.text(r"\end{center}")

    pdf.build()


# ==================================================
@timer
def create_group_info_msg(info):
    pdf = PDFviaLaTeX("MSG", dir=h_dir, pt=10, style="narrowest", english=True)

    rc = "ID"
    col = ["tag", "BNS", "OG", "crystal", "PG", "SG", "MPG", "type"]

    row = []
    tbl = []
    hl = []
    sg_no = 1
    n = -1
    for i in info["id_set"]["MSG"]["all"]:
        g = Group(i)
        g_info = g.info
        no = r"\texttt{" + g.ID + "}"
        tag = r"\texttt{" + "MSG:" + g_info.tag + "}"
        BNS = "$" + g_info.BNS + "$"
        OG = "$" + g_info.OG + "$"
        crystal = g_info.crystal
        tp = g_info.type
        pg = "$" + Group(g_info.PG).info.schoenflies + "$"
        sg = "$" + Group(g_info.SG).info.schoenflies + "$"
        mpg = r"\texttt{" + g_info.MPG.split(":")[1] + "} $(" + Group(g_info.MPG).info.international + ")$"
        if int(g.ID.split(".")[0]) > sg_no:
            sg_no += 1
            hl.append(n)
        n += 1

        row.append(no)
        tbl.append([tag, BNS, OG, crystal, pg, sg, mpg, tp])

    pdf.text(r"\begin{center}")
    pdf.table(tbl, row, col, rc, caption="Magnetic space group information.", stretch=1.0, long=True, hl=hl, cpos="lllllllll")
    pdf.text(r"\end{center}")

    pdf.build()


# ==================================================
@timer
def create_group_info():
    info = BinaryManager("info")

    create_group_info_pg(info)
    create_group_info_sg(info)
    create_group_info_mpg(info)
    create_group_info_msg(info)


# ================================================== main
if __name__ == "__main__":
    create_group_info()
