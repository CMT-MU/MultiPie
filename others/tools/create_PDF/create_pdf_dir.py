"""
Create PDF folders.
"""

import os
import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
from multipie.util.util import timer
from others.tools.data.data_group_name_list import group_name_list

pdfdir = os.path.join(__top_dir__, "others/pdf/")


# ==================================================
@timer
def create_pdf_dir():
    os.makedirs(os.path.join(pdfdir, "misc"), exist_ok=True)
    os.makedirs(os.path.join(pdfdir, "info"), exist_ok=True)

    os.makedirs(os.path.join(pdfdir, "PG"), exist_ok=True)
    for id_s, (d, md, no, tag, name, pg, sg, mpg, msg) in group_name_list["PG"].items():
        os.makedirs(os.path.join(pdfdir, f"PG/{no:03d}-{tag}"), exist_ok=True)

    os.makedirs(os.path.join(pdfdir, "SG"), exist_ok=True)
    for id_s, (d, md, no, tag, name, pg, sg, mpg, msg) in group_name_list["SG"].items():
        os.makedirs(os.path.join(pdfdir, f"SG/{no:03d}-{tag.replace("^", "_")}"), exist_ok=True)

    os.makedirs(os.path.join(pdfdir, "MPG"), exist_ok=True)
    for id_s, (d, md, no, tag, name, pg, sg, mpg, msg) in group_name_list["MPG"].items():
        os.makedirs(os.path.join(pdfdir, f"MPG/{no}"), exist_ok=True)

    os.makedirs(os.path.join(pdfdir, "MSG"), exist_ok=True)
    for id_s, (d, md, no, tag, name, pg, sg, mpg, msg) in group_name_list["MSG"].items():
        os.makedirs(os.path.join(pdfdir, f"MSG/{no}"), exist_ok=True)


# ================================================== main
if __name__ == "__main__":
    create_pdf_dir()
