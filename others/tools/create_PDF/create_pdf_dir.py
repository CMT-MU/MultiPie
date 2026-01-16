"""
Create PDF folders.
"""

import os

from multipie.core.multipie_info import __top_dir__
from multipie.util.util import timer
from multipie.util.util_binary import BinaryManager

pdfdir = os.path.join(__top_dir__, "others/pdf/")


# ==================================================
@timer
def create_pdf_dir():
    info = BinaryManager("info")

    os.makedirs(os.path.join(pdfdir, "misc"), exist_ok=True)
    os.makedirs(os.path.join(pdfdir, "info"), exist_ok=True)

    os.makedirs(os.path.join(pdfdir, "PG"), exist_ok=True)
    for id_s in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        tag = info["tag"][id_s]
        no = int(id_s.split(":")[1])
        os.makedirs(os.path.join(pdfdir, f"PG/{no:03d}-{tag}"), exist_ok=True)

    os.makedirs(os.path.join(pdfdir, "SG"), exist_ok=True)
    for id_s in info["id_set"]["SG"]["all"]:
        tag = info["tag"][id_s]
        no = int(id_s.split(":")[1])
        os.makedirs(os.path.join(pdfdir, f"SG/{no:03d}-{tag.replace("^", "_")}"), exist_ok=True)

    os.makedirs(os.path.join(pdfdir, "MPG"), exist_ok=True)
    for id_s in info["id_set"]["MPG"]["all"]:
        no = id_s.split(":")[1]
        os.makedirs(os.path.join(pdfdir, f"MPG/{no}"), exist_ok=True)

    os.makedirs(os.path.join(pdfdir, "MSG"), exist_ok=True)
    for id_s in info["id_set"]["MSG"]["all"]:
        no = id_s.split(":")[1]
        os.makedirs(os.path.join(pdfdir, f"MSG/{no}"), exist_ok=True)


# ================================================== main
if __name__ == "__main__":
    create_pdf_dir()
