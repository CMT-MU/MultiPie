"""
Create PDF concerning all of crystallographic point and space groups.

To do
- bond wyckoff decomposed by vector harmonics.
"""

import logging

from create_PDF.create_pdf_dir import create_pdf_dir
from create_PDF.create_pdf_atomic_multipole_group import create_atomic_multipole_group
from create_PDF.create_pdf_atomic_multipole import create_atomic_multipole
from create_PDF.create_pdf_character import create_character
from create_PDF.create_pdf_compatibility_relation import create_compatibility_relation
from create_PDF.create_pdf_group_info import create_group_info
from create_PDF.create_pdf_harmonics_multipole import create_harmonics_multipole
from create_PDF.create_pdf_harmonics import create_harmonics
from create_PDF.create_pdf_representation_matrix import create_representation_matrix
from create_PDF.create_pdf_symmetry_operation import create_symmetry_operation
from create_PDF.create_pdf_wyckoff import create_wyckoff
from create_PDF.create_pdf_active_multipole import create_active_multipole
from create_PDF.create_pdf_markdown import create_markdown


# ==================================================
def create_pdf():
    logging.basicConfig(format="%(message)s", level=logging.INFO)

    create_pdf_dir()
    create_group_info()
    create_compatibility_relation()
    create_atomic_multipole()
    create_character()
    create_symmetry_operation()
    create_harmonics()
    create_atomic_multipole_group()
    create_wyckoff()
    create_representation_matrix()
    create_harmonics_multipole()
    create_active_multipole()

    create_markdown()


# ================================================== main
if __name__ == "__main__":
    create_pdf()
