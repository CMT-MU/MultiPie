"""
create matrix of atomic multipoles for lm and jm basis.

When "TeX capacity exceeded, sorry [main memory size=xxxx]" error occurs, increase memory by doing
$ kpsewhich texmf.cnf
$ code /usr/local/texlive/2021/texmf.cnf
- add "main_memory = 10000000" in texmf.cnf
- update .fmt
$ sudo fmtutil-sys --all
"""
from gcoreutils.pdf_via_latex import PDFviaLaTeX
from multipie import __top_dir__
from multipie.multipole.util.spin_orbital_basis import _standard_spinless_basis, _standard_spinful_basis
from multipie.multipole.util.atomic_orbital_util import split_orb_list_rank_block
from multipie.group.point_group import PointGroup


multipole_dir = __top_dir__ + "docs/pdf/atomic_multipole"


# ==================================================
def create_atomic_multipole_lm():
    orb = _standard_spinless_basis["lm"]
    orb = split_orb_list_rank_block(orb, spinful=False, crystal="cubic")
    atomic_braket = {
        "M_001": (orb[0], orb[0]),
        "M_002": (orb[0], orb[1]),
        "M_003": (orb[0], orb[2]),
        "M_004": (orb[0], orb[3]),
        "M_005": (orb[1], orb[1]),
        "M_006": (orb[1], orb[2]),
        "M_007": (orb[1], orb[3]),
        "M_008": (orb[2], orb[2]),
        "M_009": (orb[2], orb[3]),
        "M_010": (orb[3], orb[3]),
    }
    am = PointGroup.spherical_atomic_multipole_basis(atomic_braket, spinful=False)

    cnt = 1
    pdf = PDFviaLaTeX("am_lm", pt=8, dir=multipole_dir, landscape=True)
    pdf.title("Atomic Multipoles (spinless LM basis)")
    for M_i, lst in am.atomic_info.items():
        bra, ket = atomic_braket[M_i]
        tbl = []
        for amp_i in lst:
            tag, m = am[amp_i]
            tbl.append([tag.latex(), m.latex()])
        row = list(range(cnt, cnt + len(lst)))
        cnt += len(lst)
        bra = str(bra).replace("'", "").replace(" ", "")
        ket = str(ket).replace("'", "").replace(" ", "")
        pdf.text("bra: $" + "".join(bra) + "$\n\n\\noindent")
        pdf.text("ket: $" + "".join(ket) + "$")
        pdf.table(tbl, row, ["type", "basis"], "no.", rmath=True, tmath=True, long=True, stretch=1.6, hl=True)


# ==================================================
def create_atomic_multipole_jm():
    orb = _standard_spinful_basis["jm"]
    orb = split_orb_list_rank_block(orb, spinful=True, crystal="cubic")
    atomic_braket = {
        "M_001": (orb[0], orb[0]),
        "M_002": (orb[0], orb[1]),
        "M_003": (orb[0], orb[2]),
        "M_004": (orb[0], orb[3]),
        "M_005": (orb[1], orb[1]),
        "M_006": (orb[1], orb[2]),
        "M_007": (orb[1], orb[3]),
        "M_008": (orb[2], orb[2]),
        "M_009": (orb[2], orb[3]),
        "M_010": (orb[3], orb[3]),
    }
    am = PointGroup.spherical_atomic_multipole_basis(atomic_braket, spinful=True)

    cnt = 1
    pdf = PDFviaLaTeX("am_jm", pt=8, dir=multipole_dir, landscape=True)
    pdf.title("Atomic Multipoles (spinful JM basis)")
    for M_i, lst in am.atomic_info.items():
        bra, ket = atomic_braket[M_i]
        tbl = []
        for amp_i in lst:
            tag, m = am[amp_i]
            tbl.append([tag.latex(), m.latex()])
        row = list(range(cnt, cnt + len(lst)))
        cnt += len(lst)
        bra = str(bra).replace("'", "").replace(" ", "")
        ket = str(ket).replace("'", "").replace(" ", "")
        pdf.text("bra: $" + "".join(bra) + "$\n\n\\noindent")
        pdf.text("ket: $" + "".join(ket) + "$")
        pdf.table(tbl, row, ["type", "basis"], "no.", rmath=True, tmath=True, long=True, stretch=1.6, hl=True)


# ================================================== main
create_atomic_multipole_lm()
create_atomic_multipole_jm()
