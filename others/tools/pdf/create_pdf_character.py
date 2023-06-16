"""
create character table.
"""
import sympy as sp
from gcoreutils.pdf_via_latex import PDFviaLaTeX
from gcoreutils.convert_util import sympy_to_latex
from multipie.character.character_pg_set import CharacterPGSet
from multipie.tag.tag_group import TagGroup
from multipie import get_binary

point_group_dir = __file__[: __file__.rfind("/")] + "/../../docs/pdf/point_group"


# ==================================================
def create_compatible_relation():
    pdf = PDFviaLaTeX("comp_relation", landscape=True, dir=point_group_dir)
    pdf.title("Compatible Relation")

    comp_rel = CharacterPGSet.compatibility_relation()
    group = [("Oh", TagGroup.create("cubic_series")), ("D6h", TagGroup.create("hexagonal_series"))]

    for g, pg in group:
        g_name = "$" + TagGroup(g).latex() + "$ subgroup"
        row = []
        tbl = []
        for p in comp_rel.keys():
            if p not in pg:
                continue
            name = p.latex()
            setting = p.info()[3]
            if setting and "axis" not in setting:
                name += r"\,\, ({\rm " + setting + "})"
            irrep = comp_rel[p].latex()
            row.append(name)
            tbl.append(irrep)

        pdf.table(tbl[1:], row[1:], tbl[0], row[0], caption=g_name, rcmath=True, rmath=True, cmath=True, tmath=True)


# ==================================================
def create_point_group(ctset):
    pdf = PDFviaLaTeX("point_group_ct", landscape=True, dir=point_group_dir)
    pdf.title("32 Point Groups ($\\omega=e^{2\\pi i/3}$)")

    for tag, p in ctset.items():
        name = tag.info(latex=True)
        rc = p.latex()
        col = p.symmetry_operation().latex()
        row = p.irrep_list.latex()
        tbl = [sympy_to_latex(p.character(i)) for i in p.irrep_list]
        align = "c" + "r" * len(col)

        pdf.text(name)
        pdf.text('tag = "{\\tt ' + str(tag) + '}"\n')

        pdf.text("* character table")
        pdf.text("\\begin{quote}")
        pdf.table(
            tbl,
            row,
            col,
            rc,
            rcmath=True,
            rmath=True,
            cmath=True,
            tmath=True,
            cpos=align,
        )
        pdf.text("\\end{quote}")

        pdf.text("* polar $\\leftrightarrow$ axial conversion")
        pdf.text("\\begin{quote}")
        tbl = [[i.latex() + "\\,\\,(" + j.latex() + ")" for i, j in p.parity_conversion().items()]]
        pdf.simple_table(tbl, tmath=True)
        pdf.text("\\end{quote}")

        tbl = []
        for i, ir1 in enumerate(p.irrep_list):
            t1 = []
            for j, ir2 in enumerate(p.irrep_list):
                v = sp.latex(p.symmetric_product_decomposition((ir1, ir2), ret_ex=True)) if j >= i else ""
                t1.append(v)
            tbl.append(t1)

        cp = "c|" + "c" * len(row)
        pdf.text("* symmetric product")
        pdf.text("\\begin{quote}")
        pdf.table(tbl, row, row, rmath=True, cmath=True, tmath=True, cpos=cp)
        pdf.text("\\end{quote}")

        tbl = []
        for ir in p.irrep_list:
            v = sp.latex(p.anti_symmetric_product_decomposition(ir, ret_ex=True))
            if v == "0":
                v = "-"
            tbl.append(v)
        pdf.text("* anti-symmetric product")
        pdf.text("\\begin{quote}")
        pdf.table([tbl], [""], row, rmath=True, cmath=True, tmath=True)
        pdf.text("\\end{quote}")

        pdf.text("\\newpage")


# ================================================== main
core = get_binary()
ctset = core[CharacterPGSet]

create_compatible_relation()
create_point_group(ctset)
