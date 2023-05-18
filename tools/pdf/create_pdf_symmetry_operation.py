"""
create matrix of symmetry operations.
"""
from gcoreutils.nsarray import NSArray
from gcoreutils.pdf_via_latex import PDFviaLaTeX
from gcoreutils.string_util import circle_number
from multipie.symmetry_operation.symmetry_operation_g_set import SymmetryOperationGSet
from multipie.tag.tag_group import TagGroup
from multipie import __top_dir__, get_binary

point_group_dir = __top_dir__ + "docs/pdf/point_group"
space_group_dir = __top_dir__ + "docs/pdf/space_group"


# ==================================================
def create_symmetry_operation(soset, t_pg):
    outdir = point_group_dir if t_pg else space_group_dir
    outstr = "point" if t_pg else "space"
    titlestr = "32 Point" if t_pg else "230 Space"
    pdf = PDFviaLaTeX(outstr + "_group_so", landscape=True, dir=outdir)
    pdf.title(titlestr + " Groups")

    if t_pg:
        pgset = TagGroup.create()
    else:
        pgset = TagGroup.create(space_group=True)

    for tag in pgset:
        p = soset[tag]
        name = "\n" + tag.info(latex=True)
        generator = "$" + ",\,\,".join(p.gen.latex()) + "$"

        irop = [i.latex() for i in p.cc]
        cc = "* conjugacy class"
        if not t_pg:
            cc += " (point-group part)"
        cc += "\n"
        cc += "\\begin{quote}\n"
        for ir in irop:
            cc += "[ $" + ir[0] + "$ ] = \\quad "
            cc += "$" + "$,\,\, $".join(ir) + "$\\newline"
        cc += "\n\\end{quote}\n"

        n = 10 if t_pg else 5

        ops = []
        for i, op in enumerate(p.full):
            ops.append(circle_number(i + 1) + "\quad " + op.latex())

        pdf.text(name)
        if not t_pg:
            pdf.text('tag = "{\\tt ' + str(tag) + ", " + str(tag.pg) + '}"\n')
        else:
            pdf.text('tag = "{\\tt ' + str(tag) + '}"\n')

        pdf.text("* generator : " + generator + "\n")
        pdf.text(cc)

        ss = "* symmetry operation"
        if not t_pg:
            ss += r"\quad" + ",\quad ".join(["$+" + i + "$" for i in p.plus_set.latex()])
        pdf.text(ss)
        pdf.text("\\begin{quote}")
        pdf.simple_table(ops, cols=n, tmath=True, cpos="l" * n)
        pdf.text("\\end{quote}\n")

        if t_pg:
            if str(tag) not in ["D6h", "Th", "O", "Td", "Oh"]:
                row = p.key_list().latex()
                col = row
                tbl = [[p.product(i, j).latex() for j in p.keys()] for i in p.keys()]
            pdf.text("* product table")
            if str(tag) not in ["D6h", "Th", "O", "Td", "Oh"]:
                pdf.text("\\begin{quote}")
                pdf.table(tbl, row, col, rmath=True, cmath=True, tmath=True)
                pdf.text("\\end{quote}")
            else:
                pdf.text("\\quad\\quad omitted because of large table.")

        pdf.text("\n\\newpage")


# ==================================================
def create_symmetry_operation_detail(soset, t_pg):
    outdir = point_group_dir if t_pg else space_group_dir
    outstr = "point" if t_pg else "space"
    titlestr = "32 Point" if t_pg else "230 Space"
    pdf = PDFviaLaTeX(outstr + "_group_so_detail", landscape=True, dir=outdir)
    pdf.title(titlestr + " Groups (detail)")

    rv = NSArray.vector3d(head="Q")
    gv = NSArray.vector3d(head="G")

    if t_pg:
        pgset = TagGroup.create()
    else:
        pgset = TagGroup.create(space_group=True)

    for tag in pgset:
        p = soset[tag]
        name = tag.info(latex=True)
        if not t_pg:
            name += r"\quad" + ",\quad ".join(["$+" + i + "$" for i in p.plus_set.latex()])

        so_list = []
        for i, op in enumerate(p.full):
            so_list.append(circle_number(i + 1) + "\quad " + op.latex())

        pov = p.mat(axial=False)
        axv = p.mat(axial=True)

        epp = pov.apply(rv).latex()
        epa = axv.apply(gv).latex()

        sos = [[pov[i][0:3, :].latex(), axv[i][0:3, :].latex(), epp[i], epa[i]] for i in range(len(pov))]

        pdf.table(
            sos,
            so_list,
            ["polar vector", "axial vector", "EP (polar)", "EP (axial)"],
            "sym. op.",
            rmath=True,
            tmath=True,
            stretch=1.3,
            caption=name,
            cpos="lcccc",
            long=True,
        )
        pdf.text("\\newpage")


# ================================================== main
core = get_binary()
soset = core[SymmetryOperationGSet]

create_symmetry_operation(soset, t_pg=True)
create_symmetry_operation(soset, t_pg=False)

create_symmetry_operation_detail(soset, t_pg=True)
create_symmetry_operation_detail(soset, t_pg=False)
