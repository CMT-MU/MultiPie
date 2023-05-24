"""
create response tensor expressions.
"""
from gcoreutils.pdf_via_latex import PDFviaLaTeX
from gcoreutils.convert_util import sympy_to_latex
from multipie.tag.tag_group import TagGroup
from multipie.response_tensor.response_tensor_pg_set import ResponseTensorPGSet
from multipie import get_binary

rt_dir = __file__[: __file__.rfind("/")] + "/../../docs/pdf/response_tensor"


# ==================================================
def create_pdf_response_tensor():
    core = get_binary()
    rtset = core[ResponseTensorPGSet]

    for pg_tag in TagGroup.create():
        rt = rtset[pg_tag]
        rt_Q = rt.select(head="Q")
        rt_G = rt.select(head="G")

        pdf = PDFviaLaTeX("tensor_" + str(pg_tag), style="narrow", dir=rt_dir)
        pdf.title("Response Tensors up to 4th rank in $" + pg_tag.latex() + "$")

        pdf.text(r"\begin{center}\LARGE --- polar tensors ---\end{center}")
        for tag in rt_Q.keys():
            M = rt_Q[tag]
            pdf.equation(tag.latex() + " = " + M.latex())
            def_mul = rt_Q.definition[tag]
            eqs = [sympy_to_latex(t) + " = " + sympy_to_latex(d) for t, d in def_mul.items()]
            if eqs:
                pdf.text(r"\begin{quote}")
                pdf.equation(eqs)
                pdf.text(r"\end{quote}")

        pdf.text(r"\newpage")
        pdf.text(r"\begin{center}\LARGE --- axial tensors ---\end{center}")
        for tag in rt_G.keys():
            M = rt_G[tag]
            pdf.equation(tag.latex() + " = " + M.latex())
            def_mul = rt_G.definition[tag]
            eqs = [sympy_to_latex(t) + " = " + sympy_to_latex(d) for t, d in def_mul.items()]
            if eqs:
                pdf.text(r"\begin{quote}")
                pdf.equation(eqs)
                pdf.text(r"\end{quote}")


# ================================================== main
create_pdf_response_tensor()
