"""
Create PDF for material model.
"""

import numpy as np
import sympy as sp

from multipie.util.util import to_latex
from multipie.core.group import Group


# ==================================================
class ModelPDF:
    # ==================================================
    def fmt_number(self, nums, int_digits=2, frac_digits=5):
        if frac_digits != 0:
            width = int_digits + 1 + frac_digits
            if isinstance(nums, (list, np.ndarray, tuple)):
                return r"\texttt{[" + ",".join([f"{x:{width}.{frac_digits}f}" for x in nums]) + "]}"
            else:
                return r"\texttt{" + f"{nums:{width}.{frac_digits}f}" + "}"
        else:
            if isinstance(nums, (list, np.ndarray, tuple)):
                return r"\texttt{[" + ",".join([f"{x}" for x in nums]) + "]}"
            else:
                return r"\texttt{" + f"{nums}" + "}"

    # ==================================================
    def fmt_site_bond(self, site_bond):
        if ";" not in site_bond:
            s = f"'{site_bond}'"
        else:
            tail, head = site_bond.split(";")
            head, n, m = head.split("_")
            s = r"\texttt{" f"'{head}'" + r"}-\texttt{" + f"'{tail}'({int(n)}th,{m})" + "}"

        return s

    # ==================================================
    def text(self, txt=""):
        self.pdf.text(txt)

    # ==================================================
    def vspace(self, sp="5mm"):
        self.text(r"\vspace{" + sp + "}" + "\n")

    # ==================================================
    def hr(self, text="", width="0.5mm", space="3mm", pos="5mm"):
        s = r"\vspace{5mm}" + "\n"
        if text != "":
            s += r"\noindent \rule[" + width + "]{" + pos + "}{" + width + "}"
            s += r"\hspace{" + space + r"}\textbf{" + text + r"}\hspace{" + space + "}"
            s += r"\leaders\hrule height " + width + r"\hfill\kern0pt"
        else:
            s += r"\noindent\leaders\hrule height " + width + r"\hfill\kern0pt"
        s += "\n" + r"\vspace{3mm}"

        self.text(s + "\n")

    # ==================================================
    def number_tag(self, no):
        d = r"\fbox{" + no + r"}\quad "
        return d

    # ==================================================
    def header(self):
        dt = self.mm["created"]
        name = self.mm["model"]
        version = self.mm["version"]

        self.text(r"\setlength{\baselineskip}{16pt}")
        self.text(r"\setlength{\parindent}{0pt}")
        self.text(r"\setlength{\mathindent}{0pt}")
        self.text(r"\footnotesize")

        # title.
        self.pdf.title(r"Model for ``\texttt{" + f"{name}" + "}''")
        self.text(r"\begin{flushright}")
        self.text("Generated on " + dt + f" by MultiPie {version}")
        self.text(r"\end{flushright}")

        self.text(r"\setlength{\baselineskip}{16pt}")
        self.text(r"\footnotesize")

    # ==================================================
    def general_info(self):
        tp = self.mm["toroidal_priority"]
        mn = self.mm["max_neighbor"]
        scr = self.mm["search_cell_range"]
        scr = str(scr[0:2]) + ", " + str(scr[2:4]) + ", " + str(scr[4:6])
        basis_type = self.mm["basis_type"]
        select = self.mm["SAMB_select"]
        a_select = self.mm["atomic_select"]
        s_select = self.mm["site_select"]
        b_select = self.mm["bond_select"]
        irrep = ["$" + self.mm.group.tag_irrep(i, latex=True) + "$" for i in select["Gamma"]]
        a_irrep = ["$" + self.mm.group.tag_irrep(i, latex=True) + "$" for i in a_select["Gamma"]]
        s_irrep = ["$" + self.mm.group.tag_irrep(i, latex=True) + "$" for i in s_select["Gamma"]]
        b_irrep = ["$" + self.mm.group.tag_irrep(i, latex=True) + "$" for i in b_select["Gamma"]]

        self.hr("General Condition")
        self.text(r"\begin{itemize}")
        self.text(r"\item Basis type: \texttt{" + basis_type + "}\n")

        self.text(r"\item SAMB selection:" + "\n")
        self.text(r"\, - Type: \texttt{" + str(select["X"]).replace("'", "") + "}\n")
        self.text(r"\, - Rank: \texttt{" + str(select["l"]) + "}\n")
        self.text(r"\, - Irrep.: [" + ", ".join(irrep) + "]\n")
        self.text(r"\, - Spin (s): \texttt{" + str(select["s"]) + "}\n")

        self.text(r"\item Atomic selection:" + "\n")
        self.text(r"\, - Type: \texttt{" + str(a_select["X"]).replace("'", "") + "}\n")
        self.text(r"\, - Rank: \texttt{" + str(a_select["l"]) + "}\n")
        self.text(r"\, - Irrep.: [" + ", ".join(a_irrep) + "]\n")
        self.text(r"\, - Spin (s): \texttt{" + str(a_select["s"]) + "}\n")

        self.text(r"\item Site-cluster selection:" + "\n")
        self.text(r"\, - Rank: \texttt{" + str(s_select["l"]) + "}\n")
        self.text(r"\, - Irrep.: [" + ", ".join(s_irrep) + "]\n")

        self.text(r"\item Bond-cluster selection:" + "\n")
        self.text(r"\, - Type: \texttt{" + str(b_select["X"]).replace("'", "") + "}\n")
        self.text(r"\, - Rank: \texttt{" + str(b_select["l"]) + "}\n")
        self.text(r"\, - Irrep.: [" + ", ".join(b_irrep) + "]\n")

        self.text(r"\item Max. neighbor: \texttt{" + str(mn) + "}\n")
        self.text(r"\item Search cell range: \texttt{" + scr + "}\n")
        self.text(r"\item Toroidal priority: \texttt{" + str(tp).lower() + "}\n")
        self.text(r"\end{itemize}")

    # ==================================================
    def group_cell_info(self):
        gp = self.mm.group.latex(detail=True)
        cell = self.mm["cell"]

        a = self.fmt_number(cell["a"])
        b = self.fmt_number(cell["b"])
        c = self.fmt_number(cell["c"])
        alpha = self.fmt_number(cell["alpha"], frac_digits=1)
        beta = self.fmt_number(cell["beta"], frac_digits=1)
        gamma = self.fmt_number(cell["gamma"], frac_digits=1)

        a1, a2, a3 = self.mm["unit_vector"]
        a1 = self.fmt_number(a1)
        a2 = self.fmt_number(a2)
        a3 = self.fmt_number(a3)

        self.hr("Group and Unit Cell")
        self.text(r"\begin{itemize}")
        self.text(r"\item Group: " + gp)
        if not self.mm.group.is_point_group:
            ch_gp = self.mm.group.latex(detail=True, tp="PG")
            self.text(r"\item Associated point group: " + ch_gp)
        self.text(r"\end{itemize}")

        self.text(r"\begin{itemize}")

        self.text(r"\item Unit cell:" + "\n")
        txt = f"\\quad $a={a},\\,\\, b={b},\\,\\, c={c},\\,\\, \\alpha={alpha},\\,\\, \\beta={beta},\\,\\, \\gamma={gamma}$\n"
        self.text(txt)

        self.text(r"\item Lattice vectors (conventional cell):" + "\n")
        self.text(r"\quad $\bm{a}_1=" + a1 + "$\n")
        self.text(r"\quad $\bm{a}_2=" + a2 + "$\n")
        self.text(r"\quad $\bm{a}_3=" + a3 + "$\n")

        if not self.mm.group.is_point_group:
            ps = self.mm.group.symmetry_operation.get("plus_set", None)
            if len(ps) > 1:
                self.text(r"\item Plus sets:" + "\n")
                ps = "$+" + r"$,\,\,\, $+".join(to_latex(ps, "vector")) + "$"
                self.text(r"\quad " + ps + "\n")

        self.text(r"\end{itemize}")

    # ==================================================
    def cluster_samb_info(self):
        c_samb = self.mm["cluster_samb"]

        self.hr("Cluster SAMB")

        self.text(r"\noindent\textbullet\, Site cluster" + "\n")
        no = 1
        for wp, samb in c_samb.items():
            if "@" in wp:
                continue
            self.vspace("5mm")
            self.text(r"\noindent ** Wyckoff: \texttt{" + wp + "}\n")
            for idx, (mat, ex) in samb.items():
                tag = self.mm.group.tag_multipole(idx, latex=True, superscript="s")
                for e, m in zip(tag, mat):
                    d = self.number_tag("y" + str(no)) + sp.latex(sp.Symbol(e)) + "=" + to_latex(m, "vector")
                    self.pdf.equation(d, long=True)
                    no += 1

        self.vspace("5mm")
        self.text(r"\noindent\textbullet\, Bond cluster" + "\n")
        for wp, samb in c_samb.items():
            if "@" not in wp:
                continue
            self.vspace("5mm")
            self.text(r"\noindent ** Wyckoff: \texttt{" + wp + "}\n")
            for idx, (mat, ex) in samb.items():
                tag = self.mm.group.tag_multipole(idx, latex=True, superscript="s")
                for e, m in zip(tag, mat):
                    d = self.number_tag("y" + str(no)) + sp.latex(sp.Symbol(e)) + "=" + to_latex(m, "vector")
                    self.pdf.equation(d, long=True)
                    no += 1

    # ==================================================
    def atomic_samb_info(self):
        atomic_samb = self.mm["atomic_samb"]
        basis_type = self.mm["basis_type"]

        self.hr("Atomic SAMB")

        no = 1
        for name, samb in atomic_samb.items():
            if name.bh_rank > name.kt_rank:
                continue
            self.vspace("5mm")
            bs1 = [
                self.mm.group.tag_atomic_basis(
                    self.mm.group.atomic_basis(basis_type)[name.bh_rank][i], name.bh_rank, ket=False, latex=True
                )
                for i in name.bh_idx
            ]
            bs2 = [
                self.mm.group.tag_atomic_basis(
                    self.mm.group.atomic_basis(basis_type)[name.kt_rank][i], name.kt_rank, ket=True, latex=True
                )
                for i in name.kt_idx
            ]
            bra = "$" + r",\, ".join(bs1) + "$"
            ket = "$" + r",\, ".join(bs2) + "$"
            self.text(r"\noindent\textbullet\, bra: " + bra + "\n")
            self.text(r"\noindent\textbullet\, ket: " + ket + "\n")
            for idx, (mat, ex) in samb.items():
                tag = self.mm.group.tag_multipole(idx, latex=True, superscript="a")
                for e, m in zip(tag, mat):
                    d = self.number_tag("x" + str(no))
                    d += sp.latex(sp.Symbol(e)) + "=" + to_latex(m, "matrix")
                    self.pdf.equation(d, long=True)
                    no += 1

    # ==================================================
    def combined_samb_info(self):
        combined_samb = self.mm["combined_samb"]
        samb_num = self.mm["SAMB_number"]
        samb_num_min = self.mm["SAMB_number_min"]
        basis_type = self.mm["basis_type"]
        common_id = self.mm["common_id"]

        self.hr("SAMB")

        self.text(r"\noindent" + f"{samb_num_min} (all {samb_num}) SAMBs\n")

        for comb, samb in combined_samb.items():
            if not samb:
                continue

            no_tag, sb_tag = common_id[comb]

            bk = comb.bk_info
            wp = comb.wyckoff
            clustar_str = "b" if wp.count("@") > 0 else "s"

            bs1 = [
                self.mm.group.tag_atomic_basis(
                    self.mm.group.atomic_basis(basis_type)[bk.bh_rank][i], bk.bh_rank, ket=False, latex=True
                )
                for i in bk.bh_idx
            ]
            bs2 = [
                self.mm.group.tag_atomic_basis(
                    self.mm.group.atomic_basis(basis_type)[bk.kt_rank][i], bk.kt_rank, ket=True, latex=True
                )
                for i in bk.kt_idx
            ]

            bra = "$" + r",\, ".join(bs1) + "$"
            ket = "$" + r",\, ".join(bs2) + "$"
            cluster = (
                r"'\texttt{" + comb.head + "}' site-cluster"
                if clustar_str == "s"
                else r"'\texttt{" + comb.head + r"}'-'\texttt{" + comb.tail + "}' bond-cluster"
            )
            cluster += r" : \texttt{" + sb_tag[0] + "}"
            wyckoff = r"\texttt{" + wp + "}"

            self.vspace("5mm")
            self.text(r"\noindent\textbullet\, " + cluster + "\n")
            self.text(r"\noindent\quad * bra: " + bra + "\n")
            self.text(r"\noindent\quad * ket: " + ket + "\n")
            self.text(r"\noindent\quad * wyckoff: " + wyckoff + "\n")

            i = 0
            for idx, (cl, ex) in samb.items():
                tag = self.mm.group.tag_multipole(idx, latex=True, superscript="c")
                for t, m in zip(tag, cl):
                    zi = no_tag[0][i]
                    i += 1
                    ex = 0
                    for cg, t1, c1, t2, c2 in m:
                        t1 = self.mm.group.tag_multipole(t1, c1, latex=True, superscript="a")
                        t2 = self.mm.group.tag_multipole(t2, c2, latex=True, superscript=clustar_str)
                        ex += cg * sp.Symbol(t1, commutative=False) * sp.Symbol(t2, commutative=False)

                    d = self.number_tag(zi)
                    d += sp.latex(sp.Symbol(t)) + "=" + sp.latex(ex)
                    self.pdf.equation(d, long=True)

            if len(no_tag) > 1:
                if self.mm["pdf_ctrl"]["common_samb"]:
                    sb_head = "(" + ", ".join(sb_tag) + ")"
                    sb_str = np.array(no_tag).T.tolist()
                    sb_str = ", ".join(["(" + ", ".join(i) + ")" for i in sb_str])
                    self.text(r"\noindent * common SAMBs" + "\n")
                    self.text(r"\noindent\texttt{" + sb_head + ", " + sb_str + "}")
                else:
                    n = (len(no_tag) - 1) * len(no_tag[0])
                    self.text(r"\noindent\quad" + f"+ {n} common SAMBs" + "\n")

        self.vspace()

    # ==================================================
    def rep_site_bond_info(self):
        site = self.mm["site"]["representative"]
        bond = self.mm["bond"]["info"]

        self.hr("Site and Bond")

        cap = "Orbital of each site"
        lbl = ["site", "orbital"]
        row = []
        tbl = []
        hl = []
        no = 0
        for name, rep in site.items():
            orb = sum([[self.mm.group.tag_atomic_basis(i, l, latex=True) for i in ol] for l, ol in enumerate(rep.orbital)], [])
            if len(orb) < 17:
                orb = "$" + "$, $".join(orb) + "$"
                name = r"\texttt{" + name + "}"
                row.append(rep.no)
                tbl.append([name, orb])
                hl.append(no)
                no += 1
            else:
                orb1 = "$" + "$, $".join(orb[:16]) + "$"
                name = r"\texttt{" + name + "}"
                row.append(rep.no)
                tbl.append([name, orb1])
                orb2 = "$" + "$, $".join(orb[16:]) + "$"
                row.append("")
                tbl.append(["", orb2])
                hl.append(no + 1)
                no += 2
        self.pdf.table(tbl, row, lbl, r"\#", cpos="ccl", long=True, stretch=1.6, hl=hl[:-1], caption=cap)

        rank_dict = {0: "s", 1: "p", 2: "d", 3: "f"}
        cap = "Neighbor and bra-ket of each bond"
        lbl = ["head", "tail", "neighbor", "head (bra)", "tail (ket)"]
        row = []
        tbl = []
        for no, rep in enumerate(bond):
            tail = r"\texttt{" + rep.tail + "}"
            head = r"\texttt{" + rep.head + "}"
            nb = rep.neighbor
            if len(nb) > 9:
                nb = (r"\texttt{" + str(nb[:4])[:-1] + r",$\cdots$," + str(nb[-1]) + "]}").replace(" ", "")
            else:
                nb = (r"\texttt{" + str(nb) + "}").replace(" ", "")
            tr = r"\texttt{[" + ",".join([rank_dict[i] for i in rep.t_rank]) + "]}"
            hr = r"\texttt{[" + ",".join([rank_dict[i] for i in rep.h_rank]) + "]}"
            row.append(str(no + 1))
            tbl.append([head, tail, nb, hr, tr])
        self.pdf.table(tbl, row, lbl, r"\#", cpos="ccclll", long=True, stretch=1.6, hl=True, caption=cap)

    # ==================================================
    def site_info(self):
        site = self.mm["site"]["cell"]

        self.hr("Site in Unit Cell")

        self.text(r"\noindent Sites in (conventional) cell (no plus set), SL = sublattice" + "\n")

        lbl = ["position" + r" ($\bm{s})$", "mapping"]
        for name, cell_site_name in site.items():
            rep = self.mm["site"]["representative"][name]
            cap = (
                r"'\texttt{"
                + f"{name}"
                + r"}' (\#"
                + f"{rep.no}) site cluster "
                + r"(\texttt{"
                + f"{rep.wyckoff}"
                + r"}), \texttt{"
                + f"{rep.symmetry}"
                + "}"
            )
            row = []
            tbl = []
            for c in cell_site_name:
                if c.plus_set != 1:
                    continue
                pos = self.fmt_number(c.position)
                m = c.mapping
                if len(m) == 48:  # in case of Oh (too long).
                    m = (r"\texttt{" + str(m[:4])[:-1] + r",$\cdots$," + str(m[-1]) + "]}").replace(" ", "")
                else:
                    m = self.fmt_number(m, frac_digits=0)
                row.append(str(c.no))
                tbl.append([pos, m])

            self.pdf.table(tbl, row, lbl, "SL", cpos="cll", long=True, stretch=1.6, hl=True, caption=cap)

    # ==================================================
    def bond_info(self):
        bond = self.mm["bond"]["cell"]
        neighbor = self.mm["pdf_ctrl"]["max_neighbor"]

        self.hr("Bond in Unit Cell")

        self.text(
            r"\noindent Bonds in (conventional) cell (no plus set): tail, head = (SL, plus set), (N)D = (non)directional"
            + f" (listed up to {neighbor}th neighbor at most)"
            + "\n"
        )

        lbl = ["vector" + r" ($\bm{v})$", "center" + r" ($\bm{c})$", "mapping", "head", "tail", r"$\bm{R}$ (primitive)"]
        for name, cell_bond_name in bond.items():
            rep = self.mm["bond"]["representative"][name]
            if neighbor is not None and rep.neighbor > neighbor:
                continue

            d_str = "D" if rep.directional else "ND"
            ss, cn, mul = name.split("_")
            tail_atom, head_atom = ss.split(";")
            dist = round(rep.distance, 5)
            cap = (
                f"{rep.neighbor}-th "
                + r"'\texttt{"
                + f"{head_atom}"
                + r"}'-'\texttt{"
                + f"{tail_atom}"
                + "}'"
                + f" [{mul}] ("
                + r"\#"
                + f"{rep.no}"
                + ") bond cluster ("
                + r"\texttt{"
                + f"{rep.wyckoff}"
                + "}"
                + f"), {d_str}, "
                + r"$|\bm{v}|"
                + f"={dist}$ (cartesian)"
            )
            row = []
            tbl = []
            for i in cell_bond_name:
                if i.plus_set != 1:
                    continue
                v = self.fmt_number(i.vector)
                c = self.fmt_number(i.center)
                m = i.mapping
                if len(m) == 48:  # in case of Oh (too long).
                    m = (r"\texttt{" + str(m[:4])[:-1] + r",$\cdots$," + str(m[-1]) + "]}").replace(" ", "")
                else:
                    m = self.fmt_number(m, frac_digits=0)
                t, h = i.t_idx, i.h_idx
                t = (r"\texttt{" + str(t) + "}").replace(" ", "")
                h = (r"\texttt{" + str(h) + "}").replace(" ", "")
                R = self.fmt_number(i.R_primitive, frac_digits=0)
                row.append(str(i.no))
                tbl.append([v, c, m, h, t, R])

            self.pdf.table(tbl, row, lbl, "SL", cpos="cllllll", long=True, stretch=1.6, hl=True, caption=cap)

    # ==================================================
    def harmonics_info(self):
        if not self.mm["pdf_ctrl"]["harmonics"]:
            return

        harmonics = self.mm.get_multipole_expression()

        self.hr("Harmonics")

        cap = "Harmonics"
        lbl = ["symbol", "irrep.", "rank", "X", "multiplicity", "component", "symmetry"]
        no = 0
        row = []
        tbl = []
        hl = []
        for (X, irrep, rank, n), ex_set in harmonics.items():
            idx = (X, rank, irrep, n, -1, 0, 0, "q")
            irrep = "$" + self.mm.group.tag_irrep(irrep, latex=True) + "$"
            rank = str(rank)
            XX = "$Q$, $T$" if X == "Q" else "$G$, $M$"
            d = len(ex_set)
            n = "-" if n == -1 else str(n)
            if d == 1:
                sym = "$" + self.mm.group.tag_multipole(idx, comp=-1, latex=True) + "$"
                row.append(str(no + 1))
                tbl.append([sym, irrep, rank, XX, n, "-", "$" + sp.latex(ex_set[0]) + "$"])
                no += 1
            else:
                for i, ex in enumerate(ex_set):
                    sym = "$" + self.mm.group.tag_multipole(idx, comp=i, latex=True) + "$"
                    row.append(str(no + 1))
                    if i == 0:
                        tbl.append([sym, irrep, rank, XX, n, str(i + 1), "$" + sp.latex(ex) + "$"])
                    else:
                        tbl.append([sym, "", "", "", "", str(i + 1), "$" + sp.latex(ex) + "$"])
                    no += 1
            hl.append(no - 1)

        self.pdf.table(tbl, row, lbl, r"\#", cpos="clcccccl", long=True, stretch=1.6, hl=hl[:-1], caption=cap)

    # ==================================================
    def full_matrix_info(self):
        # full matrix.
        basis = self.mm.ket()

        self.hr("Basis in full matrix")

        cap = f"dimension = {len(basis)}"

        tbl = []
        for no, b in enumerate(basis):
            b = "$" + str(b) + "$"
            tbl.append([str(no), b])
        tbl = sum(tbl, [])
        self.pdf.table(
            tbl,
            [""],
            col=[r"\#", "orbital@atom(SL)"] * min(len(tbl) // 2, 5),
            caption=cap,
            cpos="c" + "|cl" * min(len(tbl) // 2, 5),
            stretch=1.6,
            long=True,
        )

        # atomic basis.
        if self.mm["basis_type"] == "jml":
            return

        site = self.mm["site"]["representative"]
        basis = self.mm.group.atomic_basis("lg")
        basis_form = Group.global_info()["harmonics"]["basis_function"]

        ranks = sorted(
            list(set(sum([[rank for rank, orb in enumerate(rep.orbital) if len(orb) > 0] for rep in site.values()], [])))
        )

        cap = "Atomic basis (orbital part only)"
        lbl = ["definition"]
        row = []
        tbl = []
        hl = []
        no = -1
        for rank in ranks:
            for i in basis[rank]:
                form = sp.latex(basis_form[i][0])
                sb = self.mm.group.tag_atomic_basis(i, rank, latex=True)
                row.append(sb)
                tbl.append([form])
                no += 1
            hl.append(no)

        self.pdf.table(
            tbl, row, lbl, "orbital", caption=cap, cpos="ll", hl=hl[:-1], rmath=True, tmath=True, stretch=1.6, long=True
        )

    # ==================================================
    def symmetry_operation_info(self):
        so = self.mm.group.symmetry_operation["tag"]

        self.hr("Symmetry Operation")

        cap = "Symmetry operation"

        # symmetry operation.
        tbl = []
        for no, op in enumerate(so):
            op = "$" + self.mm.group.tag_symmetry_operation(op, latex=True) + "$"
            tbl.append([str(no + 1), op])
        tbl = sum(tbl, [])
        self.pdf.table(
            tbl,
            [""],
            col=[r"\#", "SO"] * min(len(tbl), 5),
            caption=cap,
            cpos="c" + "|cl" * min(len(tbl), 5),
            stretch=1.6,
            long=True,
        )

    # ==================================================
    def __init__(self, mm, pdf):
        """
        create pdf file from material model.

        Args:
            mm (MaterialModel): material model class.
        """
        self.mm = mm
        self.pdf = pdf

        self.header()
        self.general_info()
        self.group_cell_info()
        self.symmetry_operation_info()
        self.harmonics_info()
        self.full_matrix_info()
        self.combined_samb_info()
        self.atomic_samb_info()
        self.cluster_samb_info()
        self.rep_site_bond_info()
        self.site_info()
        self.bond_info()

        self.pdf.build()
