"""
MaterialModel manages model information of cluster or crystal system.
"""
import datetime
import sympy as sp
from gcoreutils.nsarray import NSArray
from gcoreutils.pdf_via_latex import PDFviaLaTeX
from gcoreutils.convert_util import sympy_to_latex
from multipie.tag.tag_multipole import TagMultipole
from multipie.tag.tag_list import TagList
from multipie.model.material_model import MaterialModel
from multipie.model.util.atomic_matrix_util import convert_to_matrix
from multipie.multipole.util.atomic_orbital_util import to_latex
from multipie import __version__


# ==================================================
class ModelPDF:
    # ==================================================
    def vspace(self, pdf, sp="5mm"):
        pdf.text(r"\vspace{" + sp + "}" + "\n")

    # ==================================================
    def hr(self, pdf):
        pdf.text()
        pdf.text(r" \hfil \hrule height 1mm width \textwidth \hfil" + "\n")

    # ==================================================
    def ubar_name(self, name):
        if name.count("_") > 0:
            head, no = name.split("_")
            name = head + "$_{" + str(int(no)) + "}$"
        return name

    # ==================================================
    def comment(self, pdf, text):
        pdf.text("%" * 120)
        pdf.text("% " + text)
        pdf.text("%" * 120)

    # ==================================================
    def __init__(self, model_dict, samb_dict, mpm):
        """
        create pdf file from model and basis.

        Args:
            model_dict (_type_): _description_
            samb_dict (_type_): _description_
            mpm (_type_): _description_
        """
        self._mpm = mpm
        filename = mpm.filename(model_dict["info"]["model"] + "_samb", full=False)
        pdf = PDFviaLaTeX(filename, landscape=True, english=True, dir=mpm.dirname)

        info = model_dict["info"]
        name = model_dict["name"]
        data = model_dict["data"]

        cluster = samb_dict["info"]["site_cluster"] | samb_dict["info"]["bond_cluster"]
        samb_dict["info"]["cluster"] = cluster
        cluster = samb_dict["data"]["site_cluster"] | samb_dict["data"]["bond_cluster"]
        samb_dict["data"]["cluster"] = cluster

        pdf.text(r"\setlength{\baselineskip}{16pt}")
        pdf.text(r"\footnotesize")

        # main.
        self.header(pdf, info, name, data)
        self.vspace(pdf, "1cm")
        self.hr(pdf)
        pdf.text(r"\begin{itemize}")
        self.information(pdf, info, name, data)
        self.hr(pdf)
        if not info["molecule"]:
            self.lattice_info(pdf, info, name, data)
            self.hr(pdf)
        self.ket_space(pdf, info, name, data)
        self.basis_info(pdf, info, data, samb_dict)
        self.harmonics_info(pdf, info, samb_dict)
        self.group_info(pdf, info, name, data)
        pdf.text(r"\end{itemize}")

        # postscript.
        mpm.log(f"  * wrote '{filename}.tex.'", None)
        mpm.log(f"  * wrote '{filename}.pdf.'", None)

        del pdf

    # ==================================================
    def header(self, pdf, info, name, data):
        dt = str(datetime.datetime.now())
        dt = dt[: dt.rfind(":")]

        # title.
        name = info["model"]
        pdf.title(r"SAMB for ``\texttt{" + f"{name}" + "}''")
        pdf.text(r"\begin{flushright}")
        pdf.text("Generated on " + dt + f" by MultiPie {__version__}")
        pdf.text(r"\end{flushright}")

    # ==================================================
    def information(self, pdf, info, name, data):
        g_tag = self._mpm.group.tag
        pdf.text(r"\item Group: " + g_tag.info(latex=True) + "\n")
        if not info["molecule"]:
            gp = self._mpm.point_group.tag
            pdf.text(r"\item Associated point group: " + gp.info(latex=True) + "\n")

        self.vspace(pdf)
        self.comment(pdf, "generation condition.")
        pdf.text(r"\item Generation condition" + "\n")
        pdf.text(r"\quad - model type: \texttt{" + info["generate"]["model_type"] + "}\n")
        pdf.text(r"\quad - time-reversal type: \texttt{" + info["generate"]["time_reversal_type"] + "}\n")
        pdf.text(r"\quad - irrep: \texttt{[" + ", ".join(info["generate"]["irrep"]) + "]}\n")
        spinful_str = "spinful" if info["spinful"] else "spinless"
        pdf.text(r"\quad - \texttt{" + spinful_str + "}\n")

    # ==================================================
    def ket_space(self, pdf, info, name, data):
        if not info["site"]:  # no sites.
            return

        self.comment(pdf, "ket hilbert space.")
        pdf.text(r"\item Kets: dimension = " + str(len(info["ket"])))
        orb = to_latex([i.split("@")[0] for i in info["ket"]], info["spinful"], info["crystal"])
        s = [i.split("@")[1] for i in info["ket"]]
        ket = [i + "@" + j for i, j in zip(orb, s)]
        ket = [
            "$"
            + i.split("@")[0].replace("U", r"\uparrow").replace("D", r"\downarrow")
            + "$"
            + "@"
            + self.ubar_name(i.split("@")[1])
            for i in ket
        ]
        ket = sum([[str(no + 1), i] for no, i in enumerate(ket)], [])
        ncol = 5
        if len(ket) < 2 * ncol:
            ncol = len(ket) // 2
        pdf.table(
            ket,
            [""],
            ["No.", "ket"] * ncol,
            cpos="c" + "|cc" * ncol,
            caption="Hilbert space for full matrix.",
            stretch=1.3,
            long=True,
        )

        pdf.text()
        self.comment(pdf, "sites in unit cell.")
        pdf.text(r"\item Sites in (primitive) unit cell:")
        row = []
        tbl = []
        for cluster, lst in data["cluster_site"].items():
            _, wp, _, sym = info["rep_site"][name["alias"][cluster]]
            cluster_name = self.ubar_name(cluster) + f" [{str(wp)}: {sym}]"
            row.append(cluster_name)
            tbl1 = []
            for site in lst:
                pos, mp, _ = data["site"][site]
                head, idx = name["site"][site]
                label = [
                    f"{head}$_{idx}$",
                    "$" + NSArray(pos).semi_evalf().latex() + "$",
                    MaterialModel._mapping_str(mp),
                ]
                tbl1 = tbl1 + label
            tbl.append(tbl1)
        pdf.table(
            tbl,
            row,
            ["site", "position", "mapping"],
            cpos="cc|c|l",
            caption="Site-clusters.",
            stretch=1.3,
            long=True,
            hl=True,
        )

        if not info["bond"]:  # no bonds.
            return

        pdf.text()
        self.comment(pdf, "bonds in unit cell.")
        pdf.text(r"\item Bonds in (primitive) unit cell:")
        row = []
        tbl = []
        for cluster, lst in data["cluster_bond"].items():
            _, wp, _, _, sym = info["rep_bond"][name["alias"][cluster]]
            cluster_name = self.ubar_name(cluster) + f" [{str(wp)}: {sym}]"
            row.append(cluster_name)
            tbl1 = []
            for bond in lst:
                pos, mp, th, _, _ = data["bond"][bond]
                head_t, idx_t = name["site"][f"site_{th[0]+1:03d}"]
                t = head_t + "$_{" + str(idx_t) + "}$"
                head_h, idx_h = name["site"][f"site_{th[1]+1:03d}"]
                h = head_h + "$_{" + str(idx_h) + "}$"
                head, idx = name["bond"][bond]
                _, _, n, idx = head.split(":")
                no = int(bond.split("_")[1])
                label = [
                    r"b$_{" + str(no) + "}$",
                    t,
                    h,
                    str(int(n)),
                    str(int(idx.split("-")[0])),
                    "$" + NSArray(pos).semi_evalf().latex() + "$",
                    MaterialModel._mapping_str(mp),
                ]
                tbl1 = tbl1 + label
            tbl.append(tbl1)
        pdf.table(
            tbl,
            row,
            ["bond", "tail", "head", "$n$", "\#", r"$\bm{b}@\bm{c}$", "mapping"],
            caption="Bond-clusters.",
            cpos="cc|cc|c|c|c|l",
            stretch=1.3,
            long=True,
            hl=True,
        )
        self.hr(pdf)

    # ==================================================
    def basis_info(self, pdf, info, data, samb):
        if not samb["info"]["Z"]:  # no SAMBs.
            return

        # SAMB.
        self.comment(pdf, "SAMB.")
        pdf.text(r"\item SAMB:" + "\n")

        Zinfo = {}
        for k, v in samb["info"]["Z"].items():
            for z_tag in v:
                Zinfo[z_tag] = k
        Zdict = samb["data"]["Z"]

        if not info["molecule"]:
            Zkdict = samb["data"]["Zk"]
            for no, (ztag, zj, zkj) in enumerate(zip(Zdict.keys(), Zdict.values(), Zkdict.values())):
                # Z basis.
                tag = TagMultipole(zj[0]).latex()
                _, m_block, c_block = Zinfo[ztag]
                pdf.text(r"\vspace{4mm}")
                pdf.text(
                    r"\noindent \fbox{No. {"
                    + str(no + 1)
                    + r"}} $\,\,\,"
                    + tag
                    + "$"
                    + f" [{self.ubar_name(m_block)},"
                    + r"\,"
                    + f"{self.ubar_name(c_block)}]"
                )
                eq = sp.S(0)
                rdict = {}
                for c, at, ct in zj[1]:
                    a_no = str(int(at.split("_")[1]))
                    c_no = str(int(ct.split("_")[1]))
                    rdict["a" + a_no] = TagMultipole(samb["data"]["atomic"][at][0]).latex()
                    rdict["c" + c_no] = TagMultipole(samb["data"]["cluster"][ct][0]).latex()
                    eq += (
                        NSArray(c)
                        * sp.symbols(r"A_{" + a_no + "}[a" + a_no + "]")
                        * sp.symbols(r"B_{" + c_no + "}[c" + c_no + "]")
                    )
                eq = sp.latex(eq).replace("A", r"\mathbb{X}").replace("B", r"\otimes\mathbb{Y}")
                for r, s in rdict.items():
                    eq = eq.replace(r, s)
                eq = r"\hat{\mathbb{Z}}_{" + str(no + 1) + "}=" + eq
                pdf.equation(eq, long=True)
                # Z(k) basis.
                eqk = sp.S(0)
                rdict = {}
                for term in zkj[1]:
                    if len(term) == 4:
                        c, at, ut, ft = term
                    else:
                        c, at, ut = term
                        ft = ""
                    a_no = str(int(at.split("_")[1]))
                    u_no = str(int(ut.split("_")[1]))
                    rdict["a" + a_no] = TagMultipole(samb["data"]["atomic"][at][0]).latex()
                    a_symbol = sp.symbols(r"A_{" + a_no + "}[a" + a_no + "]")
                    rdict["u" + u_no] = TagMultipole(samb["data"]["uniform"][ut][0]).latex()
                    u_symbol = sp.symbols(r"B_{" + u_no + "}[u" + u_no + "]")
                    if ft:
                        f_no = str(int(ft.split("_")[1]))
                        f_symbol = sp.symbols(r"C_{" + f_no + "}[f" + f_no + "]")
                        rdict["f" + f_no] = TagMultipole(samb["data"]["structure"][ft][0]).latex()
                    else:
                        f_symbol = sp.S(1)
                    eqk += NSArray(c) * a_symbol * u_symbol * f_symbol
                eqk = (
                    sp.latex(eqk)
                    .replace("A", r"\mathbb{X}")
                    .replace("B", r"\otimes\mathbb{U}")
                    .replace("C", r"\otimes\mathbb{F}")
                )
                for r, s in rdict.items():
                    eqk = eqk.replace(r, s)
                eqk = r"\hat{\mathbb{Z}}_{" + str(no + 1) + r"}(\bm{k})=" + eqk
                pdf.equation(eqk, long=True)
        else:
            for no, (ztag, zj) in enumerate(zip(Zdict.keys(), Zdict.values())):
                # Z basis.
                tag = TagMultipole(zj[0]).latex()
                _, m_block, c_block = Zinfo[ztag]
                pdf.text(r"\vspace{4mm}")
                pdf.text(
                    r"\noindent \fbox{No. {"
                    + str(no + 1)
                    + r"}} $\,\,\,"
                    + tag
                    + "$"
                    + f" [{self.ubar_name(m_block)},"
                    + r"\,"
                    + f"{self.ubar_name(c_block)}]"
                )
                eq = sp.S(0)
                rdict = {}
                for c, at, ct in zj[1]:
                    a_no = str(int(at.split("_")[1]))
                    c_no = str(int(ct.split("_")[1]))
                    rdict["a" + a_no] = TagMultipole(samb["data"]["atomic"][at][0]).latex()
                    rdict["u" + c_no] = TagMultipole(samb["data"]["uniform"][ct][0]).latex()
                    eq += (
                        NSArray(c)
                        * sp.symbols(r"A_{" + a_no + "}[a" + a_no + "]")
                        * sp.symbols(r"B_{" + c_no + "}[u" + c_no + "]")
                    )
                eq = sp.latex(eq).replace("A", r"\mathbb{X}").replace("B", r"\otimes\mathbb{U}")
                for r, s in rdict.items():
                    eq = eq.replace(r, s)
                eq = r"\hat{\mathbb{Z}}_{" + str(no + 1) + "}=" + eq
                pdf.equation(eq, long=True)

        # atomic SAMB group.
        self.comment(pdf, "atomic SAMB group.")
        row = []
        tbl = []
        for m_group, (bra, ket) in data["atomic_braket"].items():
            m_group = self.ubar_name(m_group)
            bra = to_latex(bra, info["spinful"], info["crystal"])
            ket = to_latex(ket, info["spinful"], info["crystal"])
            bra = "$" + ", ".join(bra).replace("U", r"\uparrow").replace("D", r"\downarrow") + "$"
            ket = "$" + ", ".join(ket).replace("U", r"\uparrow").replace("D", r"\downarrow") + "$"
            row.append(m_group)
            tbl.append([bra, ket])
        pdf.table(
            tbl,
            row,
            ["bra", "ket"],
            rc="group",
            caption="Atomic SAMB group.",
            cpos="c|c|c",
            stretch=1.3,
            long=True,
        )

        # atomic SAMB.
        self.comment(pdf, "atomic SAMB.")
        row = []
        tbl = []
        hl_no = -1
        hl = []
        for m_group, m_lst in samb["info"]["atomic"].items():
            hl.append(hl_no)
            m_group = self.ubar_name(m_group)
            for amp in m_lst:
                tag, shape, mat = samb["data"]["atomic"][amp]
                tag = "$" + TagMultipole(tag).latex() + "$"
                no = str(int(amp.split("_")[1]))
                no = "\mathbb{X}_{" + no + "}"
                mat = "$" + convert_to_matrix(shape[0], shape[1], mat).latex() + "$"
                row.append(no)
                tbl.append([tag, m_group, mat])
                hl_no = hl_no + 1
        pdf.table(
            tbl,
            row,
            ["type", "group", "form"],
            rc="symbol",
            caption="Atomic SAMB.",
            cpos="c|c|c|c",
            stretch=1.3,
            long=True,
            rmath=True,
            hl=hl[1:],
        )

        # cluster SAMB.
        if not info["molecule"]:
            self.comment(pdf, "cluster SAMB.")
            row = []
            tbl = []
            hl_no = -1
            hl = []
            for c_group, c_lst in samb["info"]["cluster"].items():
                hl.append(hl_no)
                c_group = self.ubar_name(c_group)
                for cmp in c_lst:
                    tag, vec = samb["data"]["cluster"][cmp]
                    tag = "$" + TagMultipole(tag).latex() + "$"
                    no = str(int(cmp.split("_")[1]))
                    no = "\mathbb{Y}_{" + no + "}"
                    vec = "$" + NSArray(vec).latex() + "$"
                    row.append(no)
                    tbl.append([tag, c_group, vec])
                    hl_no = hl_no + 1
            pdf.table(
                tbl,
                row,
                ["type", "cluster", "form"],
                rc="symbol",
                caption="Cluster SAMB.",
                cpos="c|c|c|c",
                stretch=1.3,
                long=True,
                rmath=True,
                hl=hl[1:],
            )

        # uniform SAMB.
        self.comment(pdf, "uniform SAMB.")
        row = []
        tbl = []
        hl_no = -1
        hl = []
        for c_group, u_lst in samb["info"]["uniform"].items():
            hl.append(hl_no)
            c_group = self.ubar_name(c_group)
            for ump in u_lst:
                tag, shape, mat = samb["data"]["uniform"][ump]
                tag = "$" + TagMultipole(tag).latex() + "$"
                no = str(int(ump.split("_")[1]))
                no = "\mathbb{U}_{" + no + "}"
                mat = "$" + convert_to_matrix(shape[0], shape[1], mat).latex() + "$"
                row.append(no)
                tbl.append([tag, c_group, mat])
                hl_no = hl_no + 1
        pdf.table(
            tbl,
            row,
            ["type", "cluster", "form"],
            rc="symbol",
            caption="Uniform SAMB.",
            cpos="c|c|c|c",
            stretch=1.3,
            long=True,
            rmath=True,
            hl=hl[1:],
        )

        if not info["molecule"]:
            if not samb["info"]["structure"]:  # no structure SAMBs.
                return
            # structure SAMB.
            self.comment(pdf, "structure SAMB.")
            row = []
            tbl = []
            hl_no = -1
            hl = []
            for c_group, s_lst in samb["info"]["structure"].items():
                hl.append(hl_no)
                c_group = self.ubar_name(c_group)
                for smp in s_lst:
                    tag, val = samb["data"]["structure"][smp]
                    tag = "$" + TagMultipole(tag).latex() + "$"
                    no = str(int(smp.split("_")[1]))
                    no = "\mathbb{F}_{" + no + "}"
                    val = "$" + NSArray(val).latex() + "$"
                    row.append(no)
                    tbl.append([tag, c_group, val])
                    hl_no = hl_no + 1
            pdf.table(
                tbl,
                row,
                ["type", "cluster", "form"],
                rc="symbol",
                caption="Structure SAMB.",
                cpos="c|c|c|c",
                stretch=1.3,
                long=True,
                rmath=True,
                hl=hl[1:],
            )
        self.hr(pdf)

    # ==================================================
    def harmonics_info(self, pdf, info, samb):
        if len(samb["info"]["harmonics"]["Q"]) == 0 and len(samb["info"]["harmonics"]["G"]) == 0:  # no harmonics.
            return

        self.comment(pdf, "harmonics.")
        group = self._mpm.point_group

        polar = TagList.from_str(TagMultipole, samb["info"]["harmonics"]["Q"])
        axial = TagList.from_str(TagMultipole, samb["info"]["harmonics"]["G"])

        for hset, name in [(polar, "Polar"), (axial, "Axial")]:
            if not hset:
                continue
            row = []
            tbl = []
            hl_no = -1
            hl = []
            old_rank = -1
            for no, h_tag in enumerate(hset):
                head = h_tag.head
                rank = h_tag.rank
                mul = h_tag.mul
                comp = h_tag.comp
                irrep = h_tag.tag_irrep()
                dim = irrep.dim
                if mul < 1:
                    mul = "-"
                if dim == 1:
                    comp = "-"
                harm = group.harmonics[h_tag]
                harm_ltx = harm.latex().replace("h,", "")
                tbl.append([harm_ltx, str(rank), irrep.latex(), mul, comp, harm.expression(v=NSArray.vector3d(head)).latex()])
                row.append(str(no + 1))
                if rank != old_rank:
                    old_rank = rank
                    hl.append(hl_no)
                hl_no = hl_no + 1
            pdf.table(
                tbl,
                row,
                ["symbol", "rank", "irrep.", "mul.", "comp.", "form"],
                "No.",
                caption=f"{name} harmonics.",
                stretch=1.3,
                rmath=True,
                tmath=True,
                long=True,
                hl=hl[1:],
            )
        self.hr(pdf)

    # ==================================================
    def lattice_info(self, pdf, info, name, data):
        self.comment(pdf, "unit cell.")
        pdf.text(r"\item Unit cell:" + "\n")
        cell = info["cell"]
        a, b, c, alpha, beta, gamma = (
            float(cell["a"]),
            float(cell["b"]),
            float(cell["c"]),
            float(cell["alpha"]),
            float(cell["beta"]),
            float(cell["gamma"]),
        )
        txt = f"\quad $a={a},\,\, b={b},\,\, c={c},\,\, \\alpha={alpha},\,\, \\beta={beta},\,\, \\gamma={gamma}$\n"
        pdf.text(txt)

        pdf.text(r"\item Lattice vectors:" + "\n")
        a1, a2, a3 = NSArray(info["a1"]), NSArray(info["a2"]), NSArray(info["a3"])
        pdf.text(r"\quad $\bm{a}_1=" + a1.evalf().latex() + "$\n")
        pdf.text(r"\quad $\bm{a}_2=" + a2.evalf().latex() + "$\n")
        pdf.text(r"\quad $\bm{a}_3=" + a3.evalf().latex() + "$\n")
        if len(data["plus_set"]) > 1:
            pdf.text(r"\item Plus sets:" + "\n")
            for i in data["plus_set"]:
                pdf.text(r"\quad $+" + NSArray(i).latex() + "$\n")

        self.comment(pdf, "k path.")
        kpath = info["k_path"].replace("Γ", "$\Gamma$")
        tbl = sum([[k.replace("Γ", "$\Gamma$"), "$" + NSArray(v).latex() + "$"] for k, v in info["k_point"].items()], [])
        ncol = min(3, len(info["k_point"]))
        pdf.table(
            tbl,
            [""],
            ["symbol", "position"] * ncol,
            caption="High-symmetry line: " + kpath + ".",
            cpos="c" + "|cc" * ncol,
            stretch=1.3,
            long=True,
        )

    # ==================================================
    def group_info(self, pdf, info, name, data):
        group = self._mpm.group
        so = group.symmetry_operation

        self.comment(pdf, "group information.")
        pdf.text(r"\item Group info.: Generator = " + "$" + ",\,\,".join(so.gen.latex()) + "$" + "\n")

        # conjugacy class.
        irop = [i.latex() for i in so.cc]
        row = []
        tbl = []
        cap = "Conjugacy class"
        if not group.tag.is_point_group():
            cap += " (point-group part)"
        cap += "."
        for ir in irop:
            row.append("$" + ir[0] + "$")
            tbl.append(["$" + "$,\,\, $".join(ir) + "$"])
        pdf.table(
            tbl,
            row,
            ["symmetry operations"],
            rc="rep. SO",
            caption=cap,
            cpos="c|l",
            stretch=1.3,
            long=True,
            hl=True,
        )

        # symmetry operation.
        tbl = []
        for i, op in enumerate(so.full):
            tbl.append([str(i + 1), "$" + op.latex() + "$"])
        tbl = sum(tbl, [])
        pdf.table(
            tbl,
            [""],
            col=["No.", "SO"] * 5,
            caption="Symmetry operations.",
            cpos="c" + "|cc" * 5,
            stretch=1.3,
            long=True,
        )

        # character.
        if not group.tag.is_point_group():
            ch = group.pg.character
        else:
            ch = group.character
        col = ch.symmetry_operation().latex()
        row = ch.irrep_list.latex()
        tbl = [sympy_to_latex(ch.character(i)) for i in ch.irrep_list]
        cap = "Character table"
        if not group.tag.is_point_group():
            cap += " (point-group part)"
        cap += "."
        align = "c|" + "r" * len(col)
        pdf.table(tbl, row, col, caption=cap, rcmath=True, rmath=True, cmath=True, tmath=True, cpos=align, long=True)

        # parity conversion.
        tbl = [i.latex() + "\\,\\,(" + j.latex() + ")" for i, j in ch.parity_conversion().items()]
        pdf.table(tbl, [""], ["$\leftrightarrow$"] * min(len(tbl), 5), caption="Parity conversion.", tmath=True, long=True)

        # symmetric product.
        tbl = []
        for i, ir1 in enumerate(ch.irrep_list):
            t1 = []
            for j, ir2 in enumerate(ch.irrep_list):
                v = sp.latex(ch.symmetric_product_decomposition((ir1, ir2), ret_ex=True)) if j >= i else ""
                t1.append(v)
            tbl.append(t1)
        cp = "c|" + "c" * len(row)
        pdf.table(
            tbl,
            row,
            row,
            rmath=True,
            cmath=True,
            tmath=True,
            cpos=cp,
            caption="Symmetric product, $[\Gamma\otimes\Gamma']_+.$",
            long=True,
        )

        # anti-symmetric product.
        tbl = []
        for ir in ch.irrep_list:
            v = sp.latex(ch.anti_symmetric_product_decomposition(ir, ret_ex=True))
            if v == "0":
                v = "-"
            tbl.append(v)
        pdf.table(
            [tbl],
            [""],
            row,
            rmath=True,
            cmath=True,
            tmath=True,
            caption="Anti-symmetric product, $[\Gamma\otimes\Gamma]_-$.",
            long=True,
        )

        # virtual cluster.
        self.hr(pdf)

        pdf.text("{\n" + r"\scriptsize")

        if info["molecule"]:
            vc = group.virtual_cluster
        else:
            vc = group.pg.virtual_cluster
        sites = sum([[str(no + 1), i.latex()] for no, i in enumerate(vc.site)], [])
        col = ["No.", "position"] * min(4, len(sites) // 2)
        pdf.table(sites, [""], col, caption="Virtual-cluster sites.", tmath=True, stretch=1.7, long=True)

        col = [str(i + 1) for i in range(min(10, len(sites) // 2))]
        row = []
        tbl = []
        for info, basis in vc.items():
            row.append("$" + info.latex().replace("s,", "") + "$")
            basis = NSArray(basis.tolist(), "scalar")
            tbl.append(basis.latex())
        pdf.table(tbl, row, col, "symbol", hl=True, caption="Virtual-cluster basis.", tmath=True, stretch=1.7, long=True)

        pdf.text("}")
