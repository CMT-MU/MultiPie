"""
Material model analyzer class.

This module provides material model analyzer.
"""

import os
import numpy as np
import sympy as sp
import multiprocessing
from itertools import product
from collections import defaultdict

from multipie import __version__, SAMBType, UniqueSAMBType
from multipie.core.group import Group
from multipie.util.util_binary import BinaryManager
from multipie.util.util import deep_update, time_stamp, check_latex, check_qtdraw, read_dict, write_dict
from multipie.util.util_pdf_latex import PDFviaLaTeX
from multipie.util.util_crystal import get_cell_info, create_igrid, convert_to_primitive
from multipie.util.util_material_model import (
    get_basis_type,
    get_bond,
    get_tail_head,
    get_neighbor_info,
    create_site_grid,
    create_site_so,
    create_wyckoff_dict,
    create_braket_dict,
    create_full_matrix_info,
    parse_representative_site,
    parse_representative_bond,
    write_site_dict,
    write_bond_dict,
    write_site_grid,
    create_qtdraw,
    parse_samb_select,
    create_atomic_samb_qtdraw,
    create_cluster_samb_qtdraw,
)
from multipie.util.util_material_model_pdf import ModelPDF
from multipie.core.default_model import _default_model
from multipie.util.util_plot import plot_site, plot_bond, plot_site_samb, plot_bond_samb, plot_harmonics


# ==================================================
class MaterialModel(BinaryManager):
    # ==================================================
    def __init__(self, topdir=None, verbose=False):
        """
        Material model analyzer.

        Args:
            topdir (str, optional): top directory. [default: cwd]
            verbose (bool, optional): verbose comment ?
        """
        if topdir is None:
            topdir = os.getcwd()
        super().__init__(topdir=topdir, verbose=verbose)
        self._num_proc = multiprocessing.cpu_count()
        self._jl_verbose = 10 if verbose else 0

    # ==================================================
    @property
    def group(self):
        """
        Group class.

        Returns:
            - (Group) -- group class.
        """
        return self._group

    # ==================================================
    def load(self, name):
        """
        Load model.

        Args:
            name (str): model name.
        """
        self.set_subdir(name)
        self.load_binary(name + ".pkl")
        self._group = Group(self["group"])

    # ==================================================
    def save(self):
        """
        Save model data in subdir.
        """
        self.set_subdir(self["model"])
        self.save_binary(self["model"], detail=False)

    # ==================================================
    def save_view(self):
        """
        Save model view QtDraw file.
        """
        if check_qtdraw() and self["qtdraw_prop"]["create"]:
            from qtdraw import create_qtdraw_file

            cwd = os.getcwd()
            path = self.get_cwd()
            os.chdir(path)
            filename = self["model"]

            create_qtdraw_file(
                filename=f"{filename}.qtdw",
                callback=lambda qtdraw: create_qtdraw(
                    qtdraw, self.group, filename, self["cell_info"], self["site"], self["bond"], self["qtdraw_prop"]
                ),
            )

            os.chdir(cwd)
            if self.verbose:
                print(f"save qtdraw to '{path}/{filename}.qtdw'.")

    # ==================================================
    def save_pdf(self):
        """
        Save model as PDF file.
        """
        if check_latex() and self["pdf_ctrl"]["create"]:
            cwd = os.getcwd()
            path = self.get_cwd()
            os.chdir(path)
            filename = self["model"]

            pdf = PDFviaLaTeX(filename, landscape=True, english=True, dir=self.get_cwd())
            ModelPDF(self, pdf)

            os.chdir(cwd)
            if self.verbose:
                print(f"save PDF to '{path}/{filename}.pdf'.")

    # ==================================================
    def save_atomic_samb(self):
        """
        Save atomic SAMB QtDraw file.
        """
        if check_qtdraw() and self["qtdraw_prop"]["create"]:
            from qtdraw import create_qtdraw_file

            cwd = os.getcwd()
            path = self.get_cwd()
            os.chdir(path)
            os.makedirs("samb", exist_ok=True)
            os.chdir("samb")

            name = self["model"] + "_atomic_samb"
            create_qtdraw_file(filename=f"{name}.qtdw", callback=lambda qtdraw: create_atomic_samb_qtdraw(qtdraw, self, name))

            os.chdir(cwd)

    # ==================================================
    def save_cluster_samb(self, site_bond):
        """
        Save cluster SAMB QtDraw file.

        Args:
            site_bond (str or list): (list of) site or bond.
        """
        if check_qtdraw() and self["qtdraw_prop"]["create"]:
            from qtdraw import create_qtdraw_file

            cwd = os.getcwd()
            path = self.get_cwd()
            os.chdir(path)
            os.makedirs("samb", exist_ok=True)
            os.chdir("samb")

            if type(site_bond) == str:
                site_bond = [site_bond]
            name = self["model"] + "_" + "-".join(site_bond)

            create_qtdraw_file(
                filename=f"{name}.qtdw",
                callback=lambda qtdraw: create_cluster_samb_qtdraw(qtdraw, self, site_bond, name),
            )

            os.chdir(cwd)

    # ==================================================
    def save_site_bond(self, site_bond):
        """
        Save site/bond cluster QtDraw file.

        Args:
            site_bond (str or list): (list of) site or bond.
        """

        def draw(qtdraw):
            qtdraw.clear_data()
            qtdraw.set_model(name)
            qtdraw.set_crystal(self["crystal"])
            if self.group.is_point_group:
                qtdraw.set_cell("off")
            else:
                qtdraw.set_cell("single")

            for sb in site_bond:
                self.plot_site_bond(qtdraw, sb)

        if check_qtdraw() and self["qtdraw_prop"]["create"]:
            from qtdraw import create_qtdraw_file

            cwd = os.getcwd()
            path = self.get_cwd()
            os.chdir(path)
            os.makedirs("samb", exist_ok=True)
            os.chdir("samb")

            if type(site_bond) == str:
                site_bond = [site_bond]

            name = self["model"] + "_" + "-".join(site_bond) + "_def"
            create_qtdraw_file(filename=f"{name}.qtdw", callback=draw)
            os.chdir(cwd)

    # ==================================================
    def save_samb_qtdraw(self):
        """
        Save SAMB QtDraw file.
        """
        site_bond = [
            s for s in self["wyckoff"].keys() if s.count(";") == 0 or int(s.split("_")[1]) < self["qtdraw_prop"]["max_neighbor"]
        ]
        site = [s for s in site_bond if s.count(";") == 0]
        bond = [s for s in site_bond if s.count(";") > 0]

        self.save_atomic_samb()
        self.save_site_bond(site)
        self.save_cluster_samb(site)
        for b in bond:
            self.save_site_bond(b)
            self.save_cluster_samb(b)
        if self.verbose:
            print(f"save SAMB QtDraw files in '{self.get_cwd()}/samb'.")

    # ==================================================
    def save_samb_matrix(self, select, parameter=None):
        """
        Save SAMB matrix (and hr).

        Args:
            select (dict): select dict (see, select_combined_samb).
            parameter (dict, optional): parameter dict, dict[z#, value].
        """
        select_regularized = self.select_combined_samb(**select)[-1]
        combined_samb_matrix = self.get_combined_samb_matrix(**select)

        if parameter is None:
            H = None
        else:
            H = self.get_hr(parameter, combined_samb_matrix)

        # fileter selected id.
        cluster_rep = {}
        for k, lst in self["cluster_representative"].items():
            r_lst = []
            for ls, info in lst:
                f_lst = [i for i in ls if i in combined_samb_matrix.keys()]
                if f_lst:
                    r_lst.append((f_lst, info))
            if r_lst:
                cluster_rep[k] = r_lst

        # write matrix.
        cwd = os.getcwd()
        path = self.get_cwd()
        os.chdir(path)
        filename = self["model"] + "_matrix.py"
        var = self["model"]
        comment = "Selected SAMB matrix.\n"
        ket = [orbital + "@" + atom + f"({sl})" for atom, sl, rank, orbital in self["full_matrix"]["ket"]]
        site_dict = {
            k + "_" + str(vi.sublattice): vi.position_primitive.tolist()
            for k, v in self["site"]["cell"].items()
            for vi in v
            if vi.plus_set == 1
        }
        site = [site_dict[atom + "_" + str(sl)] for atom, sl, rank, orbital in self["full_matrix"]["ket"]]
        ket_site = dict(zip(ket, site))
        d = {
            "model": self["model"],
            "pkl": f"{self["model"]}.pkl ({self["created"]})",
            "select": select_regularized,
            "dimension": len(ket_site),
            "ket_site": ket_site,
            "matrix": {z: {k: str(v).replace(" ", "") for k, v in elm.items()} for z, elm in combined_samb_matrix.items()},
            "cluster_representative": cluster_rep,
        }
        write_dict(d, filename, var, comment)
        if self.verbose:
            print(f"save matrix to '{path}/{filename}'.")

        # write hr.dat.
        if H is not None:
            filename = self["model"] + "_hr.dat"
            with open(filename, mode="w", encoding="utf-8") as f:
                print(f"# SAMB matrix from {var}.pkl ({self["created"]})", file=f)
                print("# select", file=f)
                for k, v in select_regularized.items():
                    print(f"#   {k}: {str(v).replace(" ", "")}", file=f)
                print(f"# basis ({len(ket_site)})", file=f)
                for no, (b, p) in enumerate(ket_site.items()):
                    print(f"#   {no:2d} {b}: {str(p).replace(" ", "")}", file=f)
                for z, v in parameter.items():
                    print(f"# {z:<4} = {v}", file=f)
                print("#", file=f)
                print("# n1   n2   n3    m    n    re                        im", file=f)
                for (n1, n2, n3, m, n), v in H.items():
                    v = complex(v)
                    r, i = v.real, v.imag
                    s = f"{n1: 4d} {n2: 4d} {n3: 4d} {m: 4d} {n: 4d}    {r: .15e}    {i: .15e}"
                    print(s, file=f)
            if self.verbose:
                print(f"save hr to '{path}/{filename}'.")
        os.chdir(cwd)

    # ==================================================
    def analyze(self, model_in):
        """
        Set model.

        Args:
            model_in (dict or str): model input dict or file name.
        """
        if type(model_in) == str:
            if self.verbose:
                print(f"read '{model_in}'.")
            model_in = read_dict(model_in, self.topdir)

        self.clear()

        # set model based on default model.
        model = {}
        deep_update(model, _default_model)
        deep_update(model, model_in)

        # get basic data.
        cell = model["cell"]
        site_data = model["site"]
        bond_data = model["bond"]
        spinful = model["spinful"]
        max_neighbor = model["max_neighbor"]
        search_cell_range = model["search_cell_range"]
        toroidal_priority = model["toroidal_priority"]
        samb_select = model["SAMB_select"]
        atomic_select = model["atomic_select"]
        site_select = model["site_select"]
        bond_select = model["bond_select"]

        # set group data.
        group = Group(model["group"])
        crystal = group.info.crystal
        irrep = list(group.character["table"].keys())
        samb_select = parse_samb_select(samb_select, irrep)
        atomic_select = parse_samb_select(atomic_select, irrep)
        site_select = parse_samb_select(site_select, irrep)
        bond_select = parse_samb_select(bond_select, irrep)

        basis_type = get_basis_type(site_data, spinful)
        basis_list = ["jml", "lgs", "lg"]
        basis_info = {k: group.atomic_basis(k) for k in basis_list}
        basis_info_type = group.atomic_basis(basis_type)

        # set cell data.
        cell_info = get_cell_info(crystal, cell)
        G = cell_info["G"][0:3, 0:3]

        # set grid data.
        if group.is_point_group:
            igrid = None
        else:
            igrid = create_igrid(search_cell_range)

        # set site and bond.
        site_dict = parse_representative_site(group, site_data, basis_type, basis_info)
        site_grid = create_site_grid(site_dict, igrid)
        site_so = create_site_so(group, site_dict)
        bond_dict = parse_representative_bond(group, G, site_grid, site_so, site_dict, bond_data, max_neighbor, self.verbose)
        A = cell_info["A"][0:3, 0:3].T
        lattice = group.info.lattice
        Ap = convert_to_primitive(lattice, A, shift=False)

        # site_bond => wyckoff, braket.
        wyckoff_dict = create_wyckoff_dict(site_dict["representative"], bond_dict["representative"])
        braket_dict = create_braket_dict(site_dict["representative"], bond_dict["representative"], basis_info_type)

        # information for full matrix.
        full_mat_info = create_full_matrix_info(site_dict)

        # save information.
        name = model["model"]
        self._group = group
        self["model"] = name
        self["group"] = str(group)
        self["crystal"] = crystal
        self["unit_vector"] = A.tolist()
        self["unit_vector_primitive"] = Ap.tolist()
        self["basis_type"] = basis_type
        self["cell"] = cell_info["cell"]
        self["cell_info"] = cell_info
        self["toroidal_priority"] = toroidal_priority
        self["SAMB_select"] = samb_select
        self["atomic_select"] = atomic_select
        self["site_select"] = site_select
        self["bond_select"] = bond_select
        self["max_neighbor"] = model["max_neighbor"]
        self["search_cell_range"] = model["search_cell_range"]
        self["site"] = site_dict
        self["site_grid"] = site_grid
        self["bond"] = bond_dict
        self["wyckoff"] = wyckoff_dict
        self["braket"] = braket_dict
        self["full_matrix"] = full_mat_info
        self["qtdraw_prop"] = model["qtdraw"]
        self["pdf_ctrl"] = model["pdf"]

        # SAMB.
        atomic_samb, atomic_id = self.get_atomic_samb(atomic_select)
        cluster_samb, cluster_id = self.get_cluster_samb(site_select, bond_select)
        combined_samb, combined_id, combined_min_num, combined_num, common_id, cluster_rep = self.get_combined_samb(
            atomic_samb, cluster_samb, samb_select, toroidal_priority
        )

        # save samb.
        self["atomic_samb"] = atomic_samb
        self["atomic_id"] = atomic_id
        self["cluster_samb"] = cluster_samb
        self["cluster_id"] = cluster_id
        self["combined_samb"] = combined_samb
        self["combined_id"] = combined_id
        self["common_id"] = common_id
        self["cluster_representative"] = cluster_rep
        self["SAMB_number_min"] = combined_min_num
        self["SAMB_number"] = combined_num

        self["version"] = __version__
        self["created"] = time_stamp()

        comment = f"Model: {self["model"]}\n"
        comment += f"* Group: " + self.group.name() + "\n"
        comment += f"* SAMB selection: {self["SAMB_select"]}\n"
        comment += f"* atomic selection: {self["atomic_select"]}\n"
        comment += f"* site-cluster selection: {self["site_select"]}\n"
        comment += f"* bond-cluster selection: {self["bond_select"]}\n"
        comment += f"  {combined_min_num} (all {combined_num}) basis set"
        self.add_comment(comment)
        if self.verbose:
            print(f"{combined_min_num} (all {combined_num}) SAMBs are created.")

    # ==================================================
    def clear(self):
        """
        Clear all data.
        """
        super().clear()
        self._group = None

    # ==================================================
    def get_atomic_samb(self, select):
        """
        Get all atomic SAMB.

        Args:
            select (dict): select dict.

        Returns:
            - (dict) -- atomic SAMB, dict[BraketInfoType, SAMB Dict].
            - (dict) -- atomic id, dict["x#", (SAMB index, comp)].
        """
        basis_type = self["basis_type"]

        samb = {}
        for lst in self["braket"].values():
            for braket_info in lst:
                if braket_info not in samb:
                    samb[braket_info] = self.group.atomic_samb(
                        basis_type, (braket_info.bh_rank, braket_info.kt_rank), (braket_info.bh_idx, braket_info.kt_idx)
                    ).select(**select)

        no = 1
        atomic_id = {}
        for i in samb.values():
            for idx, (mat, ex) in i.items():
                for c in range(len(ex)):
                    atomic_id[f"x{no}"] = (idx, c)
                    no += 1

        return samb, atomic_id

    # ==================================================
    def get_cluster_samb(self, site_select, bond_select):
        """
        Get all cluster SAMB.

        Args:
            site_select (dict): site select dict.
            bond_select (dict): bond select dict.

        Returns:
            - (dict) -- cluster SAMB, dict[wyckoff, SAMB Dict].
            - (dict) -- cluster id, dict["y#", (wyckoff, SAMB index, comp)].
        """
        wp_lst = sorted(list(set([lst.wyckoff for lst in self["site"]["representative"].values()])), key=lambda i: int(i[:-1]))
        site_samb = {
            wp: self.group.cluster_samb(wp).select(**site_select).sort("Gamma", "l", "k", ("X", ["Q", "G", "T", "M"]), "n")
            for wp in wp_lst
        }

        wp_lst = sorted(
            list(set([lst.wyckoff for lst in self["bond"]["representative"].values()])),
            key=lambda i: int(i.split("@")[0][:-1]),
        )
        bond_samb = {
            wp: self.group.cluster_samb(wp, "bond")
            .select(**bond_select)
            .sort("Gamma", "l", "k", ("X", ["Q", "G", "T", "M"]), "n")
            for wp in wp_lst
        }

        samb = site_samb | bond_samb

        no = 1
        cluster_id = {}
        for wp, i in samb.items():
            for idx, (mat, ex) in i.items():
                for c in range(len(ex)):
                    cluster_id[f"y{no}"] = (wp, idx, c)
                    no += 1

        return samb, cluster_id

    # ==================================================
    def get_combined_samb(self, atomic_samb, cluster_samb, select, toroidal_priority=False):
        """
        Get all combined SAMB.

        Args:
            atomic_samb (dict): all atomic SAMB, Dict[braket_info, SAMB].
            cluster_samb (dict): alll cluster SAMB, Dict["site/bond", SAMB].
            select (dict): SAMB select.
            toroidal_priority (bool, optional): use (G,T) prior to (Q,M) in creation ?

        Returns:
            - (dict) -- combined SAMB (minimal), dict[SAMBType, SAMB].
            - (dict) -- id to SAMB info, dict[str, (tag, UniqueSAMBType)].
            - (int) -- no. of minimal SAMBs.
            - (int) -- no. of all SAMBs.
            - (dict) -- comomon id, dict[SAMBType, [id]].
            - (dict) -- cluster name id, dict[site/bond_name, [id]].
        """
        lst = {}
        for site_bond, wp in self["wyckoff"].items():
            tail, head = get_tail_head(site_bond)
            neighbor, n = get_neighbor_info(site_bond)
            braket_info_lst = self["braket"][site_bond]
            for braket_info in braket_info_lst:
                comb = SAMBType(head, tail, wp, braket_info)
                lst[comb] = lst.get(comb, []) + [(neighbor, n)]

        # sort (neighbor, n) for each comb
        lst = {k: sorted(v) for k, v in lst.items()}

        dic_minimal = {}
        for comb in lst.keys():
            bra_orb = self["site"]["representative"][comb.head].orbital[comb.bk_info.bh_rank]
            ket_orb = self["site"]["representative"][comb.tail].orbital[comb.bk_info.kt_rank]
            a_samb = atomic_samb[comb.bk_info]
            c_samb = cluster_samb[comb.wyckoff]

            if comb.head != comb.tail and bra_orb != ket_orb:
                c_samb = c_samb.select(X=["Q"])

            dic_minimal[comb] = self.group.combined_samb(a_samb.named_keys(), c_samb.named_keys(), toroidal_priority, **select)

        min_no = sum(sum(len(i[0]) for i in samb.values()) for samb in dic_minimal.values())
        combined_id = {}
        common_id = {comb: [[] for _ in nn_lst] for comb, nn_lst in lst.items()}

        # Flatten all (neighbor, n, comb) tuples and keep the local index i
        site_global_items = [
            (comb.tail, comb.head, neighbor, n, comb, i)
            for comb, nn_lst in lst.items()
            for i, (neighbor, n) in enumerate(nn_lst)
            if neighbor == 0 and n == -1
        ]
        site_global_items.sort(key=lambda t: (t[0], t[1], t[2], t[3]))
        bond_global_items = [
            (comb.tail, comb.head, neighbor, n, comb, i)
            for comb, nn_lst in lst.items()
            for i, (neighbor, n) in enumerate(nn_lst)
            if neighbor > 0
        ]
        bond_global_items.sort(key=lambda t: (t[0], t[1], t[2], t[3]))

        global_items = site_global_items + bond_global_items

        no = 1
        for Gamma in self["SAMB_select"]["Gamma"]:
            for _, _, neighbor, n, comb, i in global_items:
                samb_info = UniqueSAMBType(comb, neighbor, n)
                if neighbor == 0:
                    sb_tag = comb.tail
                else:
                    sb_tag = get_bond(comb.tail, comb.head, neighbor, n)
                for idx in dic_minimal[comb].select(**{"Gamma": Gamma}).keys():
                    for tag in self.group.tag_multipole(idx, latex=True, superscript="c"):
                        zi = f"z{no}"
                        combined_id[zi] = (tag, samb_info)
                        common_id[comb][i].append((zi, sb_tag))
                        no += 1
        common_id = {info: ([[i[0] for i in bk] for bk in lst], [bk[0][1] for bk in lst]) for info, lst in common_id.items()}

        # cluster-name id.
        cluster_rep = {}
        for info, lst in common_id.items():
            tail, head = info.tail, info.head
            bk_info = info.bk_info
            bra_rank = bk_info.bh_rank
            ket_rank = bk_info.kt_rank
            for i, tag in zip(*lst):
                tail, head = get_tail_head(tag)
                if tag.count(";") > 0:  # bond.
                    rep = self["bond"]["cell"][tag][0]
                    head_sl, tail_sl = rep.h_idx[0], rep.t_idx[0]
                    R = rep.R_primitive
                else:  # site.
                    rep = self["site"]["cell"][tag][0]
                    head_sl = tail_sl = rep.sublattice
                    R = (0, 0, 0)
                bra_start, bra_dim = self["full_matrix"]["index"][(head, head_sl, bra_rank)]
                ket_start, ket_dim = self["full_matrix"]["index"][(tail, tail_sl, ket_rank)]
                r = ((bra_start, bra_dim), (ket_start, ket_dim), R)
                cluster_rep[tag] = cluster_rep.get(tag, []) + [(i, r)]

        return dic_minimal, combined_id, min_no, no - 1, common_id, cluster_rep

    # ==================================================
    def select_combined_samb(self, **kwargs):
        """
        Select combined SAMB.

        Args:
            kwargs (dict, optional): select conditions for multipoles with keywords, "site/bond/X/l/Gamma/s".

        Returns:
            - (dict) -- selected combined IDs.
            - (dict) -- select dict for combined SAMB.
            - (dict) -- regularized select dict.

        Note:
            - site = [(site, *[orbital_rank])]. (* omittable).
            - bond = [(site1;site2, *rank1;rank2, *[neighbor])] or [(*site1;site2, *rank1;rank2, [neighbor])].
            - X = Q/G/T/M.
            - l = [0,1,2,3,4,5,6,7,8,9,10,11].
            - Gamma = [irreps.].
            - s = [0,1].
        """
        # regularize kwargs.
        kwargs = {k: v if isinstance(v, list) else [v] for k, v in kwargs.items()}
        inp = {}
        for k, v in kwargs.items():
            if k == "site":
                inp["site"] = []
                for i in v:
                    if isinstance(i, str):  # site.
                        inp["site"].append((i, None))
                    else:
                        a, b = i  # site, rank/[rank]
                        if not isinstance(b, list):
                            b = [b]
                        inp["site"].append((a, b))
            elif (
                k == "bond"
            ):  # name, neighbor, [neighbor], (name, neighbor), (name, [neighbor]), (name, rank), (name, rank, neighbor), (name, rank, [neighbor]).
                inp["bond"] = []
                for i in v:
                    if isinstance(i, str):  # tail;head.
                        inp["bond"].append((i, None, None))
                    elif isinstance(i, int):  # neighbor.
                        inp["bond"].append((None, None, [i]))
                    elif isinstance(i, list) and all(isinstance(x, int) for x in i):  # neighbor.
                        inp["bond"].append((None, None, i))
                    elif len(i) == 2:
                        a, b = i
                        if isinstance(b, list):  # tail;head, neighbor.
                            inp["bond"].append((a, None, b))
                        elif isinstance(b, str):  # tail;head, t_rank;h_rank.
                            inp["bond"].append((a, b, None))
                        else:
                            inp["bond"].append((a, None, [b]))  # tail;head, neighbor.
                    elif len(i) == 3:
                        a, b, c = i
                        if isinstance(c, list):
                            inp["bond"].append((a, b, c))  # tail;head, t_rank;h_rank, neighbor.
                        else:
                            inp["bond"].append((a, b, [c]))  # tail;head, t_rank;h_rank, neighbor.
            else:
                inp[k] = v

        # default filter.
        default = {}
        all_site = [(k, [no for no, i in enumerate(v.orbital) if len(i) > 0]) for k, v in self["site"]["representative"].items()]
        default["site"] = all_site
        all_bond = [(v.tail, v.head, v.neighbor, v.t_rank, v.h_rank) for v in self["bond"]["representative"].values()]
        default["bond"] = all_bond
        default.update(self["SAMB_select"])

        # final filter.
        dic = {"site": [], "bond": [], "X": [], "l": [], "Gamma": [], "s": []}
        for key, val in inp.items():
            if key in ["X", "l", "Gamma", "s"]:
                dic[key] = dic.get(key, []) + [i for i in default[key] if i in val]
            elif key == "site":
                for name, orb in val:
                    if orb is None:
                        d = [(s, o) for s, o in default[key] if s == name]
                    else:
                        d = [(s, sorted(list(set(o) & set(orb)))) for s, o in default[key] if s == name]
                    dic[key] = dic.get(key, []) + [i for i in d if len(i[1]) > 0]
                d = list({tuple(tuple(x) if isinstance(x, list) else x for x in t) for t in dic[key]})
                d = [tuple(list(x) if isinstance(x, tuple) and all(isinstance(i, int) for i in x) else x for x in t) for t in d]
                dic[key] = sorted(d)
            elif key == "bond":
                for name, rank, neighbor in val:
                    if rank is None and neighbor is None:  # name only.
                        tail, head = name.split(";")
                        d = [(h, t, n, hr, tr) for t, h, n, tr, hr in default[key] if t == tail and h == head]
                        d += [(h, t, n, hr, tr) for t, h, n, tr, hr in default[key] if h == tail and t == head]
                    elif name is None and rank is None:  # neighbor only.
                        d = [(h, t, n, hr, tr) for t, h, n, tr, hr in default[key] if n in neighbor]
                    elif rank is None:  # name, neighbor.
                        tail, head = name.split(";")
                        d = [(h, t, n, hr, tr) for t, h, n, tr, hr in default[key] if t == tail and h == head and n in neighbor]
                        d += [(h, t, n, hr, tr) for t, h, n, tr, hr in default[key] if h == tail and t == head and n in neighbor]
                    elif neighbor is None:  # name, rank.
                        tail, head = name.split(";")
                        t_rank, h_rank = rank.split(";")
                        t_rank, h_rank = {int(t_rank)}, {int(h_rank)}
                        d = [
                            (h, t, n, sorted(list(set(hr) & h_rank)), sorted(list(set(tr) & t_rank)))
                            for t, h, n, tr, hr in default[key]
                            if t == tail and h == head
                        ]
                        d += [
                            (h, t, n, sorted(list(set(hr) & h_rank)), sorted(list(set(tr) & t_rank)))
                            for t, h, n, tr, hr in default[key]
                            if h == tail and t == head
                        ]
                    else:
                        tail, head = name.split(";")
                        t_rank, h_rank = rank.split(";")
                        t_rank, h_rank = set(int(t_rank)), set(int(h_rank))
                        d = [
                            (h, t, n, sorted(list(set(hr) & h_rank)), sorted(list(set(tr) & t_rank)))
                            for t, h, n, tr, hr in default[key]
                            if t == tail and h == head and n in neighbor
                        ]
                        d += [
                            (h, t, n, sorted(list(set(hr) & h_rank)), sorted(list(set(tr) & t_rank)))
                            for t, h, n, tr, hr in default[key]
                            if h == tail and t == head and n in neighbor
                        ]
                    dic[key] = dic.get(key, []) + d
                d = list({tuple(tuple(x) if isinstance(x, list) else x for x in t) for t in dic[key]})
                d = [tuple(list(x) if isinstance(x, tuple) and all(isinstance(i, int) for i in x) else x for x in t) for t in d]
                dic[key] = sorted(d)
        for key, val in default.items():
            if key not in kwargs.keys():
                dic[key] = val

        # create select_dict for SAMB and combined_id.
        select_dict = {k: dic[k] for k in ["X", "l", "Gamma", "s"]}

        combined_id = {}
        for zi, (tag, samb_info) in self["combined_id"].items():
            samb_type = samb_info.samb_type
            head, tail = samb_type.head, samb_type.tail
            bk_info = samb_type.bk_info
            neighbor = samb_info.neighbor
            is_site = not ("@" in samb_type.wyckoff)
            for site, rank in dic["site"]:
                if is_site and head == site and bk_info.bh_rank in rank:
                    for idx in self["combined_samb"][samb_type].select(**select_dict).keys():
                        for tag_ in self.group.tag_multipole(idx, latex=True, superscript="c"):
                            if tag == tag_:
                                combined_id[zi] = (tag, samb_info)
            for h, t, nr, h_rank, t_rank in dic["bond"]:
                if (
                    not is_site
                    and tail == t
                    and head == h
                    and bk_info.kt_rank in t_rank
                    and bk_info.bh_rank in h_rank
                    and neighbor == nr
                ) or (
                    not is_site
                    and tail == h
                    and head == t
                    and bk_info.kt_rank in h_rank
                    and bk_info.bh_rank in t_rank
                    and neighbor == nr
                ):
                    for idx in self["combined_samb"][samb_type].select(**select_dict).keys():
                        for tag_ in self.group.tag_multipole(idx, latex=True, superscript="c"):
                            if tag == tag_:
                                combined_id[zi] = (tag, samb_info)

        return combined_id, select_dict, dic

    # ==================================================
    def get_combined_samb_matrix(self, fmt="sympy", digit=14, **kwargs):
        """
        Get combined SAMBs in matrix form (real-space).

        Args:
            fmt (str, optional): sympy/value.
            digit (int, optional): digit for value output.
            kwargs (dict, optional): select conditions for multipoles with keywords (see, select_combined_samb).

        Returns:
            - (dict) -- combined SAMB in matrix form, dict[zj, dict[ (n1, n2, n3, m, n), matrix element] ].

        Note:
            - R = (n1,n2,n3) and m and n are lattice indices, bra and ket indexes, respectively.
        """
        # check format key.
        if fmt not in ["sympy", "value"]:
            raise KeyError(f"unknown format = {fmt} is given.")

        # check fileter key.
        filter_key = ["site", "bond", "neighbor", "n", "bk_info"]
        for k in kwargs.keys():
            if k not in filter_key + list(self["SAMB_select"].keys()):
                raise KeyError(f"invalid key, {k}.")

        # filter combined id.
        combined_id, select_dict, _ = self.select_combined_samb(**kwargs)

        cluster_samb = self["cluster_samb"]
        atomic_samb = self["atomic_samb"]
        combined_samb = self["combined_samb"]

        X_cache = {}
        matrix = {}

        def _X(a_samb, t1, c1):
            key = (id(a_samb), t1, c1)
            m = X_cache.get(key)
            if m is None:
                m = sp.Matrix(a_samb[t1][0][c1])
                X_cache[key] = m
            return m

        def _hermite_conj(head, head_sl, bra_rank, tail, tail_sl, ket_rank):
            hermite_conj = False
            if (head, head_sl, bra_rank) not in self["full_matrix"]["index"]:
                hermite_conj = True
            if (tail, tail_sl, ket_rank) not in self["full_matrix"]["index"]:
                hermite_conj = True
            return hermite_conj

        def _format_val(v):
            v = sp.expand(v)
            if fmt == "value":
                v = complex(v)
                return round(v.real, digit) + round(v.imag, digit) * 1j
            return v

        for zi, (tag, samb_info) in combined_id.items():
            samb_type = samb_info.samb_type

            tail, head = samb_type.tail, samb_type.head
            wp = samb_type.wyckoff

            bk_info = samb_type.bk_info
            bra_rank, bra_idx = bk_info.bh_rank, bk_info.bh_idx
            ket_rank, ket_idx = bk_info.kt_rank, bk_info.kt_idx

            bn, bm = samb_info.neighbor, samb_info.n

            is_site = "@" not in wp
            if is_site:
                head_sl = tail_sl = 1
            else:
                name = get_bond(tail, head, bn, bm)
                bond = self["bond"]["cell"][name][0]
                head_sl, tail_sl = bond.h_idx[0], bond.t_idx[0]

            hermite_conj = _hermite_conj(head, head_sl, bra_rank, tail, tail_sl, ket_rank)
            if hermite_conj:
                bra_idx, bra_rank, ket_idx, ket_rank = ket_idx, ket_rank, bra_idx, bra_rank

            a_samb = atomic_samb[bk_info]
            c_samb = cluster_samb[wp]

            for idx, (cl, _) in combined_samb[samb_type].select(**select_dict).items():
                for tag_, m in zip(self.group.tag_multipole(idx, latex=True, superscript="c"), cl):
                    if tag_ != tag:
                        continue

                    d = defaultdict(lambda: sp.S(0) if fmt == "sympy" else 0.0)
                    for cg, t1, c1, t2, c2 in m:
                        X = _X(a_samb, t1, c1)
                        Xd = X.adjoint()
                        Y = c_samb[t2][0][c2]
                        for i, y in enumerate(Y):
                            if is_site:
                                head_sl = tail_sl = i + 1
                                n1 = n2 = n3 = 0
                            else:
                                bond = self["bond"]["cell"][name][i]
                                n1, n2, n3 = bond.R_primitive
                                head_sl, tail_sl = bond.h_idx[0], bond.t_idx[0]

                            bra_start, bra_dim = self["full_matrix"]["index"][(head, head_sl, bra_rank)]
                            ket_start, ket_dim = self["full_matrix"]["index"][(tail, tail_sl, ket_rank)]

                            for r, c in product(range(bra_dim), range(ket_dim)):
                                val = cg * y * (Xd[r, c] if hermite_conj else X[r, c])
                                key = (n1, n2, n3, bra_start + r, ket_start + c)
                                d[key] += val

                            if not is_site:
                                if head == tail and (bra_rank, bra_idx) != (ket_rank, ket_idx):
                                    if head_sl == tail_sl:
                                        bra2_start = ket_start
                                        ket2_start = bra_start
                                    else:
                                        ket_list = [o for o in self["full_matrix"]["ket"] if o[0] == head and o[1] == 1]
                                        bra_orb = [o for o in ket_list if o[2] == bra_rank][0]
                                        ket_orb = [o for o in ket_list if o[2] == ket_rank][0]
                                        bra_orb_idx = ket_list.index(bra_orb)
                                        ket_orb_idx = ket_list.index(ket_orb)

                                        bra2_start = bra_start + (ket_orb_idx - bra_orb_idx)
                                        ket2_start = ket_start - (ket_orb_idx - bra_orb_idx)

                                    for r, c in product(range(ket_dim), range(bra_dim)):
                                        val = cg * y * Xd[r, c]
                                        key = (n1, n2, n3, bra2_start + r, ket2_start + c)
                                        d[key] += val

                    # hermite conjugate
                    for Rmn, val in list(d.items()):
                        n1_, n2_, n3_, m_, n_ = Rmn
                        d[(-n1_, -n2_, -n3_, n_, m_)] += sp.conjugate(val)

                    norm = sp.sqrt(sp.expand(sum([v * sp.conjugate(v) for v in d.values()], sp.S(0))))

                    matrix[zi] = {Rmn: _format_val(v / norm) for Rmn, v in d.items() if not v.is_zero}

        matrix = dict(sorted(matrix.items(), key=lambda x: int(x[0][1:])))

        return matrix

    # ==================================================
    def get_hr(self, coeff, combined_samb_matrix=None, fmt="sympy", digit=14, **kwargs):
        """
        Get Hamiltonian matrix (real-space).

        Args:
            coeff (dict): coefficient of SAMBs, dict[zj, float/sympy].
            combined_samb_matrix (dict, optional): combined SAMBs in matrix form (real-space), { zj: {(n1, n2, n3, m, n): matrix element} }.
            fmt (str, optional): sympy/value.
            digit (int, optional): digit for value output.
            kwargs (dict, optional): select conditions for multipoles with keywords (see, select_combined_samb).

        Returns:
            - (dict) -- Hamiltonian matrix (real-space), dict[(n1, n2, n3, m, n), matrix element].

        Note:
            - R = (n1,n2,n3) and m and n are a lattie vector, bra and ket indexes, respectively.
        """
        if combined_samb_matrix is None:
            combined_samb_matrix = self.get_combined_samb_matrix(fmt, digit, **kwargs)

        Hamiltonian = defaultdict(lambda: sp.S(0) if fmt == "sympy" else 0.0)
        for zj, cj in coeff.items():
            if zj not in combined_samb_matrix.keys():
                raise Exception(f"parameter {zj} is missing.")
            d = combined_samb_matrix[zj]
            for Rmn, Zj in d.items():
                Hamiltonian[Rmn] += cj * Zj

        return Hamiltonian

    # ==================================================
    def get_multipole_expression(self):
        """
        Get all of harmonics expression.

        Returns:
            - (dict) -- harmonics expression, dict[(X,Gamma,l,n), [expression]].
        """
        amp = [v[:4] for v in sum([list(samb.keys()) for samb in self["atomic_samb"].values()], [])]
        scmp = [v[:4] for v in sum([list(samb.keys()) for samb in self["cluster_samb"].values()], [])]
        cmp = [v[:4] for v in sum([list(samb.keys()) for samb in self["combined_samb"].values()], [])]
        lst = sorted(list(set(amp + scmp + cmp)), key=lambda i: (i[2], i[1], i[0], i[3]))

        harmonics = self.group.harmonics
        harmonics_lst = {}
        for X, l, Gamma, n in lst:
            X = X.replace("T", "Q").replace("M", "G")
            idx = (X, l, Gamma, n, -1, 0, 0, "q")
            if idx in harmonics.keys():
                harmonics_lst[(X, Gamma, l, n)] = harmonics[idx][0]
            # else: # should be never reached.
            #    raise Exception(f"cannot find idx={idx1[1:]}.")

        return harmonics_lst

    # ==================================================
    def ket(self):
        """
        Get ket string list.

        Returns:
            - (list) -- ket string in LaTeX.
        """
        lst = []
        for atom, sublattice, rank, orbital in self["full_matrix"]["ket"]:
            orb = self.group.tag_atomic_basis(orbital, rank, latex=True)
            orb += "@" + r"{\rm " + atom + "}(" + str(sublattice) + ")"
            lst.append(orb)

        return lst

    # ==================================================
    def cell_site_primitive(self):
        """
        Get site position in primitive unit cell.

        Returns:
            - (dict) -- site dict, dict[name, [position]].
        """
        lst = {
            name: [i.position_primitive.tolist() for i in val if i.plus_set == 1] for name, val in self["site"]["cell"].items()
        }
        return lst

    # ==================================================
    def _write_site(self):
        """
        Write site info.
        """
        write_site_dict(self["site"])

    # ==================================================
    def _write_bond(self):
        """
        Write bond info.
        """
        write_bond_dict(self["bond"])

    # ==================================================
    def _write_grid(self):
        """
        Write grid info.
        """
        write_site_grid(self["site_grid"])

    # ==================================================
    def plot_site_bond(self, qtdraw, site_bond, rep=True):
        """
        Plot definition of site or bond.

        Args:
            qtdraw (QtDraw): QtDraw object.
            site_bond (str): site or bond, see "wyckoff".
            rep (bool, optional): highlight representative ?
        """
        if site_bond not in self["wyckoff"].keys():
            raise KeyError(f"unknown site_bond, {site_bond}.")
        if site_bond.count(";") > 0:
            bonds = np.array([np.concat([i.vector, i.center]) for i in self["bond"]["cell"][site_bond]])
            plot_bond(qtdraw, site_bond, bonds)
        else:
            sites = np.array([i.position for i in self["site"]["cell"][site_bond]])
            plot_site(qtdraw, site_bond, sites, rep)

    # ==================================================
    def plot_cluster_samb(self, qtdraw, site_bond, cluster_id, label=True):
        """
        Plot cluster SAMB.

        Args:
            qtdraw (QtDraw): QtDraw object.
            site_bond (str): site or bond, see "wyckoff".
            cluster_id (str): cluster id.
            label (boo, optional): display label ?
        """
        if site_bond not in self["wyckoff"].keys():
            raise KeyError(f"unknown site_bond, {site_bond}.")
        if cluster_id not in self["cluster_id"].keys():
            raise KeyError(f"unknown cluster id, {cluster_id}.")

        wp, idx, comp = self["cluster_id"][cluster_id]
        if self["wyckoff"][site_bond] != wp:
            raise KeyError(f"invalid wyckoff, {wp}.")

        samb = self["cluster_samb"][wp][idx][0][comp]

        if label:
            qtdraw.add_text2d(f"idx = ({",".join(map(str,idx[:4]))},{comp})")
        if site_bond.count(";") > 0:
            bonds = np.array([np.concat([i.vector, i.center]) for i in self["bond"]["cell"][site_bond]])
            sym = idx[0] not in ["T", "M"]
            if not sym:
                samb = np.vectorize(sp.im)(samb)
            plot_bond_samb(qtdraw, cluster_id + "@" + site_bond, bonds, samb, sym, cluster_id)
        else:
            sites = np.array([i.position for i in self["site"]["cell"][site_bond]])
            plot_site_samb(qtdraw, cluster_id + "@" + site_bond, sites, samb, cluster_id)

    # ==================================================
    def plot_atomic_samb(self, qtdraw, atomic_id, site_bond=None, label=True):
        """
        Plot atomic SAMB.

        Args:
            qtdraw (QtDraw): QtDraw object.
            atomic_id (str): atomic id.
            site_bond (str, optional): site or bond, see "wyckoff".
            label (boo, optional): display label ?
        """
        conv_dict = {"Q": "Q", "G": "G", "T": "Q", "M": "G"}
        if site_bond is not None and site_bond not in self["wyckoff"].keys():
            raise KeyError(f"unknown site_bond, {site_bond}.")
        if atomic_id not in self["atomic_id"].keys():
            raise KeyError(f"unknown atomic id, {atomic_id}.")

        idx, comp = self["atomic_id"][atomic_id]

        if label:
            qtdraw.add_text2d(f"idx = ({",".join(map(str,idx[:4]))},{comp})")
        if site_bond is None:
            point = [[0, 0, 0]]
        else:
            if site_bond.count(";") > 0:
                point = np.array([np.concat([i.vector, i.center]) for i in self["bond"]["cell"][site_bond]])
            else:
                point = np.array([i.position for i in self["site"]["cell"][site_bond]])

        idx = tuple([conv_dict[idx[0]]] + list(idx[1:]))
        samb = self.group.harmonics[idx][0][comp]
        name = atomic_id if site_bond is None else atomic_id + "@" + site_bond
        plot_harmonics(qtdraw, name, samb, point, atomic_id)
