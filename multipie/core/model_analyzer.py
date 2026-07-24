"""
Model analyzer class.

This module provides model analyzer.
"""

import os
import logging
import numpy as np
import seekpath
from multipie.core.material_model import MaterialModel
from multipie.core.default_control import default_control
from multipie.util.util_model_analyzer import (
    grid_path,
    fourier_r_to_k,
    fourier_k_to_r,
    output_dispersion,
    create_gnuplot_cmd,
    plot_save_dispersion,
    create_all_local_operator,
    create_local_operator,
    create_k_multipole,
    create_k_matrix,
    add_local_parameter,
    convert_zj_atomic_var,
    fermi_dirac,
)
from multipie.util.util_wannier import (
    read_win,
    read_nnkp,
    merge_wannier_info,
    read_hr,
    read_mmn,
    read_spn,
    read_uHu,
    read_uIu,
    build_ket_wannier,
    sort_ket_list,
    sort_ket_matrix_dict,
    decompose_operator_by_SAMB,
)
from multipie.util.util import read_dict, str_to_sympy, write_dict

_k_matrix_comment = """Selected SAMB matrix in momentum representation.
- model (str): model name.
- source (str): source binary.
- date (str): binary created date.
- dimension (int): matrix size.
- ket_site (list): ket info., [ket_name].
- index (dict): ket index, dict[(site,sublattice,rank), (top_index,size)].
- cluster_vector (dict): cluster vector, dict[site/bond name, dict[kb, expression] ].
- k_multipole (dict): momentum multipole in terms of p_n=k.b_n, dict[wyckoff, dict[idx, (k_multipole, symmetry)] ].
- k_matrix (dict): momentum matrix, dict[tag, (site/bond_name, wyckoff, dict[(m,n), value]) ].
"""

_zj_var_comment = """Correspondence between zj and atomic variable.
- correspondence for each bond cluster, dict[bond_name, dict[zj, expression in terms of atomic variables] ].
- only for SAMB with identity irrep.
"""

_param_comment = """Parameter dict (sorted by descending absolute value).
- finite parameter, dict[zj, value].
"""


# ==================================================
class ModelAnalyzer(dict):
    # ==================================================
    def __init__(self, N1=50, N2=50, N3=50, topdir=None, verbose=False):
        """
        Model analyzer.

        Args:
            N1 (int, optional): number of divisions in a1.
            N2 (int, optional): number of divisions in a2.
            N3 (int, optional): number of divisions in a3.
            topdir (str, optional): top directory. [default: cwd]
            verbose (bool, optional): verbose comment ?
        """
        if topdir is None:
            topdir = os.getcwd()

        self._topdir = topdir
        self._verbose = verbose
        self._mm = MaterialModel(topdir, verbose=verbose)
        self._local = create_all_local_operator()

        self._mode = default_control["mode"]
        self._samb = default_control["samb"]
        self._wannier = default_control["wannier"]
        self._output = default_control["output"]

        self._name = None
        self._parameter = None
        self._HR = None
        self._basis_type = None
        self._basis = None

        self.set_grid_size(N1, N2, N3)
        self["samb"] = {}
        self["wannier"] = {}
        self["output"] = {}

        os.chdir(self._topdir)

    # ==================================================
    @property
    def samb(self):
        """
        SAMB control.

        Returns:
            - (dict) -- SAMB control.
        """
        return self._samb

    # ==================================================
    @property
    def wannier(self):
        """
        Wannier control.

        Returns:
            - (dict) -- Wannier control.
        """
        return self._wannier

    # ==================================================
    @property
    def output(self):
        """
        Output for physical quantities control.

        Returns:
            - (dict) -- output for physical quantities control.
        """
        return self._output

    # ==================================================
    @property
    def model(self):
        """
        Material Model.

        Returns:
            - (MaterialModel) -- matrial model.
        """
        return self._mm

    # ==================================================
    @property
    def mode(self):
        """
        Analysis mode.

        Returns:
            - (str) -- mode, "samb/wannier/symcw".
        """
        return self._mode

    # ==================================================
    @property
    def name(self):
        """
        Model name.

        Returns:
            - (str) -- model name.
        """
        return self._name

    # ==================================================
    @property
    def parameter(self):
        """
        Parameter of SAMBs.

        Returns:
            - (dict) -- parameter dict.
        """
        return self._parameter

    # ==================================================
    @property
    def basis_type(self):
        """
        Atomic basis type.

        Returns:
            - (str) -- basis type, "lg/lgs/jml".
        """
        return self._basis_type

    # ==================================================
    @property
    def basis(self):
        """
        Full-matrix basis.

        Returns:
            - (list) -- list of basis, (atom, sublattice, rank, component, tag).
        """
        return self._basis

    # ==================================================
    @property
    def HR(self):
        """
        Hamiltonian matrix in real space.

        Returns:
            - (dict) -- H(R).
        """
        return self._HR

    # ==================================================
    def set_grid_size(self, N1, N2, N3):
        """
        Set grid size.

        Args:
            N1 (int): number of divisions in a1.
            N2 (int): number of divisions in a2.
            N3 (int): number of divisions in a3.
        """
        self["grid"] = [N1, N2, N3]

    # ==================================================
    def set_primitive_cell(self, A):
        """
        Set primitive cell info.

        Args:
            A (ndarray): [a1, a2, a3] (3x3).

        :meta private:
        """
        B = 2 * np.pi * np.linalg.inv(A).T
        self["A"] = A  # primitive cell.
        self["B"] = B  # reciprocal cell.
        self["unit_cell_volume"] = float(np.dot(A[0], np.cross(A[1], A[2])))  # volume of primitive cell.

    # ==================================================
    def analyze(self, control):
        """
        Analyze model with control file.

        Args:
            control (str or dict): control file (.py) or model name.
        """
        # read control.
        if type(control) == str:  # read control file or from dict.
            if control.endswith(".py"):
                file = os.path.join(self._topdir, control)
                control = read_dict(file)
            else:  # w/o control and model_name is given.
                self.model.load(control)
                matrix_info = self.model.get_samb_matrix({})
                self.model.save_samb_matrix(matrix_info)
                return

        # set property dict.
        self._samb |= control.get("samb", {})
        self._wannier |= control.get("wannier", {})
        self._output |= control.get("output", {})

        # update mode.
        mode = control.get("mode", None)
        if mode:
            self._mode = mode

        # set SAMB.
        if mode in ["samb", "symcw"]:
            self.set_samb()  # create SAMBs, and H(R) if zj are provided.

        # set wannier.
        if mode in ["wannier", "symcw"]:
            hr_dict = self.set_wannier()  # create H(R), and zj in case of "symcw".
            if mode == "symcw":
                matrix_info = self["samb"]["matrix_info"]
                Zr_dict = matrix_info["matrix"]
                parameter = decompose_operator_by_SAMB(hr_dict, Zr_dict)
                self._HR = self.model.get_hr(parameter, Zr_dict)
                self.model.save_samb_hr(matrix_info, parameter, self.HR)
                self._parameter = parameter
            else:
                self._HR = hr_dict

        # create z file.
        if self.parameter:
            z_file = self.name + "_z.py"
            comment = _param_comment + f"- by using '{mode}' mode.\n"
            dic = {tag: float(v) for tag, v in self.parameter.items()}
            dic = dict(sorted(dic.items(), key=lambda item: abs(item[1]), reverse=True))
            write_dict(dic, z_file, comment=comment, w_dir=self.name)
            if self._verbose:
                print(f"save z to '{self._topdir}/{z_file}'.")

        # compute physical quanties and output data.
        self.compute_physical_quantity()

    # ==================================================
    def set_samb(self):
        """
        Calculate SAMB related quantities.

        :meta private:
        """
        name = self.samb["model"]
        self._name = name

        # read model.
        self.model.load(name)
        self._basis_type = self.model["basis_type"]
        self._basis = self.model["full_matrix"]["ket"]

        # set primitive cell.
        self.set_primitive_cell(np.array(self.model["unit_vector_primitive"]))

        # set selected SAMBs.
        select = self.samb.get("select", {})
        matrix_info = self.model.get_samb_matrix(select)
        self["samb"]["matrix_info"] = matrix_info

        # create var file.
        IR = next(iter(self.model.group.character["table"].keys()))  # identity irrep.
        conv_dict = convert_zj_atomic_var(matrix_info, self.model["combined_cluster"], self.model["combined_id"], IR)
        conv_dict = {name: {zj: str(ex).replace(" ", "") for zj, ex in dic.items()} for name, dic in conv_dict.items()}
        var_file = name + "_var.py"
        write_dict(conv_dict, var_file, comment=_zj_var_comment, w_dir=name)
        if self._verbose:
            print(f"save var to '{self._topdir}/{var_file}'.")

        # create k-multipole file.
        if self.samb.get("k_multipole", False):
            k_multipole = self.set_k_multipole(matrix_info)
            k_file = name + "_k.py"
            write_dict(k_multipole, k_file, comment=_k_matrix_comment, w_dir=name)
            if self._verbose:
                print(f"save k-multipole to '{self._topdir}/{k_file}'.")
            self["samb"]["k_multipole"] = k_multipole
        else:
            self["samb"]["k_multipole"] = {}

        # create SAMB qtdraw.
        if self.samb.get("samb_figure", False):
            self.model.save_samb_qtdraw()

        parameter = self.samb.get("parameter", {})
        if type(parameter) == str:  # when parameter is str, read z file.
            z_file = parameter
            parameter = read_dict(z_file, name)
            parameter = {tag: float(str_to_sympy(v, rational=False)) if type(v) == str else v for tag, v in parameter.items()}
            if self._verbose:
                print(f"load parameter from '{self._topdir}/{z_file}'.")

        # determine local weight if NG_sum_rule is True.
        if self.samb.get("NG_sum_rule", False) and parameter:
            parameter = add_local_parameter(matrix_info, parameter)

        # output matrix.py and hr.dat.
        if parameter:
            self._HR = self.model.get_hr(parameter, matrix_info["matrix"])
            self.model.save_samb_hr(matrix_info, parameter, self.HR)
        else:
            self._HR = None
        self.model.save_samb_matrix(matrix_info)

        self._parameter = parameter

    # ==================================================
    def set_wannier(self):
        """
        Set data for wannier-based input.

        :meta private:
        """
        topdir = os.path.join(self._topdir, self.name, "wannier")
        seedname = self.wannier.get("seedname", None)

        # read seedname.win
        win = read_win(topdir, seedname)
        # read seedname.nnkp
        nnkp = read_nnkp(topdir, seedname)
        # read seedname_hr.dat
        hr_dict, irvec, _ = read_hr(topdir, self._wannier.get("hr_file", None))
        # read seedname.mmn
        # Mkb = read_mmn(topdir, seedname)
        # read seedname.spn
        # Sk = read_spn(topdir, seedname)
        # read seedname.uHu
        # uHu = read_uHu(topdir, seedname)
        # read seedname.uIu
        # uIu = read_uIu(topdir, seedname)

        # check common values and merge.
        wannier_info = merge_wannier_info(win, nnkp, seedname)

        atoms_list = list(wannier_info["atoms_frac"].values())
        atoms_frac = np.array([atoms_list[i] for i in wannier_info["nw2n"]])

        atoms_list = list(wannier_info["atoms_cart"].values())
        atoms_cart = np.array([atoms_list[i] for i in wannier_info["nw2n"]])

        if True:  # MultiPie dependent part.
            # sort wannier basis as those of MultiPie.
            ket_wannier = self.wannier.get("ket_wannier", [])
            if not ket_wannier:
                site_dict = {
                    (k, vi.sublattice): vi.position_primitive.tolist()
                    for k, v in self._mm["site"]["cell"].items()
                    for vi in v
                    if vi.plus_set == 1
                }
                ket_wannier = build_ket_wannier(nnkp, site_dict, rtol=1e-4, atol=1e-4)

            # convert ket_wannier to ket_multipie.
            ket_multipie = self.model["full_matrix"]["ket"]
            atoms_frac = sort_ket_list(atoms_frac, ket_wannier, ket_multipie)
            atoms_cart = sort_ket_list(atoms_cart, ket_wannier, ket_multipie)
            hr_dict = sort_ket_matrix_dict(hr_dict, ket_wannier, ket_multipie)

        info = {
            # Wannier info.
            "A": wannier_info["A"],
            "B": wannier_info["B"],
            # "ket": ket_multipie,  # sorted as MultiPie ket.
            "num_wann": wannier_info["num_wann"],
            "atoms_frac": atoms_frac,
            "atoms_cart": atoms_cart,
            "spinors": wannier_info["spinors"],
            "fermi_energy": wannier_info["fermi_energy"],
            # DFT info.
            "num_bands": wannier_info["num_bands"],
            "mp_grid": wannier_info["mp_grid"],
            "num_k": wannier_info["num_k"],
            "num_b": wannier_info["num_b"],
            "kpoints": wannier_info["kpoints"],
            "nnkpts": wannier_info["nnkpts"],
            "bvec_cart": wannier_info["bvec_cart"],
            "bvec_crys": wannier_info["bvec_crys"],
            "wb": wannier_info["wb"],
            "wk": wannier_info["wk"],
            "bveck": wannier_info["bveck"],
            "kb2k": wannier_info["kb2k"],
        }

        self["wannier"]["info"] = info
        # self["wannier"]["z_j_exp"] = z_j_exp
        # self["wannier"]["mmn"] = mmn
        # self["wannier"]["spn"] = spn
        # self["wannier"]["uHu"] = uHu
        # self["wannier"]["uIu"] = uIu

        ### physical qunatity.
        # nk = np.array([np.diag(fermi_dirac(eki - win["fermi_energy"], T=0.0)) for eki in Ek], dtype=float)
        # nk = Uk.transpose(0, 2, 1).conjugate() @ nk @ Uk
        # nr_dict = fourier_k_to_r(nk, win["kpoints"], irvec, s=False)
        # nr_dict = sort_ket_matrix_dict(nr_dict, ket_wannier, ket_multipie)
        # z_j_exp = decompose_operator_by_SAMB(nr_dict, Zr_dict)

        return hr_dict

    # ==================================================
    def compute_physical_quantity(self):
        """
        Compute physical quantities by parsing the control file.

        Args:
            name (str): model name.

        :meta private:
        """
        if self._HR is None:
            if self._verbose:
                print("set H(R) first before calculating physical quantities.")
            return

        cwd = os.getcwd()

        outdir = os.path.join(self._topdir, self.name, self.output["dir"])
        os.makedirs(outdir, exist_ok=True)
        os.chdir(outdir)

        self.set_eigen_system()
        self.compute_dispersion()
        self.compute_dos()

        os.chdir(cwd)

    # ==================================================
    def set_eigen_system(self):
        """
        Set eigen system by checking control/output if E and/or U is required.

        :meta private:
        """
        pass

    # ==================================================
    def local_operator(self, name):
        """
        Create local operator.

        Args:
            name (str): operator name, "Sx/Sy/Sz/Lx/Ly/Lz/Qu/Qv/Qyz/Qzx/Qxy".

        Returns:
            - (ndarray) -- operator matrix (dim x dim).

        :meta private:
        """
        return create_local_operator(self.basis, name, self._local, self.basis_type == "lgs")

    # ==================================================
    def compute_dispersion(self):
        """
        Compute dispersion.

        :meta private:
        """
        # check if dispersion can be computed.
        if "dispersion" not in self.output:
            return
        k_path = self.output["dispersion"].get("k_path", None)
        if k_path is None or self.model.group.group_type != "SG":
            return

        name = self.name

        # get k_point and k_path.
        k_point, k_path = self.get_kpath(k_path)
        k_point_path, k_linear, k_dis_pos = grid_path(k_point, k_path, self["grid"][0], self["B"])

        # get local operator list.
        if self._basis_type == "jml" or self._basis_type is None:
            op_lst = []
        else:
            op_lst = self.output["dispersion"].get("local", [])

        # get info.
        tb_gauge = self.output["fourier"]["tb_gauge"]
        atom = np.asarray(list(self.model.get_ket_site().values()), dtype=float)

        # set H(R) and H(k).
        HR = {((n1, n2, n3), m, n): complex(v) for (n1, n2, n3, m, n), v in self._HR.items()}
        Hk = fourier_r_to_k(HR, atom, k_point_path, tb_gauge)

        # set eigen system.
        Ek, Uk = np.linalg.eigh(Hk)
        power = self.output["dispersion"].get("power", None)
        if power is not None:
            Ek = np.power(Ek, power)

        # set local operators.
        Ok = [np.einsum("kmi,mn,kni->ki", Uk.conj(), self.local_operator(name), Uk).real for name in op_lst]

        # output dispersion data, plot, and gnuplot.
        fname = name + "_dispersion.txt"
        colormap = len(Ok) > 0
        if Ok:
            output_dispersion(fname, k_linear, Ek, Ok, op_lst)
        else:
            output_dispersion(fname, k_linear, Ek)
        plot_save_dispersion(fname, k_dis_pos, colormap)
        create_gnuplot_cmd(fname, k_dis_pos, np.max(k_linear), np.max(Ek), np.min(Ek), colormap)
        if self._verbose:
            print(f"save dispersion files into '{self._topdir}/{name}/{self.output["dir"]}.")

        # save dispersion info.
        self["output"]["dispersion"] = {
            "k_path": k_path,
            "k_point": k_point,
            "e_max": np.max(Ek),
            "e_min": np.min(Ek),
            "ef": 0.0,
        }

    # ==================================================
    def compute_dos(self):
        """
        Compute DOS.

        :meta private:
        """
        if not self.output.get("dos", False):
            return

        print("compute and output dos.")

    # ==================================================
    def get_kpath(self, k_path):
        """
        Get k path.

        Args:
            k_path (str): k path.
        Returns:
            - (dict) -- k point dict.
            - (str) -- k path.

        :meta private:
        """
        if k_path == "":  # create default path.
            A = self["A"]
            gp = next(reversed(self.model.group.wyckoff["site"].values()))  # general point.
            positions = gp["reference"].astype(float)  # fractional, conventional, plus set.
            numbers = np.full(len(positions), 1, dtype=int)

            structure = (A, positions, numbers)
            info = seekpath.get_path(structure)

            if info["spacegroup_number"] != int(self.model.group.ID):
                logging.exception("obtained SG is different with given group.")
                raise

            k_point = info["point_coords"]
            k_point["Γ"] = k_point["GAMMA"]
            del k_point["GAMMA"]

            path = info["path"]
            k_path = path[0][0] + "-" + path[0][1]
            for (a, b), (c, d) in zip(path, path[1:]):
                if b == c:
                    k_path += "-" + d
                else:
                    k_path += "|" + c + "-" + d
            k_path = k_path.replace("GAMMA", "Γ")
        else:
            k_point = self.output["dispersion"].get("k_point", {})
            k_point = {k: str_to_sympy(v).astype(float) for k, v, in k_point.items()}

        return k_point, k_path

    # ==================================================
    def set_k_multipole(self, matrix_info):
        """
        Set momentum multipole.

        Args:
            matrix_info (dict): matrix info.

        Returns:
            - (dict) -- k-multipole dict.

        Notes:
            - only tight-binding gauge is supported.

        :meta private:
        """
        if not self.samb.get("k_multipole", False):
            return {}

        combined_id = self.model["combined_id"]

        k_multipole, cluster_vec = create_k_multipole(self.model["cluster_samb"], self.model["cluster_vector"])
        cluster_vec = {sb: {str(kb): str(v).replace(" ", "") for kb, v in lst.items()} for sb, lst in cluster_vec.items()}
        k_matrix = create_k_matrix(matrix_info["matrix"], matrix_info["cluster"], matrix_info["vector"])
        k_matrix = {
            tag: (matrix_info["cluster"][tag], combined_id[tag][1].samb_type.wyckoff, mat) for tag, mat in k_matrix.items()
        }

        # convert to str for output.
        k_multipole = {
            wp: {idx: (str(samb.tolist()).replace(" ", ""), str(sym.tolist()).replace(" ", "")) for idx, (samb, sym) in v.items()}
            for wp, v in k_multipole.items()
        }
        k_matrix = {
            tag: (cn, wp, {Rmn: str(v).replace(" ", "") for Rmn, v in mat.items()}) for tag, (cn, wp, mat) in k_matrix.items()
        }

        k_multipole = {
            "model": matrix_info["model"],
            "source": matrix_info["source"],
            "date": matrix_info["date"],
            "dimension": matrix_info["dimension"],
            "ket_site": list(matrix_info["ket_site"].keys()),
            "index": matrix_info["index"],
            "cluster_vector": cluster_vec,
            "k_multipole": k_multipole,
            "k_matrix": k_matrix,
        }

        return k_multipole
