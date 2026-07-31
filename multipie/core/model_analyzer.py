"""
Model analyzer class.

This module provides model analyzer.
"""

import os
import logging
import numpy as np
import copy
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
    add_local_parameter_sym,
    convert_zj_atomic_var,
    fermi_dirac,
    convert_orbital_to_detail,
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
    create_ket_wannier_multipie,
)
from multipie.util.util import read_dict, str_to_sympy, write_dict, deep_update

_matrix_comment = """Selected SAMB matrix.
- model (str): model name.
- source (str): source binary.
- created (str): binary created date.
- select (dict): select condition used.
- dimension (int): matrix size.
- ket_site (dict): ket info., dict[ket_name, position (fractional, primitive)].
- index (dict): ket index, dict[(site,sublattice,rank), (top_index,size)].
- vector (dict): primitive bond vector, dict[cluster name, [primitive bond vector]].
- cluster (dict): cluster name, dict[SAMB ID, cluster name].
- matrix (dict): matrix, dict[zi, dict[(R,row,column), (value, bond_no)] ] (R=n1,n2,n3, primitive).
"""

_k_matrix_comment = """Selected SAMB matrix in momentum representation.
- model (str): model name.
- source (str): source binary.
- created (str): binary created date.
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
- atomic variable of (m,n) component at bond 1 is given by g(m,n) +i h(m,n).
"""

_param_comment = """Parameter dict (sorted by descending absolute value).
- finite parameter, dict[zj, value].
"""


# ==================================================
class ModelAnalyzer(dict):
    # ==================================================
    def __init__(self, topdir=None, verbose=False):
        """
        Model analyzer.

        Args:
            topdir (str, optional): top directory. [default: cwd]
            verbose (bool, optional): verbose comment ?
        """
        if topdir is None:
            topdir = os.getcwd()

        self._topdir = topdir
        self._verbose = verbose
        self._mm = MaterialModel(topdir, verbose=verbose)
        self._local = create_all_local_operator()
        os.chdir(self._topdir)

        self.reset()

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
            - (list) -- basis, list of "orbital@atom(sublattice)".
        """
        return self._basis

    # ==================================================
    @property
    def HR(self):
        """
        Real-space hamiltonian, H(R).

        Returns:
            - (dict) -- H(R), dict[(n1,n2,n3,m,n), val or (val, bond_no)].
        """
        return self._HR

    # ==================================================
    def write_dict(self, dic, ext, comment="", w_dir=None):
        """
        Write dict.

        Args:
            dic (dict): dict to write.
            ext (str): filename is name_ + ext.
            comment (str, optional): comment.
            w_dir (str, optional): directory (topdir/name/w_dir) to write, if None, topdir/name is used.
        """
        name = self["info"]["name"]
        path = os.path.join(self._topdir, name)
        if w_dir is not None:
            path = os.path.join(path, w_dir)
        os.makedirs(path, exist_ok=True)
        filename = name + "_" + ext
        write_dict(dic, filename, comment=comment, w_dir=path)

        if self._verbose:
            ext = ext[: ext.rfind(".")]
            print(f"save {ext} to '{path}/{filename}'.")

    # ==================================================
    def reset(self, control=None):
        if control is None:
            control = {}

        self["info"] = {}
        self["samb"] = {}
        self["wannier"] = {}
        self["output"] = {}

        self._samb = copy.deepcopy(default_control["samb"])
        deep_update(self._samb, control.get("samb", {}))
        self._wannier = copy.deepcopy(default_control["wannier"])
        deep_update(self._wannier, control.get("wannier", {}))
        self._output = copy.deepcopy(default_control["output"])
        deep_update(self._output, control.get("output", {}))

        self.set_mode(control.get("mode", default_control["mode"]))
        self.set_grid(*control.get("grid", default_control["grid"]))
        self.set_name(self.samb["model"])
        self.set_parameter(None)
        self.set_basis_type(None)
        self.set_basis(None)
        self.set_HR(None)

    # ==================================================
    def set_mode(self, mode):
        """
        Set analysis mode.

        Args:
            mode (str): analysis mode, "samb/wannier/symcw".
        """
        self["info"]["mode"] = mode

    # ==================================================
    def set_name(self, name):
        """
        Set model name.

        Args:
            name (str): model name.
        """
        self["info"]["name"] = name

    # ==================================================
    def set_grid(self, N1, N2, N3):
        """
        Set k-grid size.

        Args:
            N1 (int): number of divisions in b1.
            N2 (int): number of divisions in b2.
            N3 (int): number of divisions in b3.
        """
        self["info"]["grid"] = [N1, N2, N3]

    # ==================================================
    def set_parameter(self, parameter):
        """
        Set parameter, zj.

        Args:
            parameter (dict): parameter dict.
        """
        self._parameter = parameter

    # ==================================================
    def set_basis_type(self, basis_type):
        """
        Set atomic basis type.

        Args:
            basis_type (str): basis type, "lg/lgs/jml".
        """
        self._basis_type = basis_type

    # ==================================================
    def set_basis(self, basis):
        """
        Set full-matrix basis.

        Args:
            basis (list): basis, list of "orbital@atom(sublattice)".
        """
        self._basis = basis

    # ==================================================
    def set_HR(self, HR):
        """
        Set real-space hamiltonian, H(R).

        Args:
            HR (dict): H(R), dict[(n1,n2,n3,m,n), val or (val, bond_no)].
        """
        self._HR = HR

    # ==================================================
    def set_primitive_cell(self, A):
        """
        Set primitive cell info., A, B, volume.

        Args:
            A (ndarray): translational vectors of primitive cell, [a1, a2, a3] (3x3).

        :meta private:
        """
        A = np.asarray(A, dtype=float)
        B = 2 * np.pi * np.linalg.inv(A).T
        self["info"]["A"] = A.tolist()  # primitive cell.
        self["info"]["B"] = B.tolist()  # reciprocal cell.
        self["info"]["volume"] = float(np.dot(A[0], np.cross(A[1], A[2])))  # volume of primitive cell.

    # ==================================================
    def read_controle(self, control):
        if control.endswith(".py"):
            control = read_dict(control)
            return False
        else:
            return True

    # ==================================================
    def write_samb_matrix(self, matrix_info):
        """
        Save SAMB matrix.

        Args:
            matrix_info (dict): matrix info.
        """
        # convert sympy to str.
        mi = matrix_info.copy()
        matrix = matrix_info["matrix"]
        mi["matrix"] = {z: {k: (str(v[0]).replace(" ", ""), v[1]) for k, v in elm.items()} for z, elm in matrix.items()}

        self.write_dict(mi, "matrix.py", _matrix_comment, "info")

    # ==================================================
    def write_samb_hr(self, matrix_info, parameter):
        """
        Save SAMB matrix in hr format.

        Args:
            matrix_info (dict): matrix info.
            parameter (dict): parameter dict, dict[z#, value].
            HR (dict, optional): H(R) matrix. if None, H(R) is generated.
            verbose (bool, optional): verboser ?
        """
        if self.HR is None:
            return
        HR = self.HR

        # write hr.
        name = self["info"]["name"]
        filename = os.path.join(self._topdir, name, "info", name + "_hr.dat")
        with open(filename, mode="w", encoding="utf-8") as f:
            print(f"# SAMB matrix from {matrix_info["source"]} ({matrix_info["created"]})", file=f)
            print("# select", file=f)
            for k, v in matrix_info["select"].items():
                print(f"#   {k}: {str(v).replace(" ", "")}", file=f)
            print(f"# basis ({matrix_info["dimension"]})", file=f)
            for no, (b, p) in enumerate(matrix_info["ket_site"].items()):
                print(f"#   {no:2d} {b}: [{p[0]: .6f}, {p[1]: .6f}, {p[2]: .6f}]", file=f)
            for z, v in parameter.items():
                print(f"# {z:<4} = {v}", file=f)
            print("#", file=f)
            print("# n1   n2   n3    m    n    re                        im", file=f)
            for (n1, n2, n3, m, n), v in HR.items():
                v = complex(v)
                r, i = v.real, v.imag
                s = f"{n1: 4d} {n2: 4d} {n3: 4d} {m: 4d} {n: 4d}    {r: .15e}    {i: .15e}"
                print(s, file=f)

        if self._verbose:
            print(f"save hr to '{filename}'.")

    # ==================================================
    def write_info(self):
        if "samb" in self.keys():
            dic = {k: v for k, v in self["samb"].items() if k not in ["matrix_info", "var", "k_multipole"]}
            self.write_dict(dic, "info_samb.py", w_dir="info")
        if "wannier" in self.keys():
            dic = {k: v for k, v in self["wannier"].items() if k not in ["kpoints", "nnkpts", "wk", "bveck", "kb2k"]}
            self.write_dict(dic, "info_wannier.py", w_dir="info")
        if "output" in self.keys():
            dic = {k: v for k, v in self["output"].items() if k not in []}
            self.write_dict(dic, "info_output.py", w_dir="info")
        dic = {k: v for k, v in self.items() if k not in ["samb", "wannier", "output"]}
        self.write_dict(dic, "info.py", w_dir="info")

    # ==================================================
    def write_z_file(self):
        parameter = self.parameter
        if not parameter:
            return

        dic = {tag: float(v) for tag, v in parameter.items()}
        dic = dict(sorted(dic.items(), key=lambda item: abs(item[1]), reverse=True))

        comment = _param_comment + f"- by using '{self["info"]["mode"]}' mode.\n"
        self.write_dict(dic, "z.py", comment, "info")

    # ==================================================
    def set_var(self, matrix_info):
        IR = next(iter(self.model.group.character["table"].keys()))  # identity irrep.
        conv_dict = convert_zj_atomic_var(matrix_info, self.model["combined_cluster"], self.model["combined_id"], IR)
        return conv_dict

    # ==================================================
    def write_var_file(self):
        d = self["samb"]["var"]
        d = {name: {zj: str(ex).replace(" ", "") for zj, ex in dic.items()} for name, dic in d.items()}

        self.write_dict(d, "var.py", _zj_var_comment, "info")

    # ==================================================
    def write_k_multipole(self):
        d = self["samb"]["k_multipole"]
        if d:
            self.write_dict(d, "k.py", _k_matrix_comment, "info")

    # ==================================================
    def analyze(self, control):
        """
        Analyze model with control file.

        Args:
            control (str or dict): control file (.py) or model name.
        """
        # read control.
        if type(control) == str:  # read control file.
            if self.read_controle(control):  # if error.
                return

        self.reset(control)
        mode = self["info"]["mode"]

        # set SAMB.
        if mode in ["samb", "symcw"]:
            self.set_samb()  # create SAMBs, and H(R) if zj are provided.

        # set wannier.
        if mode in ["wannier", "symcw"]:
            HR = self.set_wannier()  # create H(R), and zj in case of "symcw".
            self.set_HR(HR)
            if mode == "symcw":
                matrix_info = self["samb"]["matrix_info"]
                Zr_dict = matrix_info["matrix"]
                parameter = decompose_operator_by_SAMB(HR, Zr_dict)
                HR = self.model.get_hr(parameter, Zr_dict)
                self.set_HR(HR)
                self.write_samb_hr(matrix_info, parameter)  # overwrite HR by MultiPie.
                self.set_parameter(parameter)

        # create z file.
        self.write_z_file()

        # compute physical quanties and output data.
        self.compute_physical_quantity()

        # output info.
        self.write_info()

    # ==================================================
    def read_parameter(self, filename=None):
        name = self["info"]["name"]
        if filename is None:
            filename = f"info/{name}_z.py"
        filename = os.path.join(self._topdir, name, filename)
        parameter = read_dict(filename)
        parameter = {tag: float(str_to_sympy(v, rational=False)) if type(v) == str else v for tag, v in parameter.items()}
        if self._verbose:
            print(f"load parameter from '{filename}'.")

    # ==================================================
    def set_samb(self):
        """
        Calculate SAMB related quantities.

        :meta private:
        """
        name = self.samb["model"]
        self.set_name(name)

        # read and set model.
        self.model.load(name)
        self.set_basis_type(self.model["basis_type"])
        self.set_basis(self.model["full_matrix"]["ket"])
        self.set_primitive_cell(self.model["unit_vector_primitive"])

        # set selected SAMBs.
        self["samb"]["matrix_info"] = self.model.get_samb_matrix(self.samb["select"])

        # create var file.
        self["samb"]["var"] = self.set_var(self["samb"]["matrix_info"])
        self.write_var_file()

        # create k-multipole file.
        if self.samb["k_multipole"]:
            self["samb"]["k_multipole"] = self.set_k_multipole(self["samb"]["matrix_info"])
            self.write_k_multipole()

        # create SAMB qtdraw.
        if self.samb["samb_figure"]:
            self.model.save_samb_qtdraw()

        parameter = self.samb["parameter"]
        if type(parameter) == str:  # when parameter is str, read z file.
            if parameter:
                z_file = "info/" + parameter
            else:
                z_file = None  # user default.
            parameter = self.read_parameter(z_file)

        # determine local weight if NG_sum_rule is True.
        ng = self.samb["NG_sum_rule"]
        if ng and parameter:
            parameter = add_local_parameter(self["samb"]["matrix_info"], parameter, self.model["full_matrix"]["ket"])
        if ng:
            parameter_sym = add_local_parameter_sym(self["samb"]["matrix_info"], self.model["full_matrix"]["ket"])
            self["samb"]["NG_sum_rule"] = parameter_sym
        else:
            self["samb"]["NG_sum_rule"] = None

        # output matrix.py and hr.dat.
        if parameter:
            self.set_HR(self.model.get_hr(parameter, self["samb"]["matrix_info"]["matrix"]))
            self.write_samb_hr(self["samb"]["matrix_info"], parameter)
        self.write_samb_matrix(self["samb"]["matrix_info"])

        self.set_parameter(parameter)
        self["info"]["fermi_energy"] = 0.0

    # ==================================================
    def set_wannier(self):
        """
        Set data for wannier-based input.

        :meta private:
        """
        seedname = self.wannier["seedname"]
        wannier_dir = os.path.join(self._topdir, seedname, self.wannier["dir"])

        # read seedname.win
        win = read_win(seedname, wannier_dir)
        # read seedname.nnkp
        nnkp = read_nnkp(seedname, wannier_dir)
        # check common values and merge.
        wannier_info = merge_wannier_info(win, nnkp, seedname)

        wannier_ket_info = {
            "A": win["A"],
            "atoms_frac": win["atoms_frac"],
            "nw2n": nnkp["nw2n"],
            "nw2l": nnkp["nw2l"],
            "nw2m": nnkp["nw2m"],
            "nw2r": nnkp["nw2r"],
            "nw2s": nnkp["nw2s"],
        }
        wannier2multipie, ket_wannier, ket_multipie = create_ket_wannier_multipie(wannier_info)
        ket = self.wannier.get("ket_wannier", [])
        if ket:
            wannier2multipie = [ket.index(i) for i in ket_multipie]

        nw2n = nnkp["nw2n"]
        atoms_list = list(wannier_info["atoms_frac"].values())
        atoms_frac = [atoms_list[nw2n[i]] for i in wannier2multipie]
        atoms_list = list(wannier_info["atoms_cart"].values())
        atoms_cart = [atoms_list[nw2n[i]] for i in wannier2multipie]

        if self.wannier["read_KS"]:
            # read KS Ek and Uk, and convert to MultiPie standard order(*) of ket by changing indices of Uk.
            # create H(R), and various matrix elements in real space by Fourier transformation with DFT k-grid.
            # (*) use wannier2multipie.
            #
            # read seedname.mmn
            # Mkb = read_mmn(seedname, wannier_dir)
            # read seedname.spn
            # Sk = read_spn(seedname, wannier_dir)
            # read seedname.uHu
            # uHu = read_uHu(seedname, wannier_dir)
            # read seedname.uIu
            # uIu = read_uIu(seedname, wannier_dir)
            pass
        else:
            hr_file = seedname + "_hr.dat"
            # read seedname_hr.dat
            hr_dict, irvec, ndegen = read_hr(hr_file, wannier_dir)
            idx = {i: no for no, i in enumerate(wannier2multipie)}
            hr_dict = {(n1, n2, n3, idx[a], idx[b]): complex(v) for (n1, n2, n3, a, b), v in hr_dict.items()}

        info = {
            # Wannier info.
            "A": wannier_info["A"],
            "B": wannier_info["B"],
            "ket_wannier": ket_wannier,
            "ket_multipie": ket_multipie,
            "wannier_to_multipie": wannier2multipie,  # multipie = wannier[wannier2multipie].
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

        self["wannier"] = info
        self["info"]["fermi_energy"] = info["fermi_energy"]

        ### physical qunatity.
        # nk = np.array([np.diag(fermi_dirac(eki - win["fermi_energy"], T=0.0)) for eki in Ek], dtype=float)
        # nk = Uk.transpose(0, 2, 1).conjugate() @ nk @ Uk
        # nr_dict = fourier_k_to_r(nk, win["kpoints"], irvec, s=False)
        # nr_dict = sort_ket_matrix_dict(nr_dict, ket_wannier, ket_multipie)
        # z_j_exp = decompose_operator_by_SAMB(nr_dict, Zr_dict)
        # self["wannier"]["z_j_exp"] = z_j_exp
        # self["wannier"]["mmn"] = mmn
        # self["wannier"]["spn"] = spn
        # self["wannier"]["uHu"] = uHu
        # self["wannier"]["uIu"] = uIu

        return hr_dict

    # ==================================================
    def compute_physical_quantity(self):
        """
        Compute physical quantities by parsing the control file.

        Args:
            name (str): model name.

        :meta private:
        """
        if self.HR is None and self._verbose:
            print("set H(R) first before calculating physical quantities.")
            return

        name = self["info"]["name"]
        path = os.path.join(self._topdir, name, self.output["dir"])
        os.makedirs(path, exist_ok=True)

        self.set_eigen_system()
        self.compute_dispersion()
        self.compute_dos()

    # ==================================================
    def set_eigen_system(self):
        """
        Set eigen system by checking control/output if E and/or U is required.

        :meta private:
        """
        pass

    # ==================================================
    def local_operator(self, tag):
        """
        Create local operator.

        Args:
            tag (str): operator name, "Sx/Sy/Sz/Lx/Ly/Lz/Qu/Qv/Qyz/Qzx/Qxy".

        Returns:
            - (ndarray) -- operator matrix (dim x dim).

        :meta private:
        """
        spinful = self.basis_type == "lgs"
        return create_local_operator(self.basis, tag, self._local, spinful)

    # ==================================================
    def write_dispersion(self, Ek, Ok, op_lst, k_linear, k_dis_pos):
        name = self["info"]["name"]
        path = os.path.join(self._topdir, name, self.output["dir"])
        fname = os.path.join(path, name + "_dispersion.txt")
        colormap = len(Ok) > 0
        ef = self["info"]["fermi_energy"]
        if Ok:
            output_dispersion(fname, k_linear, ef, Ek, Ok, op_lst)
        else:
            output_dispersion(fname, k_linear, ef, Ek)
        plot_save_dispersion(fname, k_dis_pos, ef, colormap)
        create_gnuplot_cmd(fname, k_dis_pos, np.max(k_linear), np.max(Ek), np.min(Ek), ef, colormap)
        if self._verbose:
            print(f"save dispersion files into '{path}'.")

    # ==================================================
    def compute_dispersion(self):
        """
        Compute dispersion.

        :meta private:
        """
        # check if dispersion can be computed.
        if "dispersion" not in self.output:
            return
        k_path = self.output["dispersion"]["k_path"]
        if k_path is None or self.model.group.group_type != "SG":
            return

        # get k_point and k_path.
        k_point, k_path = self.get_kpath(k_path)
        N1 = self["info"]["grid"][0]
        B = np.asarray(self["info"]["B"])
        k_point_path, k_linear, k_dis_pos = grid_path(k_point, k_path, N1, B)

        # get local operator list.
        if self.basis_type == "jml" or self.basis_type is None:
            op_lst = []
        else:
            op_lst = self.output["dispersion"]["local"]

        # get info.
        tb_gauge = self.output["fourier"]["tb_gauge"]
        atom = np.asarray(list(self.model.get_ket_site().values()), dtype=float)

        # set H(R) and H(k).
        HR = {((n1, n2, n3), m, n): complex(v) for (n1, n2, n3, m, n), v in self.HR.items()}
        Hk = fourier_r_to_k(HR, atom, k_point_path, tb_gauge)

        # set eigen system.
        Ek, Uk = np.linalg.eigh(Hk)
        power = self.output["dispersion"]["power"]
        if power is not None:
            Ek = np.power(Ek, power)

        # set local operators.
        Ok = [np.einsum("kmi,mn,kni->ki", Uk.conj(), self.local_operator(tag), Uk).real for tag in op_lst]

        # output dispersion data, plot, and gnuplot.
        self.write_dispersion(Ek, Ok, op_lst, k_linear, k_dis_pos)

        # save dispersion info.
        self["output"]["dispersion"] = {
            "k_path": k_path,
            "k_point": {k: str(v.tolist()).replace(" ", "") for k, v in k_point.items()},
            "e_max": float(np.max(Ek)),
            "e_min": float(np.min(Ek)),
        }

    # ==================================================
    def compute_dos(self):
        """
        Compute DOS.

        :meta private:
        """
        if not self.output["dos"]:
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
        if not self.samb["k_multipole"]:
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
            "created": matrix_info["created"],
            "dimension": matrix_info["dimension"],
            "ket_site": list(matrix_info["ket_site"].keys()),
            "index": matrix_info["index"],
            "cluster_vector": cluster_vec,
            "k_multipole": k_multipole,
            "k_matrix": k_matrix,
        }

        return k_multipole
