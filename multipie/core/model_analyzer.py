"""
Model analyzer class.

This module provides model analyzer.
"""

import os
import numpy as np
import seekpath
from multipie.core.material_model import MaterialModel
from multipie.core.default_control import default_control
from multipie.util.util_model_analyzer import (
    fourier_r_to_k,
    grid_path,
    output_linear_dispersion_eig,
    create_all_local_operator,
    create_local_operator,
    create_k_multipole,
    create_k_matrix,
)
from multipie.util.util import read_dict, str_to_sympy, write_dict

_k_matrix_comment = """Selected SAMB matrix in momentum representation.
- dimension (int): matrix size.
- ket_site (list): ket info., [ket_name].
- index (dict): ket index, dict[(site,sublattice,rank), (top_index,size)].
- cluster (dict): cluster name, dict[SAMB ID, cluster name].
- cluster_vector (dict): cluster vector, dict[site/bond name, dict[kb, expression] ].
- k_multipole (dict): momentum multipole in terms of k.b_n, dict[wyckoff, dict[idx, (k_multipole, symmetry)] ].
- k_matrix (dict): momentum matrix, dict[tag, dict[(n1,n2,n3,m,n), value] ].
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

        self._samb = default_control["samb"]
        self._wannier = default_control["wannier"]
        self._output = default_control["output"]

        self._topdir = topdir
        self._verbose = verbose
        self._HR = None
        self._name = ""
        self._mm = MaterialModel(topdir, verbose=verbose)
        self._local = create_all_local_operator()
        self.set_grid_size(N1, N2, N3)

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
    def name(self):
        """
        Model name.

        Returns:
            - (str) -- model name.
        """
        return self._name

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
        self["mp_grid"] = [N1, N2, N3]

    # ==================================================
    def analyze(self, control):
        """
        Analyze model with control file.

        Args:
            control (str or dict): control file (.py) or model name.
        """
        if type(control) == str:  # read control file or from dict.
            if control.endswith(".py"):
                file = os.path.join(self._topdir, control)
                control = read_dict(file)
            else:  # w/o control and model_name is given.
                self.model.load(control)
                matrix_info = self.model.get_samb_matrix({})
                self.model.save_samb_matrix(matrix_info)
                return

        self._samb |= control.get("samb", {})
        self._wannier |= control.get("wannier", {})
        self._output |= control.get("output", {})

        # exec. SAMB control.
        if self.samb.get("model", None) is not None:
            self._name = self.samb["model"]
            self.set_samb()

        # exec. wannier control.
        if self._wannier.get("cw", None) is not None:
            # self._name = self.wannier["model"]
            self.set_from_wannier()

        # compute physical quanties and output data.
        self.compute_physical_quantity()

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
    def local_operator(self, name):
        """
        Create local operator.

        Args:
            name (str): operator name, "Sx/Sy/Sz/Lx/Ly/Lz/Qu/Qv/Qyz/Qzx/Qxy".

        Returns:
            - (ndarray) -- operator matrix (dim x dim).

        :meta private:
        """
        ket = self.model["full_matrix"]["ket"]
        basis_type = self.model["basis_type"]
        return create_local_operator(ket, name, self._local, basis_type == "lgs")

    # ==================================================
    def set_samb(self):
        """
        Calculate SAMB related quantities.

        :meta private:
        """
        self.model.load(self.name)
        self.set_primitive_cell(np.array(self.model["unit_vector_primitive"]))
        select = self.samb.get("select", {})

        parameter = self.samb.get("parameter", {})
        if type(parameter) == str:  # when parameter is str, which means filename of z_j dict.
            z_file = os.path.join(self._topdir, self.name, parameter)
            parameter = read_dict(z_file)
            self._samb["parameter"] = parameter

        if self.samb.get("samb_figure", False):
            self.model.save_samb_qtdraw()

        matrix_info = self.model.get_samb_matrix(select)
        if parameter:
            self._HR = self.model.get_hr(parameter, matrix_info["matrix"])
            self.model.save_samb_hr(matrix_info, parameter, self._HR)
        self.model.save_samb_matrix(matrix_info)

        self.set_k_multipole(matrix_info)

    # ==================================================
    def set_from_wannier(self):
        """
        Set data for wannier-based input.

        :meta private:
        """
        pass

    # ==================================================
    def compute_physical_quantity(self):
        """
        Compute physical quantities by parsing the control file.

        :meta private:
        """
        if self._HR is None:
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
    def compute_dispersion(self):
        """
        Compute dispersion.

        :meta private:
        """
        if "dispersion" not in self.output:
            return

        k_path = self.output["dispersion"].get("k_path", None)
        if k_path is None or self.model.group.group_type != "SG":
            return

        k_point, k_path = self.get_kpath(k_path)
        k_point_path, k_linear, k_dis_pos = grid_path(k_point, k_path, self["mp_grid"][0], self["B"])
        self._output["dispersion"]["k_path"] = k_path
        self._output["dispersion"]["k_point"] = k_point

        tb_gauge = self.output["fourier"]["tb_gauge"]
        atom = np.asarray(list(self.model.get_ket_site().values()), dtype=float)
        basis_type = self.model["basis_type"]
        if basis_type == "jml":
            op_lst = []
        else:
            op_lst = self.output["dispersion"].get("local", [])

        HR = {((n1, n2, n3), m, n): complex(v) for (n1, n2, n3, m, n), v in self._HR.items()}
        Hk = fourier_r_to_k(HR, atom, k_point_path, tb_gauge)

        Ek, Uk = np.linalg.eigh(Hk)
        Ok = [np.einsum("kmi,mn,kni->ki", Uk.conj(), self.local_operator(name), Uk).real for name in op_lst]
        fname = self.name + "_dispersion.txt"
        if Ok:
            output_linear_dispersion_eig(fname, k_linear, Ek, Ok, k_dis_pos=k_dis_pos, colormap=True)
        else:
            output_linear_dispersion_eig(fname, k_linear, Ek, k_dis_pos=k_dis_pos)

        if self._verbose:
            outdir = os.path.join(self._topdir, self.name, self.output["dir"])
            print(f"save dispersion to '{outdir}/{fname}'.")

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
            A = self.model["unit_vector_primitive"]
            d = next(reversed(self.model.group.wyckoff["site"].values()))  # general point.
            positions = d["reference"].astype(float)  # fractional, conventional, plus set.
            numbers = np.full(len(positions), 1, dtype=int)

            structure = (A, positions, numbers)
            info = seekpath.get_path(structure)

            if info["spacegroup_number"] != int(self.model.group.ID):
                print("obtained SG is different with given group.")
                raise Exception

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

        Notes:
            - only tight-binding gauge is supported.

        :meta private:
        """
        if not self.samb.get("k_multipole", False):
            return

        k_multipole, cluster_vec = create_k_multipole(self.model["cluster_samb"], self.model["cluster_vector"])
        cluster_vec = {sb: {str(kb): str(v).replace(" ", "") for kb, v in lst.items()} for sb, lst in cluster_vec.items()}
        k_matrix = create_k_matrix(matrix_info["matrix"], matrix_info["cluster"], matrix_info["vector"])

        # convert to str for output.
        k_multipole = {
            wp: {idx: (str(samb.tolist()).replace(" ", ""), str(sym.tolist()).replace(" ", "")) for idx, (samb, sym) in v.items()}
            for wp, v in k_multipole.items()
        }
        k_matrix = {tag: {Rmn: str(v).replace(" ", "") for Rmn, v in mat.items()} for tag, mat in k_matrix.items()}

        # output.
        outdir = os.path.join(self._topdir, self.name, self.output["dir"])
        fname = self.name + "_k.py"
        write_dict(
            {
                "dimension": matrix_info["dimension"],
                "ket_site": list(matrix_info["ket_site"].keys()),
                "index": matrix_info["index"],
                "cluster": matrix_info["cluster"],
                "cluster_vector": cluster_vec,
                "k_multipole": k_multipole,
                "k_matrix": k_matrix,
            },
            fname,
            comment=_k_matrix_comment,
            w_dir=outdir,
        )

        if self._verbose:
            print(f"save k-multipole to '{outdir}/{fname}'.")
