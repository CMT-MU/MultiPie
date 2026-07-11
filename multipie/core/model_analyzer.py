"""
Model analyzer class.

This module provides model analyzer.
"""

import os
import numpy as np
from multipie.core.material_model import MaterialModel
from multipie.core.default_control import default_control
from multipie.util.util_model_analyzer import (
    fourier_r_to_k,
    grid_path,
    output_linear_dispersion_eig,
    create_all_local_operator,
    create_local_operator,
)
from multipie.util.util import read_dict, str_to_sympy


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
    def exec(self, control):
        """
        Execute analyzer with control file.

        Args:
            control (str or dict): control file (.py) or model name.
        """
        if type(control) == str:  # control file.
            if control.endswith(".py"):
                file = os.path.join(self._topdir, control)
                control = read_dict(file)
            else:  # model name.
                self.model.load(control)
                self.model.save_samb_matrix({})
                return

        self._samb |= control.get("samb", {})
        self._wannier |= control.get("wannier", {})
        self._output |= control.get("output", {})

        if self.samb.get("model", None) is not None:
            self._name = self.samb["model"]
            self.model.load(self.name)
            self.set_primitive_cell(np.array(self.model["unit_vector_primitive"]))
            select = self.samb.get("select", {})
            parameter = self.samb.get("parameter", {})
            if self.samb.get("samb_figure", False):
                self.model.save_samb_qtdraw()
            if parameter:
                md = self.model.save_samb_matrix(select)
                self._HR = self.model.save_samb_hr(md, parameter)
            else:
                self.model.save_samb_matrix(select)
                return

        if self._wannier.get("cw", None) is not None:
            # self._name = self.wannier["model"]
            self.set_from_wannier()

        self.compute_physical_quantity()

    # ==================================================
    def set_primitive_cell(self, A):
        """
        Set primitive cell info.

        Args:
            A (ndarray): [a1, a2, a3] (3x3).
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
        """
        ket = self.model["full_matrix"]["ket"]
        basis_type = self.model["basis_type"]
        return create_local_operator(ket, name, self._local, basis_type == "lgs")

    # ==================================================
    def set_from_wannier(self):
        pass

    # ==================================================
    def compute_physical_quantity(self):
        if self._HR is None:
            print("set H(R) first before calculating physical quantities.")
            return

        cwd = os.getcwd()

        outdir = os.path.join(self._topdir, self.name, self.output["dir"])
        os.makedirs(outdir, exist_ok=True)
        os.chdir(outdir)

        self.set_eigen_system()

        if "dispersion" in self.output:
            self.compute_dispersion()

        if "dos" in self.output:
            self.compute_dos()

        os.chdir(cwd)

    # ==================================================
    def set_eigen_system(self):
        # check output to determine for calculating E only or both E, U.
        # and then compute E (and U).
        pass

    # ==================================================
    def compute_dispersion(self):
        print("compute and output dispersion.")
        tb_gauge = self.output["fourier"]["tb_gauge"]
        k_point = self.output["dispersion"]["k_point"]
        k_path = self.output["dispersion"]["k_path"]
        k_point = {k: str_to_sympy(v).astype(float) for k, v, in k_point.items()}
        k_point_path, k_linear, k_dis_pos = grid_path(k_point, k_path, self["mp_grid"][0], self["B"])
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
        if Ok:
            output_linear_dispersion_eig(".", self.name + "_dispersion", k_linear, Ek, Ok, k_dis_pos=k_dis_pos, colormap=True)
        else:
            output_linear_dispersion_eig(".", self.name + "_dispersion", k_linear, Ek, k_dis_pos=k_dis_pos)

    # ==================================================
    def compute_dos(self):
        print("compute and output dos.")
