"""
Model analyzer class.

This module provides model analyzer.
"""

import os
import re
import logging
import importlib.util
import numpy as np
import copy
import functools
import seekpath
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import gaussian_filter1d
from multipie.core.material_model import MaterialModel
from multipie.util.util_crystal import convert_to_primitive
from multipie.core.default_control import default_control
from multipie.util.util_model_analyzer import (
    grid_path,
    grid_all,
    fourier_r_to_k,
    output_dispersion,
    create_gnuplot_cmd,
    plot_save_dispersion,
    plot_save_weight,
    create_gnuplot_weight,
    output_dos,
    plot_save_dos,
    create_gnuplot_dos,
    create_all_local_operator,
    create_local_operator,
    create_k_multipole,
    create_k_matrix,
    add_local_parameter,
    add_local_parameter_sym,
    convert_zj_atomic_var,
)
from multipie.util.util_wannier import (
    read_win,
    read_nnkp,
    merge_wannier_info,
    read_hr,
    read_wsvec,
    apply_ws_degeneracy,
    decompose_operator_by_SAMB,
    create_ket_wannier_multipie,
    map_wannier_to_model,
    convert_hr_to_model,
    model_primitive_vector,
)
from multipie.util.util import read_dict, str_to_sympy, write_dict, deep_update
from multipie.util.util import check_dict_keys

logger = logging.getLogger(__name__)

_matrix_comment = """Selected SAMB matrix.
- model (str): model name.
- source (str): source binary.
- created (str): binary created date.
- select (dict): select condition used.
- dimension (int): matrix size.
- ket (list): ket name of full matrix, [name].
- ket_pos (list): ket position (fractional, primitive), [pos].
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
- ket (list): ket name of full matrix, [name].
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
def import_libtetrabz():
    """
    Import libtetrabz (optional dependency for the tetrahedron method).

    Returns:
        - (module) -- libtetrabz.
    """
    try:
        import libtetrabz
    except ImportError as e:
        raise ImportError(
            "libtetrabz is required for 'output/dos_method' = 'tetrahedron'. install it by 'pip install libtetrabz' "
            "(or 'pip install multipie[tetrahedron]')."
        ) from e

    return libtetrabz


# ==================================================
def _k_chunk(dim, budget=64e6):
    """
    Number of k points in a chunk, so that H(k) and its eigen system for the chunk fit in the memory budget.

    Args:
        dim (int): matrix size of H(k).
        budget (float, optional): memory budget in bytes.

    Returns:
        - (int) -- number of k points.
    """
    return max(1, int(budget // (48 * dim * dim)))  # complex H(k), U(k), and work space.


# ==================================================
def _join_kpath(path):
    """
    Join k path segments of seekpath.

    Args:
        path (list): list of (start, end).

    Returns:
        - (str) -- k path, e.g., "Γ-X-M|X-R".
    """
    k_path = path[0][0] + "-" + path[0][1]
    for (a, b), (c, d) in zip(path, path[1:]):
        if b == c:
            k_path += "-" + d
        else:
            k_path += "|" + c + "-" + d

    return k_path.replace("GAMMA", "Γ")


# ==================================================
def _kpath_structure(group, A):
    """
    Structure for seekpath, which has the symmetry of the group in its primitive cell.

    Args:
        group (Group): space group.
        A (array-like): primitive lattice vectors, [a1,a2,a3] (Cartesian, rows).

    Returns:
        - (tuple) -- (A, positions (fractional), numbers) in the primitive cell.

    Note:
        - two orbits of generic points with different species are used, since a single orbit can have higher symmetry than the group, e.g., inversion for a polar group.
    """
    positions, numbers = [], []
    for no, pos in enumerate(["[0.1213,0.2347,0.3571]", "[0.3119,0.1723,0.0791]"]):
        _, sites = group.find_wyckoff_site(pos)
        sites = np.asarray(convert_to_primitive(group.info.lattice, sites, shift=True), dtype=float)
        sites = np.unique(np.round(sites, 8) % 1, axis=0)
        positions += sites.tolist()
        numbers += [no + 1] * len(sites)

    return (np.asarray(A, dtype=float).tolist(), positions, numbers)


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

        # absolute path, so that the results do not depend on the current directory.
        self._topdir = os.path.abspath(topdir)
        self._verbose = verbose
        self._mm = MaterialModel(self._topdir, verbose=verbose)
        self._local = create_all_local_operator()

        self.reset()

    # ==================================================
    @property
    def samb(self):
        """
        SAMB control.

        Returns:
            - (dict) -- SAMB control.

        :meta private:
        """
        return self._samb

    # ==================================================
    @property
    def wannier(self):
        """
        Wannier control.

        Returns:
            - (dict) -- Wannier control.

        :meta private:
        """
        return self._wannier

    # ==================================================
    @property
    def output(self):
        """
        Output for physical quantities control.

        Returns:
            - (dict) -- output for physical quantities control.

        :meta private:
        """
        return self._output

    # ==================================================
    @property
    def model(self):
        """
        Material Model.

        Returns:
            - (MaterialModel) -- matrial model.

        :meta private:
        """
        return self._mm

    # ==================================================
    @property
    def parameter(self):
        """
        Parameter of SAMBs.

        Returns:
            - (dict) -- parameter dict.

        :meta private:
        """
        return self._parameter

    # ==================================================
    @property
    def basis_type(self):
        """
        Atomic basis type.

        Returns:
            - (str) -- basis type, "lg/lgs/jml".

        :meta private:
        """
        return self._basis_type

    # ==================================================
    @property
    def basis(self):
        """
        Full-matrix basis.

        Returns:
            - (list) -- basis, list of "orbital@atom(sublattice)".

        :meta private:
        """
        return self._basis

    # ==================================================
    @property
    def HR(self):
        """
        Real-space hamiltonian, H(R).

        Returns:
            - (dict) -- H(R), dict[(n1,n2,n3,m,n), val or (val, bond_no)].

        :meta private:
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

        :meta private:
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
    def write_samb_matrix(self, matrix_info):
        """
        Save SAMB matrix.

        Args:
            matrix_info (dict): matrix info.

        :meta private:
        """
        # convert sympy to str (the same expressions appear many times, so that each is converted once).
        to_str = functools.cache(lambda v: str(v).replace(" ", ""))
        mi = matrix_info.copy()
        matrix = matrix_info["matrix"]
        mi["matrix"] = {z: {k: (to_str(v[0]), v[1]) for k, v in elm.items()} for z, elm in matrix.items()}

        self.write_dict(mi, "matrix.py", _matrix_comment, "info")

    # ==================================================
    def write_samb_hr(self, matrix_info, parameter, HR):
        """
        Save SAMB matrix in hr format.

        Args:
            matrix_info (dict): matrix info.
            parameter (dict): parameter dict, dict[z#, value].
            HR (dict, optional): H(R) matrix. if None, H(R) is generated.

        :meta private:
        """
        if HR is None:
            return

        # write hr.
        name = self["info"]["name"]
        filename = os.path.join(self._topdir, name, "info", name + "_hr.dat")
        ket = matrix_info["ket"]
        pos = matrix_info["ket_pos"]
        with open(filename, mode="w", encoding="utf-8") as f:
            print(f"# SAMB matrix from {matrix_info['source']} ({matrix_info['created']})", file=f)
            print("# select", file=f)
            for k, v in matrix_info["select"].items():
                print(f"#   {k}: {str(v).replace(' ', '')}", file=f)
            print(f"# basis ({matrix_info['dimension']})", file=f)
            for no, (b, p) in enumerate(zip(ket, pos)):
                p = [float(i) for i in p]
                print(f"#   {no:2d} {b}: [{p[0]: .6f}, {p[1]: .6f}, {p[2]: .6f}]", file=f)
            for z, v in parameter.items():
                print(f"# {z:<4} = {v}", file=f)
            print("#", file=f)
            print("# n1   n2   n3    m    n    re                        im", file=f)
            for (n1, n2, n3, m, n), v in HR.items():
                n1, n2, n3, m, n = int(n1), int(n2), int(n3), int(m), int(n)
                v = complex(v)
                r, i = v.real, v.imag
                s = f"{n1: 4d} {n2: 4d} {n3: 4d} {m: 4d} {n: 4d}    {r: .15e}    {i: .15e}"
                print(s, file=f)

        if self._verbose:
            print(f"save hr to '{filename}'.")

    # ==================================================
    def write_info(self):
        """
        Save info. (samb, wannier, output, and info).

        :meta private:
        """
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
    def write_z_file(self, parameter):
        """
        Save z file.

        Args:
            parameter (dict): parameter dict.

        :meta private:
        """
        if not parameter:
            return

        dic = {tag: float(v) for tag, v in parameter.items()}
        dic = dict(sorted(dic.items(), key=lambda item: abs(item[1]), reverse=True))

        comment = _param_comment + f"- by using '{self['info']['mode']}' mode.\n"
        self.write_dict(dic, "z.py", comment, "info")

    # ==================================================
    def write_var_file(self, var):
        """
        Save var file.

        Args:
            var (dict): var dict.

        :meta private:
        """
        d = {name: {zj: str(ex).replace(" ", "") for zj, ex in dic.items()} for name, dic in var.items()}

        self.write_dict(d, "var.py", _zj_var_comment, "info")

    # ==================================================
    def write_k_multipole(self, k_multipole):
        """
        Save k multipole.

        Args:
            k_multipole (dict): k-multipole dict.

        :meta private:
        """
        self.write_dict(k_multipole, "k.py", _k_matrix_comment, "info")

    # ==================================================
    def write_dispersion(self, Ek, Ok, op_lst, k_linear, k_dis_pos):
        """
        Save dispersion.

        Args:
            Ek (ndarray): energy eigen values.
            Ok (ndarray): expectation values of local operators.
            op_lst (list_): operator name list.
            k_linear (ndarray): linear k position along high-symmetry line.
            k_dis_pos (dict): k discrete point.

        :meta private:
        """
        name = self["info"]["name"]
        ef = self["info"]["fermi_energy"]

        path = os.path.join(self._topdir, name, self.output["dir"])
        fname = os.path.join(path, name + "_dispersion.txt")
        colormap = len(Ok) > 0
        if Ok:
            output_dispersion(fname, k_linear, ef, Ek, Ok, op_lst)
        else:
            output_dispersion(fname, k_linear, ef, Ek)
        plot_save_dispersion(fname, k_dis_pos, ef, colormap)
        create_gnuplot_cmd(fname, k_dis_pos, np.max(k_linear), np.nanmax(Ek), np.nanmin(Ek), ef, colormap)
        if self._verbose:
            print(f"save dispersion files into '{path}'.")

    # ==================================================
    def read_control(self, control):
        """
        Read control file.

        Args:
            control (str): control file name (relative to topdir).

        Returns:
            - (dict) -- control dict.

        :meta private:
        """
        if not control.endswith(".py"):
            raise ValueError(f"control file must be '.py' file, '{control}' is given.")

        return read_dict(control, self._topdir)

    # ==================================================
    def read_parameter(self, filename=None):
        """
        Read parameter file.

        Args:
            filename (str, optional): '.py' file name relative to 'topdir/name/info'. for empty str, use default, 'topdir/name/info/name_z.py'.

        Returns:
            - (dict) -- parameter dict.

        :meta private:
        """
        name = self["info"]["name"]
        if not filename:
            filename = f"{name}_z.py"
        path = os.path.join(self._topdir, name, "info")
        if not filename.endswith(".py"):
            raise ValueError(f"parameter file must be '.py' file, '{filename}' is given.")
        if os.path.isabs(filename):
            raise ValueError(f"parameter file must be relative to '{path}', '{filename}' is given.")

        filename = os.path.join(path, filename)
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"parameter file '{filename}' is not found (samb/parameter is relative to '{path}').")
        parameter = read_dict(filename)
        if self._verbose:
            print(f"load parameter from '{filename}'.")

        return parameter

    # ==================================================
    def reset(self, control=None):
        """
        Reset all data, and overwrite from control.

        Args:
            control (dict, optional): control dict.

        :meta private:
        """
        if control is None:
            control = {}

        # wannier/read_KS has been removed: H(R) is always read from seedname_hr.dat.
        if "read_KS" in control.get("wannier", {}):
            if control["wannier"]["read_KS"]:
                raise ValueError("'wannier/read_KS' has been removed. H(R) is read from seedname_hr.dat.")
            logger.warning("'wannier/read_KS' has been removed and is ignored. remove it from the control.")
            control = copy.deepcopy(control)
            del control["wannier"]["read_KS"]

        # check unknown keys, e.g., typo.
        check_dict_keys(
            control,
            default_control,
            name="control",
            allowed={"samb/select": ["site", "bond", "X", "l", "Gamma", "s"]},
            free=["samb/parameter", "output/dispersion/k_point"],
        )

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
        grid = control.get("grid", default_control["grid"])
        if not isinstance(grid, (list, tuple)) or len(grid) != 3:
            raise ValueError(f"'grid' must be three positive integers, '{grid}' is given.")
        self.set_grid(*grid)
        # in wannier mode without a model, the name is seedname.
        self.set_name(self.samb["model"] if self.samb["model"] is not None else self.wannier["seedname"])
        self._use_model = False
        self._win_structure = None
        self.set_parameter(None)
        self.set_basis_type(None)
        self.set_basis(None)
        self.set_HR(None)
        self._k_path_info = None  # k path info. of the dispersion, used for callback.
        self._eigen_path = None  # (Ek, Uk) along k path.
        self._eigenvalue_full = None  # (kv, Ek) on full k grid.

    # ==================================================
    def set_mode(self, mode):
        """
        Set analysis mode.

        Args:
            mode (str): analysis mode, "samb/wannier/symcw".

        :meta private:
        """
        self["info"]["mode"] = mode

    # ==================================================
    def set_name(self, name):
        """
        Set model name.

        Args:
            name (str): model name.

        :meta private:
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

        :meta private:
        """
        if any(isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 1 for n in (N1, N2, N3)):
            raise ValueError(f"'grid' must be three positive integers, '{(N1, N2, N3)}' is given.")
        self["info"]["grid"] = [int(N1), int(N2), int(N3)]

    # ==================================================
    def set_parameter(self, parameter):
        """
        Set parameter, zj.

        Args:
            parameter (dict): parameter dict.

        :meta private:
        """
        self._parameter = parameter

    # ==================================================
    def set_basis_type(self, basis_type):
        """
        Set atomic basis type.

        Args:
            basis_type (str): basis type, "lg/lgs/jml".

        :meta private:
        """
        self._basis_type = basis_type

    # ==================================================
    def set_basis(self, basis):
        """
        Set full-matrix basis.

        Args:
            basis (list): basis, list of "orbital@atom(sublattice)".

        :meta private:
        """
        self._basis = basis

    # ==================================================
    def set_HR(self, HR):
        """
        Set real-space hamiltonian, H(R).

        Args:
            HR (dict): H(R), dict[(n1,n2,n3,m,n), val or (val, bond_no)].

        :meta private:
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
    def set_fermi_energy(self, ef):
        """
        Set Fermi energy.

        Args:
            ef (float): Fermi energy.

        :meta private:
        """
        self["info"]["fermi_energy"] = ef

    # ==================================================
    def get_var(self, matrix_info):
        """
        Get var dict.

        Args:
            matrix_info (dict): matrix info dict.

        Returns:
            - (dict) -- var dict.

        :meta private:
        """
        IR = next(iter(self.model.group.character["table"].keys()))  # identity irrep.
        conv_dict = convert_zj_atomic_var(matrix_info, self.model["combined_cluster"], self.model["combined_id"], IR)
        return conv_dict

    # ==================================================
    def get_local_operator(self, tag):
        """
        Create local operator.

        Args:
            tag (str): operator name, "Sx/Sy/Sz/Lx/Ly/Lz/Qu/Qv/Qyz/Qxz/Qxy".

        Returns:
            - (ndarray) -- operator matrix (dim x dim).

        :meta private:
        """
        spinful = self.basis_type == "lgs"
        return create_local_operator(self.basis, tag, self._local, spinful)

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
        if k_path == "":  # create default path for the primitive cell, info/A.
            if self._use_model:
                structure = _kpath_structure(self.model.group, self["info"]["A"])
            else:
                structure = self._win_structure

            info = seekpath.get_path_orig_cell(structure)
            if self._use_model and info["spacegroup_number"] != int(self.model.group.ID):
                raise ValueError(
                    f"space group found by seekpath (No. {info['spacegroup_number']}) differs from that of the model, "
                    f"{self.model.group}. Give k_path and k_point explicitly."
                )
            k_point = {k: [float(i) for i in v] for k, v in info["point_coords"].items()}
            k_point["Γ"] = k_point.pop("GAMMA")
            k_path = _join_kpath(info["path"])
        else:
            k_point = self.output["dispersion"].get("k_point", {})
            k_point = {k: str_to_sympy(v).astype(float) for k, v, in k_point.items()}

        return k_point, k_path

    # ==================================================
    def get_k_multipole(self, matrix_info):
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
            "ket": matrix_info["ket"],
            "index": matrix_info["index"],
            "cluster_vector": cluster_vec,
            "k_multipole": k_multipole,
            "k_matrix": k_matrix,
        }

        return k_multipole

    # ==================================================
    def analyze(self, control):
        """
        Analyze model with control file.

        Args:
            control (str or dict): control file (.py) or control dict.
        """
        # read control.
        if isinstance(control, str):  # read control file.
            control = self.read_control(control)

        self.reset(control)
        mode = self["info"]["mode"]

        # execute SAMB mode.
        if mode in ["samb", "symcw"]:
            self.exec_samb()  # create SAMBs, and H(R) if zj are provided.

        # execute wannier mode.
        if mode in ["wannier", "symcw"]:
            if mode == "wannier" and self.samb["model"] is not None:
                self.load_model(self.samb["model"])  # use lattice, k path, and ket of the model.
            self.exec_wannier()  # create H(R), and zj in case of "symcw".
            if mode == "symcw":
                matrix_info = self["samb"]["matrix_info"]  # created by exec_samb.
                Zr_dict = matrix_info["matrix"]
                parameter = decompose_operator_by_SAMB(self.HR, Zr_dict)
                HR = self.model.get_hr(parameter, Zr_dict)  # overwrite HR by MultiPie.
                self.write_samb_hr(matrix_info, parameter, HR)
                self.set_HR(HR)
                self.set_parameter(parameter)

        # create z file.
        self.write_z_file(self.parameter)

        # compute physical quanties and output data.
        self.compute_physical_quantity()

        # output info.
        self.write_info()

    # ==================================================
    def exec_samb(self):
        """
        Execute SAMB mode.

        :meta private:
        """
        # read and set model.
        self.load_model(self.samb["model"])

        # set selected SAMBs.
        matrix_info = self.model.get_samb_matrix(self.samb["select"])

        # create var file.
        var = self.get_var(matrix_info)
        self.write_var_file(var)

        # create k-multipole file.
        if self.samb["k_multipole"]:
            k_multipole = self.get_k_multipole(matrix_info)
            if k_multipole:
                self.write_k_multipole(k_multipole)

        # create SAMB qtdraw.
        if self.samb["samb_figure"]:
            self.model.save_samb_qtdraw()

        parameter = self.samb["parameter"]
        if isinstance(parameter, str):  # when parameter is str, read z file.
            parameter = self.read_parameter(parameter)
        # a value given as str, e.g., "1/2" or "sqrt(2)", is evaluated by SymPy.
        parameter = {tag: float(str_to_sympy(v, rational=False)) if isinstance(v, str) else v for tag, v in parameter.items()}

        # determine local weight if NG_sum_rule is True.
        ng = self.samb["NG_sum_rule"]
        if ng and parameter:
            parameter = add_local_parameter(matrix_info, parameter, self.model["full_matrix"]["ket"])
        if ng:
            parameter_sym = add_local_parameter_sym(matrix_info, self.model["full_matrix"]["ket"])
            self["samb"]["NG_sum_rule"] = parameter_sym

        # output matrix.py and hr.dat.
        if parameter:
            HR = self.model.get_hr(parameter, matrix_info["matrix"])
            self.write_samb_hr(matrix_info, parameter, HR)
            self.set_HR(HR)
        self.write_samb_matrix(matrix_info)

        self.set_parameter(parameter)
        self.set_fermi_energy(0.0)
        self["samb"]["matrix_info"] = matrix_info

    # ==================================================
    def load_model(self, name):
        """
        Load model, and set name, basis, and primitive cell of the model.

        Args:
            name (str): model name.

        :meta private:
        """
        if name is None:
            raise ValueError("no model is specified in 'samb/model'.")
        self.set_name(name)
        self.model.load(name)
        self.set_basis_type(self.model["basis_type"])
        self.set_basis(self.model["full_matrix"]["ket"])
        self.set_primitive_cell(model_primitive_vector(self.model))
        self._use_model = True

    # ==================================================
    def exec_wannier(self):
        """
        Execute wannier mode.

        :meta private:
        """
        seedname = self.wannier["seedname"]
        wannier_dir = os.path.join(self._topdir, seedname, self.wannier["dir"])

        # read seedname.win
        win = read_win(seedname, wannier_dir)
        # read seedname.nnkp
        nnkp = read_nnkp(seedname, wannier_dir)

        if self._use_model:
            # map Wannier functions onto the model (loaded by exec_samb or load_model).
            mapping = map_wannier_to_model(nnkp, win["A"], self.model, self.wannier.get("ket_wannier", []))
            w2m = mapping["w2m"]
            m2w = [w2m.index(m) for m in range(len(w2m))]
            model_ket = self.model.get_ket_site()
            ket_multipie = list(model_ket.keys())
            atoms_frac = [list(v) for v in model_ket.values()]
            atoms_cart = (np.asarray(atoms_frac, dtype=float) @ model_primitive_vector(self.model)).tolist()
        else:
            missing = [key for key in ("nw2n", "nw2l", "nw2m", "nw2r", "nw2s", "atom_pos_r") if nnkp.get(key) is None]
            if missing:
                raise ValueError(
                    f"projection information is missing in {seedname}.nnkp (e.g., auto_projections): {', '.join(missing)}."
                )
            wannier_ket_info = {
                "A": win["A"],
                "atoms_frac": win["atoms_frac"],
                "atoms_cart": win["atoms_cart"],
                "fermi_energy": win["fermi_energy"],
                "nw2n": nnkp["nw2n"],
                "nw2l": nnkp["nw2l"],
                "nw2m": nnkp["nw2m"],
                "nw2r": nnkp["nw2r"],
                "nw2s": nnkp["nw2s"],
                "atom_pos_r": nnkp["atom_pos_r"],
            }
            w2m, m2w, ket_multipie, atoms_frac, atoms_cart = create_ket_wannier_multipie(wannier_ket_info)
            w_ket = self.wannier.get("ket_wannier", [])
            if w_ket:
                if len(set(ket_multipie)) != len(ket_multipie):
                    raise ValueError(f"ket names are not unique, {ket_multipie}, so ket_wannier cannot be used.")
                if sorted(w_ket) != sorted(ket_multipie):
                    raise ValueError(f"ket_wannier must be a permutation of {ket_multipie}.")
                m2w = [w_ket.index(m) for m in ket_multipie]
                w2m = [no for no, i in sorted(enumerate(m2w), key=lambda x: x[1])]
                # positions (projection centres) in the given correspondence.
                centres = np.asarray(nnkp["atom_pos_r"], dtype=float)
                atoms_frac = [centres[nnkp["nw2n"][w]].tolist() for w in m2w]
                atoms_cart = (np.asarray(atoms_frac) @ np.asarray(win["A"], dtype=float)).tolist()
            mapping = None

            # lattice and structure of seedname.win (without model).
            self.set_primitive_cell(win["A"])
            species = {}
            numbers = [species.setdefault(atom, len(species) + 1) for atom, _ in win["atoms_frac"].keys()]
            self._win_structure = (win["A"], list(win["atoms_frac"].values()), numbers)

        hr_file = seedname + "_hr.dat"
        hr_dict, irvec, ndegen = read_hr(hr_file, wannier_dir)
        # real-space hoppings with Wigner-Seitz degeneracy (and use_ws_distance if seedname_wsvec.dat exists).
        hr_dict = apply_ws_degeneracy(hr_dict, irvec, ndegen, read_wsvec(seedname, wannier_dir))
        if mapping is None:
            # convert from wannier index to multipie index.
            HR = {(n1, n2, n3, w2m[w1], w2m[w2]): (complex(v), None) for (n1, n2, n3, w1, w2), v in hr_dict.items()}
        else:
            # convert to multipie index and primitive lattice of the model through bond vectors.
            HR = convert_hr_to_model(hr_dict, mapping)

        info = {
            "ket": ket_multipie,
            "atoms_frac": atoms_frac,
            "atoms_cart": atoms_cart,
            "wannier_to_multipie": w2m,
            "multipie_to_wannier": m2w,
        }
        if mapping is not None:
            info["lattice_transformation"] = mapping["U"].tolist()
            info["origin_shift"] = mapping["t"].tolist()

        self["wannier"] = info
        self.set_fermi_energy(win["fermi_energy"])

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

        self.set_HR(HR)

    # ==================================================
    def compute_physical_quantity(self):
        """
        Compute physical quantities by parsing the control file.

        :meta private:
        """
        if self.HR is None:
            if self._verbose:
                print("set H(R) first before calculating physical quantities.")
            return

        name = self["info"]["name"]
        path = os.path.join(self._topdir, name, self.output["dir"])
        os.makedirs(path, exist_ok=True)

        # check before the computation: callbacks are loaded, and the dispersion exists if it is required.
        modules = {}  # a file given for both callbacks is loaded once.
        callback = {key: self.load_callback(key, modules) for key in ["callback_k_path", "callback_k_full"]}
        for key in ["callback_k_path", "save_k_path"]:
            if self.output[key] and not self.has_dispersion():
                raise ValueError(
                    f"'output/{key}' requires the dispersion, which is not computed "
                    "('output/dispersion/k_path' is None, or the group of the model is not a space group)."
                )

        # eigen system on the full k grid is kept only when it is required (otherwise, it is computed in chunks).
        need_full = callback["callback_k_full"] is not None or self.output["save_k_full"]
        eigen_full = self.get_eigen_full() if need_full else None

        # fermi energy from the number of occupied states.
        if self.output["n_electron"] is not None:
            self.set_fermi_energy(self.compute_fermi_energy(eigen_full))

        disp = self.compute_dispersion()
        self["output"]["dispersion"] = disp

        self.compute_dos(eigen_full)
        self.save_eigen(eigen_full)

        # user-defined calculations.
        if callback["callback_k_path"] is not None:
            Ek, Uk = self._eigen_path
            info = self.get_callback_info(self._k_path_info["k_grid"])
            info["k_linear"] = self._k_path_info["k_linear"]
            info["k_dis_pos"] = self._k_path_info["k_dis_pos"]
            callback["callback_k_path"](Ek, Uk, info)
        if callback["callback_k_full"] is not None:
            kv, Ek, Uk = eigen_full
            info = self.get_callback_info(kv)
            info["grid"] = list(self["info"]["grid"])
            callback["callback_k_full"](Ek, Uk, info)

    # ==================================================
    def has_dispersion(self):
        """
        Is the dispersion computed ?

        Returns:
            - (bool) -- True if the dispersion is computed.

        :meta private:
        """
        if "dispersion" not in self.output or self.output["dispersion"]["k_path"] is None:
            return False

        return not (self._use_model and self.model.group.group_type != "SG")

    # ==================================================
    def get_ket(self):
        """
        Get ket names in the order of the matrix of H(k).

        Returns:
            - (list) -- ket names, "orbital@atom(sublattice)".

        :meta private:
        """
        if self._use_model:
            return list(self.model.get_ket_site().keys())
        else:
            return list(self["wannier"]["ket"])

    # ==================================================
    def get_hk_data(self):
        """
        Get data to compute H(k).

        Returns:
            - (dict) -- H(R), dict[((n1,n2,n3),m,n), value].
            - (ndarray) -- atomic positions (fractional).
            - (bool) -- TB gauge ?

        :meta private:
        """
        tb_gauge = self.output["fourier"]["tb_gauge"]
        if self._use_model:
            atom = np.asarray(list(self.model.get_ket_site().values()), dtype=float)
        else:
            atom = np.asarray(self["wannier"]["atoms_frac"], dtype=float)

        HR = {
            ((int(n1), int(n2), int(n3)), int(m), int(n)): complex(v[0] if isinstance(v, tuple) else v)
            for (n1, n2, n3, m, n), v in self.HR.items()
        }

        return HR, atom, tb_gauge

    # ==================================================
    def get_k_full(self):
        """
        Get k points on the full grid.

        Returns:
            - (ndarray) -- k points (fractional, primitive), [N1*N2*N3, 3], or Γ point only for a model of point group.

        :meta private:
        """
        if self._use_model and self.model.group.group_type != "SG":
            return np.zeros((1, 3))

        return grid_all(*self["info"]["grid"])

    # ==================================================
    def eigen_chunks(self, kv, vectors=True, chunk=None):
        """
        Eigen system of H(k) for chunks of k points.

        Args:
            kv (ndarray): k points (fractional).
            vectors (bool, optional): compute eigen vectors ? if False, Uk is None.
            chunk (int, optional): number of k points in a chunk. if None, it is determined by the matrix size.

        Yields:
            - (tuple) -- (Ek, Uk) of the chunk, Ek[k,n], Uk[k,:,n].

        :meta private:
        """
        HR, atom, tb_gauge = self.get_hk_data()
        if chunk is None:
            chunk = _k_chunk(len(atom))
        for i in range(0, len(kv), chunk):
            Hk = fourier_r_to_k(HR, atom, kv[i : i + chunk], tb_gauge)
            if vectors:
                yield np.linalg.eigh(Hk)
            else:
                yield np.linalg.eigvalsh(Hk), None

    # ==================================================
    def get_eigen_full(self):
        """
        Get eigen system on the full k grid.

        Returns:
            - (ndarray) -- k points (fractional, primitive).
            - (ndarray) -- eigen values, Ek[k,n].
            - (ndarray) -- eigen vectors, Uk[k,:,n].

        :meta private:
        """
        kv = self.get_k_full()
        Ek, Uk = zip(*self.eigen_chunks(kv))

        return kv, np.concatenate(Ek), np.concatenate(Uk)

    # ==================================================
    def get_callback_info(self, k_grid):
        """
        Get info. dict passed to callback.

        Args:
            k_grid (ndarray): k points of Ek and Uk (fractional, primitive).

        Returns:
            - (dict) -- info. dict.

        :meta private:
        """
        k_path = self._k_path_info if self._k_path_info is not None else {}
        info = {
            "k_path": k_path.get("k_path"),
            "k_point": k_path.get("k_point"),
            "k_grid": k_grid,
            "ket": self.get_ket(),
            "fermi_energy": self["info"]["fermi_energy"],
            "name": self["info"]["name"],
            "_topdir": self._topdir,
            "_verbose": self._verbose,
        }

        return info

    # ==================================================
    def load_callback(self, key, modules=None):
        """
        Load callback function.

        Args:
            key (str): "callback_k_path" or "callback_k_full", which is also the function name.
            modules (dict, optional): loaded modules, dict[path, module], to load a file once.

        Returns:
            - (function) -- callback function, or None if not given.

        :meta private:
        """
        filename = self.output[key]
        if filename is None:
            return None
        if not isinstance(filename, str) or not filename.endswith(".py"):
            raise ValueError(f"'output/{key}' must be '.py' file, '{filename}' is given.")

        path = os.path.join(self._topdir, filename)
        if not os.path.isfile(path):
            raise FileNotFoundError(f"callback file '{path}' is not found ('output/{key}' is relative to '{self._topdir}').")
        if modules is None:
            modules = {}
        path = os.path.realpath(path)
        if path not in modules:
            spec = importlib.util.spec_from_file_location(f"_multipie_callback_{len(modules)}", path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            modules[path] = module
        module = modules[path]
        func = getattr(module, key, None)
        if not callable(func):
            raise ValueError(f"function '{key}' is not defined in '{path}'.")
        if self._verbose:
            print(f"load {key} from '{path}'.")

        return func

    # ==================================================
    def compute_dispersion(self):
        """
        Compute dispersion.

        :meta private:
        """
        # check if dispersion can be computed.
        if not self.has_dispersion():
            return

        # get k_point and k_path.
        k_point, k_path = self.get_kpath(self.output["dispersion"]["k_path"])
        N1 = self["info"]["grid"][0]
        B = np.asarray(self["info"]["B"])  # [b1,b2,b3] as rows, k_cart = k_frac @ B.
        k_point_path, k_linear, k_dis_pos = grid_path(k_point, k_path, N1, B.T)

        # get local operator list.
        if self.basis_type == "jml" or self.basis_type is None:
            op_lst = []
        else:
            op_lst = self.output["dispersion"]["local"]

        # set eigen system.
        Ek, Uk = next(self.eigen_chunks(k_point_path, chunk=len(k_point_path)))
        self._eigen_path = (Ek, Uk)
        self._k_path_info = {
            "k_path": k_path,
            "k_point": {k: np.asarray(v, dtype=float).tolist() for k, v in k_point.items()},
            "k_grid": k_point_path,
            "k_linear": k_linear,
            "k_dis_pos": k_dis_pos,
        }
        power = self.output["dispersion"]["power"]
        if power is not None:
            with np.errstate(invalid="ignore", divide="ignore"):
                Ek = np.power(Ek, power)
            Ek = np.where(np.isfinite(Ek), Ek, np.nan)
            n_skip = np.count_nonzero(~np.isfinite(Ek))
            if n_skip == Ek.size:
                raise ValueError(f"no value of Ek^{power} can be computed (e.g., all Ek are negative). check 'power'.")
            if n_skip > 0 and self._verbose:
                print(f"skip {n_skip} values of Ek^{power} that cannot be computed (e.g., negative Ek).")

        # set local operators.
        Ok = [np.einsum("kmi,mn,kni->ki", Uk.conj(), self.get_local_operator(tag), Uk).real for tag in op_lst]

        # output dispersion data, plot, and gnuplot.
        self.write_dispersion(Ek, Ok, op_lst, k_linear, k_dis_pos)

        # orbital weight (fat band).
        weight = self.output["dispersion"]["weight"]
        if weight is not False and weight is not None:
            labels, group = self.get_ket_group(weight, self.get_ket(), "output/dispersion/weight")
            Wk = [np.sum(np.abs(Uk[:, g, :]) ** 2, axis=1) for g in group]
            self.write_weight(Ek, Wk, labels, k_linear, k_dis_pos)

        # save dispersion info.
        d = {
            "k_path": k_path,
            "k_point": {k: str(v).replace(" ", "") for k, v in k_point.items()},
            "e_max": float(np.nanmax(Ek)),
            "e_min": float(np.nanmin(Ek)),
        }

        return d

    # ==================================================
    def write_weight(self, Ek, Wk, labels, k_linear, k_dis_pos):
        """
        Save dispersion with orbital weight (fat band).

        Args:
            Ek (ndarray): energy eigen values.
            Wk (list): orbital weight for each group, [Ek.shape].
            labels (list): label of each group.
            k_linear (ndarray): linear k position along high-symmetry line.
            k_dis_pos (dict): k discrete point.

        :meta private:
        """
        name = self["info"]["name"]
        ef = self["info"]["fermi_energy"]

        path = os.path.join(self._topdir, name, self.output["dir"])
        fname = os.path.join(path, name + "_weight.txt")
        output_dispersion(fname, k_linear, ef, Ek, Wk, labels)
        plot_save_weight(fname, k_dis_pos, labels)
        create_gnuplot_weight(fname, k_dis_pos, np.max(k_linear), np.nanmax(Ek), np.nanmin(Ek), ef, labels)
        if self._verbose:
            print(f"save weight files into '{path}'.")

    # ==================================================
    def get_ket_group(self, spec, ket, key, spin=False):
        """
        Get groups of kets for partial DOS or orbital weight.

        Args:
            spec (str or list): "auto" (or True), or [[ket name]].
            ket (list): ket names.
            key (str): control key used in error message.
            spin (bool, optional): resolve spin, u/d, for "auto" ?

        Returns:
            - (list) -- label of each group.
            - (list) -- ket indices of each group.

        Note:
            - for "auto", kets are grouped by the site cluster and the orbital, and the spin, "u" or "d", is summed unless spin is True.
            - for [[ket name]], the label is the ket names joined by "+".

        :meta private:
        """
        if spec is True or spec == "auto":
            group = {}
            for no, k in enumerate(ket):
                orbital, site = k.rsplit("@", 1)
                m = re.fullmatch(r"\((.+),([ud])\)", orbital)  # (orbital,spin).
                if m and not spin:
                    orbital = m.group(1)
                site = re.sub(r"\(\d+\)$", "", site)  # remove sublattice.
                group.setdefault(f"{orbital}@{site}", []).append(no)
            return list(group.keys()), list(group.values())

        msg = f"'{key}' must be False, 'auto', or a list of lists of ket names, e.g., [['s@A(1)', 's@A(2)']]"
        if not isinstance(spec, (list, tuple)) or not spec or not all(isinstance(g, (list, tuple)) and g for g in spec):
            raise ValueError(f"{msg}, '{spec}' is given.")
        unknown = [k for g in spec for k in g if k not in ket]
        if unknown:
            raise ValueError(f"unknown ket(s) in '{key}': {unknown}. kets are {ket}.")

        labels = ["+".join(g) for g in spec]
        group = [[no for no, k in enumerate(ket) if k in g] for g in spec]

        return labels, group

    # ==================================================
    def get_dos_method(self):
        """
        Get DOS method (also used for fermi energy).

        Returns:
            - (str) -- "gaussian" or "tetrahedron".

        :meta private:
        """
        method = self.output["dos_method"]
        if method not in ["gaussian", "tetrahedron"]:
            raise ValueError(f"'output/dos_method' must be 'gaussian' or 'tetrahedron', '{method}' is given.")
        if method == "tetrahedron" and self._use_model and self.model.group.group_type != "SG":
            raise ValueError("'output/dos_method' = 'tetrahedron' cannot be used for a model of point group (Γ point only).")

        return method

    # ==================================================
    def get_eigenvalue_full(self, eigen_full=None):
        """
        Get eigen values on the full k grid (computed once).

        Args:
            eigen_full (tuple, optional): (kv, Ek, Uk) on the full k grid. if None, the eigen values are computed.

        Returns:
            - (ndarray) -- k points (fractional, primitive).
            - (ndarray) -- eigen values, Ek[k,n].

        :meta private:
        """
        if eigen_full is not None:
            return eigen_full[0], eigen_full[1]
        if self._eigenvalue_full is None:
            kv = self.get_k_full()
            Ek = np.concatenate([e for e, _ in self.eigen_chunks(kv, vectors=False)])
            self._eigenvalue_full = (kv, Ek)

        return self._eigenvalue_full

    # ==================================================
    def compute_fermi_energy(self, eigen_full=None):
        """
        Compute fermi energy from the number of occupied states, output/n_electron.

        Args:
            eigen_full (tuple, optional): (kv, Ek, Uk) on the full k grid. if None, the eigen values are computed.

        Returns:
            - (float) -- fermi energy.

        Note:
            - n_electron is the number of occupied states per primitive cell, where each band has one state.
            - for "gaussian" method, the states are occupied in ascending order at T=0 (fermi energy is at the middle of
              the highest occupied and the lowest unoccupied states). for "tetrahedron" method, libtetrabz is used, and
              the middle of the gap is used for an insulator.

        :meta private:
        """
        n = self.output["n_electron"]
        kv, Ek = self.get_eigenvalue_full(eigen_full)
        nk, nb = Ek.shape
        if isinstance(n, bool) or not isinstance(n, (int, float)) or not 0 <= n <= nb:
            raise ValueError(f"'output/n_electron' must be a number in [0, {nb}] (number of bands), '{n}' is given.")

        method = self.get_dos_method()
        if method == "tetrahedron":
            libtetrabz = import_libtetrabz()
            eig = Ek.reshape(*self["info"]["grid"], nb)
            ef, _, _ = libtetrabz.fermieng(np.asarray(self["info"]["B"]), eig, float(n))
            # in a gap, the bisection of libtetrabz stops at any point, so that the middle of the gap is used as in "gaussian".
            if n == int(n) and 0 < n < nb and np.max(Ek[:, : int(n)]) < np.min(Ek[:, int(n) :]):
                ef = 0.5 * (np.max(Ek[:, : int(n)]) + np.min(Ek[:, int(n) :]))
        else:
            e = np.sort(Ek.ravel())
            ns = n * nk  # number of occupied states.
            ns_i = int(np.rint(ns))
            if abs(ns - ns_i) > 1e-8:  # partially occupied level.
                ef = e[int(np.floor(ns))]
            elif ns_i == 0:
                ef = e[0]
            elif ns_i == len(e):
                ef = e[-1]
            else:
                ef = 0.5 * (e[ns_i - 1] + e[ns_i])
        ef = float(ef)

        self["info"]["n_electron"] = n
        if self._verbose:
            print(f"fermi energy = {ef} for n_electron = {n} ({method}).")

        return ef

    # ==================================================
    def compute_dos(self, eigen_full=None):
        """
        Compute DOS (total and partial).

        Args:
            eigen_full (tuple, optional): (kv, Ek, Uk) on the full k grid. if None, it is computed in chunks.

        Note:
            - DOS per primitive cell (summed over spin), D(E) = (1/Nk) sum_{kn} w_kn delta(E-E_kn).
            - w_kn = 1 for total DOS, and sum_{m in group} |U_k[m,n]|^2 for partial DOS. sum of partial DOS of all kets is total DOS.
            - "gaussian": delta function is broadened by Gaussian, and energy range includes all eigen values with margin of 5 x width.
              the width must not be smaller than the energy step.
            - "tetrahedron": integrated DOS N(E) by optimized tetrahedron method of libtetrabz, and DOS as (N(E+dE/2)-N(E-dE/2))/dE.
              energy range includes all eigen values with margin of 5% of the band width.

        :meta private:
        """
        dos = self.output["dos"]
        if dos is False or dos is None:
            return

        if self._verbose:
            print("compute and output dos.")

        ket = self.get_ket()
        labels, group = self.get_ket_group(dos, ket, "output/dos", self.output["dos_spin"])
        e_num = self.output["dos_e_num"]
        if isinstance(e_num, bool) or not isinstance(e_num, int) or e_num < 2:
            raise ValueError(f"'output/dos_e_num' must be integer larger than 1, '{e_num}' is given.")
        method = self.get_dos_method()

        # projection to groups.
        P = np.zeros((len(ket), len(group)))
        for g, idx in enumerate(group):
            P[idx, g] = 1.0

        # energy range from the eigen values.
        kv, Ek = self.get_eigenvalue_full(eigen_full)
        if eigen_full is None:
            chunks = self.eigen_chunks(kv)
        else:
            chunk = _k_chunk(len(ket))
            Uk = eigen_full[2]
            chunks = ((Ek[i : i + chunk], Uk[i : i + chunk]) for i in range(0, len(kv), chunk))
        e_min, e_max = float(np.min(Ek)), float(np.max(Ek))

        if method == "gaussian":
            if e_num < 12:  # energy step must not exceed the width, with the margin of 5 x width on both sides.
                raise ValueError(f"'output/dos_e_num' must be at least 12 for 'gaussian' method, '{e_num}' is given.")
            width = self.output["dos_width"]
            e_range = e_max - e_min
            if width is None:
                width = 0.01 * e_range if e_range > 0 else 0.01
                if e_num > 11:
                    width = max(width, e_range / (e_num - 11))  # not smaller than the energy step.
            elif isinstance(width, bool) or not isinstance(width, (int, float)) or width <= 0:
                raise ValueError(f"'output/dos_width' must be positive or None, '{width}' is given.")
            width = float(width)
            energy = np.linspace(e_min - 5 * width, e_max + 5 * width, e_num)
            dE = energy[1] - energy[0]
            if width < dE * (1 - 1e-10):
                raise ValueError(
                    f"Gaussian width of DOS, {width}, is smaller than the energy step, {dE}. "
                    "increase 'output/dos_e_num' or 'output/dos_width'."
                )

            # histogram on the energy grid (linear interpolation between the two nearest points, which keeps the
            # weight and the mean energy), and Gaussian broadening.
            hist = np.zeros((e_num, len(group) + 1))
            for Ec, Uc in chunks:
                x = ((Ec - energy[0]) / dE).ravel()
                i0 = np.clip(np.floor(x).astype(int), 0, e_num - 2)
                f = x - i0
                w = np.einsum("kmn,mg->kng", np.abs(Uc) ** 2, P).reshape(-1, len(group))
                for g in range(len(group) + 1):
                    wg = 1.0 if g == 0 else w[:, g - 1]
                    hist[:, g] += np.bincount(i0, weights=wg * (1 - f), minlength=e_num)
                    hist[:, g] += np.bincount(i0 + 1, weights=wg * f, minlength=e_num)
            dos = gaussian_filter1d(hist / (len(kv) * dE), width / dE, axis=0, mode="constant", truncate=5.0)
            integrated = cumulative_trapezoid(dos[:, 0], energy, initial=0.0)
        else:
            libtetrabz = import_libtetrabz()
            width = None
            margin = 0.05 * (e_max - e_min) if e_max > e_min else 0.1
            energy = np.linspace(e_min - margin, e_max + margin, e_num)
            dE = energy[1] - energy[0]

            # weight of each group, W[k,n,g].
            W = np.concatenate([np.einsum("kmn,mg->kng", np.abs(Uc) ** 2, P) for _, Uc in chunks])
            nk, nb = Ek.shape
            eig = Ek.reshape(*self["info"]["grid"], nb)
            B = np.asarray(self["info"]["B"])

            # integrated DOS by tetrahedron method at the energy grid and its mid points, for chunks of energy (to limit
            # memory). DOS is (N(E+dE/2)-N(E-dE/2))/dE, which keeps the normalization including flat bands.
            e_fine = energy[0] - dE / 2 + np.arange(2 * e_num + 1) * (dE / 2)
            N = np.zeros((len(e_fine), len(group) + 1))
            e_chunk = max(1, int(5e6 // (nk * nb)))
            for i in range(0, len(e_fine), e_chunk):
                w = libtetrabz.intdos(B, eig, e_fine[i : i + e_chunk]).reshape(nk, nb, -1)
                N[i : i + e_chunk, 0] = w.sum(axis=(0, 1))
                N[i : i + e_chunk, 1:] = np.einsum("kne,kng->eg", w, W)
            dos = (N[2::2] - N[0:-2:2]) / dE
            integrated = N[1::2, 0]

        # columns: total, integrated total, partial.
        dos = np.column_stack([dos[:, 0], integrated, dos[:, 1:]])

        self.write_dos(energy, dos, labels, [[ket[i] for i in g] for g in group], len(kv), method, width)
        self["output"]["dos"] = {
            "method": method,
            "e_max": e_max,
            "e_min": e_min,
            "width": width,
            "column": ["E", "total", "N(E)"] + labels,
        }

    # ==================================================
    def write_dos(self, energy, dos, labels, group, nk, method, width):
        """
        Save DOS, plot, and gnuplot.

        Args:
            energy (ndarray): energy grid.
            dos (ndarray): DOS, [energy, column], columns (total, integrated total, partial, ...).
            labels (list): label of partial DOS.
            group (list): ket names of partial DOS.
            nk (int): number of k points.
            method (str): DOS method.
            width (float): Gaussian broadening (None for tetrahedron method).

        :meta private:
        """
        name = self["info"]["name"]
        ef = self["info"]["fermi_energy"]
        path = os.path.join(self._topdir, name, self.output["dir"])
        fname = os.path.join(path, name + "_dos.txt")

        comment = [
            "DOS per primitive cell (summed over spin).",
            f"shifted by fermi energy = {ef}",
            f"number of k points = {nk} (grid = {list(self['info']['grid'])})" if nk > 1 else "Γ point only",
            f"Gaussian width = {width}" if method == "gaussian" else "optimized tetrahedron method (libtetrabz)",
            "column 1: E",
            "column 2: total",
            "column 3: N(E), integrated total DOS",
        ]
        comment += [f"column {no + 4}: {label} = {' + '.join(kets)}" for no, (label, kets) in enumerate(zip(labels, group))]
        output_dos(fname, energy, ef, dos, labels, comment)
        plot_save_dos(fname, labels)
        create_gnuplot_dos(fname, energy[0], energy[-1], np.max(dos[:, 0]), ef, labels)

        if self._verbose:
            print(f"save dos files into '{path}'.")

    # ==================================================
    def save_eigen(self, eigen_full=None):
        """
        Save eigen system along k path and/or on full k grid in .npz.

        Args:
            eigen_full (tuple, optional): (kv, Ek, Uk) on the full k grid.

        :meta private:
        """
        name = self["info"]["name"]
        path = os.path.join(self._topdir, name, self.output["dir"])
        common = {"ket": np.asarray(self.get_ket()), "fermi_energy": self["info"]["fermi_energy"]}

        if self.output["save_k_path"]:
            Ek, Uk = self._eigen_path
            info = self._k_path_info
            fname = os.path.join(path, name + "_eigen_k_path.npz")
            np.savez(fname, Ek=Ek, Uk=Uk, k_grid=info["k_grid"], k_linear=info["k_linear"], **common)
            if self._verbose:
                print(f"save eigen system along k path to '{fname}'.")

        if self.output["save_k_full"]:
            kv, Ek, Uk = eigen_full
            fname = os.path.join(path, name + "_eigen_k_full.npz")
            np.savez(fname, Ek=Ek, Uk=Uk, k_grid=kv, grid=np.asarray(self["info"]["grid"]), **common)
            if self._verbose:
                print(f"save eigen system on full k grid to '{fname}'.")
