"""
Group database class.

This module provides Group database maneger.
"""

import re
import numpy as np
import sympy as sp
from itertools import product
from sympy.physics.quantum.cg import CG
from joblib import Parallel, delayed


from multipie import PGMultipoleType, SphericalMultipoleType
from multipie.util.util import str_to_sympy
from multipie.util.util_binary import BinaryManager
from multipie.util.util_dict import Dict
from multipie.util.util_gram_schmidt import gram_schmidt
from multipie.util.util_wyckoff import (
    find_wyckoff_site,
    find_wyckoff_bond,
    create_cell_site,
    create_cell_bond,
    convert_to_bond,
)
from multipie.util.util_tag import TagMultipole, TagSymmetryOperation, TagBasis, TagIrrep
from multipie.util.util_samb import orthogonalize_multipole
from multipie.util.util_atomic_multipole import atomic_multipole_matrix
from multipie.util.util_response_tensor import (
    create_active_dict,
    get_response_tensor_mp,
    P0,
    P1,
    P2,
    P3,
    P4,
    mp_string,
    simplify_tensor,
)

CHOP = 1e-6


def replace_bar(s):
    return re.sub(r"\\bar\{(\d)\}", r"-\1", s).replace("{", "").replace("}", "").replace(r"\rm ", "")


# ==================================================
class Group(dict):
    _info = BinaryManager("info")

    # ==================================================
    def __init__(self, tag, with_opt=False):
        """
        Group database manager.

        Args:
            tag (str or int): group tag.
            with_opt (bool, optional): load additional data ?

        Note:
            - PG tag is PG:1-47 or Schoenflies.
            - SG tag is SG:1-230, Schoenflies, or number (1-230).
            - MPG tag is MPG:HM_ID or HM_ID, HM_ID=(PG.no.ID).
            - MSG tag is MSG:BNS_ID or BNS_ID, BNS_ID=(SG.no).
        """
        if type(tag) == int:
            tag = f"SG:{tag}"
        if tag not in self._info["id"].keys() and tag not in self._info["tag"].keys():
            raise Exception(f"unknown tag, '{tag}'.")

        self._group_dict = {"PG": None, "SG": None, "MPG": None, "MSG": None, "opt": None}

        if tag in self._info["tag"].keys():
            self._id = tag
        else:
            self._id = self._info["id"][tag]
        self._type = self._id.split(":")[0]

        file = f"{self._type}/{self._id}"

        self._group_dict[self._type] = BinaryManager(file)

        if with_opt:
            pg_file = self._group(self._type)["info"].PG
            self._group_dict["opt"] = BinaryManager("PG/" + pg_file + "_opt")

    # ==================================================
    @classmethod
    def global_info(cls):
        """
        Global information.

        Returns:
            - (dict) -- global information.
        """
        return cls._info

    # ==================================================
    @property
    def info(self):
        """
        Group info.

        Returns:
            - (namedtuple) -- group info.
        """
        return self._group_dict[self._type]["info"]

    # ==================================================
    @property
    def group_type(self):
        """
        Group type.

        Returns:
            - (str) -- PG/SG/MPG/MSG.
        """
        return self._type

    # ==================================================
    @property
    def opt(self):
        """
        Optional dict (PG only).

        Returns:
            - (dict) -- optional dict.

        Note:
            - **harmonics_multipole**: multipolar harmonics.
            - **harmonics_irrep**: harmonics grouped by irrep.
            - **representation_matrix**: representation matrix.
        """
        if self.group_type != "PG":
            return
        return self._group_dict["opt"]

    # ==================================================
    def _load_group_binary(self, tp):
        if self._group_dict[tp] is None:
            if tp == "PG":
                id_s = self.info.PG
            elif tp == "SG":
                id_s = self.info.SG
            elif tp == "MPG":
                id_s = self.info.MPG
            elif tp == "MSG":
                id_s = self.info.MSG

            file = f"{tp}/{id_s}"
            self._group_dict[tp] = BinaryManager(file)

    # ==================================================
    def _group(self, tp=None):
        if tp is None:
            tp = self.group_type
        self._load_group_binary(tp)

        return self._group_dict[tp]

    # ==================================================
    @property
    def is_point_group(self):
        """
        Check point group.

        Returns:
            - (bool) -- is point group ?
        """
        return self.group_type in ["PG", "MPG"]

    # ==================================================
    @property
    def is_magnetic_group(self):
        """
        Check magnetic group.

        Returns:
            - (bool) -- is magnetic group ?
        """
        return self.group_type in ["MPG", "MSG"]

    # ==================================================
    @property
    def is_hexagonal_subgroup(self):
        """
        Check hexagonal subgroup.

        Returns:
            - (bool) -- is hexagonal subgroup ?
        """
        return self.info.hexagonal_g

    # ==================================================
    @property
    def symmetry_operation(self):
        """
        Get symmetry operation.

        Returns:
            - (dict) -- symmetry operation.
        """
        return self._group()["symmetry_operation"]

    # ==================================================
    @property
    def wyckoff(self):
        """
        Get wyckoff site (all types) or bond (PG,SG only).

        Returns:
            - (dict) -- wyckoff site/bond.
        """
        return self._group()["wyckoff"]

    # ==================================================
    @property
    def active_multipole(self):
        """
        Get active multipole (in accordance with MPG).

        Returns:
            - (dict) -- active multipole.
        """
        return self._group("MPG")["active_multipole"]

    # ==================================================
    @property
    def harmonics(self):
        """
        Get harmonics dict (in accordance with PG).

        Returns:
            - (dict) -- harmonics dict.

        Note:
            - Harmonics is sorted in order of (Q/G/T/M), s, k, l, Gamma, n, p.
        """
        return self._group("PG")["harmonics"]

    # ==================================================
    @property
    def ID(self):
        """
        Get ID of group.

        Returns:
            - (str) -- ID of group.
        """
        if self.group_type == "PG":
            return self.info.PG.split(":")[1]
        elif self.group_type == "SG":
            return self.info.SG.split(":")[1]
        else:
            return self.info.tag

    # ==================================================
    @property
    def character(self):
        """
        Get character dict (in accordance with PG).

        Returns:
            - (dict) -- character dict.
        """
        return self._group("PG")["character"]

    # ==================================================
    @staticmethod
    def tag_multipole(tag, comp=None, latex=False, superscript="", vector=False, internal=False):
        """
        Tag multipole.

        Args:
            tag (tuple): multipole index (harmonics key with any X).
            comp (int, optional): component.
            latex (bool, optional): LaTeX format ?
            superscript (str, optional): additional superscript in LaTeX.
            vector (bool, optional): vector symbol in LaTeX ?
            internal (bool, optional): with internal variable in LaTeX ?

        Returns:
            - (str or list) -- multipole tag, if comp is None, all components in list.
        """
        dim = {"A": 1, "B": 1, "E": 2, "T": 3}

        mp_type = "point_group"
        X, l, Gamma, n, p, s, k, x = tag
        n = dim[Gamma[0]]

        if latex:
            if comp is None:
                ret = [TagMultipole.latex(tag, mp_type, i + 1, superscript, vector, internal) for i in range(n)]
            else:
                if comp != -1:
                    comp += 1
                ret = TagMultipole.latex(tag, mp_type, comp, superscript, vector, internal)
        else:
            if comp is None:
                ret = [TagMultipole.str(TagMultipole.parse(tag, mp_type, i + 1)) for i in range(n)]
            else:
                ret = TagMultipole.str(TagMultipole.parse(tag, mp_type, comp))

        return ret

    # ==================================================
    @staticmethod
    def tag_symmetry_operation(tag, latex=False):
        """
        Tag symmetry operation (fractional coordinate).

        Args:
            tag (str): symmetry operation tag.
            latex (bool, optional): LaTeX format ?

        Returns:
            - (dict or str) -- info dict (fractional coordinate) or LaTeX.
        """
        if latex:
            d = TagSymmetryOperation.latex(tag)
        else:
            d = TagSymmetryOperation.parse(tag)

        return d

    # ==================================================
    @staticmethod
    def tag_irrep(tag, latex=False):
        """
        Tag irrep.

        Args:
            tag (str): irrep. tag.
            latex (bool, optional): LaTeX format ?

        Returns:
            - (dict or str) -- info dict or LaTeX.
        """
        if latex:
            d = TagIrrep.latex(tag)
        else:
            d = TagIrrep.parse(tag)

        return d

    # ==================================================
    @staticmethod
    def tag_atomic_basis(tag, rank, ket=True, latex=False):
        """
        Tag atomic basis.

        Args:
            tag (str): basis tag.
            rank (int): rank.
            ket (bool, optional): ket or bra ?
            latex (bool, optional): LaTeX format ?

        Returns:
            - (dict or str) -- info dict or LaTeX.
        """
        if latex:
            d = TagBasis.latex(tag, rank, ket)
        else:
            d = TagBasis.parse(tag)
            d["rank"] = rank

        return d

    # ==================================================
    def __str__(self):
        """
        Group tag.

        Returns:
            - (str) -- tag.
        """
        return self.info.tag

    # ==================================================
    def latex(self, detail=False, tp=None):
        """
        Group info. in LaTeX format.

        Args:
            detail (bool, optional): detailed info ?
            tp (str, optional); group type.

        Returns:
            - (str) -- info. in LaTeX format.
        """
        if tp is None:
            tp = self.group_type

        info = self._group(tp)["info"]
        setting = info.setting

        if detail:
            ID = self.ID
            crystal = info.crystal
            sp = r"\quad"
            if tp in ["PG", "SG"]:
                S = "$" + info.schoenflies + "$"
                I = "$" + info.international + "$"
                r = f"{tp} No. {ID}{sp}{S}{sp}{I}"
                if setting != "":
                    r += f"{sp}({setting} setting)"
                r += f"{sp}[ {crystal} ]"
            elif tp == "MPG":
                I = "$" + info.international + "$"
                t = info.type
                r = f"{tp} No. {ID}{sp}{I}"
                if setting != "":
                    r += f"{sp}({setting} setting)"
                r += f"{sp}[ Type {t}, {crystal} ]"
            else:
                I = "$" + info.BNS + "$"
                t = info.type
                r = f"{tp} No. {ID}{sp}{I}"
                if setting != "":
                    r += f"{sp}({setting} setting)"
                r += f"{sp}[ Type {t}, {crystal} ]"
            return r
        else:
            if tp in ["PG", "SG"]:
                return info.schoenflies
            elif tp == "MPG":
                return info.international
            else:
                return info.BNS

    # ==================================================
    def name(self, detail=False, tp=None):
        """
        Group info. in text format.

        Args:
            detail (bool, optional): detailed info ?
            tp (str, optional); group type.

        Returns:
            - (str) -- info. in plain text format.
        """
        if tp is None:
            tp = self.group_type

        info = self._group(tp)["info"]
        setting = info.setting

        if detail:
            ID = self.ID
            crystal = info.crystal
            sp = "  "
            if tp in ["PG", "SG"]:
                S = str(self)
                I = replace_bar(info.international)
                r = f"{tp}#{ID}:{sp}{S}{sp}({I})"
                if setting != "":
                    r += f"{sp}({setting} setting)"
                r += f"{sp}[ {crystal} ]"
            elif tp == "MPG":
                I = replace_bar(info.international)
                t = info.type
                r = f"{tp}#{ID}:{sp}{I}"
                if setting != "":
                    r += f"{sp}({setting} setting)"
                r += f"{sp}[ Type {t}, {crystal} ]"
            else:
                I = replace_bar(info.BNS)
                t = info.type
                r = f"{tp}#{ID}:{sp}{I}{sp}"
                if setting != "":
                    r += f"{sp}({setting} setting)"
                r += f"{sp}[ Type {t}, {crystal} ]"
            return r
        else:
            if tp in ["PG", "SG"]:
                return replace_bar(info.schoenflies).replace("_", "")
            elif tp == "MPG":
                return replace_bar(info.international).replace("_", "")
            else:
                return replace_bar(info.BNS).replace("_", "")

    # ==================================================
    def atomic_basis(self, basis_type):
        """
        Get atomic basis set dict.

        Args:
            basis_type (str): basis type, "jml/lgs/lg".

        Returns:
            - (dict) -- atomic basis in each rank.
        """
        harm_info = self.global_info()["harmonics"]["atomic_basis"]
        if basis_type not in ["jml", "lgs", "lg"]:
            raise Exception(f"unknown basis type, '{basis_type}'.")
        if basis_type in ["jml", "lgs"]:
            basis = harm_info["spinful"][basis_type]
        else:
            basis = harm_info["spinless"][basis_type]

        return basis

    # ==================================================
    def find_wyckoff(self, site_bond):
        """
        Find Wyckoff site (all types) or bond (PG/SG only).

        Args:
            site_bond (ndarray or str): site or bond (vector+center) to find.

        Returns:
            - (str) -- Wyckoff site or bond (return None when not found).
            - (ndarray) -- set of Wyckoff site or bond (vector+center).

        Note:
            - for SG, Wyckoff site or bond in order of (plus_set1, plus_set2, ...).
        """
        if site_bond.count("[") == 2:
            return self.find_wyckoff_bond(site_bond)
        else:
            return self.find_wyckoff_site(site_bond)

    # ==================================================
    def find_wyckoff_site(self, site):
        """
        Find Wyckoff site.

        Args:
            site (ndarray or str): site to find.

        Returns:
            - (str) -- Wyckoff site (return None when not found).
            - (ndarray) -- set of Wyckoff site.

        Note:
            - site string is "[x,y,z]".
            - for SG, Wyckoff site in order of (plus_set1, plus_set2, ...).
        """
        if self.group_type == "MSG":
            return find_wyckoff_site(self._group(), site, msg=True)
        else:
            return find_wyckoff_site(self._group(), site, msg=False)

    # ==================================================
    def find_wyckoff_bond(self, bond):
        """
        Find Wyckoff bond (PG/SG only).

        Args:
            bond (ndarray or str): bond (vector+center) to find.

        Returns:
            - (str) -- Wyckoff bond (return None when not found).
            - (ndarray) -- set of Wyckoff bond (vector+center).

        Note:
            - bond string is "[tail];[head] / [vector]@[center] / [start]:[vector]".
            - for SG, Wyckoff bond in order of (plus_set1, plus_set2, ...).
        """
        if self.group_type in ["MPG", "MSG"]:
            return None

        return find_wyckoff_bond(self._group(), bond)

    # ==================================================
    def atomic_samb(self, basis_type, rank_bra_ket, mask_bra_ket=None):
        """
        Get atomic SAMB Dict (in accordance with PG).

        Args:
            basis_type (str): basis type, "jml/lgs/lg".
            rank_bra_ket (tuple): bra-ket rank, 0-3, (bra,ket).
            mask_bra_ket (tuple, optional): mask orbital index to use, (bra:[int],ket:[int]).

        Returns:
            - (Dict) -- atomic SAMB.

        Note:
            - SAMB is sorted in order of (Q/G/T/M), s, k, Gamma, l, n.
            - mask_bra_ket (bra,ket) can be given as ([int],[int]), ([str],[str]), (ndarray,ndarray) or (str,str).
            - ([int],[int]) orbital indices of tesseral basis for given basis_type.
            - ([str],[str]) orbital names of tesseral basis for given basis_type.
            - (ndarray,ndarray) or (str,str) unitary matrix from JML/Lgs/Lg to n basis, <JML/Lgs/Lg | n>.
        """

        def check_mask(obj):
            if not isinstance(obj, tuple) or len(obj) != 2:
                return None
            bra, ket = obj
            if isinstance(bra, np.ndarray) and isinstance(ket, np.ndarray):
                return "ndarray"
            if isinstance(bra, str) and isinstance(ket, str):
                return "ndarray"
            if isinstance(bra, (list, tuple)) and isinstance(ket, (list, tuple)):
                if not bra and not ket:  # both empty list.
                    return None
                if all(isinstance(x, int) for x in bra) and all(isinstance(x, int) for x in ket):
                    return "int"
                if all(isinstance(x, str) for x in bra) and all(isinstance(x, str) for x in ket):
                    return "str"
            return None

        if basis_type not in ["jml", "lgs", "lg"]:
            raise Exception(f"unknown basis_type, '{basis_type}'.")

        a_samb = self._group("PG")["atomic_samb"]
        bra_rank, ket_rank = rank_bra_ket

        if bra_rank > ket_rank:
            dic = a_samb[basis_type][(ket_rank, bra_rank)]
            dic = Dict(dic.key_type, {tag: (mat.transpose(0, 2, 1).conjugate(), ex) for tag, (mat, ex) in dic.items()})
        else:
            dic = a_samb[basis_type][rank_bra_ket]

        if mask_bra_ket is not None:
            tp = check_mask(mask_bra_ket)
            if tp is None:
                raise Exception(f"unknown mask_bra_ket, '{mask_bra_ket}'.")

            bra_mask, ket_mask = mask_bra_ket
            if tp in ["int", "str"]:
                if tp == "str":
                    bra_set = self.atomic_basis(basis_type)[bra_rank]
                    ket_set = self.atomic_basis(basis_type)[ket_rank]
                    bra_mask = [bra_set.index(t) for t in bra_mask]
                    ket_mask = [ket_set.index(t) for t in ket_mask]
                if len(bra_mask) == 0:
                    bra_mask = list(range(len(self.atomic_basis(basis_type)[bra_rank])))
                if len(ket_mask) == 0:
                    ket_mask = list(range(len(self.atomic_basis(basis_type)[ket_rank])))
                dic = Dict(dic.key_type, {tag: (mat[:, bra_mask][:, :, ket_mask], ex) for tag, (mat, ex) in dic.items()})
            else:
                if isinstance(bra_mask, str):
                    bra_mask = str_to_sympy(bra_mask)
                if isinstance(ket_mask, str):
                    ket_mask = str_to_sympy(ket_mask)

                U_bra = bra_mask.conjugate().T
                U_ket = ket_mask

                dic1 = Dict(PGMultipoleType)
                for idx, (mat, ex) in dic.items():
                    matp = np.vectorize(sp.expand)(np.einsum("ik,akl,lj->aij", U_bra, mat, U_ket, dtype=object))
                    dic1[tuple(idx)] = (matp, ex)
                dic = dic1

            # remove zero matrices
            dic = Dict(dic.key_type, {tag: (mat, ex) for tag, (mat, ex) in dic.items() if not np.all(mat == 0)})

            dic = orthogonalize_multipole(dic, bra_rank == ket_rank)

        return dic

    # ==================================================
    def cluster_samb(self, wp, cluster_type=None):
        """
        Get cluster SAMB Dict (PG/SG only).

        Args:
            wp (str): wyckoff position.
            cluster_type (str): cluster type, "site/bond_s/bond_a/vector/bond".

        Returns:
            - (Dict) -- cluster SAMB.

        Note:
            - SAMB is sorted in order of (Q/G/T/M), Gamma, (k,) l, n, p.
        """
        if wp.count("@") == 0 or cluster_type is None:
            cluster_type = "site"

        if cluster_type not in ["site", "bond_s", "bond_a", "vector", "bond"]:
            raise Exception(f"unknown cluster_type, '{cluster_type}'.")

        if self.group_type in ["MPG", "MSG"]:
            return None

        tp = self.group_type
        c_samb = self._group()["cluster_samb"]
        self._load_group_binary(tp)

        ct = "site" if cluster_type == "site" else "bond_s"
        if wp not in c_samb[ct].keys():
            raise Exception(f"unknown tag, '{wp}'.")

        if cluster_type == "bond":
            s = c_samb["bond_s"][wp]
            a = c_samb["bond_a"][wp]
            a = Dict(
                PGMultipoleType,
                {("T" if idx[0] == "Q" else "M", *idx[1:]): (sp.I * v, ex) for idx, (v, ex) in a.items()},
            )
            dic = Dict(s.key_type, dict(s) | dict(a))
            return dic
        else:
            return c_samb[cluster_type][wp]

    # ==================================================
    def combined_samb(self, a_samb_key, c_samb_key, toroidal_priority=False, **kwargs):
        """
        Get combined SAMB Dict (PG/SG only).

        Args:
            a_samb_key (dict_keys): named key of atomic SAMB, [PGMultipoleType].
            c_samb_key (dict_keys): named key of site or bond cluster SAMB, [PGMultipoleType].
            toroidal_priority (bool, optional): use (G,T) prior to (Q,M) in creation ?
            kwargs (dict, optional): select conditions for multipoles with keywords, ["X", "l", "Gamma", "s"].

        Returns:
            - (Dict) -- combined SAMB.

        Note:
            - SAMB is sorted in order of (Q/G/T/M), s, k, Gamma, l, n, p.
        """
        if self.group_type in ["MPG", "MSG"]:
            raise None

        harmonics = self.harmonics
        keys = list(harmonics.named_keys())
        replace_map = {"Q": "T", "G": "M"}
        Z_key = keys + [tag._replace(X=replace_map[tag.X]) for tag in keys if tag.X in replace_map]

        head_list = ["G", "Q", "T", "M"] if toroidal_priority else ["Q", "G", "M", "T"]
        Z_key = Dict.sort_key(
            Dict.select_key(Z_key, **kwargs), PGMultipoleType, ("X", head_list), "s", "k", "Gamma", "l", "n", "p"
        )

        unique = lambda seq: list(dict.fromkeys(seq))
        Z_s = kwargs.pop("s", [0, 1])
        Z_Gamma_list = unique(tag.Gamma for tag in Z_key)
        X_l_s_k_p_list = unique((t.X, t.l, t.s, t.k, t.p) for t in a_samb_key if t.s in Z_s)
        Y_l_p_list = unique((t.X, t.l, t.p) for t in c_samb_key)

        t_val = {"Q": 1, "G": 1, "T": -1, "M": -1}

        dic = Dict(PGMultipoleType)

        key_of = lambda tg, i: (tg.X, tg.l, tg.Gamma, tg.n, i, tg.s, tg.k, tg.x)

        for (X, l1, s, k, p1), (Y, l2, p2) in product(X_l_s_k_p_list, Y_l_p_list):
            tag1_list = Dict.select_key(a_samb_key, X=X, l=l1, s=s, k=k, p=p1)
            tag2_list = Dict.select_key(c_samb_key, X=Y, l=l2, p=p2)
            tag1_tag2_list = [(tag1, tag2) for tag1 in tag1_list for tag2 in tag2_list]

            Z_list = [Z for Z, t in t_val.items() if t == t_val[X] * t_val[Y]]
            l_list = list(range(abs(l1 - l2), l1 + l2 + 1))
            for Gamma in Z_Gamma_list:
                tag_list, coeff_list, tag12_list = [], [], []
                for tag in Dict.select_key(Z_key, **{"X": Z_list, "l": l_list, "Gamma": Gamma}):
                    coeff, tag12 = [], []
                    for tag1, tag2 in tag1_tag2_list:
                        cg = self.cg(tag1, tag2, tag)
                        if cg is None:
                            continue
                        dim1, dim2 = TagIrrep.parse(tag1.Gamma)["dimension"], TagIrrep.parse(tag2.Gamma)["dimension"]
                        coeff += [cg[g1, g2, :] for g1 in range(dim1) for g2 in range(dim2)]
                        tag12 += [(tuple(tag1), g1, tuple(tag2), g2) for g1 in range(dim1) for g2 in range(dim2)]
                    if coeff:
                        coeff_list.append(coeff)
                        tag_list.append(tag)
                        tag12_list.append(tag12)

                # orthogonalize
                if not coeff_list:
                    continue

                coeff_list, idx = gram_schmidt(coeff_list)
                tag_list, tag12_list = [tag_list[i] for i in idx], [tag12_list[i] for i in idx]
                for tag, tag12, coeff in zip(tag_list, tag12_list, coeff_list):
                    ex = list(harmonics[tag._replace(X="Q" if tag.X == "T" else "G" if tag.X == "M" else tag.X)][0])
                    tag = tag._replace(**{"s": s, "k": k})

                    dim = TagIrrep.parse(tag.Gamma)["dimension"]
                    lst = [[(cg[g], *entry) for entry, cg in zip(tag12, coeff) if cg[g] != 0] for g in range(dim)]

                    cnt = len(dic.select(X=tag.X, l=tag.l, Gamma=tag.Gamma, n=tag.n, s=tag.s, k=tag.k, x=tag.x))
                    if cnt == 0:
                        tag = key_of(tag, -1)
                    elif cnt == 1:
                        base, tag1, tag = key_of(tag, -1), key_of(tag, 1), key_of(tag, 2)
                        dic[tag1] = dic[base]
                        del dic[base]
                    else:
                        tag = key_of(tag, cnt + 1)

                    dic[tag] = (lst, ex)

        dic = dic.sort("Gamma", ("X", ["Q", "G", "M", "T"]), "s", "k", "l", "n", "p")

        return dic

    # ==================================================
    def multipole_cluster_samb(self, X, l, site_bond):
        """
        Create atomic multipole and cluster combined SAMB (PG/SG only).

        Args:
            X (str): head of atomic multipole.
            l (int): rank of atomic multipole.
            site_bond (str): site or bond.

        Returns:
            - (Dict) -- combined SAMB.
            - (str) -- wyckoff site/bond.
            - (ndarray) -- cluster sites or center of bonds.
        """
        t_inv = {"Q": "T", "G": "M", "T": "Q", "M": "G"}
        Xr = t_inv[X] if X in ["T", "M"] else X
        a_samb_key = self.harmonics.select(X=Xr, l=l).named_keys()
        if X in ["T", "M"]:
            a_samb_key = [tag._replace(X=t_inv[tag.X]) for tag in a_samb_key]

        wp, cluster = self.find_wyckoff(site_bond)
        if "@" in wp:
            sites = cluster[:, 3:6]
            c_samb_key = self.cluster_samb(wp, "bond").named_keys()
        else:
            sites = cluster
            c_samb_key = self.cluster_samb(wp).named_keys()
        samb = self.combined_samb(a_samb_key, c_samb_key)

        return samb, wp, sites

    # ==================================================
    def combined_object(self, wp, samb_type, samb):
        """
        Create combined object for given atomic multipole (X,l) and cluster (PG/SG only).

        Args:
            wp (str): wyckoff site/bond of cluster.
            samb_type (str): combined SAMB type.
            samb (Dict): combined SAMB dict (obtained by multipole_cluster_samb).

        Returns:
            - (ndarray) -- combined object in terms of (xyz) polynomial for each cluster site.
        """
        sgn = {"Q": 1, "G": 1, "T": -1, "M": -1}
        tp = "bond" if "@" in wp else "site"
        c_samb = self.cluster_samb(wp, tp)
        d = len(next(iter(c_samb.values()))[0][0])

        obj = np.full(d, sp.S(0))
        for coeff, a_key, a_comp, c_key, c_comp in samb:
            s = 1 if sgn[samb_type] * sgn[a_key[0]] == 1 else -sp.I
            a_key = (a_key[0].replace("T", "Q").replace("M", "G"), *a_key[1:])
            a_val = self.harmonics[a_key][0][a_comp]
            c_val = c_samb[c_key][0][c_comp]
            obj += s * coeff * a_val * c_val

        return obj

    # ==================================================
    def response_tensor_all(self, X):
        """
        Create all response tensors up to rank 4 (in accordance with MPG).

        Args:
            X (str): tensor type.

        Returns:
            - (dict) -- response tensors.
        """
        c = [1, 2, 3]
        v = [(1, 1), (2, 2), (3, 3), (2, 3), (3, 1), (1, 2)]
        u = [(1, 1), (2, 2), (3, 3), (2, 3), (3, 1), (1, 2), (3, 2), (1, 3), (2, 1)]
        a = [(2, 3), (3, 1), (1, 2)]
        t = [
            (1, 1, 1),
            (2, 2, 2),
            (3, 3, 3),
            (2, 2, 3),
            (3, 3, 1),
            (1, 1, 2),
            (2, 3, 3),
            (3, 1, 1),
            (1, 2, 2),
            (1, 2, 3),
        ]

        d = {}

        # rank 0.
        m = np.zeros((1, 1), dtype=object)
        m[0, 0] = (X, self.response_tensor(X, 0))
        d[(X, 0, "")] = simplify_tensor(m)
        # rank 1.
        m = np.zeros((1, 3), dtype=object)
        for ii, i in enumerate(c):
            m[0, ii] = (mp_string([i], X, False), self.response_tensor(X, 1, i))
        d[(X, 1, "")] = simplify_tensor(m)
        # rank 2.
        ms, ma = np.zeros((3, 3), dtype=object), np.zeros((3, 3), dtype=object)
        for ii, i in enumerate(c):
            for jj, j in enumerate(c):
                ms[ii, jj] = (mp_string([i, j], X, False), self.response_tensor(X, 2, (i, j), "s"))
                ma[ii, jj] = (mp_string([i, j], X, False), self.response_tensor(X, 2, (i, j), "a"))
        d[(X, 2, "s")] = simplify_tensor(ms)
        d[(X, 2, "a")] = simplify_tensor(ma)
        # rank 3.
        m = np.zeros((6, 3), dtype=object)
        for ii, i in enumerate(v):
            for jj, j in enumerate(c):
                m[ii, jj] = (mp_string([*i, j], X, False), self.response_tensor(X, 3, (*i, j), "s"))
        d[(X, 3, "s")] = simplify_tensor(m)
        m = np.zeros((3, 3), dtype=object)
        for ii, i in enumerate(a):
            for jj, j in enumerate(c):
                m[ii, jj] = (mp_string([*i, j], X, False), self.response_tensor(X, 3, (*i, j), "a"))
        d[(X, 3, "a")] = simplify_tensor(m)
        # rank 4.
        m = np.zeros((6, 6), dtype=object)
        for ii, i in enumerate(v):
            for jj, j in enumerate(v):
                m[ii, jj] = (mp_string([*i, *j], X, False), self.response_tensor(X, 4, (*i, *j), "sss"))
        d[(X, 4, "sss")] = simplify_tensor(m)
        m = np.zeros((6, 6), dtype=object)
        for ii, i in enumerate(v):
            for jj, j in enumerate(v):
                m[ii, jj] = (mp_string([*i, *j], X, False), self.response_tensor(X, 4, (*i, *j), "ssa"))
        d[(X, 4, "ssa")] = simplify_tensor(m)
        m = np.zeros((3, 3), dtype=object)
        for ii, i in enumerate(a):
            for jj, j in enumerate(a):
                m[ii, jj] = (mp_string([*i, *j], X, False), self.response_tensor(X, 4, (*i, *j), "aas"))
        d[(X, 4, "aas")] = simplify_tensor(m)
        m = np.zeros((3, 3), dtype=object)
        for ii, i in enumerate(a):
            for jj, j in enumerate(a):
                m[ii, jj] = (mp_string([*i, *j], X, False), self.response_tensor(X, 4, (*i, *j), "aaa"))
        d[(X, 4, "aaa")] = simplify_tensor(m)
        m = np.zeros((6, 3), dtype=object)
        for ii, i in enumerate(v):
            for jj, j in enumerate(a):
                m[ii, jj] = (mp_string([*i, *j], X, False), self.response_tensor(X, 4, (*i, *j), "sa"))
        d[(X, 4, "sa")] = simplify_tensor(m)
        m = np.zeros((3, 6), dtype=object)
        for ii, i in enumerate(a):
            for jj, j in enumerate(v):
                m[ii, jj] = (mp_string([*i, *j], X, False), self.response_tensor(X, 4, (*i, *j), "as"))
        d[(X, 4, "as")] = simplify_tensor(m)
        m = np.zeros((6, 9), dtype=object)
        for ii, i in enumerate(v):
            for jj, j in enumerate(u):
                m[ii, jj] = (mp_string([*i, *j], X, False), self.response_tensor(X, 4, (*i, *j), "s"))
        d[(X, 4, "s")] = simplify_tensor(m)
        m = np.zeros((3, 9), dtype=object)
        for ii, i in enumerate(a):
            for jj, j in enumerate(u):
                m[ii, jj] = (mp_string([*i, *j], X, False), self.response_tensor(X, 4, (*i, *j), "a"))
        d[(X, 4, "a")] = simplify_tensor(m)
        m = np.zeros((10, 3), dtype=object)
        for ii, i in enumerate(t):
            for jj, j in enumerate(c):
                m[ii, jj] = (mp_string([*i, j], X, False), self.response_tensor(X, 4, (*i, j), "t"))
        d[(X, 4, "t")] = simplify_tensor(m)

        return d

    # ==================================================
    def response_tensor(self, X, rank, comp=None, opt=None):
        """
        Get response tensor (in accordance with MPG).

        Args:
            X (str): tensor type, "Q/G/T/M".
            rank (int): tensor rank, (0,1,2,3,4).
            comp (int or tuple, optional): tensor component, index=(1:x,2:y,3:z).
            opt (str, optional): tensor type option.

        Returns:
            - (sympy) -- tensor component.

        Note:
            - tensor type, None (rank 0,1), "s/a" (rank 2,3), "sss/ssa/aas/aaa/sa/as/s/a/t" (rank 4).
        """
        hexagonal = "hexagonal" if self.is_hexagonal_subgroup else "cubic"
        cartesian_mp = self.global_info()["response_tensor"]["cartesian_multipole"][hexagonal]
        active_dict = create_active_dict(self._group("MPG")["active_multipole"], cartesian_mp)

        axial_tensor = X in ["G", "M"]
        magnetic_tensor = X in ["T", "M"]

        if rank == 0:
            p = P0()
        elif rank == 1:
            p = P1(comp)
        elif rank == 2:
            p = P2(*comp, opt)
        elif rank == 3:
            p = P3(*comp, opt)
        elif rank == 4:
            p = P4(*comp, opt)
        else:
            raise KeyError(f"invalid rank, {rank}.")

        ex = get_response_tensor_mp(p, active_dict, axial_tensor, magnetic_tensor)
        return ex

    # ==================================================
    def transform_vector(self, vector, X="Q", cartesian=False):
        """
        Transform vector.

        Args:
            vector (ndarray or str): vector.
            X (str, optional): vector type, "Q/G/T/M".
            cartesian (bool, optional): vector in cartesian coordinate ?

        Returns:
            - (ndarray) -- transformed vectors in order of symmetry operations.
        """
        if type(vector) == str:
            vector = str_to_sympy(vector)

        so = "cartesian" if cartesian else "fractional"
        so = self.symmetry_operation[so][:, 0:3, 0:3]
        d = self.symmetry_operation["det"]
        t = self.symmetry_operation.get("tr_sign", [1] * len(d))

        s = so @ vector

        if X in ["G", "M"]:
            s = np.asarray([si * di for si, di in zip(s, d)])
        if X in ["T", "M"]:
            s = np.asarray([si * ti for si, ti in zip(s, t)])

        return s

    # ==================================================
    def transform_multipole(self, multipole, X="Q"):
        """
        Transform multipole.

        Args:
            multipole (str): multipole in terms of (x,y,z,r) in cartesian coordinate.
            X (str, optional): vector type, "Q/G/T/M".

        Returns:
            - (ndarray) -- transformed multipoles in order of symmetry operations.
        """
        x, y, z = sp.symbols("x y z", real=True)
        v = np.array([x, y, z])

        so = self.symmetry_operation["cartesian"][:, 0:3, 0:3].transpose(0, 2, 1)  # transpose for each matrix.
        d = self.symmetry_operation["det"]
        t = self.symmetry_operation.get("tr_sign", [1] * len(d))
        vt = so @ v

        multipole = multipole.replace("r", "sqrt(x**2+y**2+z**2)")
        multipole = sp.S(str_to_sympy(multipole))

        s = np.asarray([multipole.subs({"x": vi[0], "y": vi[1], "z": vi[2]}, simultaneous=True) for vi in vt]).reshape(-1)

        if X in ["G", "M"]:
            s = np.asarray([si * di for si, di in zip(s, d)])
        if X in ["T", "M"]:
            s = np.asarray([si * ti for si, ti in zip(s, t)])

        return s

    # ==================================================
    def create_cell_site(self, site):
        """
        Create sites in unit cell (sorted and with plus set).

        Args:
            site (ndarray or str): site.

        Returns:
            - (ndarray) -- sites (sorted) in unit cell.
            - (list) -- SO mapping (index from 1).
            - (str) -- Wyckoff site or None.

        Note:
            - site string is "[x,y,z]".
            - if group is PG/SG, return Wyckoff site as well.
        """
        sites, mapping = create_cell_site(self._group(), site)
        wp = self.find_wyckoff_site(site)[0]
        return sites, mapping, wp

    # ==================================================
    def create_cell_bond(self, bond):
        """
        Create bonds in unit cell (sorted, nondirectional and with plus set).

        Args:
            bond (ndarray or str): bond.

        Returns:
            - (ndarray) -- bonds (sorted) in unit cell.
            - (list) -- SO mapping (index from 1).
            - (str) -- Wyckoff bond or None.

        Note:
            - bond string is "[tail];[head] / [vector]@[center] / [start]:[vector]".
            - if group is PG/SG, return Wyckoff bond as well.
        """
        bonds, mapping = create_cell_bond(self._group(), bond)
        wp = self.find_wyckoff_bond(bond)[0] if self.group_type in ["PG", "SG"] else None
        return bonds, mapping, wp

    # ==================================================
    def create_cell_vector(self, vector_site, X="Q", average=True, cartesian=False):
        """
        Create trasformed vector in unit cell.

        Args:
            vector_site (str): vector # site/bond. site/bond in fractional coordinate.
            X (str, optional): vector type, "Q/G/T/M".
            average (bool, optional): average at same site ?
            cartesian (bool, optional): vector in cartesian coordinate ?

        Returns:
            - (ndarray) -- transformed vectors at each unique site, [[va, vb],[vc, vd],...].
            - (ndarray) -- transformed unique sites (sorted).
            - (list) -- SO mapping.
            - (str) -- Wyckoff site (for bond center).

        Note:
            - if average is True, e.g., [va, vb] is averaged as (va+vb)/2.
            - if group is PG/SG, return Wyckoff site as well.
            - result in order of plus_set for SG.
        """
        if type(vector_site) == str:
            vector, site_bond = vector_site.split("#")
            vector = str_to_sympy(vector)
            if "@" in site_bond or ";" in site_bond or ":" in site_bond:
                site_bond = convert_to_bond(site_bond)
                site = site_bond[3:6]
            else:
                site = str_to_sympy(site_bond)
        else:
            raise ValueError(f"invalid format, {vector_site}")

        v = self.transform_vector(vector, X, cartesian)
        sites, mp, wp = self.create_cell_site(site)
        v = np.asarray([[v[i - 1].tolist() for i in lst] for lst in mp])

        if average:
            v = np.asarray([np.mean(vl, axis=0) for vl in v])

        return v, sites, mp, wp

    # ==================================================
    def create_cell_multipole(self, multipole_site, X="Q", average=True):
        """
        Create transformed multipole in unit cell.

        Args:
            multipole_site (str): multipole # site/bond. multipole (site/bond) in cartesian (fractional) coordinate.
            X (str, optional): vector type, "Q/G/T/M".
            average (bool, optional): average at same site ?

        Returns:
            - (ndarray) -- transformed multipoles at each unique site, [[ma, mb],[mc, md],...].
            - (ndarray) -- transformed unique sites (sorted).
            - (list) -- SO mapping.
            - (str) -- Wyckoff site (for bond center).

        Note:
            - if average is True, e.g., [ma, mb] is averaged as (ma+mb)/2.
            - if group is PG/SG, return Wyckoff site as well.
            - result in order of plus_set for SG.
        """
        if type(multipole_site) == str:
            multipole, site_bond = multipole_site.split("#")
            if "@" in site_bond or ";" in site_bond or ":" in site_bond:
                site_bond = convert_to_bond(site_bond)
                site = site_bond[3:6]
            else:
                site = str_to_sympy(site_bond)
        else:
            raise ValueError(f"invalid format, {multipole_site}")

        v = self.transform_multipole(multipole, X)
        sites, mp, wp = self.create_cell_site(site)
        v = np.asarray([[sp.S(v[i - 1]) for i in lst] for lst in mp])

        if average:
            v = [np.mean(vl) for vl in v]

        return v, sites, mp, wp

    # ==================================================
    def cg(self, tag1, tag2, tag):
        """
        Clebsch-Gordan (CG) coefficient, < tag1; tag2 | tag > (PG/SG only).

        Args:
            tag1 (tuple): multipole tag1.
            tag2 (tuple): multipole tag2.
            tag (tuple): multipole tag.

        Returns:
            - (ndarray) -- CG coefficient, CG[gamma1][gamma2][gamma].
        """
        if self.group_type in ["MPG", "MSG"]:
            return None

        harmonics = self.harmonics
        t_even = {"Q": "Q", "T": "Q", "G": "G", "M": "G"}
        t_val = {"Q": 1, "G": 1, "T": -1, "M": -1}
        p_val = {"Q": 0, "T": 0, "G": 1, "M": 1}

        t1, p1, l1 = t_val[tag1[0]], p_val[tag1[0]], tag1[1]
        t2, p2, l2 = t_val[tag2[0]], p_val[tag2[0]], tag2[1]
        t, p, l = t_val[tag[0]], p_val[tag[0]], tag[1]

        tag1 = [t_even[tag1[0]], *tag1[1:4], -1, 0, 0, "q"]
        tag2 = [t_even[tag2[0]], *tag2[1:4], -1, 0, 0, "q"]
        tag = [t_even[tag[0]], *tag[1:4], -1, 0, 0, "q"]
        u1 = harmonics[tag1][1].T
        u2 = harmonics[tag2][1].T
        u = harmonics[tag][1].T

        if t1 * t2 * t != 1 or (l1 + l2 - l + p1 + p2 - p) % 2 != 0 or not (abs(l1 - l2) <= l <= l1 + l2):
            return None

        ml1 = list(range(l1, -l1 - 1, -1))
        ml2 = list(range(l2, -l2 - 1, -1))
        ml = list(range(l, -l - 1, -1))
        phase = (-sp.I) ** (l1 + l2 - l)
        m_to_idx = {m: k for k, m in enumerate(ml)}

        cg_array = np.full((len(u1), len(u2), len(u)), sp.S(0), dtype=object)

        for g1, g2, g in product(range(len(u1)), range(len(u2)), range(len(u))):
            s = sp.S(0)
            for i, j in product(range(len(ml1)), range(len(ml2))):
                m1, m2 = ml1[i], ml2[j]
                m = m1 + m2
                if m not in m_to_idx:
                    continue
                k = m_to_idx[m]
                c = CG(l1, m1, l2, m2, l, m).doit()
                if c:
                    s += phase * c * sp.conjugate(u1[g1, i] * u2[g2, j]) * u[g, k]
            if s != 0:
                cg_array[g1, g2, g] = s

        return cg_array

    # ==================================================
    def create_atomic_samb_L(self, basis_rank):
        """
        Create atomic multipoles X_{lm}^{(0)}(0) for given LM (in accordance with PG).

        Args:
            basis_rank (int): bra-ket rank (integer only).

        Returns:
            - (dict) -- atomic multipoles, Dict[SphericalMultipoleType, ndarray].
            - (list) -- basis.
        """
        if self.group_type in ["MPG", "MSG"]:
            return None

        verbose = 0  # 0, 1, 10.

        L = basis_rank
        basis = [f"({L},{m},1/2)" for m in reversed(range(-L, L + 1))]

        def proc(xlmsk):
            (X, l, m, s, k) = xlmsk
            am = sp.Matrix.zeros(len(basis))
            for i, b1 in enumerate(basis):
                for j, b2 in enumerate(basis):
                    b1 = sp.sympify(b1)
                    b2 = sp.sympify(b2)
                    v = atomic_multipole_matrix(X, l, m, s, k, b1, b2, "lms")
                    am[i, j] = sp.simplify(v)

            return xlmsk, am

        max_l = 2 * L
        Xlmsk_list = [(X, l, m, 0, 0) for X in ["Q", "M"] for l in range(max_l + 1) for m in reversed(range(-l, l + 1))]

        sub = Parallel(n_jobs=-1, verbose=verbose)([delayed(proc)(xlmsk) for xlmsk in Xlmsk_list])

        Xlmsk = {tag: am for tag, am in sub if not am.is_zero_matrix}

        Xlsk = {}
        for (X, l, m, s, k), am in Xlmsk.items():
            Xlsk[(X, l, s, k)] = Xlsk.get((X, l, s, k), []) + [am.tolist()]

        samb = Dict(SphericalMultipoleType)
        for (X, l, s, k), mat in Xlsk.items():
            samb[(X, l, s, k, "q")] = np.asarray(mat)
        samb = samb.sort(("X", ["Q", "G", "T", "M"]), "s", "k", "l")

        basis = [",".join(i.split(",")[:2]) + ")" for i in basis]

        dic1 = {}
        for (X, l, irrep, n, p, _, _, _), (ex, U, _) in self.harmonics.items():
            for (head, _, s, k, _), mat in samb.select(l=l).items():
                if (X == "Q" and head in ["Q", "T"]) or (X == "G" and head in ["G", "M"]):
                    gx = np.vectorize(sp.expand)(np.einsum("mg,mij->gij", U, mat))
                    dic1[(head, l, irrep, s, k)] = dic1.get((head, l, irrep, s, k), []) + [(gx, ex)]

        # renumber multiplicity.
        dic = Dict(PGMultipoleType)
        for (head, l, irrep, s, k), v in dic1.items():
            if len(v) == 1:
                vi = v[0]
                dic[(head, l, irrep, -1, -1, s, k, "q")] = vi
            else:
                for no, vi in enumerate(v):
                    dic[(head, l, irrep, no + 1, -1, s, k, "q")] = vi

        return dic, basis
