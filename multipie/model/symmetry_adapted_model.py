"""
SymmetryAdaptedModel manages symmetry adapted multipole basis set for cluster or crystal system.
"""
from itertools import product
import sympy as sp
from gcoreutils.nsarray import NSArray
from multipie.tag.tag_list import TagList
from multipie.multipole.util.multipole_util import matrix_to_dict
from multipie.multipole.util.structure_samb_util import cs_list
from multipie.symmetry_operation.util.symmetry_operation_util import to_primitive
from multipie.model.util.symmetry_adapted_model_util import (
    create_atomic_samb_set,
    create_cluster_samb_set,
    create_uniform_samb_set,
    create_structure_samb_set,
    create_z_samb_set,
    create_zk_samb_set,
    redefine_index,
)
from multipie import __version__


# ==================================================
header_str = """
=== SAMB (* only for crystal with fourier transform) ===
- info
    - atomic : { "M_#" : ["amp_#"] }
    - site_cluster : { "S_#" : ["smp_#"] }
    - bond_cluster : { "B_#" : ["bmp_#"] }
    - uniform : { "S_#"/"B_#" : ["ump_#"] }
    - structure* : { "B_#" : ["kmp_#"] }
    - Z : { ("M_#", "S_#"/"B_#") : ["z_#"] }
    - version : MultiPie version
    - harmonics : { head : [TagMultipole] }

- data
    - atomic : { "amp_#" : ( TagMultipole, shape, [(i, j, matrix element)] ) }
    - site_cluster : { "smp_#" : ( TagMultipole, [vector component] ) }
    - bond_cluster : { "bmp_#" : ( TagMultipole, [vector component] ) }
    - uniform* : { "ump_#" : ( TagMultipole, shape, [(i, j, matrix element)] ) }
    - structure* : { "kmp_#" : (TagMultipole, "structure factor") }
    - Z : {"z_#" : ( TagMultipole, [(coeff, "amp_#", "smp_#"/"bmp_#/ump_#")] ) }
    - Zk* : {"z_#" : ( TagMultipole, [(coeff, "amp_#", "ump_#", "kmp_#")] ) }
"""


# ==================================================
matrix_header_str = """
== SAMB in full matrix form in real space (* only for crystal) ===
- model : model name.
- molecule : molecule or crystal ?
- group : (tag, detailed str)
- dimension : dimension of full matrix
- ket : ket basis list, orbital@site
- version : MultiPie version
- k_point* : representative k points
- k_path* : high-symmetry line in k space
- cell_site : { name_idx(pset): (position, SOs) }
- A* : transform matrix, [a1,a2,a3]
- matrix : { "z_#": "matrix"}
"""

# ==================================================
matrix_header_k_str = """
== SAMB in full matrix form with fourier transform ===
- model : model name.
- molecule : molecule or crystal ?
- group : (tag, detailed str)
- dimension : dimension of full matrix
- ket : ket basis list, orbital@site
- version : MultiPie version
- k_point : representative k points
- k_path : high-symmetry line in k space
- A : transform matrix, [a1,a2,a3]
- bond : { "bond_#": "vector" }
- matrix : { "z_#": "matrix" }
"""


# ==================================================
def _run(msg, mpm, func, *args, **kwargs):
    mpm.log(f"  * {msg} ... ", None, end="")
    mpm.set_stamp()
    res = func(*args, **kwargs)
    mpm.log("done")

    return res


# ==================================================
class SymmetryAdaptedModel(dict):
    """
    symmetry adapted multipole basis set for molecule or crystal.
    """

    # Attributes:
    #     _model (dict): molecule or crystal model information.
    #     _mpm (MultiPieManager): MultiPie manager.
    #     _fourier_transform (bool): fourier transform combined multipole?
    #     _parallel (bool): use parallel code.
    #     _verbose (bool): verbose parallel info.
    #
    # ==================================================
    def __init__(self, model, mpm, head):
        """
        initialize the class.

        Args:
            model (dict): model information of cluster or crystal system.
            mpm (MultiPieManager): MultiPie manager.
            head (list or str): heads of SAMB to be created.
        """
        super().__init__()

        # model
        self._model = model
        self._mpm = mpm
        self._parallel = self._mpm.parallel
        self._verbose = self._mpm.verbose

        toroidal_priority = model["info"]["generate"]["toroidal_priority"]
        irrep = model["info"]["generate"]["irrep"]

        site = model["data"]["site"]
        cluster_site = model["data"]["cluster_site"]

        bond = model["data"]["bond"]
        cluster_bond = model["data"]["cluster_bond"]

        cluster_atomic = model["data"]["cluster_atomic"]
        atomic_braket = model["data"]["atomic_braket"]
        alias = model["name"]["alias"]

        spinful = model["info"]["spinful"]
        is_phonon = model["info"]["generate"]["model_type"] == "phonon"

        g = mpm.group
        pg = mpm.point_group

        # atomic
        func = create_atomic_samb_set
        args = [pg, atomic_braket, spinful, is_phonon, False]
        atomic_info, atomic_data = _run("atomic", self._mpm, func, *args)

        # site/bond-cluster
        rep_site_dict = {cluster_tag: site[site_list[0]][0] for cluster_tag, site_list in cluster_site.items()}
        rep_bond_dict = {cluster_tag: bond[bond_list[0]][0] for cluster_tag, bond_list in cluster_bond.items()}

        func = create_cluster_samb_set
        args = [g, rep_site_dict, rep_bond_dict, False]
        site_cluster_info, site_cluster_data, bond_cluster_info, bond_cluster_data = _run(
            "site/bond-cluster", self._mpm, func, *args
        )

        # uniform
        if self._model["info"]["molecule"] or self._model["info"]["generate"].get("fourier_transform", False):
            cluster_samb_set = {
                S_i: [site_cluster_data[smp_i] for smp_i in lst] for S_i, lst in site_cluster_info.items()
            }
            cluster_samb_set |= {
                B_i: [bond_cluster_data[bmp_i] for bmp_i in lst] for B_i, lst in bond_cluster_info.items()
            }
            braket_indexes_dict = {
                cluster_tag: [site[site_tag][2] for site_tag in site_list]
                for cluster_tag, site_list in cluster_site.items()
            }
            braket_indexes_dict |= {
                cluster_tag: [bond[bond_tag][2] for bond_tag in bond_list]
                for cluster_tag, bond_list in cluster_bond.items()
            }
            dim = len(site)

            func = create_uniform_samb_set
            args = [cluster_samb_set, braket_indexes_dict, dim]
            uniform_info, uniform_data = _run("uniform", self._mpm, func, *args)

        # Z
        x_tag_dict = {(M_i, atomic_data[amp_i][0]): amp_i for M_i, lst in atomic_info.items() for amp_i in lst}
        if self._model["info"]["molecule"]:
            y_tag_dict = {(SB_i, uniform_data[ump_i][0]): ump_i for SB_i, lst in uniform_info.items() for ump_i in lst}
        else:
            y_tag_dict = {
                (S_i, site_cluster_data[smp_i][0]): smp_i for S_i, lst in site_cluster_info.items() for smp_i in lst
            }
            y_tag_dict |= {
                (B_i, bond_cluster_data[bmp_i][0]): bmp_i for B_i, lst in bond_cluster_info.items() for bmp_i in lst
            }

        cluster_site_bond = cluster_site | cluster_bond
        braket_site_no_dict = {S_i: site[site_list[0]][2] for S_i, site_list in cluster_site.items()}
        braket_site_no_dict |= {B_i: bond[bond_list[0]][2] for B_i, bond_list in cluster_bond.items()}
        M_SB_list = []
        for SB_i in cluster_site_bond.keys():
            for _, _, M_i in cluster_atomic[braket_site_no_dict[SB_i]]:
                M_SB_list.append((M_i, SB_i))

        func = create_z_samb_set
        args = [g, x_tag_dict, y_tag_dict, M_SB_list, atomic_braket, alias, toroidal_priority, False]
        z_info, z_data = _run("Z", self._mpm, func, *args, head=head, irrep=irrep)

        # info and data
        self["info"] = {
            "atomic": atomic_info,
            "site_cluster": site_cluster_info,
            "bond_cluster": bond_cluster_info,
        }
        self["data"] = {
            "atomic": atomic_data,
            "site_cluster": site_cluster_data,
            "bond_cluster": bond_cluster_data,
        }

        if self._model["info"]["molecule"]:
            self["info"] |= {"uniform": uniform_info}
            self["data"] |= {"uniform": uniform_data}

        self["info"] |= {"Z": z_info, "version": __version__}
        self["data"] |= {"Z": z_data}

        if not self._model["info"]["molecule"] and self._model["info"]["generate"]["fourier_transform"]:
            self["info"] |= {"uniform": uniform_info}
            self["data"] |= {"uniform": uniform_data}
            self.fourier_transform()

        if self._model["info"]["option"]["minimal_samb"]:
            self._minimal_samb()

        self["info"]["harmonics"] = self._harmonics_info()

    # ==================================================
    def _harmonics_info(self):
        """
        convert combined multipoles from dictionary to numerical matrix form (real-space).

        Args:
            fmt (str, optional): sympy/value.

        Returns:
            dict: matrix form of combined multipoles, z_dict/z_full.

        Note:
            z_dict = { "z_#" : {(n1, n2, n3, a, b): matrix element}, R=(n1,n2,n3) is a lattice vector.
        """
        atomic_data = self["data"]["atomic"]
        site_cluster_data = self["data"]["site_cluster"]
        bond_cluster_data = self["data"]["bond_cluster"]
        uniform_data = self["data"].get("uniform", None)
        structure_data = self["data"].get("structure", None)

        head_list = ["Q", "G", "M", "T"]
        harmonics_info = {"Q": [], "G": []}
        for X in head_list:
            data = atomic_data | site_cluster_data | bond_cluster_data
            if uniform_data:
                data |= uniform_data
            if structure_data:
                data |= structure_data
            tag_list = TagList([tag for _, (tag, _) in data.items()]).select(head=X)
            harmonics_list = sorted(set([tag.to_harmonics() for tag in tag_list]))

            X = "G" if X == "M" else X
            X = "Q" if X == "T" else X
            harmonics_info[X] += [str(htag) for htag in harmonics_list]

        harmonics_info = {"Q": sorted(set(harmonics_info["Q"])), "G": sorted(set(harmonics_info["G"]))}

        return harmonics_info

    # ==================================================
    def _minimal_samb(self):
        """
        delete unused samb and redefine the serial number of multipoles.
        """
        molecule = self._model["info"]["molecule"]

        atomic_info = self["info"]["atomic"]
        site_cluster_info = self["info"]["site_cluster"]
        bond_cluster_info = self["info"]["bond_cluster"]
        z_info = self["info"]["Z"]
        uniform_info = self["info"].get("uniform", None)
        structure_info = self["info"].get("structure", None)

        atomic_data = self["data"]["atomic"]
        site_cluster_data = self["data"]["site_cluster"]
        bond_cluster_data = self["data"]["bond_cluster"]
        z_data = self["data"]["Z"]
        uniform_data = self["data"].get("uniform", None)
        structure_data = self["data"].get("structure", None)
        zk_data = self["data"].get("Zk", None)

        if site_cluster_data:
            site_cluster_info, site_cluster_data, z_data = redefine_index(
                site_cluster_info, site_cluster_data, z_data, molecule=molecule
            )

        if site_cluster_info:
            init_idx = int(list(site_cluster_info.values())[-1][-1].split("_")[1]) + 1
        else:
            init_idx = 1

        if bond_cluster_data:
            bond_cluster_info, bond_cluster_data, z_data = redefine_index(
                bond_cluster_info, bond_cluster_data, z_data, molecule=molecule, init_idx=init_idx
            )

        if molecule:
            atomic_info, atomic_data, z_data = redefine_index(atomic_info, atomic_data, z_data, molecule=molecule)
            uniform_info, uniform_data, z_data = redefine_index(uniform_info, uniform_data, z_data, molecule=molecule)
        else:
            if not self._model["info"]["generate"].get("fourier_transform", False):
                atomic_info, atomic_data, z_data = redefine_index(atomic_info, atomic_data, z_data)
            else:
                _, _, z_data = redefine_index(atomic_info, atomic_data, z_data)
                atomic_info, atomic_data, zk_data = redefine_index(atomic_info, atomic_data, zk_data)
                uniform_info, uniform_data, zk_data = redefine_index(uniform_info, uniform_data, zk_data)
                structure_info, structure_data, zk_data = redefine_index(structure_info, structure_data, zk_data)

        self["info"] = {
            "atomic": atomic_info,
            "site_cluster": site_cluster_info,
            "bond_cluster": bond_cluster_info,
        }
        self["data"] = {
            "atomic": atomic_data,
            "site_cluster": site_cluster_data,
            "bond_cluster": bond_cluster_data,
        }

        if self._model["info"]["molecule"]:
            self["info"] |= {"uniform": uniform_info}
            self["data"] |= {"uniform": uniform_data}

        self["info"] |= {"Z": z_info, "version": __version__}
        self["data"] |= {"Z": z_data}

        if not self._model["info"]["molecule"] and self._model["info"]["generate"]["fourier_transform"]:
            self["info"] |= {"uniform": uniform_info}
            self["data"] |= {"uniform": uniform_data}
            self["info"] |= {"structure": structure_info}
            self["data"] |= {"structure": structure_data}
            self["data"] |= {"Zk": zk_data}

        self["info"]["harmonics"] = self._harmonics_info()

    # ==================================================
    def fourier_transform(self):
        """
        convert combined multipoles from dictionary to numerical matrix form (real-space).

        Returns:
            dict: matrix form of combined multipoles, z_dict/z_full.

        Note:
            z_dict = { "z_#" : {(n1, n2, n3, a, b): matrix element}, R=(n1,n2,n3) is a lattice vector.
        """
        if self._model["info"]["molecule"]:
            return
        if "Zk" in self["data"]:
            return
        else:
            self._model["info"]["generate"]["fourier_transform"] = True

        bond_cluster_info = self["info"]["bond_cluster"]
        bond_cluster_data = self["data"]["bond_cluster"]
        uniform_data = self["data"]["uniform"]
        z_data = self["data"]["Z"]

        cluster_bond = self._model["data"]["cluster_bond"]
        bond = self._model["data"]["bond"]
        dim = len(self._model["data"]["site"])

        if len(bond) == 0:
            self._model["info"]["generate"]["fourier_transform"] = False
            return

        # structure
        bc_samb_set = {B_i: [bond_cluster_data[bmp_i] for bmp_i in lst] for B_i, lst in bond_cluster_info.items()}

        func = create_structure_samb_set
        args = [bc_samb_set, cluster_bond, bond, self._parallel]
        structure_info, structure_data = _run("structure", self._mpm, func, *args)

        # Zk
        bc_samb_set = {
            B_i: [(bmp_i, bond_cluster_data[bmp_i][1]) for bmp_i in lst] for B_i, lst in bond_cluster_info.items()
        }
        u_samb_set = [(ump_i, m) for ump_i, (_, m) in uniform_data.items()]
        k_samb_set = [(kmp_i, fk) for kmp_i, (_, fk) in structure_data.items()]

        func = create_zk_samb_set
        args = [z_data, bc_samb_set, u_samb_set, k_samb_set, cluster_bond, bond, dim, self._parallel]
        zk_data = _run("Zk", self._mpm, func, *args)

        self["info"]["structure"] = structure_info
        self["data"]["structure"] = structure_data
        self["data"]["Zk"] = zk_data

        if self._model["info"]["option"]["minimal_samb"]:
            self._minimal_samb()

        self["info"]["harmonics"] = self._harmonics_info()

    # ==================================================
    def create_dict(self):
        """
        create dictionary containing information and data of atomic, site-cluster, bond-cluster,
        uniform, (structure), Z(combined), Zk(fourier transform of Z).

        Returns:
            dict: information and data of multipole bases.
        """
        info = self["info"]

        atomic_data = {amp_i: (str(tag), *matrix_to_dict(m)) for amp_i, (tag, m) in self["data"]["atomic"].items()}
        site_cluster_data = {smp_i: (str(tag), str(v)) for smp_i, (tag, v) in self["data"]["site_cluster"].items()}
        bond_cluster_data = {bmp_i: (str(tag), str(v)) for bmp_i, (tag, v) in self["data"]["bond_cluster"].items()}
        z_data = {z_i: (str(tag), [tuple(map(str, v)) for v in lst]) for z_i, (tag, lst) in self["data"]["Z"].items()}

        data = {
            "atomic": atomic_data,
            "site_cluster": site_cluster_data,
            "bond_cluster": bond_cluster_data,
        }

        if self._model["info"]["molecule"]:
            uniform_data = {
                ump_i: (str(tag), *matrix_to_dict(m)) for ump_i, (tag, m) in self["data"]["uniform"].items()
            }
            data |= {"uniform": uniform_data}

        data |= {"Z": z_data}

        if not self._model["info"]["molecule"] and self._model["info"]["generate"]["fourier_transform"]:
            uniform_data = {
                ump_i: (str(tag), *matrix_to_dict(m)) for ump_i, (tag, m) in self["data"]["uniform"].items()
            }
            structure_data = {
                kmp_i: (str(tag), str(fk).replace("'", "")) for kmp_i, (tag, fk) in self["data"]["structure"].items()
            }
            zk_data = {
                z_i: (str(tag), [tuple(map(str, v)) for v in lst]) for z_i, (tag, lst) in self["data"]["Zk"].items()
            }
            data |= {"uniform": uniform_data}
            data |= {"structure": structure_data}
            data |= {"Zk": zk_data}

        return {"info": info, "data": data}

    # ==================================================
    def create_matrix(self, fmt="sympy"):
        """
        convert combined multipoles from dictionary to numerical matrix form (real-space).

        Args:
            fmt (str, optional): sympy/value.

        Returns:
            dict: matrix form of combined multipoles, z_dict/z_full.

        Note:
            z_dict = { "z_#" : {(n1, n2, n3, a, b): matrix element}, R=(n1,n2,n3) is a lattice vector.
        """
        if fmt not in ["sympy", "value"]:
            raise KeyError(f"unknown format = {fmt} is given.")

        site = self._model["data"]["site"]
        cluster_site = self._model["data"]["cluster_site"]

        bond = self._model["data"]["bond"]
        cluster_bond = self._model["data"]["cluster_bond"]

        n_sgn_dict = {}
        for B_i, bond_list in cluster_bond.items():
            n_sgn_list = []
            for bond_n in bond_list:
                vec = bond[bond_n][3]
                if "bond" in vec:
                    n = int(vec.split("_")[1])
                else:
                    n = int(bond_n.split("_")[1])

                if vec[0] == "-":
                    n_sgn_list.append((n, -1))
                else:
                    n_sgn_list.append((n, +1))

            n_sgn_dict[B_i] = n_sgn_list

        cluster_atomic = self._model["data"]["cluster_atomic"]
        atomic_braket = self._model["data"]["atomic_braket"]
        alias = self._model["name"]["alias"]
        dim_full = self._model["info"]["dimension"]

        atomic_data = self["data"]["atomic"]
        site_cluster_data = self["data"]["site_cluster"]
        bond_cluster_data = self["data"]["bond_cluster"]
        z_info, z_data = self["info"]["Z"], self["data"]["Z"]

        z_samb = {}
        z_samb["model"] = self._model["info"]["model"]
        z_samb["molecule"] = self._model["info"]["molecule"]
        z_samb["group"] = self._model["info"]["group"]
        z_samb["dimension"] = self._model["info"]["dimension"]
        z_samb["ket"] = self._model["info"]["ket"]
        z_samb["cell_site"] = self._model["info"]["cell_site"]
        z_samb["version"] = __version__

        if not z_samb["molecule"]:
            z_samb["k_point"] = self._model["info"]["k_point"]
            z_samb["k_path"] = self._model["info"]["k_path"]
            z_samb["A"] = self._model["detail"]["A"]

        z_full = {}
        for (_, M_i, SB_i), z_i_list in z_info.items():
            bra_orb_list, ket_orb_list = atomic_braket[M_i]
            dim_r, dim_c = len(bra_orb_list), len(ket_orb_list)
            for z_i in z_i_list:
                lst1 = z_data[z_i][1]

                # site × atomic
                if SB_i[0] == "S":
                    site_i_list = cluster_site[SB_i]
                    Z = sp.Matrix.zeros(dim_full, dim_full)
                    for c, amp_i, smp_i in lst1:
                        X = atomic_data[amp_i][1]
                        if self._model["info"]["molecule"]:
                            smp_i = smp_i.replace("u", "s")
                        v = list(site_cluster_data[smp_i][1])
                        for i in range(len(v)):
                            vi = v[i]
                            if vi == 0:
                                continue
                            bra_idx, ket_idx = site[site_i_list[i]][2]
                            for bra_no, ket_no, M_i_ in cluster_atomic[(bra_idx, ket_idx)]:
                                if M_i_ != M_i:
                                    continue
                                Z[bra_no : bra_no + dim_r, ket_no : ket_no + dim_c] += c * vi * sp.Matrix(X)

                    if bra_orb_list != ket_orb_list:
                        Z = (Z + Z.adjoint()) / sp.sqrt(2)

                    Z = Z / Z.norm()
                    z_full[z_i] = {
                        (0, 0, 0, a, b): Z[a, b] for a in range(dim_full) for b in range(dim_full) if Z[a, b] != 0
                    }

                # bond × atomic
                else:
                    S1, S2, _, _ = alias[SB_i].split(":")

                    z_full[z_i] = {}
                    n_sgn_list = n_sgn_dict[SB_i]
                    bond_i_list = cluster_bond[SB_i]
                    for v in lst1:
                        c, amp_i, bmp_i = v
                        X = atomic_data[amp_i][1]
                        if self._model["info"]["molecule"]:
                            bmp_i = bmp_i.replace("u", "b")
                        v = list(bond_cluster_data[bmp_i][1])
                        for i in range(len(v)):
                            vi = v[i]
                            if vi == 0:
                                continue

                            n, sgn = n_sgn_list[i]
                            bra_idx, ket_idx = bond[bond_i_list[i]][2]
                            bv = sgn * NSArray(bond[f"bond_{n:03d}"][3], style="vector")
                            if not self._model["info"]["molecule"]:
                                lattice = self._model["info"]["group"][1].split("/")[1].replace(" ", "")[0]
                                bv = to_primitive(lattice, bv)

                            s1 = NSArray(list(site.values())[bra_idx][0], style="vector")
                            s2 = NSArray(list(site.values())[ket_idx][0], style="vector")
                            n1, n2, n3 = -(bv - (s1 - s2))
                            n1, n2, n3 = round(n1), round(n2), round(n3)

                            for bra_no, ket_no, M_i_ in cluster_atomic[(bra_idx, ket_idx)]:
                                if M_i_ != M_i:
                                    continue
                                Z = sp.Matrix.zeros(dim_full, dim_full)
                                Z[bra_no : bra_no + dim_r, ket_no : ket_no + dim_c] += c * vi * sp.Matrix(X)

                                if S1 == S2 and bra_orb_list != ket_orb_list:
                                    if bra_idx == ket_idx:
                                        Z[ket_no : ket_no + dim_c, bra_no : bra_no + dim_r] += (
                                            c * vi * sp.Matrix(X).adjoint()
                                        )
                                    else:
                                        bra_no += dim_r
                                        ket_no -= dim_r
                                        Z[bra_no : bra_no + dim_c, ket_no : ket_no + dim_r] += (
                                            c * vi * sp.Matrix(X).adjoint()
                                        )
                                    Z /= sp.sqrt(2)

                                for a in range(dim_full):
                                    for b in range(dim_full):
                                        z = sp.expand(Z[a, b])
                                        if z == sp.S(0) or z == 0:
                                            continue

                                        if (n1, n2, n3, a, b) not in z_full[z_i]:
                                            z_full[z_i][(n1, n2, n3, a, b)] = z / sp.sqrt(2)
                                            z_full[z_i][(-n1, -n2, -n3, b, a)] = sp.conjugate(z) / sp.sqrt(2)
                                        else:
                                            z_full[z_i][(n1, n2, n3, a, b)] += z / sp.sqrt(2)
                                            z_full[z_i][(-n1, -n2, -n3, b, a)] += sp.conjugate(z) / sp.sqrt(2)

        if fmt == "value":
            DIGIT = 14
            z_full = {
                z_i: {
                    k: str(round(complex(v).real, DIGIT) + round(complex(v).imag, DIGIT) * 1j)
                    for k, v in d.items()
                    if sp.expand(v) != 0
                }
                for z_i, d in z_full.items()
            }
        else:
            z_full = {z_i: {k: str(v) for k, v in d.items() if sp.expand(v) != 0} for z_i, d in z_full.items()}

        z_samb["matrix"] = z_full

        return z_samb

    # ==================================================
    def create_matrix_k(self, full=False, fmt="sympy"):
        """
        convert combined multipoles from dictionary to numerical matrix form (momentum-space).

        Args:
            full (bool, optional): full matrix format ?
            fmt (str, optional): sympy/value.

        Returns:
            dict: matrix form of combined multipoles, z_dict/z_full.

        Note:
            z_dict = { "z_#" : [(i, j, ""/c###/s###, v)] }
            z_full = { "z_#" : NSArray(matrix) }
        """
        if fmt not in ["sympy", "value"]:
            raise KeyError(f"unknown format = {fmt} is given.")

        molecule = self._model["info"]["molecule"]
        cluster_atomic = self._model["data"]["cluster_atomic"]
        atomic_braket = self._model["data"]["atomic_braket"]
        alias = self._model["name"]["alias"]
        dim_full = self._model["info"]["dimension"]

        atomic_data = self["data"]["atomic"]
        uniform_data = self["data"]["uniform"]
        if molecule:
            z_info, z_data = self["info"]["Z"], self["data"]["Z"]
        else:
            structure_data = self["data"]["structure"]
            z_info, z_data = self["info"]["Z"], self["data"]["Zk"]

        z_samb = {}
        z_samb["model"] = self._model["info"]["model"]
        z_samb["molecule"] = self._model["info"]["molecule"]
        z_samb["group"] = self._model["info"]["group"]
        z_samb["dimension"] = self._model["info"]["dimension"]
        z_samb["ket"] = self._model["info"]["ket"]
        z_samb["version"] = __version__

        if not molecule:
            z_samb["k_point"] = self._model["info"]["k_point"]
            z_samb["k_path"] = self._model["info"]["k_path"]
            z_samb["A"] = self._model["detail"]["A"]
            bond = {}
            for name, ex in self._model["data"]["bond"].items():
                ex = ex[0].split("@")[0]
                bond[name] = ex
            z_samb["bond"] = bond

        z_full = {}
        for (_, M_i, SB_i), z_i_list in z_info.items():
            if SB_i[0] == "S":
                S1 = S2 = alias[SB_i]
            else:
                S1, S2, _, _ = alias[SB_i].split(":")

            bra_orb_list, ket_orb_list = atomic_braket[M_i]
            dim_r, dim_c = len(bra_orb_list), len(ket_orb_list)
            for z_i in z_i_list:
                lst1 = z_data[z_i][1]
                Z = sp.Matrix.zeros(dim_full, dim_full)
                for v in lst1:
                    if len(v) == 3:
                        c, amp_i, ump_i = v
                        F = 1
                    else:
                        c, amp_i, ump_i, kmp_i = v
                        F = structure_data[kmp_i][1]
                    X = atomic_data[amp_i][1]
                    U = uniform_data[ump_i][1]
                    for (bra_idx, ket_idx), lst2 in cluster_atomic.items():
                        u = U[bra_idx, ket_idx]
                        if u != 0:
                            for bra_no, ket_no, M_i_ in lst2:
                                if M_i_ != M_i:
                                    continue

                                if S1 == S2 and bra_orb_list != ket_orb_list:
                                    Z[bra_no : bra_no + dim_r, ket_no : ket_no + dim_c] += (
                                        c * u * F * sp.Matrix(X) / sp.sqrt(2)
                                    )
                                    if bra_idx == ket_idx:
                                        Z[ket_no : ket_no + dim_c, bra_no : bra_no + dim_r] += (
                                            c * u * F * sp.Matrix(X).adjoint()
                                        ) / sp.sqrt(2)
                                    else:
                                        bra_no += dim_r
                                        ket_no -= dim_r
                                        Z[bra_no : bra_no + dim_c, ket_no : ket_no + dim_r] += (
                                            c * u * F * sp.Matrix(X).adjoint()
                                        ) / sp.sqrt(2)
                                else:
                                    Z[bra_no : bra_no + dim_r, ket_no : ket_no + dim_c] += c * u * F * sp.Matrix(X)

                if not Z.is_hermitian:
                    Z = Z + Z.adjoint()
                z_full[z_i] = sp.expand(Z)

        DIGIT = 14
        if full:
            if fmt == "sympy":
                z_full = {z_i: str(mat.tolist()) for z_i, mat in z_full.items()}
            else:
                z_full_val = {}
                for z_i, Z in z_full.items():
                    Z_val = sp.Matrix.zeros(dim_full, dim_full)
                    for i, j in product(range(dim_full), range(dim_full)):
                        vij = Z[i, j]
                        if vij != 0:
                            if len(z_data[z_i][1][0]) == 3:
                                vij = complex(vij)
                                vij = round(vij.real, DIGIT) + round(vij.imag, DIGIT) * 1j
                                Z_val[i, j] += vij
                            else:
                                d = {csn: vij.coeff(csn) for csn in cs_list(vij, real=True)}
                                for csn, vij in d.items():
                                    if vij != sp.S(0):
                                        vij = complex(vij)
                                        vij = round(vij.real, DIGIT) + round(vij.imag, DIGIT) * 1j
                                        Z_val[i, j] += vij * csn
                    z_full_val[z_i] = Z_val
                z_full = {z_i: str(NSArray(mat.tolist(), "matrix")) for z_i, mat in z_full_val.items()}
            z_samb["matrix"] = z_full
            return z_samb
        else:
            z_dict = {}
            for z_i, Z in z_full.items():
                mat = []
                for i, j in product(range(dim_full), range(dim_full)):
                    vij = Z[i, j]
                    if vij != 0:
                        if len(z_data[z_i][1][0]) == 3:
                            mat.append((i, j, "", str(vij)))
                        else:
                            d = {str(csn): vij.coeff(csn) for csn in cs_list(vij, real=True)}
                            for csn_str, vij in d.items():
                                if vij != sp.S(0):
                                    if fmt == "value":
                                        vij = complex(vij)
                                        vij = round(vij.real, DIGIT) + round(vij.imag, DIGIT) * 1j
                                    mat.append((i, j, csn_str, str(vij)))
                z_dict[z_i] = mat

            z_samb["matrix"] = z_dict

            return z_samb

    # ==================================================
    @classmethod
    def _header(cls):
        return header_str

    # ==================================================
    @classmethod
    def _matrix_header_k(cls):
        return matrix_header_k_str

    # ==================================================
    @classmethod
    def _matrix_header(cls):
        return matrix_header_str
