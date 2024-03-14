"""
MaterialModel manages model information of cluster or crystal system.
"""

import numpy as np
from gcoreutils.nsarray import NSArray
from gcoreutils.crystal_util import cell_transform_matrix
from multipie.multipole.util.atomic_orbital_util import (
    parse_orb_list,
    sort_orb_list,
    split_orb_list_rank_block,
    rank,
)
from multipie import __version__

try:
    from qtdraw.qt_draw import QtDraw
    from qtdraw.multipie.setting import default_style

    _qtdraw_available = True
except ImportError:
    _qtdraw_available = False


_default_max_neighbor = 10
_default_search_cell_range = (-2, 3, -2, 3, -2, 3)


# ==================================================
input_str = """
=== Input for SAMB construction (* optional [default]) ===
- model : model name (str).
- group : group name (Schoenflies notation) or group number (space group) (str/int).
- site : site position, orbital info. (dict) { name: ("position", orbital or orbital list) }.
- bond* : bond (list) [ ("tail", "head", (list of) neighbors) ], [[]].
- spinful* : spinful basis ? (bool), [False].
- cell* : unit-cell constants (dict) { "a", "b", "c", "alpha", "beta", "gamma" }, [a=b=c=1,alpha=beta=gamma=90].
- option*
  - view* : view point (int list), [None].
  - view_mode* : mode for QtDraw file (str) ("standard"/"arrow"/"debug"), ["standard"].
  - output* : output folder (str), [model name].
  - minimal_samb* (bool) : minimal output in _samb.pdf ? [True].
  - binary_output* (bool) : output matrix data in binary format ? [False].
- generate*
  - model_type* : model type (str), ("tight_binding"/"phonon"), ["tight_binding"].
  - time_reversal_type* : time-reversal type (str), ("electric"/"magnetic"/"both"), ["electric"].
  - irrep* : irrep. (str list), [identity irrep.] (empty list is for all irreps.), [None].
  - fourier_transform* : create fourier transformed SAMB ? [False].
  - toroidal_priority* : create toroidal multipoles (G,T) in high priority ? [False].
- k_point* : k-point (dict) {name: "position"}, [{ "Γ": "[0,0,0]", "X": "[1/2,0,0]" }].
- k_path* : k-path (str) (concatenate by "-" or "\|"), ["Γ-X"].
- detail*
  - cell_range* : search range for bonds, [(-2, 3, -2, 3, -2, 3)].
  - max_neighbor* : max. of neighbors to search, [10].
"""

# ==================================================
header_str = """
=== Molecule or Crystal Model (* only for crystal) ===
- info
    - model : model name.
    - molecule : molecule or crystal ?
    - group : (tag, detailed str)
    - crystal : crystal class
    - cell* : {a, b, c, alpha, beta, gamma}
    - volume* : unit cell volume
    - a1* : unit cell a1 vector
    - a2* : unit cell a2 vector
    - a3* : unit cell a3 vector
    - option :
        - view : view index
        - view_mode : QtDraw mode, standard/arrow/debug
        - output : output directory.
        - minimal_samb : output minimal SAMB ?
        - binary_output : output matrix data in binary format ?
    - generate :
        - model_type : tight_binding/phonon
        - time_reversal_type : electric/magnetic/both
        - irrep : irrep list
        - fourier_transform* : create fourier transformed SAMB ?
        - toroidal_priority : create toroidal multipoles (G,T) in high priority ?
    - k_point* : representative k points
    - k_path* : high-symmetry line in k space
    - dimension : dimension of full matrix
    - spinful : spinful or not
    - orbital : list of orbitals (U,D: up/down spin)
    - ket : ket basis list, orbital@site
    - ket_site : list of sites
    - site : input for "site" { name: (position, orbitals) }
    - rep_site : representative site { name: (position, wp, orbitals, site-symmetry) }
    - cell_site : { name_idx(pset): (position, SOs) }
    - bond : input for "bond" [ (tail, head, neighbors) ]
    - rep_bond : representative bond { name: (vector@center, wp, directional, neighbor, site-symmetry) }
    - cell_bond : { name_idx(pset): (vector@center, SOs) }

- name
    - alias : { cluster_tag: name or name: cluster_tag }
    - site : { site_tag: (name, no) }
    - site_name : { position : (site_tag, pset) }
    - bond : { bond_tag: (tail:head:n:multiplicity, no) }
    - bond_name : { vector@center : (bond_tag, pset) }

- data
    - plus_set* : [ plus_set list ]
    - cluster_site : { cluster_tag: site_list }
    - cluster_bond : { cluster_tag: bond_list }
    - site : { site_tag: (position, SO, (bra_site_no, ket_site_no)) }
    - bond : { bond_tag: (vector@center, SO, (bra_site_no, ket_site_no), vector, tail;head) }
    - cluster_atomic : { (bra_site_no, ket_site_no): [(bra_no, ket_no, matrix_tag)] }
    - atomic_braket : { matrix_tag : (bra_list, ket_list) }

- detail
    - rep_bond_all : { tail_head: [rep_bond in 0-9th neighbors] }
    - cell_range* : search range for bonds
    - max_neighbor : max. of neighbors to search
    - A* : transform matrix, [a1,a2,a3]
    - version : MultiPie version
"""


# ==================================================
class MaterialModel(dict):
    """
    molecule or crystal model information.
    """

    # Attributes:
    #     _mpm (MultiPieManager): MultiPieManager.
    #     _cell (dict): cell info.
    #     _volume (float): cell volume.
    #     _A (NSArray): unit vectors in each column (3x3).
    #     _G (NSArray): metric tensor (3x3).
    #     _A_norm (NSArray): normalized unit vectors.
    # ==================================================
    def __init__(self, dic, mpm):
        """
        initialize the class.

        Args:
            dic (dict): model info. dict.
            mpm (MultiPieManager): MultiPie manager.
        """
        info = {
            "model": "default",
            "molecule": False,
            "group": (1, ""),
            "crystal": "monoclinic",
            "cell": {
                "a": 1.0,
                "b": 1.0,
                "c": 1.0,
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 90.0,
            },
            "volume": 1.0,
            "a1": "[1.0, 0.0, 0.0]",
            "a2": "[0.0, 1.0, 0.0]",
            "a3": "[0.0, 0.0, 1.0]",
            "option": {
                "view": [0, 0, 1],
                "view_mode": "standard",
                "output": "material_model",
                "minimal_samb": True,
                "binary_output": False,
            },
            "generate": {
                "model_type": "tight_binding",
                "time_reversal_type": "electric",
                "irrep": ["A"],
                "fourier_transform": False,
                "toroidal_priority": False,
            },
            "k_point": {},
            "k_path": "",
            "dimension": 1,
            "spinful": True,
            "orbital": [],
            "ket": [],
            "ket_site": [],
            "site": {},
            "rep_site": {},
            "cell_site": {},
            "bond": [],
            "rep_bond": {},
            "cell_bond": {},
        }
        name = {
            "alias": {},
            "site": {},
            "site_name": {},
            "bond": {},
            "bond_name": {},
        }
        data = {
            "plus_set": [],
            "cluster_site": {},
            "cluster_bond": {},
            "site": {},
            "bond": {},
            "cluster_atomic": {},
            "atomic_braket": {},
        }
        detail = {
            "rep_bond_all": {},
            "cell_range": _default_search_cell_range,
            "max_neighbor": _default_max_neighbor,
            "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]",
            "version": __version__,
        }

        self._mpm = mpm

        info["molecule"] = self._mpm.molecule

        info["spinful"] = dic["spinful"]
        info["option"] = dic["option"]
        info["model"] = dic["model"]
        info["group"] = (str(self._mpm.group), self._group_str(info["molecule"]))
        info["crystal"] = self._mpm.group.tag.crystal

        cell_info = cell_transform_matrix(cell=dic["cell"], crystal=info["crystal"], translation=False)
        self._cell = cell_info[0]
        self._volume = cell_info[1]
        self._A = cell_info[2]
        self._G = cell_info[3]
        self._A_norm = cell_info[4]
        if not info["molecule"]:
            info["cell"] = self._cell
            info["volume"] = self._volume
            info["a1"] = str(self._A[:, 0])
            info["a2"] = str(self._A[:, 1])
            info["a3"] = str(self._A[:, 2])

        info["generate"] = dic["generate"]
        irrep = info["generate"]["irrep"]
        if irrep is None:
            irrep = self._mpm.point_group.identity_irrep
        elif len(irrep) == 0:
            irrep = [str(i) for i in self._mpm.point_group.character.irrep_list]
        if type(irrep) == str:
            irrep = [irrep]
        irrep = sorted(list(set([i.replace("a", "").replace("b", "") for i in irrep])))
        info["generate"]["irrep"] = irrep

        if info["molecule"]:
            del (
                info["cell"],
                info["volume"],
                info["a1"],
                info["a2"],
                info["a3"],
                info["k_point"],
                info["k_path"],
                info["generate"]["fourier_transform"],
            )
            del data["plus_set"]
            del detail["cell_range"], detail["A"]
        else:
            info["k_point"] = dic["k_point"]
            info["k_path"] = dic["k_path"]
            detail["A"] = str(self.A)
            data["plus_set"] = self._mpm.group.symmetry_operation.plus_set.str()

        if "detail" in dic.keys():
            if "cell_range" in dic["detail"]:
                detail["cell_range"] = dic["detail"]["cell_range"]
            if "max_neighbor" in dic["detail"]:
                detail["max_neighbor"] = dic["detail"]["max_neighbor"]

        if "site" in dic.keys():
            info["site"] = dic["site"]
            self._set_rep_site(info)
            orbital = self._set_site(info, name, data)
            ket_dict = self._set_orbital(info, name, orbital)

        if "bond" in dic.keys():
            info["bond"] = dic["bond"]
            self._set_rep_bond(info, name, detail)
            self._set_bond(info, name, data)

        self._set_matrix(name, data, orbital, ket_dict)

        self["info"] = info
        self["name"] = name
        self["data"] = data
        self["detail"] = detail

    # ==================================================
    @property
    def model(self):
        return self["info"]["model"]

    # ==================================================
    def create_view(self, scale=True, mode=None):
        """
        create QtDraw file.

        Args:
            scale (bool, optional): objects are scaled by means of volume.
            mode (str, optional): view mode, standard/arrow/debug. if None, option in the model is used.

        Returns:
            dict: QtDraw dict.
        """
        if not _qtdraw_available:
            raise Exception("QtDraw is not installed.")

        info = self["info"]
        name = self["name"]
        data = self["data"]
        molecule = info["molecule"]

        if mode is None:
            mode = info["option"]["view_mode"]

        if molecule:
            n_pset = 1
        else:
            n_pset = len(data["plus_set"])

        qtdraw = QtDraw(
            model=info["model"],
            cell=self._cell,
            view=info["option"]["view"],
            axis_type="abc",
            cluster=molecule,
            clip=False,
            background=True,
        )
        qtdraw.set_crystal(info["crystal"])
        if scale:
            scale = self._volume ** (1 / 3)
        else:
            scale = 1.0

        cluster_site_n = len(default_style["site"])
        cluster_bond_n = len(default_style["bond"])

        # plot sites.
        for pos, (site_name, pset) in name["site_name"].items():
            head, idx = name["site"][site_name]
            cluster_no = int(name["alias"][head].split("_")[1])
            prop_no = cluster_no - 1 if cluster_no < cluster_site_n else cluster_site_n - 1
            color, size, opacity = default_style["site"][prop_no]
            if n_pset > 1:
                cluster_name = f"S{cluster_no}({pset})"
            else:
                cluster_name = f"S{cluster_no}"
            mp = data["site"][site_name][1]
            if mode == "debug":
                site_no = int(site_name.split("_")[1])
                label = f"s{site_no}"
            else:
                label = f"{head}_{idx}: " + self._mapping_str(mp)
            qtdraw.plot_site(
                pos,
                name=cluster_name,
                label=label,
                color=color,
                size=size * scale,
                opacity=opacity,
            )
            if mode != "standard" and idx == 1 and pset == 1:
                qtdraw.plot_site(
                    pos,
                    name="Z",
                    label=label,
                    color="yellow",
                    size=size * scale * 1.2,
                    opacity=0.3,
                )

        # plot bonds.
        for bond, (bond_name, pset) in name["bond_name"].items():
            head, idx = name["bond"][bond_name]
            cluster_no = int(name["alias"][head].split("_")[1])
            prop_no = cluster_no - 1 if cluster_no < cluster_bond_n else cluster_bond_n - 1
            (c1, c2), width, opacity = default_style["bond"][prop_no]
            if n_pset > 1:
                cluster_name = f"B{cluster_no}({pset})"
            else:
                cluster_name = f"B{cluster_no}"
            _, mp, ij, _, _ = data["bond"][bond_name]

            bond_alias_name = name["bond"][bond_name][0]
            bond_no = int(bond_name.split("_")[1])
            _, _, nd, n, _ = info["rep_bond"][bond_alias_name]
            color1, color2 = (c1, c1) if nd == "ND" else (c1, c2)
            if mode == "debug":
                label = f"b{bond_no} ({ij[0]+1},{ij[1]+1})"
            else:
                label = f"b{bond_no}/{n}th: " + self._mapping_str(mp)
            v, c = bond.split("@")
            opa = opacity * 0.5 if mode != "standard" else opacity
            qtdraw.plot_bond(
                c,
                v,
                color=color1,
                color2=color2,
                width=width * scale,
                opacity=opa,
                name=cluster_name,
                label=label,
            )
            if mode != "standard":
                v = NSArray(v).transform(self.A)
                v_len = v.norm() * 0.25
                acolor = "red" if idx == 1 and pset == 1 else "black"
                qtdraw.plot_vector(
                    c,
                    v,
                    length=v_len,
                    width=width * scale * 0.3,
                    color=acolor,
                    opacity=opacity,
                    name=cluster_name,
                    label=label,
                )

        dic = qtdraw.create_dict(self.model)
        qtdraw.close()
        del qtdraw

        return dic

    # ==================================================
    def _set_rep_site(self, info):
        """
        set representative site.
        """
        spinful = info["spinful"]
        crystal = info["crystal"]

        rep_site = {}
        for name, (pos, orb_list) in info["site"].items():
            wp = str(self._mpm.group.find_wyckoff_position(str(pos)))
            pos = str(NSArray(pos))
            orb_list = split_orb_list_rank_block(
                sort_orb_list(parse_orb_list(orb_list, spinful, crystal), spinful, crystal),
                spinful,
                crystal,
            )
            sym = self._mpm.group.wyckoff.site_symmetry(wp)
            if not info["molecule"]:
                pos = str(self._mpm.group.shift_home_unit_cell(pos))
            rep_site[name] = (pos, wp, orb_list, sym)

        info["rep_site"].update(rep_site)

    # ==================================================
    def _set_site(self, info, name, data):
        """
        set sites in cluster or cell.
        """
        alias = {}
        cluster_site = {}
        cell_site = {}
        data_site = {}
        name_site = {}
        site_to_name = {}
        orbital = {}

        if not info["molecule"]:
            plus_set = NSArray.from_str(data["plus_set"])
            n_pset = len(plus_set)
        else:
            n_pset = 1

        site_no = 0
        for no, (site, (pos, _, orb_list, _)) in enumerate(info["rep_site"].items()):
            if info["molecule"]:
                s_map = self._mpm.group.site_mapping(pos)
                basic_num = len(s_map)
            else:
                s_map = self._mpm.group.site_mapping(pos, plus_set=True)
                basic_num = len(s_map) // n_pset
            position = list(s_map.keys())
            prop = list(s_map.values())
            cluster = f"S_{no+1:03d}"
            alias[cluster] = site
            alias[site] = cluster
            cluster_site[cluster] = []
            for rsite_no in range(basic_num):
                s, mp = (
                    position[rsite_no],
                    prop[rsite_no],
                )  # as sorted by pset group in s_map.
                s_no = site_no + rsite_no
                sname = f"site_{s_no+1:03d}"
                cluster_site[cluster].append(sname)
                ij = (s_no, s_no)
                data_site[sname] = (s, mp, ij)
                if n_pset > 1:
                    for pset in range(n_pset):
                        ss = position[pset * basic_num + rsite_no]
                        cell_site[f"{site}_{rsite_no+1}({pset+1})"] = (
                            ss,
                            self._mapping_str(mp),
                        )
                else:
                    cell_site[f"{site}_{rsite_no+1}"] = (s, self._mapping_str(mp))
                name_site[sname] = (site, rsite_no + 1)
                orbital[sname] = orb_list
            assert data_site[f"site_{site_no+1:03d}"][1][0] == 0, f"first SO is not identity operation, {mp}"
            for rsite_no in range(basic_num):
                s_no = site_no + rsite_no
                sname = f"site_{s_no+1:03d}"
                for pset in range(n_pset):
                    s = position[pset * basic_num + rsite_no]
                    site_to_name[s] = (sname, pset + 1)
            site_no = site_no + basic_num

        # check orbitals
        d = {}
        spinful = info["spinful"]
        crystal = info["crystal"]
        for sname, orb_list in orbital.items():
            for orbs in orb_list:
                o = orbs[0]
                l = rank(o, spinful, crystal)
                if l in d:
                    d[l] += [(sname, orbs)]
                else:
                    d[l] = [(sname, orbs)]
        for l, lst in d.items():
            sname, orbs = lst[0]
            if len(lst) > 1:
                for sname_, orbs_ in lst[1:]:
                    if orbs != orbs_:
                        s1 = name_site[sname][0]
                        s2 = name_site[sname_][0]
                        raise Exception(f"invalid orbitals are given, {s1} {orbs} != {s2} {orbs_}.")

        info["cell_site"].update(cell_site)

        name["alias"].update(alias),
        name["site"].update(name_site)
        name["site_name"].update(site_to_name)

        data["cluster_site"].update(cluster_site)
        data["site"].update(data_site)

        return orbital

    # ==================================================
    def _set_orbital(self, info, name, orbital):
        ket_site = [f"{head}_{idx}" for head, idx in name["site"].values()]

        unique_orb = set()
        ket = []
        ket_dict = {}
        ket_no = 0
        for no, site in enumerate(name["site"].keys()):
            s = ket_site[no]
            for orb in sum(orbital[site], []):
                ket.append(f"{orb}@{s}")
                ket_dict[(orb, no)] = ket_no
                ket_no = ket_no + 1
                unique_orb.add(orb)

        dimension = len(ket)
        orbital = sort_orb_list(list(unique_orb), info["spinful"], info["crystal"])

        info["dimension"] = dimension
        info["ket"] = ket
        info["ket_site"] = ket_site
        info["orbital"] = orbital

        return ket_dict

    # ==================================================
    def _set_rep_bond(self, info, name, detail):
        """
        set representative bond.
        """
        rep_bond_all = {}
        for tail, head, neighbor in info["bond"]:
            if type(neighbor) is not list:
                neighbor = [neighbor]
            if tail not in info["rep_site"].keys():
                raise Exception(f"{tail} is not found in sites.")
            if head not in info["rep_site"].keys():
                raise Exception(f"{head} is not found in sites.")

            # create rep. bond.
            tail_pos = info["rep_site"][tail][0]
            head_pos = info["rep_site"][head][0]
            rep_bond = self._representative_bond(tail_pos, head_pos, info, detail)
            rep_bond_list = [{}]
            for n, lst in enumerate(rep_bond):
                n_rep_bond = {}
                for i, b in enumerate(lst):
                    bb = NSArray(b)
                    v, c = bb.convert_bond("bond")
                    wp = str(self._mpm.group.find_wyckoff_position(str(c)))
                    _, nd = self._mpm.group.bond_mapping(b)
                    nd = "ND" if nd else "D"
                    b = f"{v}@{c}"
                    b = self._rep_bond_site_order(b, info, name)
                    bond_name = f"{tail}:{head}:{n+1}:{i+1}"
                    sym = self._mpm.group.wyckoff.site_symmetry(wp)
                    n_rep_bond[bond_name] = (b, wp, nd, n + 1, sym)
                rep_bond_list.append(n_rep_bond)

            bond_th = f"{tail}_{head}"
            rep_bond_all[bond_th] = rep_bond_list

            # set rep. bond for required neighbors.
            rep_bond_selected = {}
            for n in neighbor:
                b = rep_bond_all[bond_th][n]
                for k, v in b.items():
                    rep_bond_selected[k] = v

            info["rep_bond"].update(rep_bond_selected)

        detail["rep_bond_all"].update(rep_bond_all)

    # ==================================================
    def _rep_bond_site_order(self, bond, info, name):
        """
        get rep_bond in order of site numbers.
        """
        site_to_name = name["site_name"]

        if info["molecule"]:
            b_map, _ = self._mpm.group.bond_mapping(bond)
        else:
            b_map, _ = self._mpm.group.bond_mapping(bond, plus_set=True)

        lst = []
        for b in b_map.keys():
            bb = NSArray(b)
            t, h = bb.convert_bond("bond_th")
            v, c = bb.convert_bond("bond")
            if info["molecule"]:
                tail_name = site_to_name[str(t)][0]
                head_name = site_to_name[str(h)][0]
            else:
                tail_name = site_to_name[str(t.shift())][0]
                head_name = site_to_name[str(h.shift())][0]
            t_no, h_no = (
                int(tail_name.split("_")[1]) - 1,
                int(head_name.split("_")[1]) - 1,
            )
            lst.append((t_no, h_no, b))
        lst = sorted(lst)

        return lst[0][2]

    # ==================================================
    def _set_bond(self, info, name, data):
        """
        set bonds in cell.
        """
        # set bond in cluster or cell for required neighbors.
        alias = {}
        cluster_bond = {}
        cell_bond = {}
        data_bond = {}
        name_bond = {}
        bond_to_name = {}

        site_to_name = name["site_name"]

        if not info["molecule"]:
            plus_set = NSArray.from_str(data["plus_set"])
            n_pset = len(plus_set)
        else:
            n_pset = 1

        bond_no = 1  # no=0 is for null vector.
        for no, (rbname, (bond, _, _, _, _)) in enumerate(info["rep_bond"].items()):
            if info["molecule"]:
                b_map, _ = self._mpm.group.bond_mapping(bond)
                basic_num = len(b_map)
            else:
                b_map, _ = self._mpm.group.bond_mapping(bond, plus_set=True)
                basic_num = len(b_map) // n_pset
            bonds = list(b_map.keys())
            prop = list(b_map.values())
            cluster = f"B_{no+1:03d}"
            alias[cluster] = rbname
            alias[rbname] = cluster
            cluster_bond[cluster] = []
            b_vector = {}
            for rbond_no in range(basic_num):
                b, mp = (
                    bonds[rbond_no],
                    prop[rbond_no],
                )  # as sorted by pset group in s_map.
                b_no = bond_no + rbond_no
                bname = f"bond_{b_no:03d}"
                cluster_bond[cluster].append(bname)
                bb = NSArray(b)
                t, h = bb.convert_bond("bond_th")
                v, c = bb.convert_bond("bond")
                if info["molecule"]:
                    tail_name = site_to_name[str(t)][0]
                    head_name = site_to_name[str(h)][0]
                else:
                    tail_name = site_to_name[str(t.shift())][0]
                    head_name = site_to_name[str(h.shift())][0]
                ht_no = (
                    int(head_name.split("_")[1]) - 1,
                    int(tail_name.split("_")[1]) - 1,
                )
                assert ht_no[0] >= ht_no[1], f"site indices are wrong, {ht_no}"
                vc = f"{v}@{c}"
                data_bond[bname] = (vc, mp, ht_no, str(v), b)
                br = str(NSArray(vc).regular_direction()).split("@")[0]
                b_vector[br] = b_vector.get(br, []) + [bname]
                if n_pset > 1:
                    for pset in range(n_pset):
                        v, c = NSArray(bonds[pset * basic_num + rbond_no]).convert_bond("bond")
                        vc1 = f"{v}@{c}"
                        cell_bond[f"{rbname}_{rbond_no+1}({pset+1})"] = (
                            vc1,
                            self._mapping_str(mp),
                        )
                else:
                    cell_bond[f"{rbname}_{rbond_no+1}"] = (vc, self._mapping_str(mp))
                name_bond[bname] = (rbname, rbond_no + 1)
            assert data_bond[f"bond_{bond_no:03d}"][1][0] == 0, f"first SO is not identity operation, {mp}"
            for rbond_no in range(basic_num):
                b_no = bond_no + rbond_no
                bname = f"bond_{b_no:03d}"
                for pset in range(n_pset):
                    v, c = NSArray(bonds[pset * basic_num + rbond_no]).convert_bond("bond")
                    b = f"{v}@{c}"
                    bond_to_name[b] = (bname, pset + 1)
            bond_no = bond_no + basic_num

            # replace same bond vector.
            for lst in b_vector.values():
                top = lst[0]
                v0 = data_bond[top][3]
                for i in lst[1:]:
                    b, mp, ht_no, v, th = data_bond[i]
                    nv = top if v0 == v else "-" + top
                    data_bond[i] = (b, mp, ht_no, nv, th)
            data["bond"].update(data_bond)

        info["cell_bond"].update(cell_bond)

        name["alias"].update(alias)
        name["bond"].update(name_bond)
        name["bond_name"].update(bond_to_name)

        data["cluster_bond"].update(cluster_bond)
        data["bond"].update(data_bond)

    # ==================================================
    def _set_matrix(self, name, data, orbital, ket_dict):
        cluster_atomic = {}
        atomic_braket = {}

        ij_list = [v[2] for v in data["site"].values()] + [v[2] for v in data["bond"].values()]

        for ij in ij_list:
            i, j = ij
            site_i = f"site_{i+1:03d}"
            site_j = f"site_{j+1:03d}"
            Si = name["site"][list(data["site"].keys())[i]][0]
            Sj = name["site"][list(data["site"].keys())[j]][0]
            if Si == Sj:
                orbs_i = orbital[site_i]
                orbs_j = orbs_i
                braket = sum(
                    [[(orb, orbs_j[j]) for j in range(i, len(orbs_i))] for i, orb in enumerate(orbs_i)],
                    [],
                )
            else:
                orbs_i = orbital[site_i]
                orbs_j = orbital[site_j]
                braket = sum([[(orbi, orbj) for orbj in orbs_j] for orbi in orbs_i], [])
            cluster_atomic[ij] = [(f"{bra}:{ket}", (bra[0], i), (ket[0], j)) for bra, ket in braket]
            for bra, ket in braket:
                atomic_braket[f"{bra}:{ket}"] = (bra, ket)

        rs = {k: f"M_{i+1:03d}" for i, k in enumerate(atomic_braket.keys())}
        d = {rs[k]: v for k, v in atomic_braket.items()}
        atomic_braket = d

        dca = {}
        for k, lst in cluster_atomic.items():
            nlst = [(ket_dict[k1], ket_dict[k2], rs[i]) for i, k1, k2 in lst]
            dca[k] = nlst
        cluster_atomic = dca

        data["cluster_atomic"].update(cluster_atomic)
        data["atomic_braket"].update(atomic_braket)

    # ==================================================
    @classmethod
    def _mapping_str(cls, mp):
        """
        modify SO mapping (+1).

        Args:
            mp (list): mapping list.

        Returns:
            str: modified mapping list.
        """
        m = [i - 1 if i < 0 else i + 1 for i in mp]
        return str(m).replace(" ", "")

    # ==================================================
    def _group_str(self, molecule):
        """
        create group string.

        Args:
            molecule (bool): molecule or crystal ?

        Returns:
            str: info. for group.
        """
        no, _, IS, setting = self._mpm.group.tag.info()
        if setting:
            setting = "(" + setting + ")"
        SS = str(self._mpm.group.tag)
        group_str = "point" if molecule else "space"
        ret = f"{group_str} group No. {no} : {SS} / {IS}"
        if setting:
            ret += f" {setting}"
        if not molecule:
            ret += f" : PG {self._mpm.point_group}"

        return ret

    # ==================================================
    def _site_grid(self, site, info, detail, repeat=False):
        """
        site grid.

        Args:
            site (str): representative site.
            info (dict): info. dict.
            detail (dict): detail dict.

        Returns:
            NSArray: all sites for given grid.
        """
        if info["molecule"]:
            site = self._mpm.group.transform_site(site, remove_duplicate=True)
            return site

        site = self._mpm.group.transform_site(site, shift=True, remove_duplicate=True, plus_set=True)

        if repeat:
            offset = detail["cell_range"][::2]
            N = [i - j for i, j in zip(detail["cell_range"][1::2], offset)]
            ns = N[0] * N[1] * N[2]

            igrid = NSArray.igrid(N, offset)
            igrid = np.array([np.tile(v, (len(site), 1)) for v in igrid])
            igrid = igrid.reshape(ns * len(site), 3)
            all_site = NSArray(np.tile(site.numpy(), (ns, 1)) + igrid, "vector")
        else:
            all_site = site

        return all_site

    # ==================================================
    def _representative_bond(self, tail, head, info, detail):
        """
        generate representative bonds.

        Args:
            tail (str): tail position.
            head (str): head position.
            info (dict): info dict.
        Returns:
            list: list of representative bonds in each neighbor.
        """
        tail = NSArray("{" + tail + "}")
        head = self._site_grid(head, info, detail, repeat=True)

        all_bond = NSArray.distance(tail, head, self._G)
        if info["molecule"]:
            all_bond = list(all_bond.values())[1:]
        else:
            all_bond = list(all_bond.values())[1 : detail["max_neighbor"] + 1]
        rep_bond = []
        for lst in all_bond:
            th = [f"{tail[0]};{head[i]}" for _, i in lst]
            bs = NSArray.from_str(th)
            rb = self._mpm.group.find_rep_bond(bs)
            rep_bond.append(rb.str())

        return rep_bond

    # ==================================================
    @property
    def A(self):
        return self._A

    # ==================================================
    @property
    def G(self):
        return self._G

    # ==================================================
    @property
    def A_norm(self):
        return self._A_norm

    # ==================================================
    @classmethod
    def regularize(cls, d):
        """
        regularize minimal model dict.

        Args:
            d (dict): dict of model.

        Returns:
            dict: regularized model dict.
        """
        # default setting.
        if "cell" not in d.keys():
            d["cell"] = None

        model = d["model"]

        if "site" not in d.keys():
            d["site"] = {}
        if "bond" not in d.keys():
            d["bond"] = []

        if "option" not in d.keys():
            d["option"] = {
                "view": None,
                "view_mode": "standard",
                "output": model,
                "minimal_samb": True,
                "binary_output": False,
            }
        else:
            if "view" not in d["option"].keys():
                d["option"]["view"] = None
            if "view_mode" not in d["option"].keys():
                d["option"]["view_mode"] = "standard"
            if "output" not in d["option"].keys():
                d["option"]["output"] = model
            if "minimal_samb" not in d["option"].keys():
                d["option"]["minimal_samb"] = True
            if "binary_output" not in d["option"].keys():
                d["option"]["binary_output"] = False

        if "generate" not in d.keys():
            d["generate"] = {
                "model_type": "tight_binding",
                "time_reversal_type": "electric",
                "irrep": None,
                "fourier_transform": False,
                "toroidal_priority": False,
            }
        else:
            if "model_type" not in d["generate"].keys():
                d["generate"]["model_type"] = "tight_binding"
            if "time_reversal_type" not in d["generate"].keys():
                d["generate"]["time_reversal_type"] = "electric"
            if "irrep" not in d["generate"].keys():
                d["generate"]["irrep"] = None
            if "fourier_transform" not in d["generate"].keys():
                d["generate"]["fourier_transform"] = False
            if "toroidal_priority" not in d["generate"].keys():
                d["generate"]["toroidal_priority"] = False

        if "detail" not in d.keys():
            d["detail"] = {
                "rep_bond_all": {},
                "cell_range": _default_search_cell_range,
                "max_neighbor": _default_max_neighbor,
                "A": "[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]",
                "version": __version__,
            }
        else:
            if "cell_range" not in d["detail"].keys():
                d["detail"]["cell_range"] = _default_search_cell_range
            if "max_neighbor" not in d["detail"].keys():
                d["detail"]["max_neighbor"] = _default_max_neighbor

        if "spinful" not in d.keys():
            d["spinful"] = False
        else:
            if d["generate"]["model_type"] == "phonon":
                d["spinful"] = False

        if "k_point" not in d.keys():
            d["k_point"] = {"Γ": "[0, 0, 0]", "X": "[1/2, 0, 0]"}
        if "k_path" not in d.keys():
            d["k_path"] = "Γ-X"

        return d

    # ==================================================
    @classmethod
    def _header(cls):
        return header_str
