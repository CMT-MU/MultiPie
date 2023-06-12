"""
SpaceGroup manages space group.
"""
import sympy as sp
import numpy as np
from gcoreutils.nsarray import NSArray
from multipie.data.data_no_space_group import _data_no_space_group
from multipie.tag.tag_group import TagGroup
from multipie.tag.tag_list import TagList
from multipie.tag.tag_multipole import TagMultipole
from multipie.group.point_group import PointGroup
from multipie.symmetry_operation.util.symmetry_operation_util import to_conventional, to_primitive
from multipie.symmetry_operation.symmetry_operation_g import SymmetryOperationG
from multipie.symmetry_operation.symmetry_operation_g_set import SymmetryOperationGSet
from multipie.wyckoff.wyckoff_g_set import WyckoffGSet
from multipie.wyckoff.wyckoff_g import WyckoffG
from multipie.multipole.util.z_samb_util import create_z_samb
from multipie import get_binary


# ==================================================
class SpaceGroup:
    """
    space group.

    Attributes:
        tag (TagGroup): space-group tag.
        symmetry_operation (SymmetryOperationG): the class to manage symmetry operations.
        wyckoff (WyckoffG): the class to manage Wyckoff positions.
        pg (PointGroup): associated point group.
        core (BinaryManager): core for binary data.
    """

    # ==================================================
    def __init__(self, sg_tag=None, core=None, verbose=False):
        """
        initialize the class.

        Args:
            sg_tag (TagGroup or str or int, optional): tag of space group.
            core (BinaryManager, optional): core for binary data.
            verbose (bool, optional): verbose access to binary data ?

        Notes:
            - if sg_tag is None, "C1^1" is used.
        """
        if sg_tag is None:
            sg_tag = "C1^1"
        if type(sg_tag) == int:
            sg_tag = _data_no_space_group[sg_tag]
        if type(sg_tag) == str:
            sg_tag = TagGroup(sg_tag)
        assert sg_tag in TagGroup.create(space_group=True), f"{sg_tag} is not a space-group tag."

        if core is None:
            core = get_binary(verbose=verbose)

        self.tag = sg_tag
        """space-group tag."""

        self.symmetry_operation: SymmetryOperationG = core[SymmetryOperationGSet][sg_tag]
        """the class to manage symmetry operations."""

        self.wyckoff: WyckoffG = core[WyckoffGSet][sg_tag]
        """the class to manage Wyckoff positions."""

        self.pg = PointGroup(self.tag.pg, core=core, verbose=verbose)
        """associated point group."""

        self.core = core
        """binary data."""

    # ==================================================
    def __str__(self):
        return str(self.tag)

    # ==================================================
    def __repr__(self):
        return repr(self.tag)

    # ==================================================
    def latex(self):
        return self.tag.latex()

    # ==================================================
    def transform_matrix_site(self, site, cc_only=False):
        """
        transform matrix for site basis (reduced coordinate).

        Args:
            site (str or NSArray): site, "[x,y,z]".
            cc_only (bool, optional): conjugacy-class operations only ?

        Returns: tuple
            - NSArray: symmetry operation matrices.
            - NSArray: site basis.
        """
        if not isinstance(site, (str, NSArray)):
            raise KeyError(f"{type(site)} is not accepted for site.")

        s = self.symmetry_operation._equivalent_vector(site, primitive=True, shift=True, remove_duplicate=True)
        op = self.symmetry_operation.mat(cc_only=cc_only, primitive=True)

        tm = []
        for ii in range(len(op)):
            opi = op[ii]
            tmg = np.zeros((len(s), len(s))).astype(int).tolist()
            for i in range(len(s)):
                si = s[i].transform(opi).shift()
                tmg[s.index(si)][i] = 1
            tm.append(tmg)

        tm = NSArray(tm, "matrix")
        s = to_conventional(self.symmetry_operation.lattice, s).shift()

        return tm, s

    # ==================================================
    def transform_matrix_bond(self, bond, nondirectional=False, cc_only=False):
        """
        transform matrix for bond basis (reduced coordinate).

        Args:
            bond (str or NSArray): bond, "vector@center", "tail;head", "start:vector".
            nondirectional (bool, optional): ignore directional property ?
            cc_only (bool, optional): conjugacy-class operations only ?

        Returns: tuple
            - NSArray: symmetry operation matrices.
            - NSArray: bond basis.
        """
        if not isinstance(bond, (str, NSArray)):
            raise KeyError(f"{type(bond)} is not accepted for bond.")

        b = self.symmetry_operation._equivalent_bond(
            bond, primitive=True, nondirectional=nondirectional, shift=True, remove_duplicate=True
        )
        op = self.symmetry_operation.mat(cc_only=cc_only, primitive=True)

        tm = []
        for ii in range(len(op)):
            opi = op[ii]
            tmg = np.zeros((len(b), len(b))).astype(int).tolist()
            for i in range(len(b)):
                bi = b[i].transform(opi).shift()
                if nondirectional:
                    j = b.index(bi)
                    if j is None:
                        j = b.index(bi.reverse_direction())
                else:
                    j = b.index(bi)
                tmg[j][i] = 1
            tm.append(tmg)

        tm = NSArray(tm, "matrix")
        b = to_conventional(self.symmetry_operation.lattice, b).shift()

        return tm, b

    # ==================================================
    def transform_site(self, site, shift=False, remove_duplicate=False, plus_set=False):
        """
        set of transformed sites (reduced coordinate).

        Args:
            site (str or NSArray): site, "[x,y,z]".
            shift (bool, optional): transformed sites are moved into unit cell ?
            remove_duplicate (bool, optional): remove the same sites ?
            plus_set (bool, optional): add partial translations ?

        Returns:
            NSArray: transformed sites (in order of symmetry operations, and plus_set).

        Notes:
            - when remove_duplicate and shift are True, return values in the same order as site_mapping keys.
        """
        if not isinstance(site, (str, NSArray)):
            raise KeyError(f"{type(site)} is not accepted for site.")

        if remove_duplicate and shift:
            s = NSArray.from_str(self.site_mapping(site, plus_set).keys())
        else:
            s = self.symmetry_operation._equivalent_vector(
                site, plus_set=plus_set, shift=shift, remove_duplicate=remove_duplicate
            )

        return s

    # ==================================================
    def transform_bond(self, bond, nondirectional=False, shift=False, remove_duplicate=False, plus_set=False):
        """
        set of transformed bonds (reduced coordinate).

        Args:
            bond (str or NSArray): bond, "vector@center", "tail;head", "start:vector".
            nondirectional (bool, optional): ignore directional property ?
            shift (bool, optional): transformed bond starting points are moved into unit cell ?
            remove_duplicate (bool, optional): remove the same bonds ?
            plus_set (bool, optional): add partial translations ?

        Returns:
            NSArray: transformed bonds (in order of symmetry operations, and plus_set).

        Notes:
            - when remove_duplicate, nondirectional and shift are True, return values in the same order as bond_mapping keys.
        """
        if not isinstance(bond, (str, NSArray)):
            raise KeyError(f"{type(bond)} is not accepted for bond.")

        if remove_duplicate and nondirectional and shift:
            b = NSArray.from_str(self.bond_mapping(bond, plus_set=plus_set)[0].keys())
        else:
            b = self.symmetry_operation._equivalent_bond(
                bond, plus_set=plus_set, nondirectional=nondirectional, shift=shift, remove_duplicate=remove_duplicate
            )

        return b

    # ==================================================
    def shift_home_unit_cell(self, site):
        """
        shift to home unit cell.

        Args:
            site (str or NSArray): site to be shifted.

        Returns:
            NSArray: shifted site.
        """
        if type(site) == str:
            site = NSArray(site)
        lattice = self.symmetry_operation.lattice
        site = to_conventional(lattice, to_primitive(lattice, site).shift()).shift()

        return site

    # ==================================================
    def find_rep_site(self, site_list):
        """
        find representative sites from site_list.

        Args:
            site_list (NSArray): site list.

        Returns:
            NSArray: representative sites.
        """

        def remove_equivalent_site(sl, rep_site):
            if not sl:
                return rep_site
            else:
                s = list(sl)[0]
                r = set(self.transform_site(s, shift=True, remove_duplicate=True).str())
                rep_site.append(s)
                return remove_equivalent_site(sl - r, rep_site)

        assert len(site_list) > 0, "empty site list is given."

        lattice = self.symmetry_operation.lattice
        site_list = set(to_conventional(lattice, to_primitive(lattice, site_list).shift().remove_duplicate()).shift().str())

        rep_site = []
        rep_site = remove_equivalent_site(site_list, rep_site)
        rep_site = NSArray.from_str(rep_site)

        return rep_site

    # ==================================================
    def find_rep_bond(self, bond_list):
        """
        find representative bonds from bond_list.

        Args:
            bond_list (NSArray): bond list.

        Returns:
            NSArray: representative bonds.
        """

        def remove_equivalent_bond(bl, rep_bond):
            if not bl:
                return rep_bond
            else:
                b = list(bl)[0]
                r = self.transform_bond(b, shift=True, remove_duplicate=True)
                r = NSArray.concat([r, r.reverse_direction()])
                r = set(r.str())
                rep_bond.append(b)
                return remove_equivalent_bond(bl - r, rep_bond)

        assert len(bond_list) > 0, "empty bond list is given."

        # remove equivalent bonds.
        tail = bond_list[0].convert_bond("bond_th")[0]
        lattice = self.symmetry_operation.lattice
        bond_list = set(to_conventional(lattice, to_primitive(lattice, bond_list).shift().remove_duplicate()).shift().str())

        # find representative bonds.
        rep_bond = []
        rep_bond = remove_equivalent_bond(bond_list, rep_bond)
        rep_bond = NSArray.from_str(rep_bond)

        # set bonds with given tail.
        tail = to_primitive(lattice, tail).shift()
        rep_bond1 = []
        for bond in rep_bond:
            bonds = self.symmetry_operation._equivalent_bond(bond, remove_duplicate=True, primitive=True, shift=True)
            t, h = bonds.convert_bond("bond_th")
            idx = t.shift().index(tail)
            if idx is None:
                idx = h.shift().index(tail)
                assert idx is not None, "not reach here."
                b = bonds[idx].reverse_direction()
            else:
                b = bonds[idx]
            v, c = b.shift().convert_bond("bond")
            b = NSArray(f"{v}@{c}")
            b = to_conventional(lattice, b)
            rep_bond1.append(str(b))
        rep_bond = NSArray.from_str(rep_bond1)

        return rep_bond

    # ==================================================
    def site_mapping(self, site, plus_set=False):
        """
        mapping between site and symmetry operations (reduced coordinate).

        Args:
            site (str or NSArray): site, "[x,y,z]".
            plus_set (bool, optional): add partial translations ?

        Returns:
            dict: mapping from site string to symmetry-operation IDs.

        Notes:
            - sites are sorted by 1st SO.
            - if plus_set is True, return [set(t0),set(t1),...]
        """
        if not isinstance(site, (str, NSArray)):
            raise KeyError(f"{type(site)} is not accepted for vector/site.")

        s = self.symmetry_operation._equivalent_vector(site, primitive=True, shift=True)

        lst = {}
        for i in range(len(s)):
            si = str(s[i])
            lst[si] = lst.get(si, []) + [i]

        # sort in order of 1st component of SO.
        lst = dict(sorted(lst.items(), key=lambda i: i[1][0]))

        # convert to conventional cell.
        s = NSArray.from_str(lst.keys())
        s = to_conventional(self.symmetry_operation.lattice, s, plus_set=True).shift()

        if plus_set:
            mpp = list(lst.values()) * len(self.symmetry_operation.plus_set)
            lst = {str(s[i]): mpp[i] for i in range(len(s))}
        else:
            mpp = list(lst.values())
            lst = {str(s[i]): mpp[i] for i in range(len(lst))}

        return lst

    # ==================================================
    def bond_mapping(self, bond, plus_set=False):
        """
        mapping between bond and symmetry operations (reduced coordinate).

        Args:
            bond (str or NSArray): bond, "vector@center", "tail;head", "start:vector".
            plus_set (bool, optional): add partial translations ?

        Returns: tuple
            - dict: mapping from bond to symmetry-operations IDs. bond = "tail;head"
            - bool: nondirectional ?

        Notes:
            - if bond direction is reversed, symmetry operation ID with minus sign.
            - bonds are sorted by 1st SO.
            - if plus_set is True, return [set(t0),set(t1),...]
        """
        if not isinstance(bond, (str, NSArray)):
            raise KeyError(f"{type(bond)} is not accepted for bond.")
        if type(bond) == str:
            bond = NSArray(bond)

        b_d = self.symmetry_operation._equivalent_bond(
            bond, primitive=True, shift=True, nondirectional=False, remove_duplicate=True
        )
        b_all = self.symmetry_operation._equivalent_bond(
            bond, primitive=True, shift=True, nondirectional=False, remove_duplicate=False
        )
        b_nd = self.symmetry_operation._equivalent_bond(
            bond, primitive=True, shift=True, nondirectional=True, remove_duplicate=True
        )

        nd = len(b_nd) != len(b_d)
        b = b_nd if nd else b_d

        tail0 = bond.convert_bond("bond_th")[0]
        tail0 = NSArray.from_str(self.site_mapping(tail0).keys())
        tail0 = to_primitive(self.symmetry_operation.lattice, tail0).shift()
        head0 = bond.convert_bond("bond_th")[1]
        head0 = NSArray.from_str(self.site_mapping(head0).keys())
        head0 = to_primitive(self.symmetry_operation.lattice, head0).shift()

        # set bond direction with the same tail as given one.
        if str(tail0.sort()) == str(head0.sort()):
            for i in range(len(b)):
                t, h = b[i].convert_bond("bond_th")
                t_idx = tail0.index(t.shift())
                h_idx = tail0.index(h.shift())
                if t_idx > h_idx:
                    b[i] = b[i].reverse_direction()

        # create mapping in primitive cell.
        lst1 = {}
        for i in range(len(b_all)):
            bi = b_all[i]
            bib = str(bi.reverse_direction())
            bi = str(bi)
            if b.index(bi) is not None:
                lst1[bi] = lst1.get(bi, []) + [i + 1]
            elif b.index(bib) is not None:
                lst1[bib] = lst1.get(bib, []) + [-(i + 1)]
            else:
                assert False, "not reach here."

        # reverse SO if identity SO idx are negative.
        if sum(list(lst1.values()), []).count(-1) > 0:
            lst = {i: [-k for k in j] for i, j in lst1.items()}
        else:
            lst = lst1

        # reset SO number staring from zero.
        lsts = {}
        for i, j in lst.items():
            lsts[i] = [k - 1 if k > 0 else k + 1 for k in j]

        # sort in order of 1st component of SO.
        lsts = dict(sorted(lsts.items(), key=lambda i: abs(i[1][0])))

        # convert to conventional cell.
        b = NSArray.from_str(lsts.keys())
        b = to_conventional(self.symmetry_operation.lattice, b, plus_set=plus_set).shift()

        if plus_set:
            mp = list(lsts.values()) * len(self.symmetry_operation.plus_set)
            lsts = {str(b[i]): mp[i] for i in range(len(b))}

        return lsts, nd

    # ==================================================
    def find_wyckoff_position(self, site):
        """
        find Wyckoff position for a given site (reduced coordinate).

        Args:
            site (str or NSArray): site, "[x,y,z]".

        Returns:
            str: Wyckoff position.
        """
        if not isinstance(site, (str, NSArray)):
            raise KeyError(f"{type(site)} is not accepted for site.")

        wpv = self.symmetry_operation._equivalent_vector(site, plus_set=True, shift=True, remove_duplicate=True)
        n = len(wpv)
        tags = [i for i in self.wyckoff.keys() if i.n == n]
        if len(tags) == 1:
            return tags[0]

        wps = [self.wyckoff.position(tag) for tag in tags]
        wpv0 = self.symmetry_operation._equivalent_vector(site, primitive=True, shift=True, remove_duplicate=True)
        w0 = to_conventional(self.symmetry_operation.lattice, wpv0[0:1])[0]
        sht = [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]
        for tag, wp in zip(tags, wps):
            for i in wp:
                for s in sht:
                    s = NSArray(s, "vector")
                    if str(i) == str(w0 + s):
                        return tag
                    if sp.solve(list(i - w0 + s), self.wyckoff.v.tolist(), manual=True):
                        return tag

        raise Exception(f"{str(site)} cannot be found.")

    # ==================================================
    def site_cluster_samb(self, site):
        """
        create site-cluster multipole basis set.

        Args:
            site (str or NSArray): representative site, "[x,y,z]".

        Returns: tuple,
            - dict: site-cluster SAMB, {TagMultipole: NSArray(vector)}.
            - NSArray: sites, NSArray(vector).
        """
        if not isinstance(site, (str, NSArray)):
            raise KeyError(f"{type(site)} is not accepted for site.")

        site_map = self.site_mapping(site)
        vc = self.pg.virtual_cluster

        info = vc.key_list()
        bs = []
        n = len(site_map.keys())
        for basis in vc.values():
            s_basis = NSArray(str([0 for _ in range(n)]))
            for si, nos in enumerate(site_map.values()):
                for no in nos:
                    s_basis[si] += basis[no]
            bs.append(s_basis)

        bs = NSArray.concat(bs)
        bs, idx = NSArray.orthogonalize(bs)

        sc_samb = {info[i].replace(m_type="s"): bs[i].expand().simplify() for i in idx}
        site = NSArray.from_str(site_map.keys())

        return sc_samb, site

    # ==================================================
    def bond_cluster_samb(self, bond):
        """
        create bond-cluster multipole basis set.

        Args:
            bond (str or NSArray): bond, "vector@center", "tail;head", "start:vector".

        Returns: tuple,
            - dict: bond-cluster SAMB, {TagMultipole: NSArray(vector)}.
            - NSArray: bonds, NSArray(bond).
        """
        if not isinstance(bond, (str, NSArray)):
            raise KeyError(f"{type(bond)} is not accepted for bond.")

        bond_map, _ = self.bond_mapping(bond)
        vc = self.pg.virtual_cluster

        info = vc.key_list()
        n = len(bond_map.keys())
        bs_s = []
        bs_a = []
        for basis in vc.values():
            bs_basis = NSArray(str([0 for _ in range(n)]))
            ba_basis = NSArray(str([0 for _ in range(n)]))
            for bo, nos in enumerate(bond_map.values()):
                for no in nos:
                    bs_basis[bo] += basis[abs(no)]
                    if no < 0:
                        ba_basis[bo] -= basis[abs(no)]
                    else:
                        ba_basis[bo] += basis[abs(no)]
            bs_s.append(bs_basis)
            bs_a.append(ba_basis)

        bs_s = NSArray.concat(bs_s)
        bs_a = NSArray.concat(bs_a)
        bs_s, idx_s = NSArray.orthogonalize(bs_s)
        bs_a, idx_a = NSArray.orthogonalize(bs_a)

        sg_basis_s = {info[i].replace(m_type="b"): bs_s[i].expand().simplify() for i in idx_s}
        sg_basis_a = {info[i].replace(m_type="b", head="T"): sp.I * bs_a[i].expand().simplify() for i in idx_a}

        bc_samb = sg_basis_s | sg_basis_a
        bond = NSArray.from_str(bond_map.keys())

        return bc_samb, bond

    # ==================================================
    def z_samb(self, x_tag_list, y_tag_list, toroidal_priority=False, **kwargs):
        """
        create combined multipole basis set.

        Args:
            x_tag_list (str/TagMultipole/[TagMultipole]):  multipole/harmonics tag list.
            y_tag_list (str/TagMultipole/[TagMultipole]):  multipole/harmonics tag list.
            toroidal_priority (bool, optional): create toroidal multipoles (G,T) in priority? else prioritize conventional multipoles (Q,M).
            kwargs (dict, optional): select condition for multipoles, keywords in TagMultipole except for head.

        Returns:
            dict: {(TagMultipole, TagMultipole(atomic), TagMultipole(site/bond)): [(coefficient, TagMultipole(atomic), TagMultipole(site/bond)] }.
        """
        if type(x_tag_list) == str:
            x_tag_list = TagList([TagMultipole(x_tag_list)])
        if type(x_tag_list) in (TagMultipole, list):
            x_tag_list = TagList(x_tag_list)

        if type(y_tag_list) == str:
            y_tag_list = TagList([TagMultipole(y_tag_list)])
        if type(y_tag_list) in (TagMultipole, list):
            y_tag_list = TagList(y_tag_list)

        cg = self.pg.clebsch_gordan
        hs = self.pg.harmonics

        z_samb = create_z_samb(cg, hs, x_tag_list, y_tag_list, toroidal_priority, **kwargs)

        return z_samb
