"""
PointGroup manages point group.
"""
import sympy as sp
import numpy as np
from gcoreutils.nsarray import NSArray
from multipie.tag.tag_group import TagGroup
from multipie.tag.tag_wyckoff import TagWyckoff
from multipie.tag.tag_list import TagList
from multipie.tag.tag_multipole import TagMultipole
from multipie.symmetry_operation.util.symmetry_operation_util import to_cartesian
from multipie.symmetry_operation.symmetry_operation_g_set import SymmetryOperationGSet
from multipie.symmetry_operation.symmetry_operation_g import SymmetryOperationG
from multipie.wyckoff.wyckoff_g_set import WyckoffGSet
from multipie.wyckoff.wyckoff_g import WyckoffG
from multipie.harmonics.harmonics_pg_set import HarmonicsPGSet
from multipie.harmonics.harmonics_pg import HarmonicsPG
from multipie.character.character_pg_set import CharacterPGSet
from multipie.character.character_pg import CharacterPG
from multipie.virtual_cluster.virtual_cluster_pg_set import VirtualClusterPGSet
from multipie.virtual_cluster.virtual_cluster_pg import VirtualClusterPG
from multipie.clebsch_gordan.clebsch_gordan_pg_set import ClebschGordanPGSet
from multipie.clebsch_gordan.clebsch_gordan_pg import ClebschGordanPG
from multipie.response_tensor.response_tensor_pg_set import ResponseTensorPGSet
from multipie.response_tensor.response_tensor_pg import ResponseTensorPG
from multipie.multipole.base.base_atomic_multipole_dataset import BaseAtomicMultipoleDataset
from multipie.multipole.util.atomic_samb_util import create_atomic_samb
from multipie.multipole.util.atomic_orbital_util import basis_type
from multipie.multipole.util.z_samb_util import create_z_samb
from multipie import get_binary


# ==================================================
class PointGroup:
    """
    point group.

    Attributes:
        tag (TagGroup): point-group tag.
        symmetry_operation (SymmetryOperationG): the class to manage symmetry operations.
        wyckoff (WyckoffG): the class to manage Wyckoff positions.
        harmonics (HarmonicsPG): the class to manage harmonics.
        character (CharacterPG): the class to manage character table.
        virtual_cluster (VirtualClusterPG): the class to manage virtual cluster.
        clebsch_gordan (ClebschGordanPG): the class to manage Clebsch-Gordan coefficients.
        response (ResponseTensorPG) : the class to manage response tensors.
        mp_dataset (BaseAtomicMultipoleDataset) : the dict to manage atomic multipole matrices.
        core (BinaryManager): core for binary data.
    """

    # ==================================================
    def __init__(self, pg_tag=None, core=None, verbose=False):
        """
        initialize the class.

        Args:
            pg_tag (TagGroup or str, optional): point-group tag.
            core (BinaryManager, optional): core for binary data.
            verbose (bool, optional): verbose access to binary data ?

        Notes:
            - if pg_tag is None, "C1" is used.
        """
        if pg_tag is None:
            pg_tag = "C1"
        if type(pg_tag) == str:
            pg_tag = TagGroup(pg_tag)
        assert pg_tag in TagGroup.create(), f"{pg_tag} is not a point-group tag."

        if core is None:
            core = get_binary(verbose=verbose)

        self.tag = pg_tag
        """point-group tag."""

        self.harmonics: HarmonicsPG = core[HarmonicsPGSet][pg_tag]
        """the class to manage harmonics."""

        self.character: CharacterPG = core[CharacterPGSet][pg_tag]
        """the class to manage character table."""

        self.symmetry_operation: SymmetryOperationG = core[SymmetryOperationGSet][pg_tag]
        """the class to manage symmetry operations."""

        self.wyckoff: WyckoffG = core[WyckoffGSet][pg_tag]
        """the class to manage Wyckoff positions."""

        self.virtual_cluster: VirtualClusterPG = core[VirtualClusterPGSet][pg_tag]
        """the class to manage virtual cluster."""

        self.clebsch_gordan: ClebschGordanPG = core[ClebschGordanPGSet][pg_tag]
        """the class to manage Clebsch-Gordan coefficients."""

        self.response: ResponseTensorPG = core[ResponseTensorPGSet][pg_tag]
        """the class to manage response tensors."""

        self.mp_dataset: BaseAtomicMultipoleDataset = core[BaseAtomicMultipoleDataset]
        """the class to manage atomic multipole matrices."""

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
    @property
    def identity_irrep(self):
        """
        identity irreducible representation.

        Returns:
            str: identity irreducible representation.
        """
        return str(self.character.irrep_list[0])

    # ==================================================
    def irrep_decomposition_table(self):
        """
        product irreducible decomposition table.

        Returns:
            dict: (irrep1,irrep2): {irrep:n}, { (str, str): [(int,str)] }.
        """
        dic = {}
        for irrep1 in self.character.key_list():
            for irrep2 in self.character.key_list():
                d = {str(ir): n for n, ir in self.character.symmetric_product_decomposition((irrep1, irrep2), False)}
                if irrep1 == irrep2:
                    for n, ir in self.character.anti_symmetric_product_decomposition(irrep1, False):
                        d[str(ir)] = d.get(str(ir), 0) + n
                dic[(str(irrep1), str(irrep2))] = d

        return dic

    # ==================================================
    def transform_matrix_site(self, site, cc_only=False):
        """
        trasform matrix for site basis (reduced coordinate).

        Args:
            site (str or NSArray): site, "[x,y,z]".
            cc_only (bool, optional): conjugacy-class operations only ?

        Returns: tuple
            - NSArray: symmetry operation matrices.
            - NSArray: site basis.
        """
        if not isinstance(site, (str, NSArray)):
            raise KeyError(f"{type(site)} is not accepted for site.")

        s = self.symmetry_operation._equivalent_vector(site, remove_duplicate=True)
        op = self.symmetry_operation.mat(cc_only=cc_only)

        tm = []
        for ii in range(len(op)):
            opi = op[ii]
            tmg = np.zeros((len(s), len(s))).astype(int).tolist()
            for i in range(len(s)):
                si = s[i].transform(opi)
                tmg[s.index(si)][i] = 1
            tm.append(tmg)

        tm = NSArray(tm, "matrix")

        return tm, s

    # ==================================================
    def transform_matrix_bond(self, bond, nondirectional=False, cc_only=False):
        """
        trasform matrix for bond basis (reduced coordinate).

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

        b = self.symmetry_operation._equivalent_bond(bond, nondirectional=nondirectional, remove_duplicate=True)
        op = self.symmetry_operation.mat(cc_only=cc_only)

        tm = []
        for ii in range(len(op)):
            opi = op[ii]
            tmg = np.zeros((len(b), len(b))).astype(int).tolist()
            for i in range(len(b)):
                bi = b[i].transform(opi)
                if nondirectional:
                    j = b.index(bi)
                    if j is None:
                        j = b.index(bi.reverse_direction())
                else:
                    j = b.index(bi)
                tmg[j][i] = 1
            tm.append(tmg)

        tm = NSArray(tm, "matrix")

        return tm, b

    # ==================================================
    def transform_matrix_vector(self, axial=False, cc_only=False):
        """
        trasform matrix for vector basis (cartesian coordinate).

        Args:
            axial (bool, optional); axial vector ?
            cc_only (bool, optional): conjugacy-class operations only ?

        Returns: tuple
            - NSArray: symmetry operation matrices.
            - NSArray: vector basis.
        """
        tm = self.symmetry_operation.mat(cc_only=cc_only, axial=axial, cartesian=True)
        head = "G" if axial else "Q"
        v = NSArray.vector3d(head=head)

        return tm, v

    # ==================================================
    def transform_matrix_orbital(self, axial=False, cc_only=False, **kwargs):
        """
        transform matrix for harmonics basis (cartesian coordinate).

        Args:
            axial (bool, optional): axial harmonics ?
            cc_only (bool, optional): conjugacy-class operations only ?
            kwargs (dict, optional): select conditions for harmonics, (rank/irrep/mul/comp).

        Returns: tuple
            - NSArray: symmetry operation matrices.
            - NSArray: harmonics basis.
        """
        head = "G" if axial else "Q"
        orbitals = self.harmonics.select(head=head, **kwargs)
        v = NSArray.vector3d(head=head)
        t_pos = self.symmetry_operation._equivalent_vector(v, cc_only=cc_only, axial=False, cartesian=True)

        hset = NSArray([i.expression(v=v) for i in orbitals], "scalar")

        if axial:
            p_det = self.symmetry_operation.mat(axial=False, cc_only=cc_only).det()
            tm = self.symmetry_operation._transform_matrix_orbital_basis(v, t_pos, hset, p_det)
        else:
            tm = self.symmetry_operation._transform_matrix_orbital_basis(v, t_pos, hset)

        return tm, hset

    # ==================================================
    def transform_site(self, site, remove_duplicate=False):
        """
        set of transformed sites (reduced coordinate).

        Args:
            site (str or NSArray): site, "[x,y,z]".
            remove_duplicate (bool, optional): remove the same sites ?

        Returns:
            NSArray: transformed sites (in order of symmetry operations).

        Notes:
            - when remove_duplicate is True, return values in the same order as site_mapping keys.
        """
        if not isinstance(site, (str, NSArray)):
            raise KeyError(f"{type(site)} is not accepted for site.")

        if remove_duplicate:
            s = NSArray.from_str(self.site_mapping(site).keys())
        else:
            s = self.symmetry_operation._equivalent_vector(site, remove_duplicate=False)

        return s

    # ==================================================
    def transform_bond(self, bond, nondirectional=False, remove_duplicate=False):
        """
        set of transformed bonds (reduced coordinate).

        Args:
            bond (str or NSArray): bond, "vector@center", "tail;head", "start:vector".
            nondirectional (bool, optional): ignore directional property ?
            remove_duplicate (bool, optional): remove the same bonds ?

        Returns:
            NSArray: transformed bonds (in order of symmetry operations).

        Notes:
            - when remove_duplicate and nondirectional are True, return values in the same order as bond_mapping keys.
        """
        if not isinstance(bond, (str, NSArray)):
            raise KeyError(f"{type(bond)} is not accepted for bond.")

        if remove_duplicate and nondirectional:
            b = NSArray.from_str(self.bond_mapping(bond)[0].keys())
        else:
            b = self.symmetry_operation._equivalent_bond(bond, nondirectional=nondirectional, remove_duplicate=remove_duplicate)

        return b

    # ==================================================
    def transform_vector(self, vector, axial=False):
        """
        set of transformed vectors (cartesian coordinate).

        Args:
            vector (str or NSArray): vector, "[x,y,z]".
            axial (bool, optional); axial vector ?

        Returns:
            NSArray: transformed vectors (in order of symmetry operations).
        """
        if not isinstance(vector, (str, NSArray)):
            raise KeyError(f"{type(vector)} is not accepted for vector.")

        if type(vector) == str:
            vector = NSArray(vector)
        op = self.symmetry_operation.mat(cc_only=False, axial=axial, cartesian=True)
        vector = op.apply(vector)

        return vector

    # ==================================================
    def transform_orbital(self, orbital, axial=False):
        """
        set of transformed orbitals (cartesian coordinate).

        Args:
            orbital (str): (xyz)-polynomial expression.
            axial (bool, optional): axial orbital ?

        Returns:
            NSArray: transformed orbitals (in order of symmetry operations).
        """
        v = NSArray.vector3d()
        tm = self.symmetry_operation.mat(cc_only=False, axial=False, cartesian=True).inverse()
        tr = tm.apply(v)

        h = NSArray(orbital).subs({"x": v[0], "y": v[1], "z": v[2]})
        if axial:
            so = self.symmetry_operation.mat(axial=False, cc_only=False).det()
            t_orb = [(so[i] * self.symmetry_operation._transform_variable(h, v, tr[i])).tolist() for i in range(len(tr))]
        else:
            t_orb = [self.symmetry_operation._transform_variable(h, v, tr[i]).tolist() for i in range(len(tr))]

        return t_orb

    # ==================================================
    def irrep_decomposition(self, point_group, rank, axial=False):
        """
        irrep. decomposition of harmonics in terms of given point group.

        Args:
            point_group (str): point-group tag to be decompsed.
            rank (int): rank of harmonics.
            axial (bool, optional): axia harmonics ?

        Returns:
            list: decomposition info., [(harm_tag, [(coeff,harmonics basis of this group)].
        """
        head = "G" if axial else "Q"

        hset = self.harmonics.select(rank=rank, head=head)  # self

        pg = PointGroup(point_group, core=self.core)
        other = pg.harmonics.select(rank=rank, head=head)  # other

        decomp = []
        for h in hset:
            Uself = h.u_matrix()
            c = []
            for o in other:
                Uother = o.u_matrix()
                ci = NSArray.dot(Uother, Uself)
                if ci != 0:
                    c.append((ci, str(o)))
            decomp.append((str(h), c))

        return decomp

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
                r = set(self.transform_site(s, remove_duplicate=True).str())
                rep_site.append(s)
                return remove_equivalent_site(sl - r, rep_site)

        assert len(site_list) > 0, "empty site list is given."

        site_list = set(site_list.str())

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
                r = self.transform_bond(b, remove_duplicate=True)
                r = NSArray.concat([r, r.reverse_direction()])
                r = set(r.str())
                rep_bond.append(b)
                return remove_equivalent_bond(bl - r, rep_bond)

        assert len(bond_list) > 0, "empty bond list is given."

        # remove equivalent bonds.
        tail = bond_list[0].convert_bond("bond_th")[0]
        bond_list = set(bond_list.str())

        # find representative bonds.
        rep_bond = []
        rep_bond = remove_equivalent_bond(bond_list, rep_bond)
        rep_bond = NSArray.from_str(rep_bond)

        # set bonds with given tail.
        rep_bond1 = []
        for bond in rep_bond:
            bonds = self.transform_bond(bond, remove_duplicate=True)
            t, h = bonds.convert_bond("bond_th")
            idx = t.index(tail)
            if idx is None:
                idx = h.index(tail)
                assert idx is not None, "not reach here."
                b = bonds[idx].reverse_direction()
            else:
                b = bonds[idx]
            v, c = b.convert_bond("bond")
            rep_bond1.append(f"{v}@{c}")
        rep_bond = NSArray.from_str(rep_bond1)

        return rep_bond

    # ==================================================
    def site_mapping(self, site):
        """
        mapping between site and symmetry operations (reduced coordinate).

        Args:
            site (str or NSArray): site, "[x,y,z]".

        Returns:
            dict: mapping from site string to symmetry-operation IDs.

        Notes:
            - sites are sorted by 1st SO.
        """
        if not isinstance(site, (str, NSArray)):
            raise KeyError(f"{type(site)} is not accepted for site.")

        s = self.symmetry_operation._equivalent_vector(site)

        lst = {}
        for i in range(len(s)):
            si = str(s[i])
            lst[si] = lst.get(si, []) + [i]

        # sort in order of 1st component of SO.
        lst = dict(sorted(lst.items(), key=lambda i: i[1][0]))

        return lst

    # ==================================================
    def bond_mapping(self, bond):
        """
        mapping between bond and symmetry operations (reduced coordinate).

        Args:
            bond (str or NSArray): bond, "vector@center", "tail;head", "start:vector".

        Returns: tuple
            - dict: mapping from bond to symmetry-operation IDs. bond = "tail;head".
            - bool: nondirectional ?

        Notes:
            - if bond direction is reversed, symmetry operation ID with minus sign.
            - bonds are sorted by 1st SO.
        """
        if not isinstance(bond, (str, NSArray)):
            raise KeyError(f"{type(bond)} is not accepted for bond.")
        if type(bond) == str:
            bond = NSArray(bond)

        b_d = self.symmetry_operation._equivalent_bond(bond, nondirectional=False, remove_duplicate=True)
        b_all = self.symmetry_operation._equivalent_bond(bond, nondirectional=False, remove_duplicate=False)
        b_nd = self.symmetry_operation._equivalent_bond(bond, nondirectional=True, remove_duplicate=True)

        nd = len(b_nd) != len(b_d)
        b = b_nd if nd else b_d

        tail0 = bond.convert_bond("bond_th")[0]
        tail0 = NSArray.from_str(self.site_mapping(tail0).keys())
        head0 = bond.convert_bond("bond_th")[1]
        head0 = NSArray.from_str(self.site_mapping(head0).keys())

        # set bond direction with the same tail as given one.
        if str(tail0.sort()) == str(head0.sort()):
            for i in range(len(b)):
                t, h = b[i].convert_bond("bond_th")
                t_idx = tail0.index(t)
                h_idx = tail0.index(h)
                if t_idx > h_idx:
                    b[i] = b[i].reverse_direction()

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

        wpv = self.symmetry_operation._equivalent_vector(site, remove_duplicate=True)
        n = len(wpv)
        tags = [i for i in self.wyckoff.keys() if i.n == n]
        if len(tags) == 1:
            return tags[0]

        wps = [self.wyckoff.position(tag) for tag in tags]
        w0 = wpv[0]
        for tag, wp in zip(tags, wps):
            for i in wp:
                if str(i) == str(w0):
                    return tag
                if sp.solve(list(i - w0), self.wyckoff.v.tolist(), manual=True):
                    return tag

        raise Exception(f"{str(site)} cannot be found.")

    # ==================================================
    def virtual_cluster_basis(self, wyckoff=None, ortho=True, create=False):
        """
        virtual-cluster basis set.

        Args:
            wyckoff (str, optional): Wyckoff position.
            ortho (bool, optional): orthogonalize ?
            create (bool, optional): create from scratch ?

        Returns: tuple
            - TagMultipole or NSArray: virtual-cluster basis set, (multipole tag, basis vector).
            - NSArray: site basis (reduced coordinate).
        """
        if wyckoff is not None and TagWyckoff(wyckoff) in self.wyckoff.keys():
            wp = wyckoff
        else:
            wp = self.wyckoff.key_list()[-1]

        if not create:
            sites = self.virtual_cluster._site[str(wp)]
            vc_basis = self.virtual_cluster._basis[str(wp)]
            return vc_basis, sites

        sites = self.wyckoff.position(wp, default=True)
        sites_c = to_cartesian(self.symmetry_operation.crystal, sites)

        basis = []
        info = []
        for rank in range(12):
            hs = self.harmonics.select(rank=rank, head="Q")
            for h in hs:
                fp = NSArray([h.expression(v=sites_c[p]).tolist() for p in range(len(sites_c))], "vector").expand()
                if np.any(fp != sp.S(0)):
                    basis.append(fp)
                    info.append(h.tag)

        basis = NSArray.concat(basis)

        if ortho:
            basis, idx = NSArray.orthogonalize(basis, len(sites))
            basis = basis[idx]
            info = [info[i] for i in idx]

        vc_basis = {i: basis[no].expand() for no, i in enumerate(info)}

        return vc_basis, sites

    # ==================================================
    def site_cluster_samb(self, site):
        """
        create site-cluster multipole basis set.

        Args:
            site (str or NSArray): representative site, "[x,y,z]".

        Returns: tuple.
            - dict: site-cluster SAMB, {TagMultipole: NSArray(vector)}.
            - NSArray: sites, NSArray(vector).
        """
        if not isinstance(site, (str, NSArray)):
            raise KeyError(f"{type(site)} is not accepted for site.")

        site_map = self.site_mapping(site)
        vc = self.virtual_cluster

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

        Returns: tuple.
            - dict: bond-cluster SAMB, {TagMultipole: NSArray(vector)}.
            - NSArray: bonds, NSArray(bond).
        """
        if not isinstance(bond, (str, NSArray)):
            raise KeyError(f"{type(bond)} is not accepted for bond.")

        bond_map, _ = self.bond_mapping(bond)
        vc = self.virtual_cluster

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
    def atomic_samb(self, bra_list, ket_list, spinful=False):
        """
        create atomic multipole basis set.

        Args:
            bra_list ([str]): atomic basis set for bra. see the examples in :class:`.AtomicBasisSet`
            ket_list ([str]): atomic basis set for ket. see the examples in :class:`.AtomicBasisSet`
            spinful (bool, optional): spinful ?

        Returns:
            dict: atomic SAMB, {TagMultipole: NSArray(matrix)}.
        """
        b_type = basis_type(bra_list[0], self.tag.crystal)

        bam = self.mp_dataset[b_type]
        hs = self.harmonics

        a_samb = create_atomic_samb(bra_list, ket_list, spinful, self.tag.crystal, bam, hs, ortho=True)

        return a_samb

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

        cg = self.clebsch_gordan
        hs = self.harmonics

        z_samb = create_z_samb(cg, hs, x_tag_list, y_tag_list, toroidal_priority, **kwargs)

        return z_samb

    # ==================================================
    @classmethod
    def spherical_atomic_multipole_basis(cls, bra_list, ket_list, spinful=False, crystal="cubic", core=None, verbose=False):
        """
        create spherical atomic multipole basis set.

        Args:
            bra_list ([str]): atomic basis set for bra. see the examples in :class:`.AtomicBasisSet`
            ket_list ([str]): atomic basis set for ket. see the examples in :class:`.AtomicBasisSet`
            spinful (bool, optional): spinful ?
            crystal (str, optional): seven crystal systems, triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic.
            core (BinaryManager, optional): core for binary data.
            verbose (bool, optional): verbose parallel info and verbose access to binary data ?

        Returns:
            dict: spherical atomic SAMB, {TagMultipole: NSArray(matrix)}.
        """
        if core is None:
            core = get_binary(verbose=verbose)

        b_type = basis_type(bra_list[0], crystal)

        bam = core[BaseAtomicMultipoleDataset][b_type]
        hs = None

        a_samb = create_atomic_samb(bra_list, ket_list, spinful, crystal, bam, hs, ortho=False)

        return a_samb
