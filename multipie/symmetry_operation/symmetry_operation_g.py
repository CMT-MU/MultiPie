"""
SymmetryOperationG manages point/space-group symmetry operations.
"""
import sympy as sp
from gcoreutils.nsarray import NSArray
from multipie.tag.tag_list import TagList
from multipie.tag.tag_group import TagGroup
from multipie.tag.tag_symmetry_operation import TagSymmetryOperation
from multipie.symmetry_operation.util.symmetry_operation_util import latticeT, to_cartesian, to_primitive, to_conventional
from multipie.symmetry_operation.symmetry_operation import SymmetryOperation
from multipie.data.data_symmetry_operation_pg import _data_symmetry_operation_pg
from multipie.data.data_symmetry_operation_sg import _data_symmetry_operation_sg


# ==================================================
class SymmetryOperationG(
    dict
):  # dict of (symmetry operation tag, symmetry operation), {TagSymmetryOperation: SymmetryOperation}.
    """
    point/space-group symmetry operations.

    Attributes:
        tag (TagGroup): point-group tag.
        gen ([TagSymmetryOperation]): generators.
        cc ([[TagSymmetryOperation]]): symmetry operations categorized by conjugacy class.
        cc1 ([TagSymmetryOperation]): representative symmetry operations in conjugacy class.
        full ([TagSymmetryOperation]): all symmetry operations.
        cc_num ([int]): the number of symmetry operations in conjugacy class.
        crystal (str): crystal type, (triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic).
        lattice (str): lattice type, (A/B/C/P/I/F/R/null).
        translation (NSArray): partial translation vectors (reduced coordinate).
    """

    #
    #     __idx_cc1 ([int]): index of top of conjugacy class.
    #     __product ({(TagSymmetryOperation,TagSymmetryOperation): TagSymmetryOperation}): product table (point-group only).
    #     __inv ({TagSymmetryOperation:TagSymmetryOperation}): inverse operation.
    #
    # ==================================================
    def __init__(self, g_tag):
        """
        initialize the class.

        Args:
            g_tag (TagGroup or str): tag of point/space group.
        """
        g_tag = TagGroup(str(g_tag))
        self.tag = g_tag
        """point-group tag."""

        if self.tag.is_point_group():
            gen, cc = _data_symmetry_operation_pg[str(g_tag)]
        else:
            gen, cc = _data_symmetry_operation_sg[str(g_tag)]

        self.gen = TagList.from_str(TagSymmetryOperation, gen)
        """generators."""

        self.cc = [TagList.from_str(TagSymmetryOperation, i) for i in cc]
        """symmetry operations categorized by conjugacy class."""

        self.cc1 = TagList.from_str(TagSymmetryOperation, [i[0] for i in cc])
        """representative symmetry operations in conjugacy class."""

        self.cc_num = [len(i) for i in self.cc]
        """the number of symmetry operations in conjugacy class."""

        self.full = TagList.from_str(TagSymmetryOperation, sum(cc, []))
        """all symmetry operations."""

        self.__idx_cc1 = [self.full.index(i) for i in self.cc1]

        self.crystal = self.tag.crystal
        """crystal type, (triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic)."""

        self.lattice = self.tag.lattice
        """lattice type, (A/B/C/P/I/F/R/null)."""

        if self.lattice != "":
            self.translation = latticeT[self.lattice]
            """partial translation vectors (reduced coordinate)."""
        else:
            self.translation = None

        for s in self.full:
            self[s] = SymmetryOperation(str(s), self.crystal)

        # product table.
        self.__product = {}
        self.__inv = {}
        if self.tag.is_point_group():
            so = self.mat()
            so_inv = so.inverse()
            t = self.key_list()
            for i in range(len(so)):
                soi = so[i]
                self.__inv[t[i]] = t[so.index(so_inv[i])]
                for j in range(len(so)):
                    soj = so[j]
                    self.__product[(t[i], t[j])] = t[so.index(soi @ soj)]

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
    def plus_set(self):
        """
        additional partial translation vectors in space group.

        Returns:
            NSArray: translation vectors.
        """
        assert not self.tag.is_point_group(), "plus_set is only for space group."
        return self.translation

    # ==================================================
    def mat(self, axial=False, cc_only=False, primitive=False, cartesian=False):
        """
        matrix for symmetry operatoins.

        Args:
            axial (bool, optional): matrix for axial vector ?
            cc_only (bool, optional): conjugacy-class operations only ?
            primitive (bool, optional): for primitive cell ? (space group only)
            cartesian (bool, optional): in cartesian coordinate ?

        Returns:
            list: list of operations (3x3: point group) or (4x4: space group), [NSArray].
        """
        if cc_only:
            idx = TagList([self.full[i] for i in self.__idx_cc1])
        else:
            idx = self.full
        if axial:
            ops = NSArray([self[i].m_axial.tolist() for i in idx], "matrix")
        else:
            ops = NSArray([self[i].m_polar.tolist() for i in idx], "matrix")
        if primitive:
            ops = to_primitive(self.lattice, ops)
        if cartesian:
            ops = to_cartesian(self.crystal, ops, self.tag.is_point_group())

        return ops

    # ==================================================
    def __getitem__(self, item):
        if type(item) == str:
            return self.get(TagSymmetryOperation(item))
        else:
            return self.get(item)

    # ==================================================
    def key_list(self):
        return TagList(self.keys())

    # ==================================================
    def inverse(self, i):
        """
        inverse of symmetry operation (point group only).

        Args:
            i (TagSymmetryOperation): symmetry operation.

        Returns:
            TagSymmetryOperation: inverse of symmetry operation.
        """
        assert self.tag.is_point_group(), "product table is available only for point group."

        if type(i) == str:
            i = TagSymmetryOperation(i)
        return self.__inv[i]

    # ==================================================
    def product(self, i, j):
        """
        product of symmetry operations (point group only).

        Args:
            i (TagSymmetryOperation): symmetry operation (left).
            j (TagSymmetryOperation): symmetry operation (right).

        Returns:
            TagSymmetryOperation: product of symmetry operations.
        """
        assert self.tag.is_point_group(), "product table is available only for point group."

        if type(i) == str:
            i = TagSymmetryOperation(i)
        if type(j) == str:
            j = TagSymmetryOperation(j)
        return self.__product[(i, j)]

    # ==================================================
    def _transform_matrix(self, phi_p, phi, vs):
        """
        transform matrix for polynomial basis (phi_i' = T_ij phi_j).

        Args:
            phi_p (NSArray): transformed polynomials.
            phi (NSArray): polynomials to be transformed.
            vs (NSArray): list of variables in polynomials.

        Returns:
            NSArray: transform matrix (phi' = M.phi).
        """
        cv = NSArray(list(sp.symbols(f"coeffvar0:{phi.size}")), "vector")
        eqs = [phi_p[i].tolist() - phi @ cv for i in range(phi_p.size)]
        eqs = [i.expand() for i in eqs]
        tm = sp.conjugate(sp.Matrix([list(sp.linsolve(sp.Poly(eq, vs.tolist()).coeffs(), cv.tolist()))[0] for eq in eqs]))

        return NSArray(tm.tolist(), "matrix")

    # ==================================================
    def _transform_variable(self, phi, v_old, v_new):
        """
        transform polynomials by substituting transformed variables.

        Args:
            phi (NSArray): polynomials to be transformed.
            v_old (NSArray): list of variables in polynomials.
            v_new (NSArray): list of transformed variables.

        Returns:
            NSArray: transformed polynomial or list of polynomials.
        """
        return phi.subs({v_old[i]: v_new[i] for i in range(v_old.size)})

    # ==================================================
    def _transform_matrix_orbital_basis(self, v, vp, phi, pdet=None):
        """
        transform matrix for polynomial basis, (phi_i' = T_ij phi_j).

        Args:
            v (NSArray): polar vector variable (cartesian coordinate).
            vp (NSArray): list of transformed polar vectors, r' (cartesian coordinate).
            phi (NSArray): list of basis (cartesian coordinate).
            pdet (NSArray, optional): list of determinant.

        Returns:
            NSArray: list of transform matrices.
        """
        tm = []
        if pdet is None:
            for i in range(len(vp)):
                phi_p = self._transform_variable(phi, v, vp[i])
                tm.append(self._transform_matrix(phi_p, phi, v).tolist())
        else:
            for i in range(len(vp)):
                phi_p = self._transform_variable(phi, v, vp[i])
                phi_p = [phi_p[j] * pdet[j] for j in range(len(phi_p))]
                tm.append(self._transform_matrix(phi_p, phi, v).tolist())

        return NSArray(tm, "matrix")

    # ==================================================
    def _equivalent_vector(
        self, v, cc_only=False, axial=False, remove_duplicate=False, cartesian=False, primitive=False, plus_set=False, shift=False
    ):
        """
        create symmetry-operation equivalent vectors.

        Args:
            v (str or NSArray): vector/site position, "[x,y,z]".
            cc_only (bool, optional): conjugacy-class operations only ?
            axial (bool, optional): axial vector ?
            remove_duplicate (bool, optional): remove the same sites ?
            cartesian (bool, optional): cartesian coordinate ? (point group only).
            primitive (bool, optional): primitive cell ? (space group only).
            plus_set (bool, optional): add partial translations ? (space group only).
            shift (bool, optional): move into unit cell ? (space group only).

        Returns:
            NSArray: equivalent vectors.

        Notes:
            - if remove_duplicate is True, return values are sorted.
        """
        if not isinstance(v, (str, NSArray)):
            raise KeyError(f"{type(v)} is not accepted for vector/site.")

        if type(v) == str:
            v = NSArray(v)

        if self.tag.is_point_group():
            ops = self.mat(axial=axial, cc_only=cc_only, cartesian=cartesian)

            ev = ops.apply(v)
            if remove_duplicate:
                ev = ev.remove_duplicate()
        else:
            ops = self.mat(axial=axial, cc_only=cc_only)
            ops = to_primitive(self.lattice, ops)
            v = to_primitive(self.lattice, v)
            ev = ops.apply(v)

            if shift:
                ev = ev.shift()

            if remove_duplicate:
                ev = ev.remove_duplicate().sort()

            if not primitive:
                ev = to_conventional(self.lattice, ev, plus_set=plus_set)
                if type(ev) == list:
                    ev = NSArray.concat(ev)
                if shift:
                    ev = ev.shift()
                if remove_duplicate:
                    ev = ev.remove_duplicate().sort()

        return ev

    # ==================================================
    def _equivalent_bond(
        self,
        bond,
        cc_only=False,
        remove_duplicate=False,
        cartesian=False,
        primitive=False,
        plus_set=False,
        shift=False,
        nondirectional=False,
    ):
        """
        create SO-equivalent bonds.

        Args:
            bond (str or NSArray): bond, "vector@center", "tail;head", "start:vector".
            cc_only (bool, optional): conjugacy-class operations only ?
            remove_duplicate (bool, optional): remove the same bonds ?
            cartesian (bool, optional): cartesian coordinate ? (point group only).
            primitive (bool, optional): primitive cell ? (space group only).
            plus_set (bool, optional): add partial translations ? (space group only).
            shift (bool, optional): move into unit cell ? (space group only).
            nondirectional (bool, optional): ignore directional property ?

        Returns:
            NSArray: equivalent bonds.

        Notes:
            - if bond is string, vector-center "[vx,vy,vz]@[cx,cy,cz]", tail-head "[tx,ty,tz];[hx,hy,hz], or start-vector "[sx,sy,sz]:[vx,vy,vz]" style is accepted.
            - if remove_duplicate is True, return values are sorted.
        """
        if not isinstance(bond, (str, NSArray)):
            raise KeyError(f"{type(bond)} is not accepted for bond.")

        if type(bond) == str:
            bond = NSArray(bond)

        t, h = bond.convert_bond("bond_th")
        if self.tag.is_point_group():
            t = self._equivalent_vector(t, cc_only=cc_only, cartesian=cartesian)
            h = self._equivalent_vector(h, cc_only=cc_only, cartesian=cartesian)
            b = NSArray.create_bond_from_pair(t, h, "bond_th")
            if nondirectional:
                b = b.regular_direction()
            if remove_duplicate:
                b = b.remove_duplicate().sort()
        else:
            t = self._equivalent_vector(t, cc_only=cc_only)
            h = self._equivalent_vector(h, cc_only=cc_only)
            t = to_primitive(self.lattice, t)
            h = to_primitive(self.lattice, h)
            b = NSArray.create_bond_from_pair(t, h, "bond_th")
            if shift:
                b = b.shift()
            if nondirectional:
                b = b.regular_direction()
            if remove_duplicate:
                b = b.remove_duplicate().sort()

            if not primitive:
                t, h = b.convert_bond("bond_th")
                t = to_conventional(self.lattice, t, plus_set=plus_set)
                h = to_conventional(self.lattice, h, plus_set=plus_set)
                b = NSArray.create_bond_from_pair(t, h, "bond_th")
                if shift:
                    b = b.shift()
                if nondirectional:
                    b = b.regular_direction()
                if remove_duplicate:
                    b = b.remove_duplicate().sort()

        return b
