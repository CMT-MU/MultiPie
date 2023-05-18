"""
CharacterPG manages point-group character table.
"""
import sympy as sp
from gcoreutils.convert_util import text_to_sympy
from multipie.tag.tag_group import TagGroup
from multipie.tag.tag_irrep import TagIrrep
from multipie.tag.tag_symmetry_operation import TagSymmetryOperation
from multipie.tag.tag_list import TagList
from multipie.data.data_character_table import _data_character
from multipie.data.data_product_decomp import _data_product_decomp
from multipie.data.data_compatibility_relation import _data_compatibility_relation


# ==================================================
class CharacterPG(dict):  # dict of (irrep, characters), {TagIrrep: [sympy]}.
    """
    point-group character table.

    Attributes:
        tag (TagGroup): point-group tag.
    """

    #
    #     __irrep (list): irrep., [TagIrrep].
    #     __symmetry_operation (list): symmetry operations of conjugacy class, [TagSymmetryOperation].
    #     __so_number (list): number of SOs in each conjugacy class, [int].
    #     __character_symbol (dict): character table in symbolic form, {TagIrrep:[sympy]}.
    #     __parity_conv (dict): parity conversion, {TagIrrep: TagIrrep}.
    #     __product_decomp_s (dict): symmetric product decomposition, {(TagIrrep,TagIrrep):[(int,TagIrrep)]}.
    #     __product_decomp_a (dict): anti-symmetric product decomposition, {TagIrrep:[(int,TagIrrep)]}.
    #
    # ==================================================
    def __init__(self, pg_tag):
        """
        initialize the class.

        Args:
            pg_tag (TagGroup or str): point-group tag.
        """
        pg_tag = TagGroup(str(pg_tag))

        # set tag.
        self.tag = pg_tag
        """point-group tag."""

        so, so_n, ct_dict = _data_character[str(pg_tag)]
        s, a = _data_product_decomp[str(pg_tag)]

        # set irrep.
        self.__irrep = TagList.from_str(TagIrrep, ct_dict.keys())

        # symmetry operation
        self.__symmetry_operation = TagList.from_str(TagSymmetryOperation, so)
        self.__so_number = so_n

        # set character table.
        w = sp.symbols(r"\omega \omega^*")
        wpv = text_to_sympy("-1/2+sqrt(3)*I/2")
        wmv = text_to_sympy("-1/2-sqrt(3)*I/2")
        self.__character_symbol = {}
        for irrep, ct in ct_dict.items():
            self[TagIrrep(irrep)] = text_to_sympy(ct, local={"wp": wpv, "wm": wmv})
            self.__character_symbol[TagIrrep(irrep)] = text_to_sympy(ct, local={"wp": w[0], "wm": w[1]})

        # set parity conversion.
        cr = TagList.from_str(TagIrrep, _data_compatibility_relation[str(pg_tag)])
        n = len(cr) // 2
        cr_conv = cr[n:] + cr[:n]
        self.__parity_conv = {i: j for i, j in zip(cr, cr_conv)}

        # set product decomposition.
        self.__product_decomp_s = {}
        for ir1, p in zip(self.__irrep, s):
            for ir2, v in zip(self.__irrep, p):
                v = [(n, TagIrrep(ir)) for n, ir in v]
                self.__product_decomp_s[(ir1, ir2)] = v
        self.__product_decomp_a = {}
        for ir, v in zip(self.__irrep, a):
            v = [(n, TagIrrep(ir)) for n, ir in v]
            self.__product_decomp_a[ir] = v

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
    def __getitem__(self, tag):
        if type(tag) == str:
            return self.get(TagIrrep(tag))
        else:
            return self.get(tag)

    # ==================================================
    def key_list(self):
        return TagList(self.keys())

    # ==================================================
    @property
    def irrep_list(self):
        """
        list of irreps.

        Returns:
            [TagIrrep]: list of irreps.
        Notes:
            - first element is identity irrep.
        """
        return self.__irrep

    # ==================================================
    def character(self, irrep=None, all_so=False, explicit=False):
        """
        character table.

        Args:
            irrep (str or TagIrrep, optional): irrep. tag ("" = identity irrep.).
            all_so (bool, optional): return characters for all symmetry operations ?
            explicit (bool, optional): use explicit expression of character ?

        Returns:
            - [sympy]: character table (symbol) for given irrep.
            or
            - {TagIrrep:[sympy]}: all irreps.
        """
        if irrep == "":
            irrep = self.irrep_list[0]
        ctbl = self if explicit else self.__character_symbol

        if irrep is None:
            if all_so:
                return {i: self._to_all(c) for i, c in ctbl.items()}
            else:
                return ctbl
        else:
            if type(irrep) == str:
                lst = list(map(str, self.irrep_list))
                assert irrep in lst, f"{irrep} is not find in {lst}"
                irrep = TagIrrep(irrep)
            if all_so:
                return self._to_all(ctbl[irrep])
            else:
                return ctbl[irrep]

    # ==================================================
    def _to_all(self, ch_cc):
        """
        convert character table for all symmetry operations.

        Args:
            ch_cc ([sympy]): characters for cc operations.

        Returns:
            [sympy]: characters for all symmetry operations.
        """
        ch = []
        for chi, n in zip(ch_cc, self.__so_number):
            ch.extend([chi] * n)

        return ch

    # ==================================================
    def symmetry_operation(self, ret_num=False):
        """
        representative symmetry operations in conjugacy class.

        Args:
            ret_num (bool, optional): return with the number of cc operations ?

        Returns:
            - [TagSymmetryOperation]: symmetry operation.
            or
            - [TagSymmetryOperation]: symmetry operations.
            - [int]: the number of symmetry operations in conjugacy class.
        """
        if ret_num:
            return self.__symmetry_operation, self.__so_number
        else:
            return self.__symmetry_operation

    # ==================================================
    def compatibility_relation(self, pg_tag):
        """
        compatibility_relation from the present point group to a given one.

        Args:
            pg_tag (TagGroup or str): point-group tag to obtain compatibility relation.

        Returns:
            {TagIrrep: [TagIrrep]}: compatibility relation from self to other.
        """
        if type(pg_tag) == str:
            other = TagGroup(pg_tag)
        else:
            other = pg_tag
        assert self.tag.subgroup == other.subgroup, "diffrent subgroup is given."
        p = self.irrep_list
        c = TagList.from_str(TagIrrep, _data_compatibility_relation[str(other)])
        return {i: j for i, j in zip(p, c)}

    # ==================================================
    def parity_conversion(self):
        """
        parity conversion.

        Returns:
            {TagIrrep: TagIrrep: correspondence from self to converted one.
        """
        return self.__parity_conv

    # ==================================================
    def symmetric_product_decomposition(self, irrep_pair, ret_ex=False):
        """
        symmetric product irreducible decomposition.

        Args:
            irrep_pair ((TagIrrep or str, TagIrrep or str)): a pair of irreps.
            ret_ex (bool, optional): return as an expression ?

        Returns:
            - [(int,TagIrrep)]: (the number of appearances, irrep.).
            or
            - sympy: decomposed expression.
        """
        ir1, ir2 = irrep_pair
        if type(ir1) == str:
            ir1 = TagIrrep(ir1)
        if type(ir2) == str:
            ir2 = TagIrrep(ir2)
        val = self.__product_decomp_s[(ir1, ir2)]

        if ret_ex:
            ex = 0
            for n, v in val:
                ex += n * v.symbol()
            val = ex

        return val

    # ==================================================
    def anti_symmetric_product_decomposition(self, irrep, ret_ex=False):
        """
        anti-symmetric product irreducible decomposition.

        Args:
            irrep (TagIrrep or str): irrep. tag.
            ret_ex (bool, optional): return as an expression ?

        Returns:
            - [(int,TagIrrep)]: (the number of appearances, irrep.).
            or
            - sympy: decomposed expression.
        """
        if type(irrep) == str:
            irrep = TagIrrep(irrep)
        val = self.__product_decomp_a[irrep]

        if ret_ex:
            ex = 0
            for n, v in val:
                ex += n * v.symbol()
            val = ex

        return val

    # ==================================================
    def irrep_decomposition(self, char):
        """
        irreducible decomposition.

        Args:
            char ([sympy]): characters to be decomposed for symmetry operations in conjugacy class, .

        Returns:
            [(int,TagIrrep)]: (the number of appearances, irrep.).
        """
        decomp = []
        g = sum(self.__so_number)
        for irrep, ir_char in self.items():
            s = 0
            for c, ic, d in zip(char, ir_char, self.__so_number):
                s += d * sp.conjugate(c) * ic
            n = sp.simplify(s) / g
            if n != 0:
                decomp.append((n, irrep))
        return decomp
