"""
ResponseTensorPG manages physical response tensors upto rank 4.
"""
import sympy as sp
from gcoreutils.nsarray import NSArray
from multipie.tag.tag_group import TagGroup
from multipie.tag.tag_list import TagList
from multipie.response_tensor.util.response_tensor_util import simplify_tensor, create_tensor
from multipie.tag.tag_response_tensor import TagResponseTensor
from multipie.const import __def_dict__


# ==================================================
class ResponseTensorPG(dict):  # dict of (response tag, matrix), { TagResponseTensor: NSArray }.
    """
    a set of response tensors upto rank 4.

    Attributes:
        tag (TagGroup): point-group tag.
        definition ({TagResponseTensor:{sympy:sympy}}): definition in terms of multipoles for each component.
    """

    # ==================================================
    def __init__(self, pg_tag):
        """
        initialize the class.

        Args:
            pg_tag (TagGroup or str): point group tag.
        """
        if pg_tag is None:
            self.tag = None
            """point-group tag."""
            self.definition = {}
        else:
            self.tag = TagGroup(str(pg_tag))
            self._initialize()

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
    def _initialize(self):
        """
        create response tensors.
        """
        # create response tensor data, { tag: [[(tensor symbol, component)]] }.
        d = {}
        for head in __def_dict__["head"]:
            for rank, comp in __def_dict__["response_head"].keys():
                tag = TagResponseTensor.create(head, rank, comp)
                d[tag] = create_tensor(self.tag, tag)

        # create active tensor, { tag: Matrix }.
        mat = {}
        for tag, M in d.items():
            M = simplify_tensor(M)
            zero = sp.zeros(len(M), len(M[0]))
            x = zero.copy()
            for i, lst in enumerate(M):
                for j, (c, m) in enumerate(lst):
                    if m == sp.S(0):
                        x[i, j] = sp.S(0)
                    else:
                        x[i, j] = c
            if not x.equals(zero):
                mat[tag] = x

        # create definition of each element, { tag: {tensor_symbol: definition} }.
        de = {}
        for tag, M1 in mat.items():
            de[tag] = []
            M2 = sp.Matrix([[m for _, m in lst] for lst in d[tag]])
            Cs = []

            for i in range(M1.rows):
                for j in range(M1.cols):
                    m = M1[i, j]
                    if m != sp.S(0):
                        Cs += [i for i in m.atoms(sp.Symbol)]
            Cs = list(set(Cs))
            for m1, m2 in zip(M1, M2):
                c_lst = [i for i in m1.atoms(sp.Symbol)]
                if all([c in Cs for c in c_lst]) and m1 != sp.S(0):
                    de[tag].append((m1, m2))
                    for c in c_lst:
                        Cs.remove(c)

        for k, v in mat.items():
            self[k] = NSArray(v.tolist(), "matrix")
        self.definition = {k: dict(v) for k, v in de.items()}
        """definition in terms of multipoles for each component."""

    # ==================================================
    def select(self, **kwargs):
        """
        select response tensors with given keywords.

        Args:
            kwargs (dict): select conditions for response tensors, (head/rank/comp).

        Returns:
            ResponseTensorPG: selected response tensors.
        """
        c_rt = ResponseTensorPG(None)
        for t in self.key_list().select(**kwargs):
            c_rt[t] = self[t]
            c_rt.definition[t] = self.definition[t]
        c_rt.tag = self.tag

        return c_rt

    # ==================================================
    def _dump(self):
        """
        dump response tensors.
        """
        for tag in self.keys():
            print(f"=== rank {tag.rank} {tag.i_type} {tag.t_type} tensor '{tag.comp}' component ===")
            print(tag.symbol(), "=", self[tag])

            def_mul = self.definition[tag]
            for t, d in def_mul.items():
                print("   ", t, ":=", d)

    # ==================================================
    def _info(self):
        """
        get information of the response tensors.
        """
        s = "* -------------------------------------------------------------------------------------------------------- *" + "\n"
        s += " The symmetry properties of the response tensors 'C' of the 0-th to 4-th orders are summarized" + "\n"
        s += " in terms of the electric, electric toroidal, magnetic, and magnetic toroidal multipoles, Q, G, M, and T." + "\n"
        s += " The basis of the matrices is represented by using the Voigt notation (1=xx, 2=yy, 3=zz, 4=yz, 5=zx, 6=xy)." + "\n"
        s += "* -------------------------------------------------------------------------------------------------------- *" + "\n"
        s += "- rank 0 tensor" + "\n"
        s += "  C = Q_{0}" + "\n"
        s += "    C = (Q_{0})" + "\n"
        s += "- rank 1 tensor" + "\n"
        s += "  C = Q_{1}" + "\n"
        s += "    C = (Q_{x}, Q_{y}, Q_{z})" + "\n"
        s += "- rank 2 tensor" + "\n"
        s += "  C12 = S12 + A12" + "\n"
        s += "    S12 = S21" + "\n"
        s += "    A12 = -A21" + "\n"
        s += "- rank 3 tensor" + "\n"
        s += "  C123 = S12;3 + A12;3" + "\n"
        s += "    S12;3 = S21;3" + "\n"
        s += "    A12;3 = -A21;3" + "\n"
        s += "- rank 4 tensor" + "\n"
        s += "  C1234 = S12;34 + Sb12;34 + A12;34 + Ab12;34 + M12;34 + Mb12;34" + "\n"
        s += "    S12;34 = S21;34 = S12;43 = S34;12" + "\n"
        s += "    Sb12;34 = Sb21;34 = Sb12;43 = -Sb34;12" + "\n"
        s += "    A12;34 = -A21;34 = -A12;43 = A34;12" + "\n"
        s += "    Ab12;34 = -Ab21;34 = -Ab12;43 = -Ab34;12" + "\n"
        s += "    M12;34 = M21;34 = -M12;43" + "\n"
        s += "    Mb12;34 = -Mb21;34 = Mb12;43" + "\n"

        print(s)

    # ==================================================
    def __getitem__(self, tag):
        if type(tag) == str:
            return self.get(TagResponseTensor(tag))
        else:
            return self.get(tag)

    # ==================================================
    def key_list(self):
        return TagList(self.keys())
