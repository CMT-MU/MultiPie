"""
TagMultipole manages tag of multipole.
"""
from dataclasses import dataclass
from multipie.tag.tag import Tag
from multipie.tag.tag_list import TagList
from multipie.tag.tag_irrep import TagIrrep
from multipie.data.data_harmonics import (
    _data_harmonics_polar,
    _data_harmonics_axial,
)
from multipie.const import __def_dict__, is_magnetic, is_axial


# ==================================================
@dataclass(frozen=True, order=True)
class TagMultipole(Tag):
    """
    tag of multipole.
    """

    head: str = "Q"  # type of multipole, Q/G/T/M.
    m_type: str = (
        "h"  # type of multipole, h(harmonics)/a(atomic)/c(cluster)/s(site-cluster)/b(bond-cluster)/u(uniform)/k(structure).
    )
    rank: int = 0  # rank or l.
    irrep: str = ""  # irrep. of point group.
    mul: int = 0  # multiplicity (>0 if multiplicity is required).
    comp: int = -1  # component or m (-l <= m <= l).
    s: int = -1  # spinless(0)/spinful(1).
    k: int = -1  # 0/(-1,0,1).

    dim: int = 1

    # ==================================================
    @classmethod
    def _from_str(cls, tag):
        head = tag[0]
        if tag[1] in ("h", "a", "c", "s", "b", "u", "k"):
            m_type = tag[1]
            tag = tag[3:]
        else:
            m_type = ""
            tag = tag[2:]
        tag = tag.replace(")", "")
        if tag.count("|") > 0:
            a, b = tuple(tag.split("|"))
            s, k = tuple(b.split(","))
            s = int(s)
            k = int(k)
        elif m_type == "a":
            a = tag
            s, k = 0, 0
        else:  # harmonics
            a = tag
            s, k = -1, -1

        rank, irrep, mul, comp = tuple(a.split(","))
        rank = int(rank)
        mul = int(mul) if mul else 0
        comp = int(comp) if comp else -1

        if irrep == "":
            dim = 2 * rank + 1
        elif irrep.count("a") > 0 or irrep.count("b") > 0:
            dim = 1
        else:
            dim = __def_dict__["irrep_dim"][irrep[0]]

        d = {
            "head": head,
            "m_type": m_type,
            "rank": rank,
            "irrep": irrep,
            "mul": mul,
            "comp": comp,
            "s": s,
            "k": k,
            "dim": dim,
        }

        return d

    # ==================================================
    def __str__(self):
        head = self.head
        m_type = self.m_type
        rank = self.rank
        irrep = self.irrep
        mul = self.mul
        comp = self.comp
        s = self.s
        k = self.k

        dim = self.dim

        b1, b2 = "(", ")"

        ss = head + m_type + b1
        ss += str(rank) + "," + irrep + ","
        if mul > 0:
            ss += str(mul)
        ss += ","
        if irrep == "" or (irrep != "" and dim > 1):
            ss += str(comp)
        if s == 1:
            ss += "|" + str(s) + "," + str(k)
        ss += b2

        return ss

    # ==================================================
    def latex(self):
        h = r"\mathbb{" + self.head + "}"
        if self.m_type == "":
            h = r"\hat{" + h + "}"

        sup = "(" + self.m_type
        if self.irrep:
            sup += "," if self.m_type != "" else ""
            sup += TagIrrep(self.irrep).latex()
        if self.mul > 0:
            sup += "," + str(self.mul)
        sup += ")"

        sub = str(self.rank)
        if self.dim > 1:
            sub += "," + str(self.m)

        end = ""
        if self.s == 1:
            end += "(" + str(self.s) + "," + str(self.k) + ")"

        s = h
        if sub:
            s += "_{" + sub + "}"
        if sup:
            s += "^{" + sup + "}"
        if end:
            s += end

        return s

    # ==================================================
    def tag_irrep(self):
        """
        tag of irrep.

        Returns:
            TagIrrep: tag of irrep. or None (for spherical case).
        """
        return TagIrrep(self.irrep) if self.irrep else None

    # ==================================================
    @property
    def i_type(self):
        """
        parity type.

        Returns:
            str: parity type, (polar/axial).
        """
        return __def_dict__["head_i"][self.head]

    # ==================================================
    @property
    def t_type(self):
        """
        time-reversal type.

        Returns:
            str: time-reversal type, (electric/magnetic).
        """
        return __def_dict__["head_t"][self.head]

    # ==================================================
    @property
    def l(self):
        """
        rank of multipole.

        Returns:
            int: rank.
        """
        return self.rank

    # ==================================================
    @property
    def m(self):
        """
        component of multipole (:math:`-l\le m \le l` for spherical version).

        Returns:
            int: component.
        """
        return self.comp

    # ==================================================
    def is_complex(self):
        return "a" in self.irrep or "b" in self.irrep or (self.irrep == "" and self.m != 0)

    # ==================================================
    def is_axial(self):
        return is_axial(self.i_type)

    # ==================================================
    def is_magnetic(self):
        return is_magnetic(self.t_type)

    # ==================================================
    def to_harmonics(self):
        """
        convert tag to corresponding harmonics.

        Returns:
            TagMultipole: converted tag.
        """
        head = __def_dict__["head_harm"][self.head]
        return self.replace(head=head, m_type="h", s=-1, k=-1)

    # ==================================================
    def to_polar(self):
        """
        convert tag to corresponding polar one.

        Returns:
            TagMultipole: converted tag.
        """
        head = __def_dict__["to_polar"][self.head]
        return self.replace(head=head)

    # ==================================================
    def reverse_t_type(self):
        """
        time-reversal tag.

        Returns:
            TagMultipole: time-reversal tag.
        """
        head = __def_dict__["head_tr"][self.head]
        return self.replace(head=head)

    # ==================================================
    def reverse_i_type(self):
        """
        space-reversal tag.

        Returns:
            TagMultipole: space-reversal tag.
        """
        head = __def_dict__["head_ir"][self.head]
        return self.replace(head=head)

    # ==================================================
    @classmethod
    def create(cls, axial=False):
        """
        harmonics tags.

        Args:
            axial (bool, optional): axial ?

        Returns:
            TagList: tags.
        """
        if axial:
            d = _data_harmonics_axial.values()
        else:
            d = _data_harmonics_polar.values()

        tags = sum([list(pg.keys()) for pg in d], [])

        tags = TagList.from_str(TagMultipole, tags)

        return tags

    # ==================================================
    @classmethod
    def create_spherical(cls, head="Q", rank=0, mul=0, comp=1, s=-1, k=-1):
        """
        create tag for spherical (atomic) multipole.

        Args:
            head (str, optional): _description_. Defaults to "Q".
            rank (int, optional): _description_. Defaults to 0.
            mul (int, optional): _description_. Defaults to 0.
            comp (int, optional): _description_. Defaults to 1.
            s (int, optional): _description_. Defaults to -1.
            k (int, optional): _description_. Defaults to -1.

        Returns:
            _type_: _description_
        """
        ret = TagMultipole("")
        ret = ret.replace(
            head=head,
            rank=rank,
            mul=mul,
            comp=comp,
            s=s,
            k=k,
            m_type="a",
            dim=2 * rank + 1,
        )
        return ret
