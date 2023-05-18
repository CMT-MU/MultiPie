"""
TagResponseTensor manages tag of response tensor.
"""
from dataclasses import dataclass
from multipie.tag.tag import Tag
from multipie.const import __def_dict__


# ==================================================
@dataclass(frozen=True, order=True)
class TagResponseTensor(Tag):
    """
    tag of response tensor.
    """

    head: str = "Q"  #
    rank: int = 0  #
    comp: str = "s"  # rank(0):(s), (1):(s), (2):(s/a), (3):(s/a), (4):(sss/ssa/aas/aaa/sa/as).

    # ==================================================
    @classmethod
    def _from_str(cls, tag):
        tensor_head, tag = tag.split("^")
        rank, tag = tag[2:].split(",")
        rank = int(rank)
        head = tag[0]
        comp = ""
        for (r, c), th in __def_dict__["response_head"].items():
            if r == rank and th == tensor_head:
                comp = c
                break
        d = {
            "head": head,
            "rank": rank,
            "comp": comp,
        }
        return d

    # ==================================================
    def __str__(self):
        return TagResponseTensor._tensor_str(self.head, self.rank, self.comp)

    # ==================================================
    def latex(self):
        ss = str(self).replace("Sb", "\\bar{S}").replace("Ab", "\\bar{A}").replace("Mb", "\\bar{M}")
        return ss

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
    def tensor_head(self):
        return __def_dict__["response_head"][(self.rank, self.comp)]

    # ==================================================
    @classmethod
    def _tensor_str(cls, head, rank, comp):
        ss = __def_dict__["response_head"][(rank, comp)] + "^{(" + str(rank) + "," + head + ")}"
        return ss

    # ==================================================
    @classmethod
    def create(cls, head, rank, comp):
        """
        create response tensor tag.

        Args:
            head (str): head, (Q/G/T/M).
            rank (int): rank, (0/1/2/3/4).
            comp (str): (0):(s), (1):(s), (2):(s/a), (3):(s/a), (4):(sss/ssa/aas/aaa/sa/as).

        Returns:
            TagResponseTensor: tag.
        """
        return TagResponseTensor(TagResponseTensor._tensor_str(head, rank, comp))
