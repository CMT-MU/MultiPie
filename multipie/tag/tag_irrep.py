"""
TagIrrep manages tag of irreducible representation.
"""
from dataclasses import dataclass
from multipie.tag.tag import Tag
from multipie.tag.tag_list import TagList
from multipie.data.data_tag_irrep import _data_tag_irrep
from multipie.const import __def_dict__


# ==================================================
@dataclass(frozen=True, order=True)
class TagIrrep(Tag):
    """
    tag of irreducible representation.
    """

    tag: str = ""  # tag of irrep.
    dim: int = 1  # dimension.
    parity: str = ""  # null/g/u.
    cmplx: str = ""  # null/a/b.

    head: str = ""
    sub: str = ""
    sup: str = ""

    # ==================================================
    @classmethod
    def _from_str(cls, tag):
        head = tag[0]
        sub = tag[1:]
        sub = sub.replace("'", "")
        np = tag.count("'")
        sup = "'" * np

        if tag.count("a"):
            c = "a"
            d = 1
        elif tag.count("b"):
            c = "b"
            d = 1
        else:
            c = ""
            d = __def_dict__["irrep_dim"][tag[0]]
        p = ""
        if tag.count("g"):
            p = "g"
        elif tag.count("u"):
            p = "u"

        d = {"tag": tag, "dim": d, "parity": p, "cmplx": c, "head": head, "sub": sub, "sup": sup}

        return d

    # ==================================================
    def __str__(self):
        s = self.head + self.sub + self.sup
        return s

    # ==================================================
    def latex(self):
        sup = self.sup.replace("'", r"\prime")
        sub = self.sub
        if sub.count("a"):
            sup += "(a)"
            sub = sub.replace("a", "")
        elif sub.count("b"):
            sup += "(b)"
            sub = sub.replace("b", "")

        s = self.head
        if sub:
            s += "_{" + sub + "}"
        if sup:
            s += "^{" + sup + "}"
        return s

    # ==================================================
    @classmethod
    def create(cls):
        """
        irrep. tags for point groups.

        Returns:
            TagList: all irrep. tags.
        """
        return TagList.from_str(TagIrrep, _data_tag_irrep)
