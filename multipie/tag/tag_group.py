"""
TagGroup manages tag of group.
"""
from dataclasses import dataclass
from multipie.tag.tag import Tag
from multipie.tag.tag_list import TagList
from multipie.data.data_tag_space_group import _data_tag_space_group, _data_tag_space_group_crystal
from multipie.data.data_tag_point_group import (
    _data_tag_point_group,
    _data_tag_point_group_complex,
    _data_tag_point_group_real,
    _data_tag_point_group_crystal,
    _data_tag_point_group_cubic_series,
    _data_tag_point_group_hexagonal_series,
)
from multipie.const import __def_dict__


# ==================================================
@dataclass(frozen=True, order=True)
class TagGroup(Tag):
    """
    tag of group.
    """

    no: int = 1  # group No.
    SS: str = "C_1"  # Schoenflies symbol in LaTeX.
    IS: str = "C_1"  # international symbol(short) in LaTeX.
    crystal: str = "triclinic"  # crystal system.
    setting: str = ""  # setting.
    pg: str = ""  # point-group tag for space group, or "" for point group.
    subgroup: str = "cubic"  # subgroup (cubic/hexagonal).

    # ==================================================
    @classmethod
    def _from_str(cls, tag):
        if tag.count("^") > 0:
            no, SS, IS, setting, crystal, pg = _data_tag_space_group[tag]
        else:
            no, SS, IS, setting, crystal = _data_tag_point_group[tag]
            pg = ""

        subgroup = __def_dict__["subgroup"][crystal]

        d = {"no": no, "SS": SS, "IS": IS, "setting": setting, "crystal": crystal, "pg": pg, "subgroup": subgroup}

        return d

    # ==================================================
    def __str__(self):
        return self.SS.replace("_", "").replace("{", "").replace("}", "")

    # ==================================================
    def latex(self):
        return self.SS

    # ==================================================
    def info(self, latex=False):
        """
        info. of point group.

        Args:
            latex (bool, optional): in LaTeX ?

        Returns:
            - int: No. of group.
            - str: Schoenflies symbol.
            - str: International (short) symbol.
            - str: setting.
            or
            - str: info. string in LaTeX.
        """
        if latex:
            sp = r"\quad"
            r = f"No. {self.no}{sp}${self.SS}${sp}${self.IS}$"
            if self.setting:
                r += f"{sp}({self.setting} setting)"
            r += f"{sp}[ {self.crystal} ]"
            return r
        else:
            return self.no, self.SS, self.IS, self.setting

    # ==================================================
    def is_point_group(self):
        """
        is point group ?

        Returns:
            bool: is point group ?
        """
        return self.pg == ""

    # ==================================================
    def is_hexagonal_subgroup(self):
        """
        is hexagonal subgroup ?

        Returns:
            bool: is hexagonal subgroup ?
        """
        return self.subgroup == "hexagonal"

    # ==================================================
    @property
    def lattice(self):
        """
        Bravais lattice, (A/B/C/P/I/F/R) or null for point group.

        Returns:
            str: lattice letter.
        """
        if self.is_point_group():
            return ""
        else:
            return self.IS[0]

    # ==================================================
    @classmethod
    def create(cls, crystal="", space_group=False):
        """
        create group tags.

        Args:
            crystal (str, optional): crystal system, or "complex/real/cubic_series/hexagonal_series" for point group.
            space_group (bool, optional): space group or point_group ?

        Returns:
            TagList: tags.
        """
        if space_group:
            if crystal in __def_dict__["crystal"]:
                return TagList.from_str(TagGroup, _data_tag_space_group_crystal[crystal])
            else:
                return TagList.from_str(TagGroup, _data_tag_space_group.keys())
        else:
            if crystal in __def_dict__["crystal"]:
                return TagList.from_str(TagGroup, _data_tag_point_group_crystal[crystal])
            elif crystal == "complex":
                return TagList.from_str(TagGroup, _data_tag_point_group_complex)
            elif crystal == "real":
                return TagList.from_str(TagGroup, _data_tag_point_group_real)
            elif crystal == "cubic_series":
                return TagList.from_str(TagGroup, _data_tag_point_group_cubic_series)
            elif crystal == "hexagonal_series":
                return TagList.from_str(TagGroup, _data_tag_point_group_hexagonal_series)
            else:
                return TagList.from_str(TagGroup, _data_tag_point_group.keys())
