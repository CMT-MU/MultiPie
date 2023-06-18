"""
WyckoffG manages point/space-group Wyckoff positions (conventional, reduced coordinate).
"""
from gcoreutils.nsarray import NSArray
from multipie.tag.tag_group import TagGroup
from multipie.tag.tag_wyckoff import TagWyckoff
from multipie.tag.tag_list import TagList
from multipie.data.data_wyckoff_pg import _data_wyckoff_pg
from multipie.data.data_wyckoff_sg import _data_wyckoff_sg


# ==================================================
class WyckoffG(dict):  # dict of (wyckoff tag, wyckoff), {TagWyckoff: NSArray}.
    """
    point/space-group Wyckoff positions.

    Attributes:
        tag (TagGroup): group tag.
        v (sp.Matrix): vector variable
    """

    #
    #     __default (dict): default position for point group, {TagWyckoff: NSArray}.
    #     __sym (dict): {wyckoff letter: site symmetry}
    #
    # ==================================================
    def __init__(self, g_tag):
        """
        initialize the class.

        Args:
            g_tag (TagGroup or str): tag of group.
        """
        g_tag = TagGroup(str(g_tag))
        self.tag = g_tag
        """group tag."""

        self.v = NSArray.vector3d()
        """(x,y,z) vector."""

        self.__sym = {}
        if self.tag.is_point_group():
            self.__default = {}
            for w, pos, d, sym in _data_wyckoff_pg[str(g_tag)]:
                v = dict(zip(["x", "y", "z"], self.v))
                self[TagWyckoff(w)] = NSArray(pos).subs(v)
                v = dict(zip(["x", "y", "z"], NSArray(d)))
                self.__default[TagWyckoff(w)] = NSArray(pos).subs(v)
                self.__sym[TagWyckoff(w)] = sym
        else:
            self.__default = None
            for w, pos, sym in _data_wyckoff_sg[str(g_tag)]:
                v = dict(zip(["x", "y", "z"], self.v))
                self[TagWyckoff(w)] = NSArray(pos).subs(v)
                self.__sym[TagWyckoff(w)] = sym

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
    def position(self, wp, v=None, default=False):
        """
        wyckoff position (reduced coordinate).

        Args:
            wp (TagWyckoff or str): Wyckoff position tag.
            v (str or NSArray, optional): vector variable.
            default (bool, optional): use default value ? (point group only).

        Returns:
            NSArray: wyckoff positions.
        """
        if not self.tag.is_point_group() and default is True:
            raise NotImplementedError("default value is unsupported for space group.")

        if type(wp) == str:
            wp = TagWyckoff(wp)

        if default:
            return self.__default[wp]

        if v is None:
            return self[wp]
        else:
            v = dict(zip(["x", "y", "z"], NSArray(v)))
            return self[wp].subs(v)

    # ==================================================
    def __getitem__(self, wp):
        if type(wp) == str:
            return self.get(TagWyckoff(wp))
        else:
            return self.get(wp)

    # ==================================================
    def key_list(self):
        return TagList(self.keys())

    # ==================================================
    def site_symmetry(self, wp):
        """
        site symmetry.

        Args:
            wp (TagWyckoff or str): Wyckoff position tag.

        Returns:
            str: site symmetry.
        """
        if type(wp) == str:
            wp = TagWyckoff(wp)

        return self.__sym[wp]
