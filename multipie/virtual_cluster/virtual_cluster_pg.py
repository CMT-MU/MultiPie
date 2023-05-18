"""
VirtualClusterPG manages point-group virtual cluster (real).
"""
from gcoreutils.nsarray import NSArray
from multipie.tag.tag_group import TagGroup
from multipie.tag.tag_multipole import TagMultipole
from multipie.tag.tag_list import TagList
from multipie.data.data_virtual_cluster_real import (
    _data_virtual_cluster_basis_real,
    _data_virtual_cluster_site_real,
)


# ==================================================
class VirtualClusterPG(dict):  # dict of (multipole tag, vc basis), {TagMultipole: NSArray} (orthonormalized version).
    """
    point/space-group virtual cluster.

    Attributes:
        tag (TagGroup): group tag.
        site (NSArray): site positions (reduced coordinate).
    """

    # ==================================================
    def __init__(self, pg_tag):
        """
        initialize the class.

        Args:
            pg_tag (TagGroup or str): point-group tag.
        """
        pg_tag = TagGroup(pg_tag)
        self.tag = pg_tag
        """point-group tag."""

        assert self.tag.is_point_group(), "it is valid only for point group."

        vc = _data_virtual_cluster_basis_real
        vc_site = _data_virtual_cluster_site_real

        self._basis = {}
        self._site = {}
        for wp, dic in vc[str(pg_tag)].items():
            self._basis[wp] = {}
            for m, b in dic.items():
                self._basis[wp][TagMultipole(m)] = NSArray(b)
        for wp, vec in vc_site[str(pg_tag)].items():
            self._site[wp] = NSArray(vec)

        wp0 = list(self._basis.keys())[0]  # assume general point is the first.
        d = {}
        for m, b in self._basis[wp0].items():
            d[m] = b
        self.update(d)
        self.site = self._site[wp0]
        """site positions (cartesian coordinate)."""

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
            return self.get(TagMultipole(tag))
        else:
            return self.get(tag)

    # ==================================================
    def get(self, wyckoff, tag):
        if type(tag) == str:
            return self._basis[wyckoff][TagMultipole(tag)]
        else:
            return self._basis[wyckoff][tag]

    # ==================================================
    def key_list(self, wyckoff=None):
        if wyckoff is None:
            return TagList(self.keys())
        else:
            return TagList(self._basis[wyckoff].keys())
