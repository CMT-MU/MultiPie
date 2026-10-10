"""
Regression tests for the group database (Group).
"""

from collections import Counter
from collections.abc import MutableMapping

import numpy as np
import pytest

from multipie.core.group import Group, _load_group_data, _load_group_opt_data


# ==================================================
def test_number_of_groups():
    count = Counter(tag.split(":")[0] for tag in Group.global_info()["tag"].keys())
    assert count == {"PG": 47, "SG": 230, "MPG": 157, "MSG": 1651}


# ==================================================
@pytest.mark.parametrize(
    "tag, group_id, n_so, irreps",
    [
        ("Oh", "PG:32", 48, ["A1g", "A2g", "Eg", "T1g", "T2g", "A1u", "A2u", "Eu", "T1u", "T2u"]),
        (
            "D6h",
            "PG:27",
            24,
            ["A1g", "A2g", "B1g", "B2g", "E1g", "E2g", "A1u", "A2u", "B1u", "B2u", "E1u", "E2u"],
        ),
        ("C3v", "PG:19", 6, ["A1", "A2", "E"]),
        ("SG:221", "SG:221", 48, ["A1g", "A2g", "Eg", "T1g", "T2g", "A1u", "A2u", "Eu", "T1u", "T2u"]),
        (
            "SG:194",
            "SG:194",
            24,
            ["A1g", "A2g", "B1g", "B2g", "E1g", "E2g", "A1u", "A2u", "B1u", "B2u", "E1u", "E2u"],
        ),
        ("MPG:1.1.1", "MPG:1.1.1", 1, ["A"]),
    ],
)
def test_group_data(tag, group_id, n_so, irreps):
    g = Group(tag)
    assert g._id == group_id
    assert len(g.symmetry_operation["tag"]) == n_so
    assert list(g.character["table"].keys()) == irreps


# ==================================================
@pytest.mark.parametrize("no", [221, np.int64(221)])
def test_space_group_number(no):
    assert Group(no)._id == "SG:221"


# ==================================================
@pytest.mark.parametrize("tag", ["XYZ", True, 231])
def test_unknown_tag(tag):
    with pytest.raises(ValueError, match="unknown tag"):
        Group(tag)


# ==================================================
@pytest.mark.parametrize(
    "tag, site, wyckoff",
    [
        ("SG:221", "[0,0,0]", "1a"),
        ("SG:221", "[1/2,1/2,1/2]", "1b"),
    ],
)
def test_find_wyckoff(tag, site, wyckoff):
    assert Group(tag).find_wyckoff(site)[0] == wyckoff


# ==================================================
def test_group_data_cache():
    # the data of a group is loaded once and shared (read only) by the Group objects of the same group.
    Group.clear_cache()
    g1, g2 = Group("Oh"), Group("PG:32")
    so = g1.symmetry_operation["cartesian"]
    assert g2.symmetry_operation["cartesian"] is so
    with pytest.raises(ValueError, match="read-only"):
        so[0, 0, 0] = 2
    Group.clear_cache()
    so2 = Group("Oh").symmetry_operation["cartesian"]
    assert so2 is not so
    assert (so2 == so).all()


# ==================================================
def test_group_data_containers_not_shared():
    # dicts and lists are copied for each Group object, so a change of one does not affect the others.
    g1, g2 = Group("C2v"), Group("C2v")
    so = g2.symmetry_operation["cartesian"]
    g1.symmetry_operation["cartesian"] = None
    g1.symmetry_operation["tag"].append("x")
    wp = next(iter(g1.wyckoff["site"]))
    del g1.wyckoff["site"][wp]
    for g in [g2, Group("C2v")]:
        assert g.symmetry_operation["cartesian"] is so
        assert "x" not in g.symmetry_operation["tag"]
        assert wp in g.wyckoff["site"]


# ==================================================
def test_group_opt_data_cache():
    # the optional data is cached in the same way (kept for two groups), with read-only arrays.
    Group.clear_cache()
    g1, g2 = Group("C1", with_opt=True), Group("C1", with_opt=True)  # loaded once (C1: the fastest to load).
    assert _load_group_opt_data.cache_info().misses == 1
    key = next(iter(g1.opt))
    g1.opt[key] = None
    assert g2.opt[key] is not None
    with pytest.raises(ValueError, match="read-only"):
        first_array(g2.opt["representation_matrix"]).flat[0] = 2
    Group.clear_cache()
    assert _load_group_data.cache_info().currsize == _load_group_opt_data.cache_info().currsize == 0


def first_array(obj):
    """
    First NumPy array in nested dicts, lists and tuples.
    """
    if isinstance(obj, np.ndarray):
        return obj
    values = obj.values() if isinstance(obj, MutableMapping) else obj if isinstance(obj, (list, tuple)) else []
    for v in values:
        a = first_array(v)
        if a is not None:
            return a
    return None


# ==================================================
def test_group_data_lazy_and_after_clear():
    # the data of a related group, loaded when used (here the PG of an SG), is copied in the same way,
    # and a Group object can be used after clear_cache.
    g1, g2 = Group(221), Group(221)
    ct = g1.character["table"]
    ct.clear()
    assert g2.character["table"] and Group(221).character["table"]
    Group.clear_cache()
    assert g2.character["table"] and len(g2.symmetry_operation["tag"]) == 48


# ==================================================
def _cg_direct(group, tag1, tag2, tag):
    # Group.cg computed term by term, as in the definition (reference for the implementation with cached coefficients).
    import sympy as sp
    from itertools import product
    from sympy.physics.quantum.cg import CG

    t_even = {"Q": "Q", "T": "Q", "G": "G", "M": "G"}
    u1, u2, u = (group.harmonics[[t_even[t[0]], *t[1:4], -1, 0, 0, "q"]][1].T for t in (tag1, tag2, tag))
    l1, l2, l = tag1[1], tag2[1], tag[1]
    ml1, ml2, ml = (list(range(i, -i - 1, -1)) for i in (l1, l2, l))
    cg = np.full((len(u1), len(u2), len(u)), sp.S(0), dtype=object)
    for g1, g2, g in product(range(len(u1)), range(len(u2)), range(len(u))):
        s = sp.S(0)
        for i, j in product(range(len(ml1)), range(len(ml2))):
            m = ml1[i] + ml2[j]
            if m in ml:
                s += (
                    (-sp.I) ** (l1 + l2 - l)
                    * CG(l1, ml1[i], l2, ml2[j], l, m).doit()
                    * sp.conjugate(u1[g1, i] * u2[g2, j])
                    * u[g, ml.index(m)]
                )
        cg[g1, g2, g] = s
    return cg


@pytest.mark.parametrize("tag", ["D4h", "Oh", "D3d"])
def test_cg(tag):
    import sympy as sp

    g = Group(tag)
    keys = [k for k in g.harmonics.named_keys() if k.X == "Q" and k.l <= 2]
    n = 0
    for tag1 in keys:
        for tag2 in keys:
            for tag3 in keys:
                cg = g.cg(tag1, tag2, tag3)
                if cg is None:
                    continue
                ref = _cg_direct(g, tag1, tag2, tag3)
                assert cg.shape == ref.shape
                assert all(sp.simplify(a - b) == 0 for a, b in zip(cg.ravel(), ref.ravel())), (tag1, tag2, tag3)
                n += 1
    assert n > 0


# ==================================================
def test_find_xyz():
    from multipie.util.util import str_to_sympy
    from multipie.util.util_wyckoff import find_xyz

    pos = str_to_sympy("[x, 2*x + 1/2, -z]")
    sol = find_xyz(pos, [0.1, 0.7, 0.3])
    assert sol.keys() == {"x", "y", "z"} and np.isclose(sol["x"], 0.1) and np.isclose(sol["z"], -0.3)
    assert find_xyz(pos, [0.1, 0.6, 0.3]) is None  # y != 2x+1/2.
    pos = str_to_sympy("[X, -X, 0]")
    assert np.isclose(find_xyz(pos, [0.25, -0.25, 0], ["X", "Y", "Z"])["X"], 0.25)
    assert find_xyz(pos, [0.25, 0.25, 0], ["X", "Y", "Z"]) is None
