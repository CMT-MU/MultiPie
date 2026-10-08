"""
Regression tests for the group database (Group).
"""

from collections import Counter

import pytest

from multipie.core.group import Group


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
def test_space_group_number():
    assert Group(221)._id == "SG:221"


# ==================================================
def test_unknown_tag():
    with pytest.raises(Exception):
        Group("XYZ")


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
