"""
Global definition.

This module provides global definitions for MultiPie.
"""

from typing import TYPE_CHECKING
from collections import namedtuple
from multipie.core.multipie_info import __version__

# ==================================================
PGInfoType = namedtuple(
    "PGInfoType",
    [
        "tag",
        "international",
        "schoenflies",
        "crystal",
        "setting",
        "PG",
        "SG",
        "MPG",
        "MSG",
        "lattice",
        "hexagonal_g",
        "SO",
        "SO_gen",
    ],
)
"""Point Group info type."""

SGInfoType = namedtuple(
    "SGInfoType",
    [
        "tag",
        "international",
        "schoenflies",
        "crystal",
        "setting",
        "PG",
        "SG",
        "MPG",
        "MSG",
        "lattice",
        "hexagonal_g",
        "SO",
        "SO_gen",
    ],
)
"""Space Group info type."""

MPGInfoType = namedtuple(
    "MPGInfoType",
    [
        "tag",
        "international",
        "schoenflies",
        "crystal",
        "setting",
        "PG",
        "SG",
        "MPG",
        "MSG",
        "lattice",
        "hexagonal_g",
        "type",
        "SO",
    ],
)
"""Magnetic Point Group info type."""

MSGInfoType = namedtuple(
    "MSGInfoType",
    [
        "tag",
        "BNS",
        "OG",
        "schoenflies",
        "crystal",
        "setting",
        "PG",
        "SG",
        "MPG",
        "MSG",
        "lattice",
        "hexagonal_g",
        "type",
        "SO",
    ],
)
"""Magnetic Space Group info type."""

# ==================================================
SphericalMultipoleType = namedtuple("SphericalMultipoleType", ["X", "l", "s", "k", "x"])
"""Spherical multipole type (type, rank, s, k, internal-type)."""

PGMultipoleType = namedtuple("PGMultipoleType", ["X", "l", "Gamma", "n", "p", "s", "k", "x"])
"""Point-group multipole type (type, rank, irrep, multiplicity, add_multiplicity, s, k, internal-type)."""

RepSiteType = namedtuple("RepSiteType", ["no", "wyckoff", "symmetry", "position", "orbital"])
"""Representative site type."""

CellSiteType = namedtuple("CellSiteType", ["no", "position", "position_primitive", "mapping", "sublattice", "plus_set"])
"""Cell site type."""

BondInfoType = namedtuple("BondInfoType", ["tail", "head", "neighbor", "t_rank", "h_rank"])
"""Bond info type."""

RepBondType = namedtuple(
    "RepBondType",
    ["no", "tail", "head", "neighbor", "wyckoff", "directional", "vector", "center", "distance", "t_rank", "h_rank"],
)
"""Representative bond type."""

CellBondType = namedtuple(
    "CellBondType",
    [
        "no",
        "vector",
        "vector_primitive",
        "center",
        "center_primitive",
        "mapping",
        "sublattice",
        "plus_set",
        "t_idx",
        "h_idx",
        "R_primitive",
    ],
)
"""Cell bond type."""

BraketInfoType = namedtuple("BraketInfoType", ["bh_rank", "bh_idx", "kt_rank", "kt_idx"])
"""Braket info type."""

SAMBType = namedtuple("SAMBType", ["head", "tail", "wyckoff", "bk_info"])
"""Combined SAMB type."""

UniqueSAMBType = namedtuple("UniqueSAMBType", ["samb_type", "neighbor", "n"])
"""Unique combined SAMB type."""

SOMatrixType = namedtuple("SOMatrixType", ["X", "s"])
"""Symmetry operation matrix type."""

IrrepListType = namedtuple("IrrepListType", ["X", "Gamma"])
"""Irrep list type."""

InternalBasisType = namedtuple("InternalBasisType", ["s", "x"])
"""Internal basis type."""

RepMatrixType = namedtuple("RepMatrixType", ["Gamma"])
"""Representation matrix type."""


__all__ = ["Group", "MaterialModel", "create_samb", "create_samb_qtdraw", "create_samb_matrix"]

if TYPE_CHECKING:
    from multipie.core.group import Group
    from multipie.core.cmd import create_samb
    from multipie.core.cmd import create_samb_qtdraw
    from multipie.core.cmd import create_samb_matrix
    from multipie.core.material_model import MaterialModel


def __getattr__(name):
    if name == "Group":
        from multipie.core.group import Group

        return Group

    if name == "create_samb":
        from multipie.core.cmd import create_samb

        return create_samb

    if name == "create_samb_qtdraw":
        from multipie.core.cmd import create_samb_qtdraw

        return create_samb_qtdraw

    if name == "create_samb_matrix":
        from multipie.core.cmd import create_samb_matrix

        return create_samb_matrix

    if name == "MaterialModel":
        from multipie.core.material_model import MaterialModel

        return MaterialModel

    raise AttributeError
