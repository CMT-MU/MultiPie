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
"""
Point Group info type.

Fields:
  - tag (str): Schoenflies in text.
  - international (str): international short symbol in LaTeX.
  - schoenflies (str): Schoenflies symbol in LaTeX.
  - crystal (str): triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
  - setting (str): setting.
  - PG (str): = id.
  - SG (str): associated SG, 1st SG in the same PGs.
  - MPG (str): associated MPG, 1st MPG with type II in the same PGs.
  - MSG (str): associated MSG, PG -> MPG -> MSG.
  - lattice (str): "0".
  - hexagonal_g (bool): trigonal or hexagonal ?
  - SO (list): symmetry operations.
  - SO_gen (list): generator of SOs.
"""

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
"""
Space Group info type.

Fields:
  - tag (str): Schoenflies in text.
  - international (str): international short symbol in LaTeX.
  - schoenflies (str): Schoenflies symbol in LaTeX.
  - crystal (str): triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
  - setting (str): setting comment.
  - PG (str): associated PG, unique.
  - SG (str): = id.
  - MPG (str): associated MPG, SG -> PG -> MPG.
  - MSG (str): associated MSG, 1st MSG with type II in the same SGs.
  - lattice (str): A, B, C, P, I, F, R.
  - hexagonal_g (bool): trigonal or hexagonal ?
  - SO (list): symmetry operations.
  - SO_gen (list): generator of SOs.
"""

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
"""
Magnetic Point Group info type.

Fields:
  - tag (str): MPG_id (PG.no.ID).
  - international (str): international short symbol in LaTeX.
  - schoenflies (str): Schoenflies symbol in LaTeX.
  - crystal (str): triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
  - setting (str): setting.
  - PG (str): associated PG, unique.
  - SG (str): associated SG, MPG -> MSG -> SG.
  - MPG (str): = id.
  - MSG (str): associated MSG, 1st MSG in the same MPGs.
  - lattice (str): "0".
  - hexagonal_g (bool): trigonal or hexagonal ?
  - type (str): type I, II, III.
  - SO (list): symmetry operations.
"""

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
"""
Magnetic Space Group info type.
  - tag (str): BNS_id (SG.no).
  - BNS (str): Belov-Neronova-Smirnova (international) notation in LaTeX.
  - OG (str): Opechowski-Guccione notation in LaTeX.
  - schoenflies (str): Schoenflies symbol in LaTeX.
  - crystal (str): triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
  - setting (str): setting.
  - PG (str): associated PG, unique.
  - SG (str): associated SG, unique.
  - MPG (str): associated MPG, unique.
  - MSG (str): = id.
  - lattice (str): A, B, C, P, I, F, R.
  - hexagonal_g (bool): trigonal or hexagonal ?
  - type (str): type I, II, III, IV.
  - SO (list): symmetry operations.
"""

# ==================================================
SphericalMultipoleType = namedtuple("SphericalMultipoleType", ["X", "l", "s", "k", "x"])
"""
Spherical multipole type.

Fields:
  - type (str): type, Q, G, T, M.
  - rank (int): rank, 0-11.
  - s (int): internal rank, 0, 1.
  - k (int): internal component, -1, 0, 1.
  - internal-type (str) internal type, q or g.
"""

PGMultipoleType = namedtuple("PGMultipoleType", ["X", "l", "Gamma", "n", "p", "s", "k", "x"])
"""Point-group multipole type.

Fields:
  - type (str): type, Q, G, T, M.
  - rank (int): rank, 0-11.
  - irrep, (str): irrep.
  - multiplicity (int): multiplicity.
  - add_multiplicity (int): additional multiplicity.
  - s (int): internal rank, 0, 1.
  - k (int): internal component, -1, 0, 1.
  - internal-type (str): internal type, q or g.
"""

RepSiteType = namedtuple("RepSiteType", ["no", "wyckoff", "symmetry", "position", "orbital"])
"""
Representative site type.

Fields:
  - no (int): ID (from 1).
  - wyckoff (str): Wyckoff position.
  - symmetry (str): site symmetry.
  - position (list): site position.
  - orbital (list): orbital list.
"""

CellSiteType = namedtuple("CellSiteType", ["no", "position", "position_primitive", "mapping", "sublattice", "plus_set"])
"""
Cell site type.

Fields:
  - no (int): ID (from 1).
  - position (list): site position.
  - position_primitive (list): position in primitive cell.
  - mapping (list): SO mapping.
  - sublattice (int): sublattice no.
  - plus_set (int): plus set no.
"""

BondInfoType = namedtuple("BondInfoType", ["tail", "head", "neighbor", "t_rank", "h_rank"])
"""
Bond info type.

Fields:
  - tail (str): tail name.
  - head (str): head name.
  - neighbor (int): neighbor.
  - t_rank (int): tail orbital rank.
  - h_rank (int): head orbital rank.
"""

RepBondType = namedtuple(
    "RepBondType",
    ["no", "tail", "head", "neighbor", "wyckoff", "directional", "vector", "center", "distance", "t_rank", "h_rank"],
)
"""
Representative bond type.

Fields:
  - no (int): ID (from 1).
  - tail (str): tail name.
  - head (str): head name.
  - neighbor (int): neighbor.
  - wyckoff (str): Wyckoff position.
  - directional (str): (non) directional bond (N)D.
  - vector (list): bond vector.
  - center (list): bond center.
  - distance (float): bond length.
  - t_rank (int): tail orbital rank.
  - h_rank (int): head orbital rank.
"""

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
"""
Cell bond type.

Feilds:
  - no (int): ID (from 1).
  - vector (list): bond vector.
  - vector_primitive (list) bond vector in primitive cell.
  - center (list): bond center.
  - center_primitive (list): bond center in primitive cell.
  - mapping (list): SO mapping.
  - sublattice (int): sublattice no.
  - plus_set (int): plus set no.
  - t_idx (tuple): tail index, (sublattice, plus set).
  - h_idx (tuple): head index, (sublattice, plus set).
  - R_primitive (list): R vector in primitive cell.
"""

BraketInfoType = namedtuple("BraketInfoType", ["bh_rank", "bh_idx", "kt_rank", "kt_idx"])
"""
Braket info type.

Fields:
  - bh_rank (int): bra-head orbital rank.
  - bh_idx (int): bra-head orbital index.
  - kt_rank (int): ket-tail orbital rank.
  - kt_idx (int): ket-tail orbital index.
"""

SAMBType = namedtuple("SAMBType", ["head", "tail", "wyckoff", "bk_info"])
"""
Combined SAMB type.

Fields:
  - head (str): head name.
  - tail (str): tail name.
  - wyckoff (str): Wyckoff position.
  - bk_info (BraketInfoType): braket info.
"""

UniqueSAMBType = namedtuple("UniqueSAMBType", ["samb_type", "neighbor", "n"])
"""
Unique combined SAMB type.

Fields:
  - samb_type(SAMTType): SAMB type.
  - neighbor (int): neighbor.
  - n (int): multiplicity for same neighbor.
"""

SOMatrixType = namedtuple("SOMatrixType", ["X", "s"])
"""
Symmetry operation matrix type.

Fields:
  - X (str): SO matrix type, Q, G, T, M.
  - s (int): So matrix rank.
"""

IrrepListType = namedtuple("IrrepListType", ["X", "Gamma"])
"""
Irrep list type.

Fields:
  - X (str): type, Q, G, T, M.
  - Gamma (str): irrep.
"""

InternalBasisType = namedtuple("InternalBasisType", ["s", "x"])
"""
Internal basis type.

Fields:
  - s (int): internal rank.
  - x (str): internal type, q or g.
"""

RepMatrixType = namedtuple("RepMatrixType", ["Gamma"])
"""
Representation matrix type.

Fields:
  - Gamma (str): representation matrix irrep.
"""


__all__ = ["Group", "MaterialModel", "ModelAnalyzer", "create_samb", "create_samb_qtdraw", "create_samb_matrix"]

if TYPE_CHECKING:
    from multipie.core.group import Group
    from multipie.core.cmd import create_samb
    from multipie.core.cmd import create_samb_qtdraw
    from multipie.core.cmd import create_samb_matrix
    from multipie.core.material_model import MaterialModel
    from multipie.core.model_analyzer import ModelAnalyzer


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

    if name == "ModelAnalyzer":
        from multipie.core.model_analyzer import ModelAnalyzer

        return ModelAnalyzer

    raise AttributeError
