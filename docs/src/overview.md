# Overview

**MultiPie** is a database and group-theoretical manipulation tool for crystallographic point groups, space groups, and magnetic point and space groups. Using these functionalities, one can construct the **symmetry-adapted multipole basis (SAMB)** for both point and space groups.

The SAMB forms a complete basis set and is therefore useful for describing arbitrary physical quantities according to their symmetry classifications. Any physical quantity can also be decomposed in terms of the SAMB.

By installing the Python package [QtDraw](https://github.com/CMT-MU/QtDraw), crystal and molecular structures, as well as polar and axial vectors, orbitals, and the SAMB, can be visualized. This visualization capability makes the tool especially helpful for learning group theory in a transparent and intuitive manner.

Information for all supported groups is summarized in [Group Information](group_doc.md).

For implementation-related notes (in Japanese) and the API, please refer to the following:
- [Technical Note](tech_note.md) (in Japanese)
- API
    - [Core module](api_core.md) for `Group` and `MaterialModel`
    - [Utility module](api_util.md) for various utities for implementation purpose.
