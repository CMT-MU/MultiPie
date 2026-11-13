# Overview

**MultiPie** is a database and group-theoretical manipulation tool for crystallographic point groups, space groups, and magnetic point and space groups. Using these functionalities, one can construct the **symmetry-adapted multipole basis (SAMB)** for both point and space groups.

The SAMB forms a complete basis set and is therefore useful for describing arbitrary physical quantities according to their symmetry classifications. Any physical quantity can also be decomposed in terms of the SAMB.

The brief manual and a tutorial of MultiPie in PDF are given at the end of the {ref}`User Guide <brief-manual>`. Moreover, information for all supported groups is summarized in [Group Information](https://cmt-mu.github.io/MultiPieGen/src/group_info.html).

The [User Guide](guide.md) explains how to read and use the harmonics and the SAMB (in English and Japanese), and [Examples of materials](materials.md) builds the models of SrVO$_3$, MoS$_2$, and Te step by step. For the users of Ver.2, [Changes from Ver.2](changes.md) summarizes what changes in Ver.3.

By installing the Python package [QtDraw](https://github.com/CMT-MU/QtDraw), crystal and molecular structures, as well as polar and axial vectors, orbitals, and the SAMB, can be visualized. This visualization capability makes the tool especially helpful for learning group theory in a transparent and intuitive manner. A practical tutorial of MultiPie via QtDraw is given by the {ref}`Tutorial via QtDraw <brief-manual>`.

For the API, please refer to the following:

- [Core module](api_core.md) for `Group`, `MaterialModel`, `ModelAnalyzer`, and the create/analyze functions (`cmd`)
- [Utility module](api_util.md) for various utilities for implementation purpose.
