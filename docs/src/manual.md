# Brief Manual

Two documents in PDF introduce MultiPie: a brief manual of the program, and a tutorial article on learning multipoles
and group theory with MultiPie and QtDraw. Both are available in English (en) and Japanese (ja).

:::{note}
Both documents were written for MultiPie Ver.2 (the brief manual for Version 2.3, the tutorial for MultiPie 2.0.7 and
QtDraw 3.1.0). The concepts and the workflow are the same in Ver.3, but some names and results differ, e.g., the f
orbitals `f2`, `f1`, `f3x`, `f3y` are renamed, `Group(tag, with_opt=True)` is removed, and the numbering and the signs of
SAMB can differ. See [Changes from Ver.2](changes.md) and the [User Guide](guide.md) for the current version.
:::

## MultiPie brief manual

[[en](brief_manual_en.pdf), [ja](brief_manual_ja.pdf)]

An overview of MultiPie in slides:

- **Overview.** The group-theoretical data (symmetry operations, Wyckoff sites and bonds, characters, harmonics, and
  active multipoles of the point, space, magnetic point, and magnetic space groups), and the model construction and
  the physical quantities computed with SAMB.
- **SAMB.** The physical degrees of freedom classified by the time-reversal parity, the spatial-inversion property, and
  the spatial anisotropy: the electric (`Q`), magnetic (`M`), magnetic toroidal (`T`), and electric toroidal (`G`)
  multipoles.
- **Example of graphene.** The model input, the Wyckoff sites and bonds, the site and bond clusters, the combined SAMB,
  the band dispersion, the matrix elements, and how to use the basis.
- **Workflow.** Model construction with `mp_create` from a model input file, and the computation of the band
  dispersion and other physical quantities with `mp_analyze` from a control file, with all the keys of the two files.
- **Classes.** Overviews of `Group`, `MaterialModel`, and `ModelAnalyzer`, with their initialization, properties, and
  main functions.

## Tutorial of MultiPie via QtDraw

[[en](tutorial_via_qtdraw_en.pdf), [ja](tutorial_via_qtdraw_ja.pdf)]

*Multipole Representation and Group Theory Learned through MultiPie and QtDraw*, by Junya Otsuki (Okayama University):
an introduction to MultiPie from the user's point of view, using the graphical interface of QtDraw.

- **Basic usage.** Setting up the environment, starting QtDraw, and the window for the functions of MultiPie.
- **Group theory.** The crystal-field splitting of orbitals by the harmonics decomposition, and the irreducible
  decomposition of the symmetric and anti-symmetric products of irreps.
- **Atomic multipoles.** The multipole matrices of the p orbitals and their normalization.
- **Cluster multipoles.** The crystal structure of URhSn (D$_{\rm 3h}$), and the classification of the electric
  quadrupoles on the U sites by cluster multipoles, e.g., a polar and a chiral cluster order.
- **iPython.** Using MultiPie from the iPython interactive interface, e.g., to generate figures of all multipoles.
