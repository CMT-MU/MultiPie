# Example for SAMB Construction

1. Prepare input file. Example input files for C3v molecule and graphene are given as follows.
    - C3v molecule
    ```python
    # C3v.py
    C3v = {
        "model": "C3v",  # name of model.
        "group": "C3v-1",  # name of point group.
        "cell": {"c": 10},  # set large enough interlayer distance.
        #
        "site": {"A": ("[-1/6,-1/6,0]", "s"), "B": ("[-2/3,0,0]", "p")},  # positions of A and B sites and their orbitals.
        "bond": [("A", "A", 1), ("A", "B", 1)],  # nearest-neighbor A-A and B-B bonds.
        #
        "spinful": False,  # spinless.
    }
    ```
    - graphene
    ```python
    # graphene.py
    graphene = {
        "model": "graphene",  #  name of model.
        "group": 191,  # No. of space group.
        "cell": {"c": 4},  # set large enough interlayer distance.
        #
        "site": {"C": ("[1/3,2/3,0]", "pz")},  # positions of C site and its orbital.
        "bond": [("C", "C", [1, 2, 3, 4, 5, 6])],  # C-C bonds up to 6th neighbors.
        #
        "spinful": False,  # spinless.
        #
        "k_point": {"Γ": "[0, 0, 0]", "M": "[1/2, 0, 0]", "K": "[1/3, 1/3, 0]"},  # def. of k points.
        "k_path": "Γ-K-M-Γ",  # high-symmetry line.
    }
    ```
2. At the folder where the input file exists, do the following to create SAMB.
See for more detail, try `create_samb --help` command.
    ```
    $ create_samb C3v graphene
    ```
1. The following files are created in `C3v` and `graphene` folders.
    - C3v_model.py, C3v_samb.py, C3v_matrix.py, C3v_samb.tex, C3v_samb.pdf, C3v_view.qtdw
    - graphene_model.py, graphene_samb.py, graphene_matrix.py, graphene_samb.tex, graphene_samb.pdf, graphene_view.qtdw

    Here, `.tex` and `.pdf` are created if [TeXLive](https://www.tug.org/texlive/) is installed, and `.qtdw` is created if [QtDraw](https://github.com/CMT-MU/QtDraw) is installed.

    Each file contains
    - `_model.py` : model information.
    - `_samb.py` : detailed information on SAMB.
    - `_matrix.py` : full-matrix form of SAMB.
    - `_samb.tex` and `_samb.pdf` : detailed information on SAMB in LaTeX and PDF format.
    - `_view.qtdw` : molecular or crystal structure file for QtDraw.

The detailed file formats are given in **File Format** section.
