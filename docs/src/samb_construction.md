# Example for SAMB Construction

1. Prepare input file. Example input files for C3v molecule and graphene are given as follows.
    - C3v molecule
    ```{literalinclude} ../examples/C3v.py
    ```
    - graphene
    ```{literalinclude} ../examples/graphene.py
    ```

1. At the folder where the input file exists, do the following to create SAMB.
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
