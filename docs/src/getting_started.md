# Getting Started

## Tutorial

As a tutorial, we describe the procedure for generating model in the case of **graphene**.

1. Download `graphene_in.py` and `graphene_ctrl.py` from [examples](https://github.com/CMT-MU/MultiPie/tree/main/docs/src/examples) directory.
2. At the same place where the above two files are, run the following:

    ```bash
    $ mp_create -v graphene_in.py
    ```

    It creates the following files under `graphene` directory:

     - `graphene.pkl` : Model information file (binary)
     - `graphene.tex`, `graphene.pdf` : Summary of the model information
     - `graphene.qtdw` : QtDraw file of the model structure

    The `.tex`/`.pdf` files are created only when `"pdf": {"create": True}` (default) and LaTeX is installed (see [Installation](install.md) for the required TeX components).
    The `.qtdw` file is created only when `"qtdraw": {"create": True}` (default) and QtDraw is installed.

3. To analyze the model, e.g., draw dispersion, run the following:

    ```bash
    $ mp_analyze -v graphene_ctrl.py
    ```

    It creates the following files under `graphene` directory:

     - `info/graphene_matrix.py`, `info/graphene_hr.dat` : Full matrix information for selected SAMBs
     - `output/` : Various output for physical quantities
     - `samb/` : QtDraw files for SAMBs. The directory is created when `samb_figure` is `True`; the files are written only when QtDraw is installed and the model was created with `"qtdraw": {"create": True}`.

4. In order to handle `graphene.pkl` interactively, use IPython interface. See in detail [analyze_model.ipynb](examples/analyze_model.ipynb)

## Model input file

In case of graphene, the model input file is

```{literalinclude} examples/graphene_in.py
```

The other setting are provided as default values, which are given as follows:

```{literalinclude} examples/default_model.py
```

You can overwrite whatever you want as in the case of graphene.

## Control file

In case of graphene, the control file is

```{literalinclude} examples/graphene_ctrl.py
```

The default values of control file are provided as follows:

```{literalinclude} examples/default_control.py
```

The typical use of generating SAMBs, you first create the SAMBs for all irreps., and then choose the necessary SAMBs, such as the symmetry-breaking terms in addition to the identity irreps., by specifying `samb/select` and/or `samb/parameter` in the control file.

## Format of input files

Model input files and control files are Python files that contain dictionaries, e.g., `graphene_in = {...}`.

- They are not executed as Python code. Each dictionary is read with `ast.literal_eval`, so only literals (strings, numbers, lists, tuples, dictionaries, `True`/`False`/`None`) are allowed. Arithmetic such as `2*2`, variables and `import` cannot be used.
- Fractions and symbolic values are given as strings, e.g., `"[1/3,2/3,0]"`. Strings representing mathematical expressions are parsed by SymPy.
- With `mp_create` (or `create_model`), every dictionary in a model input file, `name = {...}`, is treated as a separate model. `MaterialModel.analyze(filename)` expects one model dictionary per file.

```{note}
Strings representing mathematical expressions in input files are parsed by SymPy's `parse_expr`, which uses `eval` internally, and the model file `model_name.pkl` is stored with Python's `pickle`.
Use input files and `.pkl` files only from trusted sources.
```

## Output files

The description of the output files is as follows (with model_name prefixed):

- **model_name** : all files are created under `model_name`.
  - **.pdf** : model info. in PDF.
  - **.pkl** : model data (binary).
  - **.qtdw** : model structure in QtDraw.
  - **.tex** : model info. source.
  - **samb**/
    - **_atomic_samb.qtdw** : atomic SAMBs.
    - **_A_def.qtdw** : definition of A-site cluster.
    - **_A.qtdw** : A-site-cluster SAMBs.
    - **_A;B_n_m_def.qtdw** : definition of A;B (n-neighbor, m) bond cluster.
    - **_A;B_n_m.qtdw** : A;B (n-neighbor, m) bond-cluster SAMBs.
  - **info/**
    - **_hr.dat** : H[R] matrix data.
    - **_info_output.py** : info. of output.
    - **_info_samb.py** : info. of model and SAMBs.
    - **_info_wannier.py** : info. of wannier.
    - **_info.py** : global info.
    - **_k.py** : SAMBs in momentum rep.
    - **_matrix.py** : SAMB matrix data.
    - **_var.py** : relation between zj and atomic parameters at bond 1.
    - **_z.py** : zj parameters.
  - **output/**
    - **_dispersion.eps** : dispersion in EPS.
    - **_dispersion.pdf** : dispersion in PDF.
    - **_dispersion.txt** : dispersion data.
    - **plot_band.gnu** : command script to plot by gnuplot.
