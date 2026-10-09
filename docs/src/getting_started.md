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

    The `.tex`/`.pdf` files are created only when `"pdf": {"create": True}` (default) and the TeX components listed in [Installation](install.md) (pLaTeX, `ptex2pdf` and the LaTeX packages) are available:

     - If `latex` is not found, neither file is written (with `-v`, a notice is printed).
     - Otherwise the `.tex` file is written first. If `ptex2pdf` or a LaTeX package is then found missing, or the compilation fails or does not finish, a warning is printed and the `.pdf` file, from an earlier run or incomplete, is removed (unless it cannot be removed, e.g., for lack of permission).
     - The `.pkl` file and the other outputs are written in any case.

    The `.qtdw` file is created only when `"qtdraw": {"create": True}` (default) and QtDraw is installed (with `-v`, a notice is printed if QtDraw is not found).

3. To analyze the model, e.g., draw dispersion, run the following:

    ```bash
    $ mp_analyze -v graphene_ctrl.py
    ```

    It creates the following files under `graphene` directory:

     - `info/graphene_matrix.py`, `info/graphene_hr.dat` : Full matrix information for selected SAMBs
     - `output/` : Various output for physical quantities
     - `samb/` : QtDraw files for SAMBs, written when `samb_figure` is `True`, QtDraw is installed, and the model was created with `"qtdraw": {"create": True}`. If `samb_figure` is `True` but one of the others is not satisfied, a warning is printed and the directory is not created.

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
- With `mp_create` (or `create_model`), every dictionary in a model input file, `name = {...}`, is treated as a separate model, and with `mp_analyze` (or `analyze_model`), every dictionary in a control file is treated as a separate control. `MaterialModel.analyze(filename)` and `ModelAnalyzer.analyze(filename)` expect one dictionary per file.
- Unknown keys, e.g., a typo such as `"spinfull"`, are errors: `ValueError: unknown key(s) in model 'x': 'spinfull' (did you mean 'spinful'?).` The accepted keys are those in the default model and the default control above, and also
  - `cell`: `a`, `b`, `c`, `alpha`, `beta`, `gamma`,
  - `SAMB_select`, `atomic_select`, `site_select`, `bond_select`: `X`, `l`, `Gamma`, `s`,
  - `samb/select`: `site`, `bond`, `X`, `l`, `Gamma`, `s`.

  The keys of the following dictionaries are not checked against the defaults, but have their own rules:

  - site names in `site`: non-empty strings without `;`, and without `_` for sites used in `bond` (`;` and `_` are used in the names of site and bond clusters). Letters and digits, e.g., `"Fe1"`, are recommended for all sites, since the names also appear in the PDF (LaTeX) and in file names. LaTeX special characters, e.g., `#` or `%`, in the model and site names are escaped in the PDF.
  - SAMB names in `samb/parameter`: names of the generated SAMBs, e.g., `"z1"`. `samb/parameter` may also be the name of a Python file containing one dictionary of the same form, given relative to `model_name/info/`, with the extension, e.g., `"my_z.py"` for `model_name/info/my_z.py` (not `"info/my_z.py"`, `"my_z"` or an absolute path). When the parameters are not empty, `mp_analyze` writes them to `model_name/info/model_name_z.py`, which is overwritten at each run; copy it to another name to edit it.
  - k-point labels in `output/dispersion/k_point`: labels used in `k_path`, without `-`, `|` and spaces, which separate the points in `k_path`.

```{note}
Strings representing mathematical expressions in input files are parsed by SymPy's `parse_expr`, which uses `eval` internally, and the model file `model_name.pkl` is stored with Python's `pickle`.
Use input files and `.pkl` files only from trusted sources.
```

## Command-line tools

- `mp_create` and `mp_analyze` exit with status 1 when an error occurs, and with status 2 for a usage error, e.g., no input file or an unknown option. The error is printed in one message with the input file and, when applicable, the model and the site concerned; with `-v`, the full traceback is printed instead.
- When several input files are given, all of them are read before any model is created or analyzed.
- `mp_create -i` and `mp_analyze -i` print the default model and control, which can be used as a template of an input file.
- `python -m multipie.scripts.mp_create` and `python -m multipie.scripts.mp_analyze` work in the same way.
- Messages such as the elapsed time and warnings are written to stdout through the `multipie` logger. `create_model` and `analyze_model` (and so `mp_create` and `mp_analyze`) set it up unless it already has a handler; the logging configuration of the application, e.g., the root logger in Jupyter, is not changed.

## Parallel computation

The atomic multipole matrices, `create_atomic_multipole_matrix` in `multipie.util.util_atomic_multipole`, and the atomic SAMBs, `Group.create_atomic_samb_L`, are computed in parallel with joblib.
These are used to generate the database; `mp_create` and `mp_analyze` use the precomputed data and are not affected.
The number of processes is given by the environment variable `MULTIPIE_N_JOBS`:

- not set or empty: all cores (`-1`, default),
- `1`: serial,
- `n > 1`: `n` processes, `n < 0`: number of cores + 1 + n processes (at least 1), as in joblib,
- `0` or a non-integer value is an error.

```bash
$ MULTIPIE_N_JOBS=4 python my_script.py   # a script calling the functions above
```

In PowerShell on Windows, use `$env:MULTIPIE_N_JOBS = "4"` before running the script; the setting remains for the session until it is removed by `Remove-Item env:MULTIPIE_N_JOBS`.

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
