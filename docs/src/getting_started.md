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

4. In order to handle `graphene.pkl` interactively, use the IPython interface. See [analyze_model.ipynb](examples/analyze_model.ipynb) for details.

## Model input file

In case of graphene, the model input file is

```{literalinclude} examples/graphene_in.py
```

The other settings are provided as default values, which are given as follows:

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

In a typical use, you first create the SAMBs for all irreps, and then choose the necessary SAMBs, such as the symmetry-breaking terms in addition to the identity irreps., by specifying `samb/select` and/or `samb/parameter` in the control file.

### Density of states

When `output/dos` is not `False`, the density of states (DOS) per primitive cell is computed on the k grid given by `grid` (the Γ point only for a model of a point group), and written to `output/model_name_dos.txt`.

- `"auto"` (or `True`): the partial DOS for each site cluster and orbital, with the spin summed, e.g., `s@A` for `(s,u)@A(1)`, `(s,d)@A(1)`, `(s,u)@A(2)`, .... With `"dos_spin": True`, the spin of a spinful basis is resolved, e.g., `(s,u)@A` and `(s,d)@A`.
- `[[ket name]]`: for each list of ket names, the sum of the partial DOS of the kets, e.g., `[["pz@C(1)"], ["pz@C(1)", "pz@C(2)"]]`, labeled by the ket names joined by `+`. The ket names are those of the basis in `info/model_name_hr.dat` (or `info/model_name_info_wannier.py` in wannier mode).

The columns are the energy measured from the Fermi energy, as in the dispersion, the total DOS, the integrated total DOS N(E), and the partial DOS; the comment lines at the top give the content of each column. The DOS is normalized so that the total DOS integrates to the number of bands. The energy range includes all the eigenvalues, with `output/dos_e_num` points. The plot (`_dos.pdf`, `_dos.eps`) and the gnuplot script (`plot_dos.gnu`) are also written.

The method is chosen by `output/dos_method`:

- `"gaussian"` (default): the delta functions are broadened by a Gaussian with the standard deviation `output/dos_width` (1% of the band width if `None`, and not smaller than the energy step). A width smaller than the energy step is an error; increase `dos_e_num` in that case (`dos_e_num` must be at least 12).
- `"tetrahedron"`: the optimized tetrahedron method by [libtetrabz](https://mitsuaki1987.github.io/libtetrabz/), which is installed separately (see [Installation](install.md)). N(E) is computed by the tetrahedron method, and the DOS is its difference over each energy step, (N(E+dE/2)-N(E-dE/2))/dE, so that the normalization is kept also for flat bands. It gives a sharp DOS without broadening (`dos_width` is not used), but cannot be used for a model of a point group. It is slower and uses more memory than `"gaussian"`; use a moderate `grid` and `dos_e_num` for a large model.

### Fermi energy

The Fermi energy is 0 in samb mode, and that of `seedname.win` in wannier and symcw modes. When `output/n_electron` is given, the Fermi energy is determined from the eigenvalues on the full k grid so that `n_electron` states per primitive cell are occupied, and it is used for all the output (dispersion, DOS, and the values passed to the callbacks).
Each band holds one state, so that `n_electron` is the number of electrons for a spinful model, and half of it for a spinless model. With `"dos_method": "gaussian"`, the states are occupied in ascending order of energy (T = 0), and the Fermi energy is placed at the middle of the highest occupied and the lowest unoccupied states. With `"tetrahedron"`, it is computed by libtetrabz. For an insulator (integer `n_electron` with a gap), both methods place the Fermi energy at the middle of the gap.

### Orbital weight (fat band)

When `output/dispersion/weight` is given in the same way as `output/dos` (`"auto"` or `[[ket name]]`), the weight of the kets in each eigenvector, sum of |U_k[m,n]|^2 over the kets m, is written along the k path in `output/model_name_weight.txt` (columns: k, energy measured from the Fermi energy as in the dispersion, weights). `"auto"` always sums over the spin (`dos_spin` is not used); give `[[ket name]]` for spin-resolved weights, e.g., `[["(pz,u)@C(1)", "(pz,u)@C(2)"]]`. The bands are plotted with the color of the weight, one page for each weight, in `_weight.pdf` and by the gnuplot script `plot_weight.gnu`.

### Eigenvalues and eigenvectors

With `"save_k_path": True` and/or `"save_k_full": True` in `output`, the eigenvalues and eigenvectors along the k path and/or on the full k grid are saved in `output/model_name_eigen_k_path.npz` and `output/model_name_eigen_k_full.npz`, which can be read by `numpy.load`. They contain `Ek`, `Uk`, `k_grid`, `ket`, and `fermi_energy` as described for the callbacks below, and also `k_linear` (along the k path) or `grid` (on the full k grid). As for the callbacks, `save_k_path` requires the dispersion, and `save_k_full` keeps all the eigenvectors on the full grid in memory.

For example, the following `output` in the control file computes the DOS by the tetrahedron method with the Fermi energy for 4 occupied states per primitive cell, the fat band of the s orbital at A, and calls the callbacks in `my_callback.py` (see below):

```python
"output": {
    "dos": "auto",
    "dos_method": "tetrahedron",
    "n_electron": 4,
    "callback_k_path": "my_callback.py",
    "dispersion": {"k_path": "", "weight": [["s@A(1)"]]},
},
```

### Callback for custom calculations

Quantities not provided by MultiPie can be computed by user-defined functions, which receive the eigenvalues and eigenvectors of H(k).
Give a Python file, relative to the top directory (where the control file is), in `output/callback_k_path` and/or `output/callback_k_full`, and define the function of the same name in the file (one file can define both):

```python
# my_callback.py
import os
import numpy as np


def callback_k_path(Ek, Uk, info):
    # eigen system along the k path of the dispersion.
    # "output" is the default of output/dir; change it if output/dir is changed.
    np.savetxt(os.path.join(info["_topdir"], info["name"], "output", "my_band.txt"), Ek - info["fermi_energy"])


def callback_k_full(Ek, Uk, info):
    # eigen system on the full k grid.
    ...
```

- `Ek[k, n]`: eigenvalues in ascending order (not shifted by the Fermi energy, and not raised to `output/dispersion/power`).
- `Uk[k, :, n]`: eigenvector of `Ek[k, n]`, i.e., `H(k) = Uk @ diag(Ek) @ Uk^†` for each k, in the order of `info["ket"]` and in the gauge given by `output/fourier/tb_gauge`.
- `info` (dict):
  - `k_path` (str), `k_point` (dict, label to fractional coordinates): k path and k points of the dispersion (`None` if the dispersion is not computed).
  - `k_grid` (ndarray): k points of `Ek` and `Uk` (fractional, in the primitive cell), `[number of k, 3]`.
  - `ket` (list): ket names.
  - `fermi_energy` (float), `name` (str): Fermi energy, and model name (or seedname in wannier mode without model).
  - `_topdir` (str), `_verbose` (bool): top directory (absolute path), and verbose option.
  - for `callback_k_path` only, `k_linear` (ndarray) and `k_dis_pos` (dict): position along the k path, and positions of the labels, as used to plot the dispersion.
  - for `callback_k_full` only, `grid` (list): `[N1, N2, N3]`.

The functions are called after the other output, at the end of the analysis. `callback_k_path` requires the dispersion, i.e., `output/dispersion/k_path` is not `None`. For `callback_k_full`, all the eigenvectors on the full grid are kept in memory, `N1 x N2 x N3 x (number of kets)^2` complex numbers, so that a coarse `grid` may be needed for a large model.

```{note}
Unlike the input files, the callback files are executed as Python code. Use only callback files you trust.
```

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
- `wannier/read_KS` has been removed: `"read_KS": True` is an error, and `"read_KS": False` in an old control file is ignored with a warning (remove it from the control file).

```{note}
Strings representing mathematical expressions in input files are parsed by SymPy's `parse_expr`, which uses `eval` internally, and the model file `model_name.pkl` is stored with Python's `pickle`.
Use input files and `.pkl` files only from trusted sources.
```

## Command-line tools

- `mp_create` and `mp_analyze` exit with status 1 when an error occurs, and with status 2 for a usage error, e.g., no input file or an unknown option. The error is printed in one message with the input file and, when applicable, the model and the site concerned; with `-v`, the full traceback is printed instead.
- When several input files are given, all of them are read before any model is created or analyzed.
- `mp_create -i` and `mp_analyze -i` print the default model and control, which can be used as a template of an input file.
- `python -m multipie.scripts.mp_create` and `python -m multipie.scripts.mp_analyze` work in the same way.
- Messages such as the elapsed time and warnings are written to stdout through the `multipie` logger. `create_model` and `analyze_model` (and so `mp_create` and `mp_analyze`) set it up unless it already has a handler; the logging configuration of the application, e.g., the root logger in Jupyter, is not changed. To collect the messages elsewhere, add a handler to the `multipie` logger before calling them.

## Parallel computation

The atomic multipole matrices, `create_atomic_multipole_matrix` in `multipie.util.util_atomic_multipole`, and the atomic SAMBs, `Group.create_atomic_samb_L`, are computed in parallel with joblib.
These are used to generate the database.
In `mp_create` (`MaterialModel.analyze`), the orthogonalization of the atomic and combined SAMBs is also done in parallel for a large model, e.g., a spinful model with d orbitals; a small model is processed serially, since starting the processes would take longer than the computation.
The output does not depend on the number of processes.
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
    - **_A-B-..._def.qtdw** : definition of the site clusters (all sites in one file).
    - **_A-B-....qtdw** : site-cluster SAMBs (all sites in one file).
    - **_A;B_nnn_m_def.qtdw** : definition of A;B (nnn-neighbor, m) bond cluster (nnn: neighbor, 3 digits, e.g., 001).
    - **_A;B_nnn_m.qtdw** : A;B (nnn-neighbor, m) bond-cluster SAMBs.
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
    - **_dos.eps** : DOS in EPS.
    - **_dos.pdf** : DOS in PDF.
    - **_dos.txt** : DOS data (when `output/dos` is not `False`).
    - **_eigen_k_path.npz**, **_eigen_k_full.npz** : eigenvalues and eigenvectors (when `output/save_k_path` or `output/save_k_full` is `True`).
    - **_weight.pdf** : dispersion colored by the orbital weight.
    - **_weight.txt** : dispersion with the orbital weight (when `output/dispersion/weight` is not `False`).
    - **plot_dos.gnu**, **plot_weight.gnu** : command scripts to plot DOS and orbital weight by gnuplot.
    - **plot_band.gnu** : command script to plot by gnuplot.
