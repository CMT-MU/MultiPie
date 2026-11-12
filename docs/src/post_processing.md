# Post-processing tools

After a tight-binding model is built (by SAMBs, or from Wannier90 `seedname_hr.dat`), `mp_analyze` computes physical quantities from its eigenvalues and eigenvectors.
This page explains the tools, set in `output` of the control file:

| tool | keys in `output` | output files (in `model_name/output/`, prefixed with `model_name`) |
|---|---|---|
| [Density of states](#pp-density-of-states) | `dos`, `dos_e_num`, `dos_width`, `dos_method`, `dos_spin` | `_dos.txt`, `_dos.pdf`, `_dos.eps`, `plot_dos.gnu` |
| [Fermi energy from the filling](#pp-fermi-energy) | `n_electron` | (used for all output) |
| [Orbital weight (fat band)](#pp-orbital-weight-fat-band) | `dispersion/weight` | `_weight.txt`, `_weight.pdf`, `plot_weight.gnu` |
| [Eigenvalues and eigenvectors](#pp-eigenvalues-and-eigenvectors) | `save_k_path`, `save_k_full` | `_eigen_k_path.npz`, `_eigen_k_full.npz` |
| [Callbacks for custom calculations](#pp-callbacks-for-custom-calculations) | `callback_k_path`, `callback_k_full`, `callback_info` | (written by the callbacks) |
| [Local operators](#pp-local-operators) | `dispersion/local`, `callback_info` | `_dispersion.txt` (with `dispersion/local`), or passed to the callbacks |

All of them are off by default, and the dispersion (`_dispersion.txt`, `.pdf`, `.eps`, `plot_band.gnu`) is computed as before.

## Example: graphene

The files are in [docs/src/examples](https://github.com/CMT-MU/MultiPie/tree/main/docs/src/examples): the model input `graphene_in.py` (see [Getting Started](getting_started.md)), the control file `graphene_post_ctrl.py`, and the callbacks `graphene_callback.py`.

```{literalinclude} examples/graphene_post_ctrl.py
```

Put the three files in one directory, and run

```bash
$ mp_create graphene_in.py
$ mp_analyze -v graphene_post_ctrl.py
```

Then `graphene/output/` contains, in addition to the dispersion,

- `graphene_dos.txt`, `graphene_dos.pdf`, `graphene_dos.eps`, `plot_dos.gnu` : DOS.
- `graphene_weight.txt`, `graphene_weight.pdf`, `plot_weight.gnu` : orbital weight of `pz@C(1)`.
- `graphene_eigen_k_path.npz` : eigenvalues and eigenvectors along the k path.
- `my_band.txt`, `my_occupation.txt` : written by the callbacks.

The Fermi energy is 0 for `"n_electron": 1` (half filling of the two bands), and `-v` prints `fermi energy = 0.0 for n_electron = 1 (gaussian).`

(pp-density-of-states)=
## Density of states

When `output/dos` is not `False`, the density of states (DOS) per primitive cell is computed on the k grid given by `grid` (the Γ point only for a model of a point group).

- `"auto"` (or `True`): the partial DOS for each site cluster and orbital, with the spin summed, e.g., `s@A` for `(s,u)@A(1)`, `(s,d)@A(1)`, `(s,u)@A(2)`, .... With `"dos_spin": True`, the spin of a spinful basis is resolved, e.g., `(s,u)@A` and `(s,d)@A`.
- `[[ket name]]`: for each list of ket names, the sum of the partial DOS of the kets, e.g., `[["pz@C(1)"], ["pz@C(1)", "pz@C(2)"]]`, labeled by the ket names joined by `+`. The ket names are those of the basis in `info/model_name_hr.dat` (or `info/model_name_info_wannier.py` in wannier mode).

`model_name_dos.txt` has the columns: the energy measured from the Fermi energy (as in the dispersion), the total DOS, the integrated total DOS N(E), and the partial DOS. The comment lines at the top give the content of each column, e.g., for graphene,

```text
# DOS per primitive cell (summed over spin).
# shifted by fermi energy = 0.0
# number of k points = 14400 (grid = [120, 120, 1])
# Gaussian width = 0.02449489742783178
# column 1: E
# column 2: total
# column 3: N(E), integrated total DOS
# column 4: pz@C = pz@C(1) + pz@C(2)
```

The DOS is normalized so that the total DOS integrates to the number of bands. The energy range includes all the eigenvalues, with `dos_e_num` points (1000 by default).

The method is chosen by `dos_method`:

- `"gaussian"` (default): the delta functions are broadened by a Gaussian with the standard deviation `dos_width` (1% of the band width if `None`, and not smaller than the energy step). A width smaller than the energy step is an error; increase `dos_e_num` in that case (`dos_e_num` must be at least 12).
- `"tetrahedron"`: the optimized tetrahedron method by [libtetrabz](https://mitsuaki1987.github.io/libtetrabz/), which is installed separately (see [Installation](install.md)). N(E) is computed by the tetrahedron method, and the DOS is its difference over each energy step, $(N(E+dE/2)-N(E-dE/2))/dE$, so that the normalization is kept also for flat bands. It gives a sharp DOS without broadening (`dos_width` is not used), but cannot be used for a model of a point group. It is slower and uses more memory than `"gaussian"`; use a moderate `grid` and `dos_e_num` for a large model.

DOS of graphene with `"gaussian"` (left) and `"tetrahedron"` (right). The tetrahedron method resolves the van Hove singularities and the band edges sharply:

```{image} figures/post_processing/graphene_dos_gaussian.png
:width: 48%
```
```{image} figures/post_processing/graphene_dos_tetrahedron.png
:width: 48%
```

(pp-fermi-energy)=
## Fermi energy

The Fermi energy is 0 in samb mode, and that of `seedname.win` in wannier and symcw modes. When `output/n_electron` is given, the Fermi energy is determined from the eigenvalues on the full k grid so that `n_electron` states per primitive cell are occupied, and it is used for all the output (dispersion, DOS, orbital weight, saved eigen systems, and the values passed to the callbacks).

Each band holds one state, so that `n_electron` is the number of electrons for a spinful model, and half of it for a spinless model, e.g., `"n_electron": 1` for graphene at half filling.
With `"dos_method": "gaussian"`, the states are occupied in ascending order of energy (T = 0), and the Fermi energy is placed at the middle of the highest occupied and the lowest unoccupied states. With `"tetrahedron"`, it is computed by libtetrabz. For an insulator (integer `n_electron` with a gap), both methods place the Fermi energy at the middle of the gap.

With `output/dispersion/power`, a positive number $p$, the energies are raised to the power, $E_k^p$, e.g., $\omega_k = E_k^{1/2}$ for bosonic excitations. The Fermi energy is set to 0, and `output/n_electron` cannot be used together (also in wannier mode, the Fermi energy of `seedname.win` is not used). The dispersion and the orbital weight are written and plotted as $E_k^p$, and the DOS is that of $E_k^p$, i.e., $\delta(E - E_k^p)$ is used in place of $\delta(E - E_k)$, for both `"gaussian"` and `"tetrahedron"`. Where $E_k^p$ is not a real number, i.e., negative $E_k$ with non-integer $p$, $-|E_k|^p$ is used, so that the negative energies remain negative. For integer $p$, $E_k^p$ is used as it is, e.g., it is positive for even $p$. The eigenvalues in the `.npz` files and passed to the callbacks are not raised to the power.

(pp-orbital-weight-fat-band)=
## Orbital weight (fat band)

When `output/dispersion/weight` is given in the same way as `output/dos` (`"auto"` or `[[ket name]]`), the weight of the kets in each eigenvector, the sum of $|U_k[m,n]|^2$ over the kets $m$, is computed along the k path.
`model_name_weight.txt` has the columns: k, the energy measured from the Fermi energy (as in the dispersion), and the weights. The bands are plotted with the color of the weight, one page for each weight, in `_weight.pdf` and by `plot_weight.gnu`.
`"auto"` always sums over the spin (`dos_spin` is not used); give `[[ket name]]` for spin-resolved weights, e.g., `[["(pz,u)@C(1)", "(pz,u)@C(2)"]]`.

For graphene, the weight of `pz@C(1)` is 1/2 on every band, by the symmetry between the two sublattices:

```{image} figures/post_processing/graphene_weight.png
:width: 60%
```

## Example with several orbitals

For a model with several orbitals, the partial DOS and the orbital weights show the character of the bands. The following is a spinful model with s and p orbitals at A and an s orbital at B (CsCl type, space group 221, SAMBs of the identity irrep. with arbitrary parameters), computed with `"dos": "auto"`, `"dos_method": "tetrahedron"`, `"n_electron": 4`, and `"weight": "auto"`.

Dispersion and DOS (the sharp peak is the narrow band of s at B):

```{image} figures/post_processing/sp_dispersion.png
:width: 48%
```
```{image} figures/post_processing/sp_dos.png
:width: 48%
```

Orbital weights of `s@A` and `s@B` (two of the five pages of `_weight.pdf`):

```{image} figures/post_processing/sp_weight_s_A.png
:width: 48%
```
```{image} figures/post_processing/sp_weight_s_B.png
:width: 48%
```

(pp-eigenvalues-and-eigenvectors)=
## Eigenvalues and eigenvectors

With `"save_k_path": True` and/or `"save_k_full": True`, the eigenvalues and eigenvectors along the k path and/or on the full k grid are saved in `model_name_eigen_k_path.npz` and `model_name_eigen_k_full.npz`. They contain

- `Ek[k, n]`, `Uk[k, :, n]`, `k_grid`, `ket`, `fermi_energy`: as described for the callbacks below,
- `k_linear` (along the k path) or `grid` (on the full k grid).

They can be read by `numpy.load`, e.g., to analyze the eigenvectors in Jupyter:

```python
import numpy as np

z = np.load("graphene/output/graphene_eigen_k_path.npz")
Ek, Uk, ket = z["Ek"], z["Uk"], list(z["ket"])
weight = np.abs(Uk[:, ket.index("pz@C(1)"), :]) ** 2  # weight[k, n] of pz@C(1).
```

As for the callbacks, `save_k_path` requires the dispersion (`output/dispersion/k_path` is not `None`, and the model is of a space group), and `save_k_full` keeps all the eigenvectors on the full grid in memory.

(pp-callbacks-for-custom-calculations)=
## Callbacks for custom calculations

Quantities not provided by MultiPie can be computed by user-defined functions, which receive the eigenvalues and eigenvectors of H(k).
Give a Python file, relative to the top directory (where the control file is), in `output/callback_k_path` and/or `output/callback_k_full`, and define the function of the same name in the file. One file can define both, and it is loaded once.

The example `graphene_callback.py` writes the band energies along the k path, and the occupation of each ket at $T = 0$, $n_m = (1/N_k) \sum_{k,n} \theta(E_F - E_{kn}) |U_k[m,n]|^2$:

```{literalinclude} examples/graphene_callback.py
```

For graphene at half filling, `my_occupation.txt` is

```text
pz@C(1) 0.500000
pz@C(2) 0.500000
```

The arguments of the callbacks are

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
  - for `callback_k_full` only, the quantities given in `output/callback_info`, see below.

The functions are called after the other output, at the end of the analysis. `callback_k_path` requires the dispersion. For `callback_k_full`, all the eigenvectors on the full grid are kept in memory, `N1 x N2 x N3 x (number of kets)^2` complex numbers, so that a coarse `grid` may be needed for a large model.

Set `output/callback_info` to pass the expectation value of local operators in each eigenstate on the full k grid, as `info[keyword]`, `[k, n]`. The keywords are those of `output/dispersion/local`: spin `Sx`, `Sy`, `Sz`, orbital angular momentum `Lx`, `Ly`, `Lz`, and quadrupole `Qu`, `Qv`, `Qyz`, `Qxz`, `Qxy` (see [Local operators](#pp-local-operators) for their definitions). For example,

```python
"output": {
    "callback_k_full": "my_callback.py",
    "callback_info": ["Sx", "Sy", "Sz"],
},
```

gives `info["Sx"]`, `info["Sy"]`, and `info["Sz"]` in `callback_k_full(Ek, Uk, info)`. `output/callback_info` requires `output/callback_k_full` and a model in the `lg` or `lgs` basis, and the spin operators require a spinful model.

```{note}
Within degenerate bands, the eigenvectors, and hence the value of each band, are arbitrary. For example, the bands of a spinful model with both inversion and time-reversal symmetries are doubly degenerate at every k. Only the sum over the degenerate bands is meaningful. The same holds for `output/dispersion/local`.

The operators are defined along the axes of the model. With a Wannier input, the axes of the projections in `seedname.nnkp` (z-axis, x-axis, and spin-quantization axis) are not taken into account.
```

```{note}
Unlike the input files, the callback files are executed as Python code. Use only callback files you trust.
```

(pp-local-operators)=
## Local operators

The local operators of `output/dispersion/local` and `output/callback_info` are matrices $O_{mm'}$ in the atomic basis of the model, i.e., the tesseral orbitals `s`, `px`, `py`, `pz`, `dv`, `dxy`, `dxz`, `dyz`, `du`, `fc`, `fs`, `fbz`, `f3`, `fx`, `fy`, `faz`, with the spin `u`, `d` for a spinful model. The value for the eigenstate $n$ at $k$ is the expectation value $\langle O \rangle_{kn} = \sum_{mm'} U_k[m,n]^* \, O_{mm'} \, U_k[m',n]$.

- **Spin** `Sx`, `Sy`, `Sz`: $\boldsymbol{S} = \boldsymbol{\sigma}/2$ (in units of $\hbar$), where $\boldsymbol{\sigma}$ are the Pauli matrices in the order of (`u`, `d`).
- **Orbital angular momentum** `Lx`, `Ly`, `Lz`: the usual angular momentum operators (in units of $\hbar$) for the orbitals of rank $l$, satisfying $[L_x, L_y] = iL_z$ and $\boldsymbol{L}^2 = l(l+1)$.
- **Quadrupole** `Qu`, `Qv`, `Qyz`, `Qxz`, `Qxy`: the matrix elements $\langle l m | Q(\hat{\boldsymbol{r}}) | l m' \rangle$ of the angular functions

  | keyword | $Q(\hat{\boldsymbol{r}})$ |
  |---|---|
  | `Qu` | $(3z^2 - r^2)/(2r^2)$ |
  | `Qv` | $\sqrt{3}(x^2 - y^2)/(2r^2)$ |
  | `Qyz` | $\sqrt{3}\,yz/r^2$ |
  | `Qxz` | $\sqrt{3}\,zx/r^2$ |
  | `Qxy` | $\sqrt{3}\,xy/r^2$ |

  For example, $\langle p_z | Q_u | p_z \rangle = 2/5$, $\langle d_u | Q_u | d_u \rangle = 2/7$, and $\langle f_{az} | Q_u | f_{az} \rangle = 4/15$.

The operators act within each orbital rank $l$ of each site (atom and sublattice): the matrix elements between different sites, and between different ranks such as p and d, are zero. $\boldsymbol{L}$ and $Q$ are diagonal in the spin, and $\boldsymbol{S}$ is diagonal in the orbital. For the s orbital ($l = 0$), $\boldsymbol{L}$ and $Q$ are zero. The local operators are available for a model in the `lg` or `lgs` basis, not in the `jml` basis.

## Plotting by gnuplot

The plots in PDF and EPS are written by matplotlib. The gnuplot scripts write the same plots, and are convenient to adjust them:

```bash
$ cd graphene/output
$ gnuplot plot_dos.gnu      # graphene_dos.pdf and .eps
$ gnuplot plot_weight.gnu   # graphene_weight.pdf (one page for each weight)
```

The PDF/EPS files are overwritten. The plot ranges and the Fermi level are variables at the top of the scripts, e.g., in `plot_dos.gnu`,

```text
Ef = 0.0  # fermi energy (data are already shifted by Ef).
Ef_line = 0.0  # position of fermi level in plot.
xmin = -1.347219358530748
xmax = 1.347219358530748
ymin = 0.0
ymax = 1.9823887755144438
```

and `plot_weight.gnu` has `labels`, the list of the weights to plot.
