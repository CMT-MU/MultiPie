# ModelAnalyzer

## Dict data

This class has the following dict data.

- **info** (dict): General info.
  - **mode** (str): analysis mode, samb, wannier, or symcw.
  - **grid** (list): grid size, N1, N2, N3.
  - **name** (str): model and analysis name (seedname in wannier mode without model).
  - **A** (list): translational vectors for primitive cell, [a1,a2,a3] (3x3).
  - **B** (list): translational vectors for primitive reciprocal cell, [b1,b2,b3] (3x3).
  - **volume** (float): volume of primitive cell.
  - **fermi_energy** (float): Fermi energy.

- **samb** (dict): SAMB related.

- **wannier** (dict): Wannier related.
  - **ket** (list): ket name of full matrix in MultiPie standard.
  - **atoms_frac** (list): atom position in primitive cell (fractional) [in multipie order].
  - **atoms_cart** (list): atom position in primitive cell (cartesian) [in multipie order].
  - **wannier_to_multipie** (list): converting index from wannier to multipie.
  - **multipie_to_wannier** (list): converting index from multipie to wannier.
  - **lattice_transformation** (list): integer matrix U, x_model = x_wannier U + t (with model only).
  - **origin_shift** (list): origin shift t (with model only).

In symcw mode, and in wannier mode with `samb/model`, the Wannier functions are mapped onto the kets of the model by their projection centres and orbitals, and H(R) in `seedname_hr.dat` is re-indexed into the primitive cell of the model through bond vectors.
The lattice vectors in `seedname.win` may be any primitive cell of the model lattice in the same Cartesian frame (e.g., `ibrav=7` of Quantum ESPRESSO for a body-centred lattice), and the origin may be shifted.
Otherwise, the fractional coordinates in `seedname.win` are used as those of the model primitive cell with a warning.

In wannier and symcw modes, H(R) in `seedname_hr.dat` is divided by the Wigner-Seitz degeneracy `ndegen(R)`, as in the Fourier interpolation of Wannier90.
If `seedname_wsvec.dat` exists (`use_ws_distance`), each matrix element is further distributed over its shortest images R+T with the weight 1/nT.

In wannier mode, the model is optional.
Without `samb/model`, the name is `wannier/seedname`, and the lattice, the atom positions and the default k path are those of `seedname.win`.
The kets are named after the atoms on which the projection centres lie; a projection centre not on an atom (e.g., a bond centre) is named X1, X2, ....
With `samb/model`, the model is loaded, and the lattice and the default k path are those of the model, so that the dispersion can be compared directly with symcw mode.

The fractional coordinates of `output/dispersion/k_point` refer to the reciprocal cell of `info/B`, i.e., the primitive cell of the model in samb and symcw modes and in wannier mode with model, and the cell of `seedname.win` in wannier mode without model.
The same `k_point` means different k points in these cells if they differ (e.g., `ibrav=7` of Quantum ESPRESSO for a body-centred lattice).

- **output** (dict): output of physical quantities.
  - **dispersion** (dict): dispersion related.
    - **k_path** (str): k path.
    - **k_point** (dict): definition of k points.
    - **e_max** (float): max. of energy.
    - **e_min** (float): min. of energy.

## ModelAnalyzer Class

```{eval-rst}
.. automodule:: multipie.core.model_analyzer
```
