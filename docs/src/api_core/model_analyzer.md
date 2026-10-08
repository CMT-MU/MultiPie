# ModelAnalyzer

## Dict data

This class has the following dict data.

- **info** (dict): General info.
  - **mode** (str): analysis mode, samb, wannier, or symcw.
  - **grid** (list): grid size, N1, N2, N3.
  - **name** (str): model and analysis name.
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
  - **lattice_transformation** (list): integer matrix U, x_model = x_wannier U + t (symcw only).
  - **origin_shift** (list): origin shift t (symcw only).

In symcw mode, the Wannier functions are mapped onto the kets of the model by their projection centres and orbitals, and H(R) in `seedname_hr.dat` is re-indexed into the primitive cell of the model through bond vectors.
The lattice vectors in `seedname.win` may be any primitive cell of the model lattice in the same Cartesian frame (e.g., `ibrav=7` of Quantum ESPRESSO for a body-centred lattice), and the origin may be shifted.
Otherwise, the fractional coordinates in `seedname.win` are used as those of the model primitive cell with a warning.

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
