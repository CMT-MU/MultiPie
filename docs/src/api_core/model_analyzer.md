# ModelAnalyzer

## Dict data

This class has the following dict data.

- General info.
  - **grid**: (list) grid size, N1, N2, N3.
  - **A** (ndarray): translational vectors for primitive cell, [a1,a2,a3] (3x3).
  - **B** (ndarray): translational vectors for primitive reciprocal cell, [b1,b2,b3] (3x3).
  - **unit_cell_volume** (float): volume of primitive cell.

- **samb** (dict): SAMB related.
  - **parameter** (dict): used parameter for each SAMB.
  - **matrix_info** (dict): matrix information.

- **wannier** (dict): Wannier related.
  - **

- **output** (dict): output of physical quantities.
  - **dispersion** (dict): dispersion related.
    - **k_path** (str): k path.
    - **k_point** (dict): definition of k points.

## ModelAnalyzer Class

```{eval-rst}
.. automodule:: multipie.core.model_analyzer
```
