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
  - **win** (dict): input file for wannier90.x in seedname.win file.
  - **nnkp** (dict): information needed to determine the required overlap elements Mmn(k,b) and projections A_{mn}(k).
  - **HR** (dict): Wannier Hamiltonian.
  - **amn** (dict): overlap matrix elements in seedname.amn file, A_{mn}(k) = <ψ^{KS}_{m}(k)|φ_{n}(k)> (ψ^{KS}_{m}(k) = u^{KS}_{m}(k) e^{ik.r}).
  - **eig** (dict): Kohn-Sham energies in seedname.eig file, E_{m}(k).
  - **mmn** (dict): overlap matrix elements in seedname.mmn file, M_{mn}(k,b) = <u^{KS}_{m}(k)|u^{KS}_{n}(k+b)> (ψ^{KS}_{m}(k) = u^{KS}_{m}(k) e^{ik.r}).
  - **umat** (dict): unitary matrix elements in seedname_u.mat (Uopt(k)) and seedname_u_dis.mat (Udis(k)) files, U(k) = Uopt(k)@Udis(k).
  - **spn** (dict): matrix elements of Pauli spin operator in seedname.spn file, <ψ^{KS}_{m}(k)|\sigma_x,y,z|ψ^{KS}_{n}(k)>.

- **output** (dict): output of physical quantities.
  - **dispersion** (dict): dispersion related.
    - **k_path** (str): k path.
    - **k_point** (dict): definition of k points.

## ModelAnalyzer Class

```{eval-rst}
.. automodule:: multipie.core.model_analyzer
```
