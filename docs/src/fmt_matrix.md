# SAMB in full matrix

The output file format for SAMB in full-matrix form is the following. The keywords with * is only for crystal.

## SAMB in full matrix form
- **model** : model name.
- **molecule** : molecule or crystal ?
- **group** : (tag, detailed str)
- **dimension** : dimension of full matrix
- **ket** : ket basis list, orbital@site
- **version** : MultiPie version
- **k_point*** : representative k points
- **k_path*** : high-symmetry line in k space
- **cell_site** : { name_idx(pset): (position, SOs) }
- **A*** : transform matrix, [a1,a2,a3]
- **matrix** : { "z_#": "matrix" }

## SAMB in full matrix form with fourier transform
- **model** : model name.
- **molecule** : molecule or crystal ?
- **group** : (tag, detailed str)
- **dimension** : dimension of full matrix
- **ket** : ket basis list, orbital@site
- **version** : MultiPie version
- **k_point** : representative k points
- **k_path** : high-symmetry line in k space
- **A** : transform matrix, [a1,a2,a3]
- **bond*** : { "bond_#": "vector" }
- **matrix** : { "z_#": "matrix" }
