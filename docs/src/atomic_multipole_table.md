# Atomic multipole operators

A Japanese version is [here](atomic_multipole_table_ja.md).

This page lists the atomic multipole operators $X_{lm}^{(s)}(k)$ up to rank $l = 2$, given by
`Group.atomic_multipole_operator` (see {ref}`Operator expression <atomic-multipole-operator>`), following H. Kusunose,
R. Oiwa, and S. Hayami, J. Phys. Soc. Jpn. **89**, 104704 (2020).

- **Components.** The real (tesseral) components are listed: $X_{l0}$, and
  $C_{lm} = [(-1)^m X_{lm} + X_{l,-m}]/\sqrt{2}$ and $S_{lm} = [(-1)^m X_{lm} - X_{l,-m}]/(\sqrt{2}\,i)$ for $m > 0$.
  They are named as the orbitals: $(x, y, z) = (C_{11}, S_{11}, C_{10})$ for $l = 1$, and
  $(u, v, yz, xz, xy) = (C_{20}, C_{22}, S_{21}, C_{21}, S_{22})$ for $l = 2$. E.g., $Q_{2,u} = z^2 - (x^2 + y^2)/2$.
  The components are hermitian, as $X_{l,-m} = (-1)^m X_{lm}^\dagger$ (Eqs. 6 and 14).
- **Operators.** The position $x, y, z$, the orbital angular momentum $l_x, l_y, l_z$ (dimensionless), and the Pauli
  matrices $\sigma_x, \sigma_y, \sigma_z$ (the spin is $\boldsymbol{\sigma}/2$), with $[l_i, l_j] = i\epsilon_{ijk} l_k$,
  $[\sigma_i, \sigma_j] = 2i\epsilon_{ijk}\sigma_k$, and $[l_i, x_j] = i\epsilon_{ijk} x_k$; $x_i$ commute with each other,
  and $\sigma_i$ commute with $l_i$ and $x_i$.
- **Order of the operators.** The order of the definition is kept: in a product such as $l_x\, y$, $l_x$ acts after
  $y$. Only a run of consecutive position operators is sorted, and the Pauli matrices are put at the right end; the
  commutation relations are not used otherwise, so that, e.g., $l_x y$ and $y l_x$ appear as separate terms.
- **Notation.** $X_{l,\gamma}$ for $s = 0$ and $X^{(1)}_{l,\gamma}(k)$ for $s = 1$. The matrix elements in an orbital
  basis are those of the atomic SAMB before the symmetry adaptation and the normalization, with the radial factor
  $\langle r^k \rangle = 1$.

The vanishing ones (e.g., $M_{0}$ and $G_{1}$ of $s = 0$) are not listed. The list is created by `docs/src/examples/atomic_multipole/create_table.py`.

```{include} examples/atomic_multipole/table_en.md
```
