# Atomic multipole operators

A Japanese version is [here](atomic_multipole_table_ja.md).

This page lists the atomic multipole operators $X_{lm}^{(s)}(k)$ up to rank $l = 2$, given by
`Group.atomic_multipole_operator` (see {ref}`Operator expression <atomic-multipole-operator>`), following H. Kusunose,
R. Oiwa, and S. Hayami, J. Phys. Soc. Jpn. **89**, 104704 (2020).

- **Components.** The real (tesseral) components are listed: $X_{l0}$, and
  $C_{lm} = [(-1)^m X_{lm} + X_{l,-m}]/\sqrt{2}$ and $S_{lm} = [(-1)^m X_{lm} - X_{l,-m}]/(\sqrt{2}\,i)$ for $m > 0$.
  They are named as the orbitals: $(x, y, z) = (C_{11}, S_{11}, C_{10})$ for $l = 1$, and
  $(u, v, yz, xz, xy) = (C_{20}, C_{22}, S_{21}, C_{21}, S_{22})$ for $l = 2$. E.g., $Q_{u} = z^2 - (x^2 + y^2)/2$.
  The components are hermitian, as $X_{l,-m} = (-1)^m X_{lm}^\dagger$ (Eqs. 6 and 14).
- **Operators.** The position $x, y, z$, the orbital angular momentum $l_x, l_y, l_z$ (dimensionless), and the Pauli
  matrices $\sigma_x, \sigma_y, \sigma_z$ (the spin is $\boldsymbol{\sigma}/2$), with $[l_i, l_j] = i\epsilon_{ijk} l_k$,
  $[\sigma_i, \sigma_j] = 2i\epsilon_{ijk}\sigma_k$, and $[l_i, x_j] = i\epsilon_{ijk} x_k$; $x_i$ commute with each other,
  and $\sigma_i$ commute with $l_i$ and $x_i$.
- **Order of the operators.** The order of the definition is kept: in a product such as $l_x\, y$, $l_x$ acts after
  $y$. Only a run of consecutive position operators is sorted, and the Pauli matrices are put at the right end; the
  commutation relations are not used otherwise, so that, e.g., $l_x y$ and $y l_x$ appear as separate terms.
- **Notation.** $X_{\gamma}$ for $s = 0$ and $X^{(1)}_{\gamma}(k)$ for $s = 1$, with the name of the component
  $\gamma$ (e.g., $Q_{u}$ and $Q^{(1)}_{x}(0)$), and $X_{0}$ and $X^{(1)}_{0}(k)$ for rank 0. The matrix elements in an
  orbital basis are those of the atomic SAMB before the symmetry adaptation and the normalization, with the radial
  factor $\langle r^k \rangle = 1$.
- **Rank 1 in vectors.** The three components of rank 1 are written together as a vector, e.g.,
  $\boldsymbol{T} = (T_{x}, T_{y}, T_{z})$ and $\boldsymbol{Q}^{(1)}(0) = (Q^{(1)}_{x}(0), Q^{(1)}_{y}(0), Q^{(1)}_{z}(0))$,
  with $\boldsymbol{r} = (x, y, z)$, $\boldsymbol{l} = (l_x, l_y, l_z)$ and $\boldsymbol{\sigma} = (\sigma_x, \sigma_y, \sigma_z)$.
  The order of the operators in the products is kept also here, e.g., $(\boldsymbol{r} \times \boldsymbol{l})_x = y l_z - z l_y$
  and $(\boldsymbol{l} \times \boldsymbol{r})_x = l_y z - l_z y$. The scalars of rank 0 are written in the same way. The
  terms with $\boldsymbol{r} \cdot \boldsymbol{l}$ or $\boldsymbol{l} \cdot \boldsymbol{r}$ are omitted, as they vanish
  ($\boldsymbol{r} \cdot \boldsymbol{l} = \boldsymbol{l} \cdot \boldsymbol{r} = 0$ for
  $\boldsymbol{l} = -i\,\boldsymbol{r} \times \nabla$); `Group.atomic_multipole_operator` gives them, as it does not use
  the commutation relations.

The vanishing ones (e.g., $M_{0}$ and $\boldsymbol{G}$ of $s = 0$) are not listed.

```{include} examples/atomic_multipole/table_en.md
```
