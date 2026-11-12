# Examples of materials

This page builds tight-binding models of three materials with MultiPie, step by step: SrVO$_3$ (a cubic perovskite),
monolayer MoS$_2$ (a hexagonal layer without inversion), and Te (a chiral crystal with the spin-orbit coupling). They
cover the main ingredients of a model: several sites with different orbitals, bonds between different sites up to
further neighbors, a point group whose irreps carry primes (D$_{\rm 3h}$), and a spinful model in a chiral group.

The same model inputs are used by the tests (`tests/test_materials.py`), which check the properties described at the
end of each section, see [Checks by symmetry](#checks-by-symmetry). The outputs on this page are created by
`docs/src/examples/materials/create_examples.py`.

The weights $z_j$ used for the band dispersions below are chosen only for illustration; they are not fitted to any
first-principles calculation. To build a realistic model, fit them to a calculation, e.g., with the `symcw` mode of
`mp_analyze` (see [Post processing](post_processing.md)).

## Overview of the procedure

For each material:

1. Write a model input file, `<name>_in.py`. It gives the space group, the cell, the sites with their orbitals, and
   the bonds by neighbor. `"SAMB_select": {"Gamma": "IR"}` selects the SAMB of the identity irrep, which are the
   ones allowed in a Hamiltonian.
2. Create the model,

    ```bash
    mp_create SrVO3_in.py
    ```

    which writes `SrVO3/SrVO3.pkl` (the model) and `SrVO3/SrVO3.pdf` (a summary with all the SAMB and their
    definitions).
3. Choose the weights $z_j$ of the SAMB $Z_j$ in a control file, `<name>_ctrl.py`. The Hamiltonian is

    $$
    H = \sum_j z_j Z_j,
    $$

    where each $Z_j$ is a hermitian matrix normalized to 1, so that $z_j$ has the unit of energy.
4. Analyze the model,

    ```bash
    mp_analyze SrVO3_ctrl.py
    ```

    which writes the band dispersion along the k path to `SrVO3/output/`.

In the tables of SAMB below, "site/bond" gives the sites of the bra and the ket, "Wyckoff" the Wyckoff position of the
site or the bond (`12b@12j` is a bond Wyckoff position `12b` with the bond centers at `12j`), "neighbor" the order of
the neighbor of the bond (`-` for an on-site term), and "bra-ket" the orbital ranks. The symbol is the combined SAMB,
$\mathbb{X}_{l}^{(c)}(\Gamma)$, with the type `X`, the rank `l`, and the irrep $\Gamma$; how to read it is explained in
the [User Guide](guide.md). The expression of each $Z_j$ in terms of the atomic and the cluster SAMB is given in the
PDF of the model.

## SrVO$_3$: cubic perovskite

SrVO$_3$ has the cubic perovskite structure (space group $Pm\bar{3}m$, No. 221, $a = 3.84$ Å). The conduction bands
come from the V $t_{2g}$ orbitals ($d_{yz}$, $d_{xz}$, $d_{xy}$), which hybridize with the O $p$ orbitals. We take the
V site at the body center (Wyckoff `1b`) and the three O sites at the face centers (`3c`), and omit Sr.

```{literalinclude} ../../tests/data/materials/SrVO3_in.py
```

The model has 10 SAMB ([SrVO3.pdf](examples/materials/SrVO3.pdf)):

```{include} examples/materials/SrVO3_samb.md
```

The number of SAMB of each kind follows from the symmetry:

- V on-site ($z_3$): the $t_{2g}$ orbitals form the irrep $T_{2g}$ of the site symmetry O$_{\rm h}$, so one energy.
- O on-site ($z_1$, $z_2$): the site symmetry of O is D$_{\rm 4h}$, where the $p$ orbitals split into the one along the
  V-O bond and the two perpendicular to it, so two energies.
- V-O hopping ($z_8$): a $t_{2g}$ orbital has no $\sigma$ bond with $p$, so only the $pd\pi$ hopping remains.
- V-V hopping ($z_9$, $z_{10}$): $dd\pi$ and $dd\delta$ between $t_{2g}$ orbitals.
- O-O hopping ($z_4$-$z_7$): four parameters, as the O-O bond has a low symmetry.

The control gives the O on-site energy, the V-O hopping, and an O-O hopping:

```{literalinclude} examples/materials/SrVO3_ctrl.py
```

```{image} figures/materials/SrVO3_dispersion.png
:width: 600px
```

The three $t_{2g}$ bands lie above the nine O $p$ bands. With the $pd\pi$ hopping alone, a $t_{2g}$ orbital hops only in
its own plane, e.g., $d_{yz}$ only within the $yz$ plane, so one of the $t_{2g}$ bands is flat along $\Gamma$-X. At
$\Gamma$ and R, the $t_{2g}$ bands are triply degenerate ($T_{2g}$), and so are the O $p$ bands at $\Gamma$ ($2T_{1u} + T_{2u}$).

## MoS$_2$: monolayer of 2H-MoS$_2$

A monolayer of 2H-MoS$_2$ has the space group $P\bar{6}m2$ (No. 187, point group D$_{\rm 3h}$), without the inversion.
The Mo site is at `1a`, and the two S sites above and below it are at `2h`, $(1/3, 2/3, \pm z)$. A large $c$ (20 Å)
separates the layers. We take the five Mo $d$ orbitals and the three S $p$ orbitals.

```{literalinclude} ../../tests/data/materials/MoS2_in.py
```

The model has 30 SAMB ([MoS2.pdf](examples/materials/MoS2.pdf)):

```{include} examples/materials/MoS2_samb.md
```

The irreps of D$_{\rm 3h}$ are written with primes, e.g., $A_1'$ (the identity irrep) and $E''$. From Ver.3, the irreps
are in the order $A_1', A_2', E', A_1'', A_2'', E''$ (see [Changes from Ver.2](changes.md)). The site symmetry splits
the Mo $d$ orbitals into $A_1'$ ($d_{z^2}$), $E'$ ($d_{x^2-y^2}$, $d_{xy}$), and $E''$ ($d_{xz}$, $d_{yz}$), giving
the three on-site energies $z_1$-$z_3$; the S $p$ orbitals split into $p_z$ and ($p_x$, $p_y$), giving $z_4$ and $z_5$.

```{literalinclude} examples/materials/MoS2_ctrl.py
```

```{image} figures/materials/MoS2_dispersion.png
:width: 600px
```

At $\Gamma$, the bands are non-degenerate ($A_1'$, $A_2''$) or doubly degenerate ($E'$, $E''$): three single bands
and four pairs for the eleven bands.

## Te: a chiral crystal with the spin-orbit coupling

Trigonal Te has the chiral space group $P3_121$ (No. 152, point group D$_3$): the Te atoms form helical chains along
$c$, and the crystal has no mirror and no inversion. The three Te sites are at `3a`, $(x, 0, 1/3)$ with $x = 0.2636$.
We take the $p$ orbitals with the spin (`"spinful": True`), and the bonds within a chain (first neighbor) and between
the chains (second neighbor).

```{literalinclude} ../../tests/data/materials/Te_in.py
```

The model has 66 SAMB ([Te.pdf](examples/materials/Te.pdf)):

:::{admonition} SAMB of Te (66)
:class: dropdown

```{include} examples/materials/Te_samb.md
```
:::

Two features show up:

- **Spin-dependent SAMB.** In a spinful model, the atomic SAMB include the spin, written with the superscript
  $(s, k)$ with $s = 1$. For example, $z_6$, $\mathbb{Q}_{0}^{(1,1;c)}(A_1)$, is the on-site spin-orbit coupling,
  proportional to $\boldsymbol{l} \cdot \boldsymbol{s}$.
- **Electric toroidal (`G`) SAMB.** As D$_3$ has only rotations, the axial (`G`) SAMB of the identity irrep are
  allowed in the Hamiltonian together with the polar (`Q`) ones, e.g., $\mathbb{G}_{0}^{(c)}(A_1)$ ($z_{19}$, $z_{20}$),
  the electric toroidal monopole of the bonds, which represents the chirality.

The control gives the spin-orbit coupling and the hoppings within and between the chains, and computes the
expectation value of $S_z$ along the k path (`"local": ["Sz"]`), shown in color:

```{literalinclude} examples/materials/Te_ctrl.py
```

```{image} figures/materials/Te_dispersion.png
:width: 600px
```

Without the inversion, the spin-orbit coupling splits the bands except at the time-reversal invariant momenta, e.g.,
$\Gamma$ and A, where the bands are doubly degenerate (Kramers). Along $\Gamma$-A (the chain direction, $k_z$), each
pair of split bands is polarized along $z$ with opposite signs (red and blue), i.e., the spin is parallel or
anti-parallel to $\boldsymbol{k}$: the radial spin texture of a chiral crystal.

(checks-by-symmetry)=
## Checks by symmetry

`tests/test_materials.py` creates the three models and gives random weights to all the SAMB. It checks properties
that hold for any complete set of SAMB, independent of their numbering, signs, and basis:

- the number of SAMB, and that all of them are `Q` or `G` of the identity irrep;
- $H(\boldsymbol{k})$ is hermitian;
- $E_n(\boldsymbol{k}) = E_n(\boldsymbol{k} W)$ for every rotation $W$ of the space group, and
  $E_n(\boldsymbol{k}) = E_n(-\boldsymbol{k})$ by the time reversal;
- the degeneracies at $\Gamma$ follow the irreps: SrVO$_3$ four triplets, MoS$_2$ three singlets and four doublets, Te
  nine Kramers pairs;
- no degeneracy remains at a generic $\boldsymbol{k}$.

They hold also when the group database changes, as long as the SAMB span the same space of Hamiltonians.
