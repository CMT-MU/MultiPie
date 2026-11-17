---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Using the Group class

`Group` gives the data of a point group (PG), a space group (SG), a magnetic point group (MPG), or a magnetic space group
(MSG) from the group database. This page shows how to use its properties and methods, with the outputs of MultiPie Ver.3.
The [User Guide](guide.md) explains what the quantities mean. A Japanese version is [here](group/group_ja.md).

Most examples use the point group D$_{\rm 4h}$, the space group P4/mmm (No. 123), and the magnetic point groups of
D$_{\rm 4h}$.

| topic | properties and methods |
|:--|:--|
| [Creating a group](#group-create) | `Group(tag)`, `str()`, `name()`, `latex()`, `ID`, `group_type`, `info`, `is_point_group`, `is_magnetic_group`, `is_hexagonal_subgroup`, `global_info()`, `clear_cache()` |
| [Symmetry operations](#group-symmetry-operation) | `symmetry_operation`, `tag_symmetry_operation()` |
| [Characters and representation matrices](#group-character) | `character`, `representation_matrix`, `tag_irrep()` |
| [Wyckoff positions](#group-wyckoff) | `wyckoff`, `find_wyckoff()`, `find_wyckoff_site()`, `find_wyckoff_bond()`, `create_cell_site()`, `create_cell_bond()` |
| [Harmonics](#group-harmonics) | `harmonics`, `multipolar_harmonics`, `tag_multipole()` |
| [Atomic basis and atomic SAMB](#group-atomic) | `atomic_basis()`, `tag_atomic_basis()`, `tag_atomic_basis_proj()`, `atomic_samb()`, `transform_atomic_samb()`, `create_atomic_samb_L()` |
| [Cluster and combined SAMB](#group-cluster) | `cluster_samb()`, `combined_samb()`, `cg()`, `multipole_cluster_samb()`, `combined_object()` |
| [Transformation by the symmetry operations](#group-transform) | `transform_vector()`, `transform_multipole()`, `create_cell_vector()`, `create_cell_multipole()` |
| [Magnetic groups and response tensors](#group-magnetic) | `active_multipole`, `response_tensor_all()`, `response_tensor()` |

(group-create)=
## 1. Creating a group

A group is given by its tag:

| group | tag | examples |
|:--|:--|:--|
| PG | `PG:1`-`PG:47`, or the Schoenflies symbol | `"PG:15"`, `"D4h"` |
| SG | `SG:1`-`SG:230`, the number, or the Schoenflies symbol | `"SG:123"`, `123`, `"D4h^1"` |
| MPG | `MPG:` + the id `PG.no.ID`, or the id | `"MPG:15.7.59"`, `"15.7.59"` |
| MSG | `MSG:` + the BNS number, or the BNS number | `"MSG:123.345"`, `"123.345"` |

The tables of all groups, with their tags, are in
[Group Information](https://cmt-mu.github.io/MultiPieGen/src/group_info.html).

```{code-cell} ipython3
from multipie import Group
g = Group("D4h")
str(g), g.group_type, g.ID
```

```{code-cell} ipython3
g.name(), g.name(detail=True)
```

```{code-cell} ipython3
g.latex()
```

```{code-cell} ipython3
g.latex(detail=True)
```

```{code-cell} ipython3
Group("PG:15").name(), Group(123).name(), Group("SG:123").name()
```

```{code-cell} ipython3
Group("15.7.59").name(detail=True)
```

```{code-cell} ipython3
Group("123.345").name(detail=True)
```

```{code-cell} ipython3
g.is_point_group, g.is_magnetic_group, g.is_hexagonal_subgroup
```

- `name()` and `latex()` give the name in text and in LaTeX; `detail=True` adds the number, the international symbol,
  and the crystal system.
- `is_point_group` is true for PG and MPG, `is_magnetic_group` for MPG and MSG, and `is_hexagonal_subgroup` for the
  subgroups of D$_{\rm 6h}$ (trigonal and hexagonal), whose harmonics use the hexagonal basis.

`info` is a named tuple of the group information. It gives the related groups (the point group, the space group, the
gray MPG, and the MSG of the same family), and `tp` of `name()` and `latex()` names one of them:

```{code-cell} ipython3
g.info.SG, g.info.MPG, g.info.MSG
```

```{code-cell} ipython3
g.name(tp="SG"), g.latex(tp="SG")
```

`Group.global_info()` is the information common to all groups, e.g., the tables of the tags and the ids:

```{code-cell} ipython3
gi = Group.global_info()
list(gi.keys())
```

```{code-cell} ipython3
gi["id"]["D4h"], gi["tag"]["PG:15"]
```

```{code-cell} ipython3
[(t, Group(t).name()) for t in gi["tag"] if t.startswith("MPG:15.")][:4]
```

The data of a group are loaded once and kept in memory; `Group.clear_cache()` releases them.

(group-symmetry-operation)=
## 2. Symmetry operations

`symmetry_operation` is a dict of the symmetry operations:

```{code-cell} ipython3
so = g.symmetry_operation
list(so.keys())
```

```{code-cell} ipython3
so["tag"][:8]
```

```{code-cell} ipython3
so["generator"]
```

```{code-cell} ipython3
so["cartesian"][2]  # 4+[001].
```

```{code-cell} ipython3
so["det"][8:12]
```

```{code-cell} ipython3
so["product"][("4+[001]", "4+[001]")]
```

| key | content |
|:--|:--|
| `tag` | the operations, e.g., `4+[001]` (rotation by $2\pi/4$ about [001]), `-1` (inversion), `m[100]` (mirror) |
| `generator` | the generators |
| `cartesian`, `fractional` | the matrices in the cartesian and fractional coordinates (`(4, 4)` with the translation for SG and MSG) |
| `det` | the determinants (+1 for proper, -1 for improper operations) |
| `product` | the product table, `{(a, b): ab}` (PG) |
| `plus_set`, `cartesian_primitive`, `fractional_primitive` | the centering translations and the matrices of the primitive cell (SG) |
| `tr_sign` | the sign of the time reversal (MPG and MSG): -1 for the operations with the time reversal, e.g., `-1'` |

`Group.tag_symmetry_operation(tag)` gives the properties of an operation:

```{code-cell} ipython3
Group.tag_symmetry_operation("4+[001]")
```

```{code-cell} ipython3
Group.tag_symmetry_operation("-4+[001]", latex=True)
```

(group-character)=
## 3. Characters and representation matrices

`character` gives the character table and the products of the irreps (of the point group):

```{code-cell} ipython3
ch = g.character
list(ch.keys())
```

```{code-cell} ipython3
ch["conjugacy"][:4]
```

```{code-cell} ipython3
ch["table"]["Eu"]  # for each conjugacy class.
```

```{code-cell} ipython3
ch["dimension"]["Eu"], ch["polar_axial_conversion"]["Eu"]
```

```{code-cell} ipython3
ch["symmetric_product"][("Eu", "Eu")]
```

```{code-cell} ipython3
ch["anti_symmetric_product"]["Eu"]
```

- `table_full` gives the characters for each operation instead of each class.
- `polar_axial_conversion` gives the irrep of an axial quantity of the same function, e.g., $(x, y)$ is $E_u$ as a
  polar vector and $E_g$ as an axial vector.
- `symmetric_product` and `anti_symmetric_product` give the decomposition of the symmetric and anti-symmetric products,
  as `(multiplicity, irrep)`.

`representation_matrix` gives the matrices of the irreps for each operation, in the order of `rm["tag"]`:

```{code-cell} ipython3
rm = g.representation_matrix
rm["tag"][2], rm["matrix"]["Eu"][2]
```

```{code-cell} ipython3
Group.tag_irrep("Eu")
```

```{code-cell} ipython3
Group.tag_irrep("Eu", latex=True)
```

(group-wyckoff)=
## 4. Wyckoff positions

`wyckoff` gives the Wyckoff sites and, for PG and SG, the Wyckoff bonds. A bond position `2b@1o` is the bond
position `2b` (2 bonds) whose centers are at the site position `1o`.

```{code-cell} ipython3
w = g.wyckoff
list(w["site"].keys())
```

```{code-cell} ipython3
w["site"]["4b"]["symmetry"], w["site"]["4b"]["conventional"]
```

```{code-cell} ipython3
w["site"]["4b"]["bond"][:3]  # bond positions centered at 4b.
```

```{code-cell} ipython3
list(w["bond"].keys())[:4]
```

```{code-cell} ipython3
w["bond"]["2b@1o"]["conventional"]  # (vector, center).
```

`find_wyckoff_site(site)` and `find_wyckoff_bond(bond)` find the Wyckoff position of a site or a bond and give its
equivalent sites or bonds; `find_wyckoff()` calls either of them. A site is written `"[x,y,z]"`, and a bond
`"[tail];[head]"`, `"[vector]@[center]"`, or `"[start]:[vector]"`. A bond is given as (vector, center).

```{code-cell} ipython3
g.find_wyckoff_site("[1,1,0]")
```

```{code-cell} ipython3
g.find_wyckoff_bond("[0,0,0];[1,0,0]")
```

```{code-cell} ipython3
g.find_wyckoff("[1,0,0]@[0,0,0]")[0]
```

For a space group, the coordinates are fractional. `create_cell_site()` and `create_cell_bond()` give the sites and the
bonds in the unit cell, the mapping from the symmetry operations (index from 1), and the Wyckoff position:

```{code-cell} ipython3
s = Group(123)
s.find_wyckoff_site("[1/2,0,0]")[0]
```

```{code-cell} ipython3
sites, mapping, wp = s.create_cell_site("[1/2,0,0]")
sites, wp
```

```{code-cell} ipython3
bonds, mapping, wp = s.create_cell_bond("[0,0,0];[1/2,0,0]")
bonds, wp
```

(group-harmonics)=
## 5. Harmonics

`harmonics` gives the point-group harmonics (rank $\le 11$) as a `Dict` keyed by the multipole key
(see {ref}`multipole-key`). Each value is the polynomials of the components and the unitary matrix $U$ from the
tesseral harmonics. A `Dict` selects its entries by the fields of the key with `select()`, and `named_keys()` gives the
keys as named tuples:

```{code-cell} ipython3
h = g.harmonics
len(h)
```

```{code-cell} ipython3
list(h.select(X="Q", l=2).keys())
```

```{code-cell} ipython3
ex, U = h[("Q", 2, "B1g", -1, -1, 0, 0, "q")]
ex, U.T
```

```{code-cell} ipython3
[(key.X, key.Gamma) for key in h.select(l=1).named_keys()]
```

`Group.tag_multipole(key, comp=None, latex=False, superscript="")` gives the symbol of a multipole; `comp` is the
component from 0, and all the components are given when it is `None`:

```{code-cell} ipython3
Group.tag_multipole(("Q", 2, "Eg", -1, -1, 0, 0, "q"))
```

```{code-cell} ipython3
Group.tag_multipole(("Q", 2, "Eg", -1, -1, 0, 0, "q"), 1, latex=True, superscript="a")
```

`multipolar_harmonics` gives the multipolar harmonics with an internal degree of freedom (internal rank `s` = 1-3,
rank $\le 3$), as `(coefficient, expression)`; see {ref}`multipolar-harmonics`. E.g., the electric toroidal monopole
$\boldsymbol{r} \cdot \boldsymbol{G}$ made of an axial vector $\boldsymbol{G}$:

```{code-cell} ipython3
mh = g.multipolar_harmonics
coeff, ex = mh[("G", 0, "A1u", -1, -1, 1, 1, "g")]
ex
```

(group-atomic)=
## 6. Atomic basis and atomic SAMB

`atomic_basis(basis_type)` gives the orbitals of each rank for `"lg"` (real orbitals), `"lgs"` (with spin), and `"jml"`
(total angular momentum); see {ref}`atomic-orbitals`. `tag_atomic_basis()` gives the properties or the LaTeX of an
orbital, and `tag_atomic_basis_proj()` the LaTeX of the orbitals selected by a mask:

```{code-cell} ipython3
g.atomic_basis("lg")[2]
```

```{code-cell} ipython3
g.atomic_basis("lgs")[1][:4]
```

```{code-cell} ipython3
g.atomic_basis("jml")[1]
```

```{code-cell} ipython3
Group.tag_atomic_basis("dxy", 2)
```

```{code-cell} ipython3
Group.tag_atomic_basis("(3/2,1/2)", 1, latex=True)
```

```{code-cell} ipython3
g.tag_atomic_basis_proj("lg", 2, mask=["dxz", "dyz"])
```

`atomic_samb(basis_type, (L1, L2), mask)` gives the atomic SAMB of the block of the bra rank `L1` and the ket rank `L2`,
as `(matrices, expression)`; see {ref}`atomic-samb`. The mask `(bra, ket)` selects the orbitals by the names or the
indices, or gives a unitary matrix to another basis; the SAMB are then orthonormalized in the selected orbitals:

```{code-cell} ipython3
a = g.atomic_samb("lg", (2, 2))  # all d orbitals.
len(a)
```

```{code-cell} ipython3
a2 = g.atomic_samb("lg", (2, 2), (["dxz", "dyz"], ["dxz", "dyz"]))
[(key.X, key.l, key.Gamma) for key in a2.named_keys()]
```

```{code-cell} ipython3
mat, ex = a2[("M", 1, "A2g", -1, -1, 0, 0, "q")]  # l_z in (dxz, dyz).
mat, ex
```

```{code-cell} ipython3
len(g.atomic_samb("lg", (1, 2))), len(g.atomic_samb("lgs", (1, 1)))  # p-d block, and p with spin.
```

- `transform_atomic_samb(samb, U_ket, U_bra)` transforms the atomic SAMB to another basis, given by the unitary
  matrices from the tesseral basis, or by a rank for the cubic or hexagonal basis of the group.
- `create_atomic_samb_L(L)` creates the atomic multipoles of the $|L, M\rangle$ basis (rank $L$) from the definition.

(group-cluster)=
## 7. Cluster and combined SAMB

`cluster_samb(wp)` gives the site cluster SAMB of a site position, and `cluster_samb(wp, "bond")` the bond cluster
SAMB of a bond position, as `(values on the sites or bonds, expression)`; see {ref}`cluster-samb`:

```{code-cell} ipython3
c = g.cluster_samb("4b")  # sites (x,x,0), (-x,-x,0), (-x,x,0), (x,-x,0).
for key, (v, ex) in c.items(): print(key[0], key[1], key[2], v.tolist(), ex.tolist())
```

```{code-cell} ipython3
b = g.cluster_samb("2b@1o", "bond")  # bonds along x and y.
[(key.X, key.l, key.Gamma) for key in b.named_keys()]
```

`combined_samb(atomic keys, cluster keys)` combines atomic SAMB and cluster SAMB into the irreps with the
Clebsch-Gordan coefficients, `cg(key1, key2, key)`; see {ref}`combined-samb`. Each component of a combined SAMB is the
list of `(coefficient, atomic key, atomic component, cluster key, cluster component)`. The keywords `X`, `l`, `Gamma`,
and `s` select the combined SAMB:

```{code-cell} ipython3
a = g.atomic_samb("lg", (1, 1), (["px", "py"], ["px", "py"]))
z = g.combined_samb(a.named_keys(), b.named_keys())
len(z)
```

```{code-cell} ipython3
z[("Q", 0, "A1g", -1, 2, 0, 0, "q")]
```

```{code-cell} ipython3
[(key.X, key.l, key.Gamma, key.p) for key in g.combined_samb(a.named_keys(), b.named_keys(), Gamma="A1g").named_keys()]
```

```{code-cell} ipython3
g.cg(("Q", 1, "Eu", -1, -1, 0, 0, "q"), ("Q", 1, "Eu", -1, -1, 0, 0, "q"), ("Q", 2, "B1g", -1, -1, 0, 0, "q"))[:, :, 0]
```

The first $A_{1g}$ SAMB ($p = 1$) is the uniform $p_x$-$p_x$ plus $p_y$-$p_y$ hopping on all the bonds. The second
($p = 2$) is the product of the $p_x$-$p_x$ minus $p_y$-$p_y$ orbital part ($B_{1g}$) and the bond part with opposite
signs on the $x$ and $y$ bonds ($B_{1g}$): the difference of the $\sigma$-type ($p_x$-$p_x$ on the $x$ bonds) and the
$\pi$-type ($p_y$-$p_y$ on the $x$ bonds) hoppings.

`multipole_cluster_samb(X, l, site)` makes the combined SAMB of the multipoles `(X, l)` on the cluster of a site or a
bond, and `combined_object(wp, X, samb)` gives a component of it as the polynomials on the sites. E.g., electric
dipoles (`"Q", 1`) on the four sites `4c` ($(\pm 1, 0, 0)$, $(0, \pm 1, 0)$):

```{code-cell} ipython3
samb, wp, sites = g.multipole_cluster_samb("Q", 1, "[1,0,0]")
wp, sites.tolist()
```

```{code-cell} ipython3
[(key.X, key.l, key.Gamma) for key in samb.named_keys()][:4]
```

```{code-cell} ipython3
g.combined_object(wp, "Q", samb[("Q", 0, "A1g", -1, -1, 0, 0, "q")][0][0])
```

```{code-cell} ipython3
g.combined_object(wp, "G", samb[("G", 1, "A2g", -1, -1, 0, 0, "q")][0][0])
```

The dipoles point outward in $\mathbb{Q}_{0}(A_{1g})$, and around the $z$ axis in $\mathbb{G}_{1}(A_{2g})$, the
electric toroidal dipole along $z$.

(group-transform)=
## 8. Transformation by the symmetry operations

`transform_vector(vector, X)` and `transform_multipole(multipole, X)` apply all the symmetry operations, in the order
of `symmetry_operation["tag"]`. `X` gives the type: `"Q"` polar, `"G"` axial (with $\det R$), `"T"` polar and `"M"` axial
with the sign of the time reversal (`tr_sign` of a magnetic group):

```{code-cell} ipython3
so["tag"][:4]
```

```{code-cell} ipython3
g.transform_vector("[1,0,0]")[:4]
```

```{code-cell} ipython3
m = Group("15.7.59")  # 4/m'm'm'.
m.symmetry_operation["tag"][8], m.symmetry_operation["tr_sign"][8]
```

```{code-cell} ipython3
m.transform_vector("[0,0,1]", X="M")[[0, 8]]  # a magnetic moment is reversed by -1'.
```

```{code-cell} ipython3
g.transform_multipole("x*y")[:4]
```

`create_cell_vector("vector#site", X)` and `create_cell_multipole("multipole#site", X)` place the transformed vectors
or multipoles on the equivalent sites (averaged over the operations that give the same site):

```{code-cell} ipython3
v, sites, mapping, wp = g.create_cell_vector("[1,0,0]#[1,0,0]")
v, wp
```

```{code-cell} ipython3
mp, sites, mapping, wp = g.create_cell_multipole("x**2-y**2#[1,0,0]")
mp
```

(group-magnetic)=
## 9. Magnetic groups and response tensors

`active_multipole` lists the multipoles allowed in the magnetic point group (of a PG, its gray group), and
`response_tensor_all(X)` gives all the response tensors of type `X` up to rank 4, as `(tensor, {component: multipoles})`;
see {ref}`response-tensors`. In 4/m'm'm', the magnetoelectric tensor (`"M"`, rank 2) is allowed:

```{code-cell} ipython3
m = Group("15.7.59")
m.name(detail=True)
```

```{code-cell} ipython3
m.active_multipole
```

```{code-cell} ipython3
r = m.response_tensor_all("M")
list(r.keys())[:4]
```

```{code-cell} ipython3
tensor, relation = r[("M", 2, "s")]
tensor
```

```{code-cell} ipython3
relation
```

`response_tensor(X, rank, comp, opt)` gives one component (indices 1, 2, 3 for $x, y, z$), e.g., $M_{zz}$ above. In the
ferromagnetic
4/mm'm', the magnetization $M_z$ and the anomalous Hall effect (`"T"`, rank 2, anti-symmetric) are both carried by the
magnetic dipole $\mathbb{M}_z$:

```{code-cell} ipython3
m.response_tensor("M", 2, (3, 3), "s")
```

```{code-cell} ipython3
f = Group("15.6.58")
f.name(), f.active_multipole
```

```{code-cell} ipython3
f.response_tensor("M", 1, 3), f.response_tensor("T", 2, (1, 2), "a")
```
