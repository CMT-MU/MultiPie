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

```python
>>> from multipie import Group
>>> g = Group("D4h")
>>> str(g), g.group_type, g.ID
('D4h', 'PG', '15')
>>> g.name(), g.name(detail=True)
('D4h', 'PG#15:  D4h  (4/mmm)  [ tetragonal ]')
>>> g.latex()
'D_{\\rm 4h}'
>>> g.latex(detail=True)
'PG No. 15\\quad$D_{\\rm 4h}$\\quad$4/mmm$\\quad[ tetragonal ]'
>>> Group("PG:15").name(), Group(123).name(), Group("SG:123").name()
('D4h', 'D4h^1', 'D4h^1')
>>> Group("15.7.59").name(detail=True)
"MPG#15.7.59:  4/m'm'm'  [ Type III, tetragonal ]"
>>> Group("123.345").name(detail=True)
"MSG#123.345:  P4/mm'm'    [ Type III, tetragonal ]"
>>> g.is_point_group, g.is_magnetic_group, g.is_hexagonal_subgroup
(True, False, False)
```

- `name()` and `latex()` give the name in text and in LaTeX; `detail=True` adds the number, the international symbol,
  and the crystal system.
- `is_point_group` is true for PG and MPG, `is_magnetic_group` for MPG and MSG, and `is_hexagonal_subgroup` for the
  subgroups of D$_{\rm 6h}$ (trigonal and hexagonal), whose harmonics use the hexagonal basis.

`info` is a named tuple of the group information. It gives the related groups (the point group, the space group, the
gray MPG, and the MSG of the same family), and `tp` of `name()` and `latex()` names one of them:

```python
>>> g.info.SG, g.info.MPG, g.info.MSG
('SG:123', 'MPG:15.2.54', 'MSG:123.340')
>>> g.name(tp="SG"), g.latex(tp="SG")
('D4h^1', 'D_{\\rm 4h}^{1}')
```

`Group.global_info()` is the information common to all groups, e.g., the tables of the tags and the ids:

```python
>>> gi = Group.global_info()
>>> list(gi.keys())
['id_set', 'tag', 'id', 'group_name', 'character', 'root_cluster', 'harmonics', 'response_tensor']
>>> gi["id"]["D4h"], gi["tag"]["PG:15"]
('PG:15', 'D4h')
>>> [(t, Group(t).name()) for t in gi["tag"] if t.startswith("MPG:15.")][:4]
[('MPG:15.1.53', '4/mmm'), ('MPG:15.2.54', "4/mmm1'"), ('MPG:15.3.55', "4/m'mm"), ('MPG:15.4.56', "4'/mm'm")]
```

The data of a group are loaded once and kept in memory; `Group.clear_cache()` releases them.

(group-symmetry-operation)=
## 2. Symmetry operations

`symmetry_operation` is a dict of the symmetry operations:

```python
>>> so = g.symmetry_operation
>>> list(so.keys())
['tag', 'generator', 'cartesian', 'fractional', 'det', 'product', 'so_matrix', 'mapping']
>>> so["tag"][:8]
['1', '2[001]', '4+[001]', '4-[001]', '2[010]', '2[100]', '2[110]', '2[1-10]']
>>> so["generator"]
['2[001]', '4+[001]', '2[010]', '-1']
>>> so["cartesian"][2]  # 4+[001].
array([[0, -1, 0],
       [1, 0, 0],
       [0, 0, 1]], dtype=object)
>>> so["det"][8:12]
array([-1, -1, -1, -1], dtype=object)
>>> so["product"][("4+[001]", "4+[001]")]
'2[001]'
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

```python
>>> Group.tag_symmetry_operation("4+[001]")
{'mirror': False, 'inversion': False, 'n': 4, 'axis': '[0,0,1]', 'det': 1, 'partial_trans': '[0,0,0]', 'zero_trans': True, 'point_group': True, 'time_reversal': False}
>>> Group.tag_symmetry_operation("-4+[001]", latex=True)
'-4^{+}_{\\,\\,001}'
```

(group-character)=
## 3. Characters and representation matrices

`character` gives the character table and the products of the irreps (of the point group):

```python
>>> ch = g.character
>>> list(ch.keys())
['conjugacy', 'dimension', 'table', 'table_full', 'polar_axial_conversion', 'symmetric_product', 'anti_symmetric_product']
>>> ch["conjugacy"][:4]
[['1'], ['2[001]'], ['2[100]', '2[010]'], ['2[110]', '2[1-10]']]
>>> ch["table"]["Eu"]  # for each conjugacy class.
array([2, -2, 0, 0, 0, -2, 2, 0, 0, 0], dtype=object)
>>> ch["dimension"]["Eu"], ch["polar_axial_conversion"]["Eu"]
(2, 'Eg')
>>> ch["symmetric_product"][("Eu", "Eu")]
[(1, 'A1g'), (1, 'B1g'), (1, 'B2g')]
>>> ch["anti_symmetric_product"]["Eu"]
[(1, 'A2g')]
```

- `table_full` gives the characters for each operation instead of each class.
- `polar_axial_conversion` gives the irrep of an axial quantity of the same function, e.g., $(x, y)$ is $E_u$ as a
  polar vector and $E_g$ as an axial vector.
- `symmetric_product` and `anti_symmetric_product` give the decomposition of the symmetric and anti-symmetric products,
  as `(multiplicity, irrep)`.

`representation_matrix` gives the matrices of the irreps for each operation, in the order of `rm["tag"]`:

```python
>>> rm = g.representation_matrix
>>> rm["tag"][2], rm["matrix"]["Eu"][2]
('4+[001]', array([[0, -1],
       [1, 0]], dtype=object))
>>> Group.tag_irrep("Eu")
{'tag': 'Eu', 'dimension': 2, 'parity': 'odd', 'complex': '', 'head': 'E', 'subscript': 'u', 'superscript': ''}
>>> Group.tag_irrep("Eu", latex=True)
'E_{u}'
```

(group-wyckoff)=
## 4. Wyckoff positions

`wyckoff` gives the Wyckoff sites and, for PG and SG, the Wyckoff bonds. A bond position `2b@1o` is the bond
position `2b` (2 bonds) whose centers are at the site position `1o`.

```python
>>> w = g.wyckoff
>>> list(w["site"].keys())
['1o', '2a', '4b', '4c', '8d', '8e', '8f', '16g']
>>> w["site"]["4b"]["symmetry"], w["site"]["4b"]["conventional"]
('m.m2', array([[x, x, 0],
       [-x, -x, 0],
       [-x, x, 0],
       [x, -x, 0]], dtype=object))
>>> w["site"]["4b"]["bond"][:3]  # bond positions centered at 4b.
['4a@4b', '4b@4b', '4c@4b']
>>> list(w["bond"].keys())[:4]
['1a@1o', '2b@1o', '2c@1o', '4d@1o']
>>> w["bond"]["2b@1o"]["conventional"]  # (vector, center).
array([[X, 0, 0, 0, 0, 0],
       [0, X, 0, 0, 0, 0]], dtype=object)
```

`find_wyckoff_site(site)` and `find_wyckoff_bond(bond)` find the Wyckoff position of a site or a bond and give its
equivalent sites or bonds; `find_wyckoff()` calls either of them. A site is written `"[x,y,z]"`, and a bond
`"[tail];[head]"`, `"[vector]@[center]"`, or `"[start]:[vector]"`. A bond is given as (vector, center).

```python
>>> g.find_wyckoff_site("[1,1,0]")
('4b', array([[ 1.,  1.,  0.],
       [-1., -1.,  0.],
       [-1.,  1.,  0.],
       [ 1., -1.,  0.]]))
>>> g.find_wyckoff_bond("[0,0,0];[1,0,0]")
('4a@4c', array([[ 1. ,  0. ,  0. ,  0.5,  0. ,  0. ],
       [-1. ,  0. ,  0. , -0.5,  0. ,  0. ],
       [ 0. ,  1. ,  0. ,  0. ,  0.5,  0. ],
       [ 0. , -1. ,  0. ,  0. , -0.5,  0. ]]))
>>> g.find_wyckoff("[1,0,0]@[0,0,0]")[0]
'2b@1o'
```

For a space group, the coordinates are fractional. `create_cell_site()` and `create_cell_bond()` give the sites and the
bonds in the unit cell, the mapping from the symmetry operations (index from 1), and the Wyckoff position:

```python
>>> s = Group(123)
>>> s.find_wyckoff_site("[1/2,0,0]")[0]
'2f'
>>> sites, mapping, wp = s.create_cell_site("[1/2,0,0]")
>>> sites, wp
(array([[0.5, 0. , 0. ],
       [0. , 0.5, 0. ]]), '2f')
>>> bonds, mapping, wp = s.create_cell_bond("[0,0,0];[1/2,0,0]")
>>> bonds, wp
(array([[ 0.5 ,  0.  ,  0.  ,  0.25,  0.  ,  0.  ],
       [-0.5 ,  0.  ,  0.  ,  0.75,  0.  ,  0.  ],
       [ 0.  ,  0.5 ,  0.  ,  0.  ,  0.25,  0.  ],
       [ 0.  , -0.5 ,  0.  ,  0.  ,  0.75,  0.  ]]), '4a@4l')
```

(group-harmonics)=
## 5. Harmonics

`harmonics` gives the point-group harmonics (rank $\le 11$) as a `Dict` keyed by the multipole key
(see {ref}`multipole-key`). Each value is the polynomials of the components and the unitary matrix $U$ from the
tesseral harmonics. A `Dict` selects its entries by the fields of the key with `select()`, and `named_keys()` gives the
keys as named tuples:

```python
>>> h = g.harmonics
>>> len(h)
216
>>> list(h.select(X="Q", l=2).keys())
[('Q', 2, 'A1g', -1, -1, 0, 0, 'q'), ('Q', 2, 'B1g', -1, -1, 0, 0, 'q'), ('Q', 2, 'B2g', -1, -1, 0, 0, 'q'), ('Q', 2, 'Eg', -1, -1, 0, 0, 'q')]
>>> ex, U = h[("Q", 2, "B1g", -1, -1, 0, 0, "q")]
>>> ex, U.T
(array([sqrt(3)*x**2/2 - sqrt(3)*y**2/2], dtype=object), array([[sqrt(2)/2, 0, 0, 0, sqrt(2)/2]], dtype=object))
>>> [(key.X, key.Gamma) for key in h.select(l=1).named_keys()]
[('G', 'A2g'), ('G', 'Eg'), ('Q', 'A2u'), ('Q', 'Eu')]
```

`Group.tag_multipole(key, comp=None, latex=False, superscript="")` gives the symbol of a multipole; `comp` is the
component from 0, and all the components are given when it is `None`:

```python
>>> Group.tag_multipole(("Q", 2, "Eg", -1, -1, 0, 0, "q"))
['Q(2,1)(Eg)', 'Q(2,2)(Eg)']
>>> Group.tag_multipole(("Q", 2, "Eg", -1, -1, 0, 0, "q"), 1, latex=True, superscript="a")
'\\mathbb{Q}_{2,2}^{(a)}(E_{g})'
```

`multipolar_harmonics` gives the multipolar harmonics with an internal degree of freedom (internal rank `s` = 1-3,
rank $\le 3$), as `(coefficient, expression)`; see {ref}`multipolar-harmonics`. E.g., the electric toroidal monopole
$\boldsymbol{r} \cdot \boldsymbol{G}$ made of an axial vector $\boldsymbol{G}$:

```python
>>> mh = g.multipolar_harmonics
>>> coeff, ex = mh[("G", 0, "A1u", -1, -1, 1, 1, "g")]
>>> ex
array([sqrt(3)*G_{x}*x/3 + sqrt(3)*G_{y}*y/3 + sqrt(3)*G_{z}*z/3],
      dtype=object)
```

(group-atomic)=
## 6. Atomic basis and atomic SAMB

`atomic_basis(basis_type)` gives the orbitals of each rank for `"lg"` (real orbitals), `"lgs"` (with spin), and `"jml"`
(total angular momentum); see {ref}`atomic-orbitals`. `tag_atomic_basis()` gives the properties or the LaTeX of an
orbital, and `tag_atomic_basis_proj()` the LaTeX of the orbitals selected by a mask:

```python
>>> g.atomic_basis("lg")[2]
['dv', 'dxy', 'dxz', 'dyz', 'du']
>>> g.atomic_basis("lgs")[1][:4]
['(px,u)', '(px,d)', '(py,u)', '(py,d)']
>>> g.atomic_basis("jml")[1]
['(1/2,1/2)', '(1/2,-1/2)', '(3/2,3/2)', '(3/2,1/2)', '(3/2,-1/2)', '(3/2,-3/2)']
>>> Group.tag_atomic_basis("dxy", 2)
{'spinless': True, 'spherical': False, 'orbital': 'xy', 'spin': None, 'rank': 2}
>>> Group.tag_atomic_basis("(3/2,1/2)", 1, latex=True)
'\\ket{\\frac{3}{2},\\frac{1}{2};p}'
>>> g.tag_atomic_basis_proj("lg", 2, mask=["dxz", "dyz"])
[\ket{d_{xz}}, \ket{d_{yz}}]
```

`atomic_samb(basis_type, (L1, L2), mask)` gives the atomic SAMB of the block of the bra rank `L1` and the ket rank `L2`,
as `(matrices, expression)`; see {ref}`atomic-samb`. The mask `(bra, ket)` selects the orbitals by the names or the
indices, or gives a unitary matrix to another basis; the SAMB are then orthonormalized in the selected orbitals:

```python
>>> a = g.atomic_samb("lg", (2, 2))  # all d orbitals.
>>> len(a)
19
>>> a2 = g.atomic_samb("lg", (2, 2), (["dxz", "dyz"], ["dxz", "dyz"]))
>>> [(key.X, key.l, key.Gamma) for key in a2.named_keys()]
[('Q', 0, 'A1g'), ('Q', 2, 'B1g'), ('Q', 2, 'B2g'), ('M', 1, 'A2g')]
>>> mat, ex = a2[("M", 1, "A2g", -1, -1, 0, 0, "q")]  # l_z in (dxz, dyz).
>>> mat, ex
(array([[[0, -sqrt(2)*I/2],
        [sqrt(2)*I/2, 0]]], dtype=object), array([z], dtype=object))
>>> len(g.atomic_samb("lg", (1, 2))), len(g.atomic_samb("lgs", (1, 1)))  # p-d block, and p with spin.
(22, 27)
```

- `transform_atomic_samb(samb, U_ket, U_bra)` transforms the atomic SAMB to another basis, given by the unitary
  matrices from the tesseral basis, or by a rank for the cubic or hexagonal basis of the group.
- `create_atomic_samb_L(L)` creates the atomic multipoles of the $|L, M\rangle$ basis (rank $L$) from the definition.

(group-cluster)=
## 7. Cluster and combined SAMB

`cluster_samb(wp)` gives the site cluster SAMB of a site position, and `cluster_samb(wp, "bond")` the bond cluster
SAMB of a bond position, as `(values on the sites or bonds, expression)`; see {ref}`cluster-samb`:

```python
>>> c = g.cluster_samb("4b")  # sites (x,x,0), (-x,-x,0), (-x,x,0), (x,-x,0).
>>> for key, (v, ex) in c.items(): print(key[0], key[1], key[2], v.tolist(), ex.tolist())
Q 0 A1g [[1/2, 1/2, 1/2, 1/2]] [1]
Q 2 B2g [[1/2, 1/2, -1/2, -1/2]] [sqrt(3)*x*y]
Q 1 Eu [[1/2, -1/2, -1/2, 1/2], [1/2, -1/2, 1/2, -1/2]] [x, y]
>>> b = g.cluster_samb("2b@1o", "bond")  # bonds along x and y.
>>> [(key.X, key.l, key.Gamma) for key in b.named_keys()]
[('Q', 0, 'A1g'), ('Q', 2, 'B1g'), ('T', 1, 'Eu')]
```

`combined_samb(atomic keys, cluster keys)` combines atomic SAMB and cluster SAMB into the irreps with the
Clebsch-Gordan coefficients, `cg(key1, key2, key)`; see {ref}`combined-samb`. Each component of a combined SAMB is the
list of `(coefficient, atomic key, atomic component, cluster key, cluster component)`. The keywords `X`, `l`, `Gamma`,
and `s` select the combined SAMB:

```python
>>> a = g.atomic_samb("lg", (1, 1), (["px", "py"], ["px", "py"]))
>>> z = g.combined_samb(a.named_keys(), b.named_keys())
>>> len(z)
12
>>> z[("Q", 0, "A1g", -1, 2, 0, 0, "q")]
([[(1, ('Q', 2, 'B1g', -1, -1, 0, 0, 'q'), 0, ('Q', 2, 'B1g', -1, -1, 0, 0, 'q'), 0)]], [1])
>>> [(key.X, key.l, key.Gamma, key.p) for key in g.combined_samb(a.named_keys(), b.named_keys(), Gamma="A1g").named_keys()]
[('Q', 0, 'A1g', 1), ('Q', 0, 'A1g', 2)]
>>> g.cg(("Q", 1, "Eu", -1, -1, 0, 0, "q"), ("Q", 1, "Eu", -1, -1, 0, 0, "q"), ("Q", 2, "B1g", -1, -1, 0, 0, "q"))[:, :, 0]
array([[sqrt(2)/2, 0],
       [0, -sqrt(2)/2]], dtype=object)
```

The first $A_{1g}$ SAMB ($p = 1$) is the uniform $p_x$-$p_x$ plus $p_y$-$p_y$ hopping on all the bonds. The second
($p = 2$) is the product of the $p_x$-$p_x$ minus $p_y$-$p_y$ orbital part ($B_{1g}$) and the bond part with opposite
signs on the $x$ and $y$ bonds ($B_{1g}$): the difference of the $\sigma$-type ($p_x$-$p_x$ on the $x$ bonds) and the
$\pi$-type ($p_y$-$p_y$ on the $x$ bonds) hoppings.

`multipole_cluster_samb(X, l, site)` makes the combined SAMB of the multipoles `(X, l)` on the cluster of a site or a
bond, and `combined_object(wp, X, samb)` gives a component of it as the polynomials on the sites. E.g., electric
dipoles (`"Q", 1`) on the four sites `4c` ($(\pm 1, 0, 0)$, $(0, \pm 1, 0)$):

```python
>>> samb, wp, sites = g.multipole_cluster_samb("Q", 1, "[1,0,0]")
>>> wp, sites.tolist()
('4c', [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]])
>>> [(key.X, key.l, key.Gamma) for key in samb.named_keys()][:4]
[('Q', 0, 'A1g'), ('G', 1, 'A2g'), ('Q', 1, 'A2u'), ('Q', 2, 'B1g')]
>>> g.combined_object(wp, "Q", samb[("Q", 0, "A1g", -1, -1, 0, 0, "q")][0][0])
array([x/2, -x/2, y/2, -y/2], dtype=object)
>>> g.combined_object(wp, "G", samb[("G", 1, "A2g", -1, -1, 0, 0, "q")][0][0])
array([-y/2, y/2, x/2, -x/2], dtype=object)
```

The dipoles point outward in $\mathbb{Q}_{0}(A_{1g})$, and around the $z$ axis in $\mathbb{G}_{1}(A_{2g})$, the
electric toroidal dipole along $z$.

(group-transform)=
## 8. Transformation by the symmetry operations

`transform_vector(vector, X)` and `transform_multipole(multipole, X)` apply all the symmetry operations, in the order
of `symmetry_operation["tag"]`. `X` gives the type: `"Q"` polar, `"G"` axial (with $\det R$), `"T"` polar and `"M"` axial
with the sign of the time reversal (`tr_sign` of a magnetic group):

```python
>>> so["tag"][:4]
['1', '2[001]', '4+[001]', '4-[001]']
>>> g.transform_vector("[1,0,0]")[:4]
array([[1, 0, 0],
       [-1, 0, 0],
       [0, 1, 0],
       [0, -1, 0]], dtype=object)
>>> m = Group("15.7.59")  # 4/m'm'm'.
>>> m.symmetry_operation["tag"][8], m.symmetry_operation["tr_sign"][8]
("-1'", -1)
>>> m.transform_vector("[0,0,1]", X="M")[[0, 8]]  # a magnetic moment is reversed by -1'.
array([[0, 0, 1],
       [0, 0, -1]], dtype=object)
>>> g.transform_multipole("x*y")[:4]
array([x*y, x*y, -x*y, -x*y], dtype=object)
```

`create_cell_vector("vector#site", X)` and `create_cell_multipole("multipole#site", X)` place the transformed vectors
or multipoles on the equivalent sites (averaged over the operations that give the same site):

```python
>>> v, sites, mapping, wp = g.create_cell_vector("[1,0,0]#[1,0,0]")
>>> v, wp
(array([[1, 0, 0],
       [-1, 0, 0],
       [0, 1, 0],
       [0, -1, 0]], dtype=object), '4c')
>>> mp, sites, mapping, wp = g.create_cell_multipole("x**2-y**2#[1,0,0]")
>>> mp
[x**2 - y**2, x**2 - y**2, -x**2 + y**2, -x**2 + y**2]
```

(group-magnetic)=
## 9. Magnetic groups and response tensors

`active_multipole` lists the multipoles allowed in the magnetic point group (of a PG, its gray group), and
`response_tensor_all(X)` gives all the response tensors of type `X` up to rank 4, as `(tensor, {component: multipoles})`;
see {ref}`response-tensors`. In 4/m'm'm', the magnetoelectric tensor (`"M"`, rank 2) is allowed:

```python
>>> m = Group("15.7.59")
>>> m.name(detail=True)
"MPG#15.7.59:  4/m'm'm'  [ Type III, tetragonal ]"
>>> m.active_multipole
['Qs', 'Qdu', 'Qg', 'Qgu', 'Ms', 'Mdu', 'Mg', 'Mgu']
>>> r = m.response_tensor_all("M")
>>> list(r.keys())[:4]
[('M', 0, ''), ('M', 1, ''), ('M', 2, 's'), ('M', 2, 'a')]
>>> tensor, relation = r[("M", 2, "s")]
>>> tensor
array([[M_{xx}, 0, 0],
       [0, M_{xx}, 0],
       [0, 0, M_{zz}]], dtype=object)
>>> relation
{M_{xx}: \mathbb{M}_{0}^{(1)} - \mathbb{M}_{u}^{(1)}, M_{zz}: \mathbb{M}_{0}^{(1)} + 2*\mathbb{M}_{u}^{(1)}}
```

`response_tensor(X, rank, comp, opt)` gives one component (indices 1, 2, 3 for $x, y, z$), e.g., $M_{zz}$ above. In the
ferromagnetic
4/mm'm', the magnetization $M_z$ and the anomalous Hall effect (`"T"`, rank 2, anti-symmetric) are both carried by the
magnetic dipole $\mathbb{M}_z$:

```python
>>> m.response_tensor("M", 2, (3, 3), "s")
\mathbb{M}_{0}^{(1)} + 2*\mathbb{M}_{u}^{(1)}
>>> f = Group("15.6.58")
>>> f.name(), f.active_multipole
("4/mm'm'", ['Qs', 'Qdu', 'Qg', 'Qgu', 'Tgaz', 'Mpz', 'Mfaz'])
>>> f.response_tensor("M", 1, 3), f.response_tensor("T", 2, (1, 2), "a")
(\mathbb{M}_{z}^{(1)}, \mathbb{M}_{z}^{(1)})
```
