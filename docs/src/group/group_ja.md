---
orphan: true
---

# Group クラスの使い方

`Group` は、群のデータベースから、点群（PG）、空間群（SG）、磁気点群（MPG）、磁気空間群（MSG）のデータを与えます。このページでは、そのプロパティとメソッドの使い方を、MultiPie Ver.3 の出力とともに示します。量の意味は[ユーザーガイド](../guide/users_guide_ja.md)で説明しています。英語版は[こちら](../group.md)です。

例の多くは、点群 D$_{\rm 4h}$、空間群 P4/mmm（No. 123）、D$_{\rm 4h}$ の磁気点群を使います。

| 項目 | プロパティとメソッド |
|:--|:--|
| [群の作成](#ja-group-create) | `Group(tag)`, `str()`, `name()`, `latex()`, `ID`, `group_type`, `info`, `is_point_group`, `is_magnetic_group`, `is_hexagonal_subgroup`, `global_info()`, `clear_cache()` |
| [対称操作](#ja-group-symmetry-operation) | `symmetry_operation`, `tag_symmetry_operation()` |
| [指標と表現行列](#ja-group-character) | `character`, `representation_matrix`, `tag_irrep()` |
| [Wyckoff 位置](#ja-group-wyckoff) | `wyckoff`, `find_wyckoff()`, `find_wyckoff_site()`, `find_wyckoff_bond()`, `create_cell_site()`, `create_cell_bond()` |
| [調和関数](#ja-group-harmonics) | `harmonics`, `multipolar_harmonics`, `tag_multipole()` |
| [原子基底と原子 SAMB](#ja-group-atomic) | `atomic_basis()`, `tag_atomic_basis()`, `tag_atomic_basis_proj()`, `atomic_samb()`, `transform_atomic_samb()`, `create_atomic_samb_L()` |
| [クラスター SAMB と結合 SAMB](#ja-group-cluster) | `cluster_samb()`, `combined_samb()`, `cg()`, `multipole_cluster_samb()`, `combined_object()` |
| [対称操作による変換](#ja-group-transform) | `transform_vector()`, `transform_multipole()`, `create_cell_vector()`, `create_cell_multipole()` |
| [磁気群と応答テンソル](#ja-group-magnetic) | `active_multipole`, `response_tensor_all()`, `response_tensor()` |

(ja-group-create)=
## 1. 群の作成

群はタグで指定します。

| 群 | タグ | 例 |
|:--|:--|:--|
| PG | `PG:1`〜`PG:47`、または Schoenflies 記号 | `"PG:15"`, `"D4h"` |
| SG | `SG:1`〜`SG:230`、番号、または Schoenflies 記号 | `"SG:123"`, `123`, `"D4h^1"` |
| MPG | `MPG:` ＋ ID `PG.no.ID`、または ID | `"MPG:15.7.59"`, `"15.7.59"` |
| MSG | `MSG:` ＋ BNS 番号、または BNS 番号 | `"MSG:123.345"`, `"123.345"` |

すべての群の表とそのタグは[群の情報](https://cmt-mu.github.io/MultiPieGen/src/group_info.html)にあります。

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

- `name()` と `latex()` は群の名前をテキストと LaTeX で与えます。`detail=True` とすると、番号、国際記号、晶系が加わります。
- `is_point_group` は PG と MPG で、`is_magnetic_group` は MPG と MSG で真になります。`is_hexagonal_subgroup` は D$_{\rm 6h}$ の部分群（三方晶と六方晶）で真になり、その調和関数は六方晶の基底を使います。

`info` は群の情報の named tuple です。関連する群（点群、空間群、灰色 MPG、同じ系列の MSG）を与え、`name()` と `latex()` の `tp` でそのどれかを指定できます。

```python
>>> g.info.SG, g.info.MPG, g.info.MSG
('SG:123', 'MPG:15.2.54', 'MSG:123.340')
>>> g.name(tp="SG"), g.latex(tp="SG")
('D4h^1', 'D_{\\rm 4h}^{1}')
```

`Group.global_info()` はすべての群に共通の情報で、例えばタグと ID の対応表を含みます。

```python
>>> gi = Group.global_info()
>>> list(gi.keys())
['id_set', 'tag', 'id', 'group_name', 'character', 'root_cluster', 'harmonics', 'response_tensor']
>>> gi["id"]["D4h"], gi["tag"]["PG:15"]
('PG:15', 'D4h')
>>> [(t, Group(t).name()) for t in gi["tag"] if t.startswith("MPG:15.")][:4]
[('MPG:15.1.53', '4/mmm'), ('MPG:15.2.54', "4/mmm1'"), ('MPG:15.3.55', "4/m'mm"), ('MPG:15.4.56', "4'/mm'm")]
```

群のデータは一度だけ読み込まれてメモリに保持されます。`Group.clear_cache()` で解放できます。

(ja-group-symmetry-operation)=
## 2. 対称操作

`symmetry_operation` は対称操作の dict です。

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

| キー | 内容 |
|:--|:--|
| `tag` | 対称操作。例えば `4+[001]`（[001] 軸まわりの $2\pi/4$ 回転）、`-1`（空間反転）、`m[100]`（鏡映） |
| `generator` | 生成元 |
| `cartesian`, `fractional` | デカルト座標と分率座標での行列（SG と MSG では並進を含む `(4, 4)`） |
| `det` | 行列式（真性回転は +1、非真性回転は -1） |
| `product` | 積の表 `{(a, b): ab}`（PG） |
| `plus_set`, `cartesian_primitive`, `fractional_primitive` | 格子の中心化の並進と、基本単位胞での行列（SG） |
| `tr_sign` | 時間反転の符号（MPG と MSG）。`-1'` のように時間反転を含む操作で -1 |

`Group.tag_symmetry_operation(tag)` は対称操作の性質を与えます。

```python
>>> Group.tag_symmetry_operation("4+[001]")
{'mirror': False, 'inversion': False, 'n': 4, 'axis': '[0,0,1]', 'det': 1, 'partial_trans': '[0,0,0]', 'zero_trans': True, 'point_group': True, 'time_reversal': False}
>>> Group.tag_symmetry_operation("-4+[001]", latex=True)
'-4^{+}_{\\,\\,001}'
```

(ja-group-character)=
## 3. 指標と表現行列

`character` は（点群の）指標表と既約表現の積を与えます。

```python
>>> ch = g.character
>>> list(ch.keys())
['conjugacy', 'dimension', 'table', 'table_full', 'polar_axial_conversion', 'symmetric_product', 'anti_symmetric_product']
>>> ch["conjugacy"][:4]
[['1'], ['2[001]'], ['2[100]', '2[010]'], ['2[110]', '2[1-10]']]
>>> ch["table"]["Eu"]  # 共役類ごと。
array([2, -2, 0, 0, 0, -2, 2, 0, 0, 0], dtype=object)
>>> ch["dimension"]["Eu"], ch["polar_axial_conversion"]["Eu"]
(2, 'Eg')
>>> ch["symmetric_product"][("Eu", "Eu")]
[(1, 'A1g'), (1, 'B1g'), (1, 'B2g')]
>>> ch["anti_symmetric_product"]["Eu"]
[(1, 'A2g')]
```

- `table_full` は共役類ごとではなく、対称操作ごとの指標を与えます。
- `polar_axial_conversion` は、同じ関数の軸性量の既約表現を与えます。例えば $(x, y)$ は、極性ベクトルとしては $E_u$、軸性ベクトルとしては $E_g$ です。
- `symmetric_product` と `anti_symmetric_product` は、対称積と反対称積の分解を `(重複度, 既約表現)` で与えます。

`representation_matrix` は、各対称操作に対する既約表現の行列を `rm["tag"]` の順に与えます。

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

(ja-group-wyckoff)=
## 4. Wyckoff 位置

`wyckoff` は Wyckoff サイトと、PG と SG では Wyckoff ボンドを与えます。ボンドの位置 `2b@1o` は、中心がサイトの位置 `1o` にあるボンドの位置 `2b`（2 本のボンド）です。

```python
>>> w = g.wyckoff
>>> list(w["site"].keys())
['1o', '2a', '4b', '4c', '8d', '8e', '8f', '16g']
>>> w["site"]["4b"]["symmetry"], w["site"]["4b"]["conventional"]
('m.m2', array([[x, x, 0],
       [-x, -x, 0],
       [-x, x, 0],
       [x, -x, 0]], dtype=object))
>>> w["site"]["4b"]["bond"][:3]  # 中心が 4b にあるボンドの位置。
['4a@4b', '4b@4b', '4c@4b']
>>> list(w["bond"].keys())[:4]
['1a@1o', '2b@1o', '2c@1o', '4d@1o']
>>> w["bond"]["2b@1o"]["conventional"]  # (ベクトル, 中心)。
array([[X, 0, 0, 0, 0, 0],
       [0, X, 0, 0, 0, 0]], dtype=object)
```

`find_wyckoff_site(site)` と `find_wyckoff_bond(bond)` は、サイトやボンドの Wyckoff 位置を探し、等価なサイトやボンドを与えます。`find_wyckoff()` はそのどちらかを呼びます。サイトは `"[x,y,z]"`、ボンドは `"[始点];[終点]"`、`"[ベクトル]@[中心]"`、`"[始点]:[ベクトル]"` のいずれかで書きます。ボンドは (ベクトル, 中心) で与えられます。

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

空間群では座標は分率座標です。`create_cell_site()` と `create_cell_bond()` は、単位胞内のサイトとボンド、対称操作からの対応（1 始まりの番号）、Wyckoff 位置を与えます。

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

(ja-group-harmonics)=
## 5. 調和関数

`harmonics` は点群調和関数（階数 $\le 11$）を、多極子のキー（{ref}`ja-multipole-key`）をキーとする `Dict` で与えます。値は各成分の多項式と、テッセラル調和関数からのユニタリー行列 $U$ です。`Dict` は `select()` でキーのフィールドによって要素を選び、`named_keys()` でキーを named tuple として与えます。

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

`Group.tag_multipole(key, comp=None, latex=False, superscript="")` は多極子の記号を与えます。`comp` は 0 から数える成分で、`None` のときはすべての成分を与えます。

```python
>>> Group.tag_multipole(("Q", 2, "Eg", -1, -1, 0, 0, "q"))
['Q(2,1)(Eg)', 'Q(2,2)(Eg)']
>>> Group.tag_multipole(("Q", 2, "Eg", -1, -1, 0, 0, "q"), 1, latex=True, superscript="a")
'\\mathbb{Q}_{2,2}^{(a)}(E_{g})'
```

`multipolar_harmonics` は、内部自由度をもつ点群多極子調和関数（内部階数 `s` = 1〜3、階数 $\le 3$）を `(係数, 式)` で与えます（{ref}`ja-multipolar-harmonics`）。例えば、軸性ベクトル $\boldsymbol{G}$ から作る電気トロイダル単極子 $\boldsymbol{r} \cdot \boldsymbol{G}$ は次のとおりです。

```python
>>> mh = g.multipolar_harmonics
>>> coeff, ex = mh[("G", 0, "A1u", -1, -1, 1, 1, "g")]
>>> ex
array([sqrt(3)*G_{x}*x/3 + sqrt(3)*G_{y}*y/3 + sqrt(3)*G_{z}*z/3],
      dtype=object)
```

(ja-group-atomic)=
## 6. 原子基底と原子 SAMB

`atomic_basis(basis_type)` は、`"lg"`（実軌道）、`"lgs"`（スピン付き）、`"jml"`（全角運動量）について、階数ごとの軌道を与えます（{ref}`ja-atomic-orbitals`）。`tag_atomic_basis()` は軌道の性質や LaTeX を、`tag_atomic_basis_proj()` はマスクで選んだ軌道の LaTeX を与えます。

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

`atomic_samb(basis_type, (L1, L2), mask)` は、bra の階数 `L1`、ket の階数 `L2` のブロックの原子 SAMB を `(行列, 式)` で与えます（{ref}`ja-atomic-samb`）。マスク `(bra, ket)` は、軌道を名前か番号で選ぶか、別の基底へのユニタリー行列を与えます。このとき SAMB は選んだ軌道の中で正規直交化されます。

```python
>>> a = g.atomic_samb("lg", (2, 2))  # すべての d 軌道。
>>> len(a)
19
>>> a2 = g.atomic_samb("lg", (2, 2), (["dxz", "dyz"], ["dxz", "dyz"]))
>>> [(key.X, key.l, key.Gamma) for key in a2.named_keys()]
[('Q', 0, 'A1g'), ('Q', 2, 'B1g'), ('Q', 2, 'B2g'), ('M', 1, 'A2g')]
>>> mat, ex = a2[("M", 1, "A2g", -1, -1, 0, 0, "q")]  # (dxz, dyz) での l_z。
>>> mat, ex
(array([[[0, -sqrt(2)*I/2],
        [sqrt(2)*I/2, 0]]], dtype=object), array([z], dtype=object))
>>> len(g.atomic_samb("lg", (1, 2))), len(g.atomic_samb("lgs", (1, 1)))  # p-d ブロックと、スピン付きの p。
(22, 27)
```

- `transform_atomic_samb(samb, U_ket, U_bra)` は原子 SAMB を別の基底に変換します。基底は、テッセラル基底からのユニタリー行列か、群の立方晶・六方晶の基底を表す階数で与えます。
- `create_atomic_samb_L(L)` は、$|L, M\rangle$ 基底（階数 $L$）の原子多極子を定義から作ります。

(ja-group-cluster)=
## 7. クラスター SAMB と結合 SAMB

`cluster_samb(wp)` はサイトの位置のサイト・クラスター SAMB を、`cluster_samb(wp, "bond")` はボンドの位置のボンド・クラスター SAMB を、`(サイトやボンドの上の値, 式)` で与えます（{ref}`ja-cluster-samb`）。

```python
>>> c = g.cluster_samb("4b")  # サイト (x,x,0), (-x,-x,0), (-x,x,0), (x,-x,0)。
>>> for key, (v, ex) in c.items(): print(key[0], key[1], key[2], v.tolist(), ex.tolist())
Q 0 A1g [[1/2, 1/2, 1/2, 1/2]] [1]
Q 2 B2g [[1/2, 1/2, -1/2, -1/2]] [sqrt(3)*x*y]
Q 1 Eu [[1/2, -1/2, -1/2, 1/2], [1/2, -1/2, 1/2, -1/2]] [x, y]
>>> b = g.cluster_samb("2b@1o", "bond")  # x 方向と y 方向のボンド。
>>> [(key.X, key.l, key.Gamma) for key in b.named_keys()]
[('Q', 0, 'A1g'), ('Q', 2, 'B1g'), ('T', 1, 'Eu')]
```

`combined_samb(原子のキー, クラスターのキー)` は、原子 SAMB とクラスター SAMB を Clebsch-Gordan 係数 `cg(key1, key2, key)` で既約表現に結合します（{ref}`ja-combined-samb`）。結合 SAMB の各成分は `(係数, 原子のキー, 原子の成分, クラスターのキー, クラスターの成分)` のリストです。キーワード `X`、`l`、`Gamma`、`s` で結合 SAMB を選べます。

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

1 つ目の $A_{1g}$ SAMB（$p = 1$）は、すべてのボンドで一様な $p_x$-$p_x$ と $p_y$-$p_y$ の和のホッピングです。2 つ目（$p = 2$）は、$p_x$-$p_x$ と $p_y$-$p_y$ の差の軌道部分（$B_{1g}$）と、$x$ 方向と $y$ 方向のボンドで符号が逆のボンド部分（$B_{1g}$）の積で、$\sigma$ 型（$x$ 方向のボンドでの $p_x$-$p_x$）と $\pi$ 型（$x$ 方向のボンドでの $p_y$-$p_y$）のホッピングの差です。

`multipole_cluster_samb(X, l, site)` は、サイトやボンドのクラスター上の多極子 `(X, l)` の結合 SAMB を作り、`combined_object(wp, X, samb)` はその成分を各サイト上の多項式として与えます。例えば、4 つのサイト `4c`（$(\pm 1, 0, 0)$、$(0, \pm 1, 0)$）上の電気双極子（`"Q", 1`）は次のとおりです。

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

$\mathbb{Q}_{0}(A_{1g})$ では双極子は外向きで、$\mathbb{G}_{1}(A_{2g})$ では $z$ 軸のまわりを回ります。後者は $z$ 方向の電気トロイダル双極子です。

(ja-group-transform)=
## 8. 対称操作による変換

`transform_vector(vector, X)` と `transform_multipole(multipole, X)` は、すべての対称操作を `symmetry_operation["tag"]` の順に作用させます。`X` は種類を表し、`"Q"` は極性、`"G"` は軸性（$\det R$ がかかる）、`"T"` と `"M"` はそれぞれ極性と軸性に時間反転の符号（磁気群の `tr_sign`）がかかります。

```python
>>> so["tag"][:4]
['1', '2[001]', '4+[001]', '4-[001]']
>>> g.transform_vector("[1,0,0]")[:4]
array([[1, 0, 0],
       [-1, 0, 0],
       [0, 1, 0],
       [0, -1, 0]], dtype=object)
>>> m = Group("15.7.59")  # 4/m'm'm'。
>>> m.symmetry_operation["tag"][8], m.symmetry_operation["tr_sign"][8]
("-1'", -1)
>>> m.transform_vector("[0,0,1]", X="M")[[0, 8]]  # 磁気モーメントは -1' で反転する。
array([[0, 0, 1],
       [0, 0, -1]], dtype=object)
>>> g.transform_multipole("x*y")[:4]
array([x*y, x*y, -x*y, -x*y], dtype=object)
```

`create_cell_vector("ベクトル#サイト", X)` と `create_cell_multipole("多極子#サイト", X)` は、変換したベクトルや多極子を等価なサイトに配置します（同じサイトを与える対称操作について平均します）。

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

(ja-group-magnetic)=
## 9. 磁気群と応答テンソル

`active_multipole` は磁気点群（PG ではその灰色群）で許される多極子の一覧を、`response_tensor_all(X)` はタイプ `X` の階数 4 までのすべての応答テンソルを `(テンソル, {成分: 多極子})` で与えます（{ref}`ja-response-tensors`）。4/m'm'm' では電気磁気テンソル（`"M"`、階数 2）が許されます。

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

`response_tensor(X, rank, comp, opt)` は 1 つの成分を与えます（添字 1, 2, 3 が $x, y, z$）。例えば上の $M_{zz}$ です。強磁性の 4/mm'm' では、磁化 $M_z$ と異常ホール効果（`"T"`、階数 2、反対称）はどちらも磁気双極子 $\mathbb{M}_z$ が担います。

```python
>>> m.response_tensor("M", 2, (3, 3), "s")
\mathbb{M}_{0}^{(1)} + 2*\mathbb{M}_{u}^{(1)}
>>> f = Group("15.6.58")
>>> f.name(), f.active_multipole
("4/mm'm'", ['Qs', 'Qdu', 'Qg', 'Qgu', 'Tgaz', 'Mpz', 'Mfaz'])
>>> f.response_tensor("M", 1, 3), f.response_tensor("T", 2, (1, 2), "a")
(\mathbb{M}_{z}^{(1)}, \mathbb{M}_{z}^{(1)})
```
