---
orphan: true
---

# ユーザーガイド

このガイドでは、MultiPie が提供する量と、その読み方・使い方を説明します。扱うのは、多極子のキーとその記号、点群調和関数、原子軌道、原子・クラスター・結合 SAMB（対称性適合多極子基底）、それらから作る強束縛模型、応答テンソルです。英語版は[こちら](users_guide.md)です。

群のデータは、MultiPieGen が生成した群のデータベースから取ります。すべての群の表は[群の情報](../group_doc.md)にあります。コードの例では、`Group`（群のデータ）と `MaterialModel`（模型）を使います。

```python
from multipie import Group, MaterialModel
```

目次：

1. [多極子のキー](#ja-multipole-key)
2. [点群調和関数](#ja-point-group-harmonics)
3. [原子軌道](#ja-atomic-orbitals)
4. [原子 SAMB](#ja-atomic-samb)
5. [クラスター SAMB](#ja-cluster-samb)
6. [結合 SAMB と模型](#ja-combined-samb)
7. [応答テンソル](#ja-response-tensors)
8. [点群多極子調和関数](#ja-multipolar-harmonics)

(ja-multipole-key)=
## 1. 多極子のキー

点群の多極子（調和関数、原子 SAMB、クラスター SAMB、結合 SAMB）は、すべて `multipie` パッケージの名前付きタプル `PGMultipoleType(X, l, Gamma, n, p, s, k, x)` で表します。その値には、多重項のすべての成分を、既約表現の成分の順に格納しています。

| フィールド | 意味 | 値 |
|:--|:--|:--|
| `X` | タイプ。空間反転（$\mathcal{P}$）と時間反転（$\mathcal{T}$）での偶奇で決まる | `Q` 電気（$\mathcal{P} = (-1)^l$、$\mathcal{T}$ 偶）、`G` 電気トロイダル（$(-1)^{l+1}$、偶）、`T` 磁気トロイダル（$(-1)^l$、奇）、`M` 磁気（$(-1)^{l+1}$、奇） |
| `l` | ランク | 0, 1, 2, ... |
| `Gamma` | 点群の既約表現 | 指標表での名前。例：`A1g`、`Eu`、`T2g` |
| `n` | 調和関数 `(X, l, Gamma)` の多重度。同じ `(X, l, Gamma)` の調和関数のうちどれか | 1, 2, ...。1 つだけなら `-1` |
| `p` | 追加の多重度。ほかのフィールドがすべて同じものに番号を付ける | 1, 2, ...。1 つだけなら `-1` |
| `s` | 内部ランク（内部自由度のランク。例：スピン） | 0（スカラー）、1 |
| `k` | 内部成分。軌道部分のランクは `l + k` | `-s, ..., s` |
| `x` | 内部タイプ | `q`（極性）、`g`（軸性） |

- `(X, l, Gamma, n)` は、常に 1 つの点群調和関数（その xyz の多項式と $U$）を指します。`T` と `M` は、`Q` と `G` の調和関数を使います。
- スカラーの量（調和関数、クラスター SAMB）は `s = 0, k = 0, x = "q"` です。
- キーは、普通のタプルで与えることもできます（例：`("Q", 1, "T1u", -1, -1, 0, 0, "q")`）。データは、これらのキーを持つ辞書 `Dict` で、フィールドで要素を選んだり並べ替えたりできます（例：`harmonics.select(X="Q", l=2)`）。

**記号。** 多極子は、そのキーから作った LaTeX 記号で書きます。

$$
\mathbb{X}_{l,\gamma}^{(s,k;\,\mathrm{tag})}(\Gamma, n\,p),
$$

| 部分 | キーとの対応 | 書く場合 |
|:--|:--|:--|
| $\mathbb{X}$ | `X` | $\mathbb{Q}, \mathbb{G}, \mathbb{T}, \mathbb{M}$ |
| $l$ | `l` | 常に |
| $\gamma$ | 成分（1, 2, ...） | 次元が 1 より大きい既約表現の多重項の、1 つの成分を指すときだけ |
| $(s, k)$ | `s`, `k` | `s > 0` のときだけ |
| tag | データの種類 | 調和関数ではなし。原子 SAMB は `a`、サイトのクラスター SAMB は `s`、ボンドのクラスター SAMB は `b`、結合 SAMB は `c` |
| $\Gamma$ | `Gamma` | 常に。例：$T_{1u}$ |
| $n$ | `n` | `n != -1` のときだけ |
| $p$ | `p` | `p != -1` のときだけ。英字で書く（`1 -> a, 2 -> b, ...`） |

`Group.tag_multipole(key, comp=None, latex=False, superscript="")` で記号を作ります。`superscript` が tag です。

```python
>>> Group.tag_multipole(("Q", 1, "T1u", -1, -1, 0, 0, "q"), 0, latex=True)
'\\mathbb{Q}_{1,1}(T_{1u})'
```

| キー | 記号 | 意味 |
|:--|:--|:--|
| `("Q", 1, "T1u", -1, -1, 0, 0, "q")` | $\mathbb{Q}_{1}(T_{1u})$ | O$_{\rm h}$ の電気双極子、$(x, y, z)$ |
| `("M", 1, "T1g", -1, -1, 0, 0, "q")` | $\mathbb{M}_{1}^{(a)}(T_{1g})$ | 原子 SAMB：軌道角運動量 |
| `("M", 1, "T1g", -1, -1, 1, -1, "q")` | $\mathbb{M}_{1}^{(1,-1;a)}(T_{1g})$ | 原子 SAMB：スピンの磁気双極子 |
| `("T", 1, "Eu", -1, 1, 0, 0, "q")` | $\mathbb{T}_{1}^{(b)}(E_{u},a)$ | D$_{\rm 4h}$ のボンドのクラスター SAMB：$(x, y)$ の対称性を持つ虚数の飛び移り |
| `("Q", 0, "A1g", -1, 2, 0, 0, "q")` | $\mathbb{Q}_{0}^{(c)}(A_{1g},b)$ | 同じ調和関数の 2 番目（$p = 2$）の結合 SAMB |

(ja-point-group-harmonics)=
## 2. 点群調和関数

**球面調和関数と方域調和関数。** ランク $l$ の調和関数は、同次多項式

$$
O_{lm}(\boldsymbol{r}) = \sqrt{\frac{4\pi}{2l+1}}\, r^l\, Y_{lm}(\theta, \phi), \qquad m = l, l-1, \dots, -l,
$$

から作ります（例：$O_{10} = z$、$O_{1,\pm1} = \mp(x \pm iy)/\sqrt{2}$）。その実の組み合わせが方域調和関数です（例：$C_{11} = x$、$S_{11} = y$、$C_{22} = \frac{\sqrt{3}}{2}(x^2 - y^2)$、$S_{22} = \sqrt{3}xy$）。

**点群調和関数。** 点群の調和関数 `(X, l, Gamma, n)` は、既約表現 $\Gamma$ の基底となるランク $l$ の関数の組です。その成分 $\gamma$ は

$$
X_{l\Gamma n}^{(\gamma)}(\boldsymbol{r}) = \sum_{m} U_{m\gamma}\, O_{lm}(\boldsymbol{r}),
$$

で、$U$ はユニタリ行列（`(2l+1) x dim`、行は $m = l, \dots, -l$）です。`Group.harmonics` は、各キーについて、成分の多項式と $U$ を与えます。

```python
>>> g = Group("Oh")
>>> ex, U = g.harmonics[("Q", 1, "T1u", -1, -1, 0, 0, "q")]
>>> ex
array([x, y, z], dtype=object)
>>> U
array([[-sqrt(2)/2, sqrt(2)*I/2, 0],
       [0, 0, 1],
       [sqrt(2)/2, sqrt(2)*I/2, 0]], dtype=object)
```

- 調和関数は、`Q`（極性）と `G`（軸性）について、ランク 11 まで与えます。軸性の調和関数は、非固有な操作 $R$ で $\det R$ の因子が余分に付きます。そのため、同じ関数の `Q` と `G` が別の既約表現に属することがあります。例えば C$_{\rm 4v}$ では、$\mathbb{Q}_{1}(A_{1}) = z$ ですが、$\mathbb{G}_{1}(A_{2}) = z$ です。回転だけからなる群では、`Q` と `G` は同じです。
- 同じランクに同じ既約表現が複数回現れるときは、`n = 1, 2, ...` で調和関数を区別します。例：C$_{\rm 4h}$ で $\mathbb{Q}_{2}(B_{g},1) = \frac{\sqrt{3}}{2}(x^2 - y^2)$、$\mathbb{Q}_{2}(B_{g},2) = \sqrt{3}xy$。
- 空間群や磁気群の調和関数は、その点群のものです。

(ja-atomic-orbitals)=
## 3. 原子軌道

ランク $L = 0, 1, 2, 3$（s, p, d, f）の軌道が、原子 SAMB の基底です。基底は 3 種類あり、`Group.atomic_basis(basis_type)` で与えられます。

**`lg`（実の軌道、スピンなし）。** 方域調和関数です。

| $L$ | 軌道 | 関数 |
|:--|:--|:--|
| 0 | `s` | $1$ |
| 1 | `px`, `py`, `pz` | $x$, $y$, $z$ |
| 2 | `dv`, `dxy`, `dxz`, `dyz`, `du` | $\frac{\sqrt{3}}{2}(x^2-y^2)$, $\sqrt{3}xy$, $\sqrt{3}xz$, $\sqrt{3}yz$, $z^2 - \frac{x^2+y^2}{2}$ |
| 3 | `fc`, `fs`, `fbz`, `f3`, `fx`, `fy`, `faz` | $\frac{\sqrt{10}}{4}x(x^2-3y^2)$, $\frac{\sqrt{10}}{4}y(3x^2-y^2)$, $\frac{\sqrt{15}}{2}z(x^2-y^2)$, $\sqrt{15}xyz$, $\frac{\sqrt{6}}{4}x(4z^2-x^2-y^2)$, $\frac{\sqrt{6}}{4}y(4z^2-x^2-y^2)$, $\frac{1}{2}z(2z^2-3x^2-3y^2)$ |

**`lgs`（実の軌道、スピンあり）。** `lg` にスピン `u`（上向き）と `d`（下向き）を掛けたもので、スピンを先に回します：`(px,u), (px,d), (py,u), ...`。

**`jml`（全角運動量）。** $|J, M\rangle$（$J = L \mp 1/2$）を、$J$ の昇順、$M$ の降順に並べたものです。例えば $L = 1$ では `(1/2,1/2), (1/2,-1/2), (3/2,3/2), ..., (3/2,-3/2)` です。

模型の入力では、サイトの軌道をこれらの名前で与えます（例：`"d"`（5 つすべて）、`["dyz", "dxz", "dxy"]`、`jml` では `"(3/2,p)"`）。スピンを含む模型（`"spinful": True`）では `lgs`（または `jml`）を使います。

(ja-atomic-samb)=
## 4. 原子 SAMB

原子 SAMB は、原子多極子演算子の軌道基底での行列を、点群に適合させたものです。原子の電子のサイト上の自由度、つまり電荷、軌道（四極子）秩序、軌道角運動量、スピンなどを表します。

- `X = Q, G, T, M` の多極子 $X_{lm}^{(s)}(k)$（`s = 0`（軌道）と `s = 1`（スピンに依存））を、調和関数の $U$ で対称性に適合させたもので、$X_{l\Gamma n}^{(s)}(k)[\gamma] = \sum_m U_{m\gamma} X_{lm}^{(s)}(k)$ です。
- 各成分のノルムは、対角ブロック（$L_1 = L_2$）では 1、非対角ブロックでは $1/\sqrt{2}$ です。

`Group.atomic_samb(basis_type, (L1, L2))` は、各キーについて、成分の行列 `(dim, n1, n2)` と調和関数の式を与えます。

```python
>>> a = g.atomic_samb("lg", (1, 1))  # p-p block of Oh.
>>> list(a.keys())
[('Q', 0, 'A1g', -1, -1, 0, 0, 'q'), ('Q', 2, 'Eg', -1, -1, 0, 0, 'q'), ('Q', 2, 'T2g', -1, -1, 0, 0, 'q'), ('M', 1, 'T1g', -1, -1, 0, 0, 'q')]
>>> mat, ex = a[("M", 1, "T1g", -1, -1, 0, 0, "q")]
>>> mat[0]  # x component, l_x / sqrt(2) in (px, py, pz).
array([[0, 0, 0],
       [0, 0, -sqrt(2)*I/2],
       [0, sqrt(2)*I/2, 0]], dtype=object)
```

O$_{\rm h}$ の p 軌道では、$\mathbb{Q}_{0}^{(a)}(A_{1g})$ は電荷（$\propto$ 単位行列）、$\mathbb{Q}_{2}^{(a)}(E_{g})$ と $\mathbb{Q}_{2}^{(a)}(T_{2g})$ は軌道の四極子、$\mathbb{M}_{1}^{(a)}(T_{1g})$ は軌道角運動量です。スピンを含めると（`lgs`）、スピンの磁気双極子は $\mathbb{M}_{1}^{(1,-1;a)}(T_{1g})$ で、p 軌道のスピン軌道相互作用 $\boldsymbol{l} \cdot \boldsymbol{s}$ は $\mathbb{Q}_{0}^{(1,1;a)}(A_{1g})$ に比例します。

(ja-cluster-samb)=
## 5. クラスター SAMB

Wyckoff 位置のクラスター SAMB は、その $N$ 個の等価なサイト（またはボンド）上の関数で、群の既約表現の基底になります。ある量がクラスターのサイトやボンドにどう分布するかを表します。

- **サイト。** サイトの Wyckoff 位置（例：`2a`）の $N$ 個のサイト上の関数です。値は実数、タイプは `Q` で、完全な正規直交系をなします。
- **ボンド。** ボンドの Wyckoff 位置の $N$ 本のボンド上の関数です。ラベル `2b@1o` は、中心がサイトの Wyckoff 位置 `1o` にある、ボンドの Wyckoff 位置 `2b`（2 本のボンド）を表します。対称なボンドの関数（ボンドの 2 つの向きで値が同じ。実数の飛び移り）は、タイプが `Q` または `G` で、値は実数です。反対称なもの（2 つの向きで符号が逆。虚数の飛び移り）は、タイプが `T` または `M` で、値は純虚数です。

`Group.cluster_samb(wp)` と `Group.cluster_samb(wp, "bond")` は、各キーについて、サイトまたはボンド上の値 `(dim, N)` と調和関数の式を与えます。D$_{\rm 4h}$ では：

```python
>>> d = Group("D4h")
>>> d.find_wyckoff_site("[0,0,1/4]")[0]  # two sites at (0, 0, ±z).
'2a'
>>> for key, (samb, ex) in d.cluster_samb("2a").items(): print(key[:3], samb, ex)
('Q', 0, 'A1g') [[sqrt(2)/2 sqrt(2)/2]] [1]
('Q', 1, 'A2u') [[sqrt(2)/2 -sqrt(2)/2]] [z]
>>> for key, (samb, ex) in d.cluster_samb("2b@1o", "bond").items(): print(key[:5], samb.tolist())
('Q', 0, 'A1g', -1, -1) [[sqrt(2)/2, sqrt(2)/2]]
('Q', 2, 'B1g', -1, -1) [[-sqrt(2)/2, sqrt(2)/2]]
('T', 1, 'Eu', -1, 1) [[I, 0], [0, I]]
```

`2b@1o` のボンドは、原点を通る $x$ 方向と $y$ 方向のものです。$\mathbb{Q}_{0}^{(b)}(A_{1g})$ は一様な実数の飛び移り、$\mathbb{Q}_{2}^{(b)}(B_{1g})$ は 2 本のボンドで符号が逆の実数の飛び移り（$x^2 - y^2$ の対称性）、$\mathbb{T}_{1}^{(b)}(E_{u},a)$ は $(x, y)$ の対称性を持つ虚数の飛び移りで、各ボンドに沿った電流です。

(ja-combined-samb)=
## 6. 結合 SAMB と模型

強束縛模型は、結合 SAMB で書きます。結合 SAMB は、原子 SAMB（サイト上の軌道、または bra と ket の軌道）とクラスター SAMB（サイトまたはボンド）の積を、Clebsch-Gordan 係数で 1 つの既約表現に結合したものです。

$$
Z_j = \sum_{\gamma_1 \gamma_2} C^{\gamma}_{\gamma_1 \gamma_2}\, X^{(\mathrm{a})}_{\gamma_1} \otimes Y^{(\mathrm{s/b})}_{\gamma_2},
\qquad H = \sum_j z_j Z_j.
$$

`MaterialModel.analyze(model_input)`（または `mp_create`）は、模型のサイト、軌道、ボンドについて、それらをすべて作ります。恒等表現に属する、タイプ `Q` と `G`（時間反転で偶）の SAMB が、ハミルトニアンに許されるものです。`"SAMB_select": {"Gamma": "IR"}` でそれらを選びます。グラフェン（C の $p_z$ を持つ 2 つのサイト、第 2 近接までのボンド）では：

```python
>>> mm = MaterialModel()
>>> mm.analyze({"model": "graphene", "group": 191, "cell": {"c": 4}, "site": {"C": ("[1/3,2/3,0]", "pz")},
...             "bond": [("C", "C", [1, 2])], "SAMB_select": {"Gamma": "IR"}, "pdf": {"create": False}})
>>> for z, (symbol, u, key, comp) in mm["combined_id"].items(): print(z, symbol, u.samb_type.wyckoff, u.neighbor)
z1 \mathbb{Q}_{0}^{(c)}(A_{1g}) 2c 0
z2 \mathbb{Q}_{0}^{(c)}(A_{1g}) 3a@3f 1
z3 \mathbb{Q}_{0}^{(c)}(A_{1g}) 6b@6l 2
```

- `mm["combined_id"]` は、各 $z_j$ について、その記号、サイトまたはボンドのクラスターとその近接の次数（`u.samb_type` と `u.neighbor`）、キー、成分を与えます。
- `mm["combined_samb"]` は、各結合 SAMB の定義を与えます。原子とクラスターの各ブロックについての、`(coefficient, atomic key, atomic component, cluster key, cluster component)` のリストです。
- `mm.get_combined_samb_matrix(fmt="value")` は、各 $Z_j$ の実空間での行列を、キー `(n1, n2, n3, m, n)` で与えます。キーは、格子ベクトル $\boldsymbol{R} = (n_1, n_2, n_3)$ と、`mm["full_matrix"]["ket"]` の bra と ket の添字です。例えば、グラフェンの最近接の飛び移り $Z_2$ の要素は $1/\sqrt{6}$ です（6 本のボンド、規格化済み）。
- 模型の PDF（`"pdf": {"create": True}` を付けた `mp_create`）には、すべての SAMB がその定義とともに載ります。

重み $z_j$ は、`ModelAnalyzer`（または `mp_analyze`）の制御ファイル（control）`{"samb": {"parameter": {"z1": ...}}}` で与えます。各 $Z_j$ は規格化されているので、$z_j$ はエネルギーの単位を持ちます。$z_j$ の番号付けは群のデータベースに依存しますが（[Ver.2 からの変更点](../changes_ja.md)を参照）、SAMB が張る空間は依存しません。完全な例は[物質の例](../materials.md)を参照してください。

(ja-response-tensors)=
## 7. 応答テンソル

ランク 0〜4 の応答テンソルは、2 つの物理量（例えば電流と電場）を結び付けます。磁気点群について、`Group.response_tensor_all(X)` は、タイプ `X`（`Q` は極性、`G` は軸性で、時間反転で偶。`T` は極性、`M` は軸性で、奇）の、対称性で許されるテンソルを与えます。

```python
>>> m = Group("MPG:15.1.53")  # 4/mmm.
>>> m.active_multipole
['Qs', 'Qdu', 'Qg', 'Qgu', 'Ts', 'Tdu', 'Tg', 'Tgu']
>>> tensor, ex = m.response_tensor_all("Q")[("Q", 2, "s")]  # e.g., the electric conductivity.
>>> tensor
array([[Q_{xx}, 0, 0],
       [0, Q_{xx}, 0],
       [0, 0, Q_{zz}]], dtype=object)
>>> ex
{Q_{xx}: \mathbb{Q}_{0}^{(1)} - \mathbb{Q}_{u}^{(1)}, Q_{zz}: \mathbb{Q}_{0}^{(1)} + 2*\mathbb{Q}_{u}^{(1)}}
```

- キーは `(X, rank, opt)` です。`opt` は、添字の対称性で指定するテンソルの部分で、ランク 0 と 1 では `""`、ランク 2 では `s`/`a`（対称/反対称）、ランク 3 では `ss`/`sa`/`a`、ランク 4 では `sss`/`ssa`/`aas`/`aaa`/`sa`/`as`/`s`/`a`/`t` です。
- 値は、成分の行列（0 になる成分は `0`）と、各成分を active multipole で書いた式です。active multipole は、時間反転を含む群のすべての操作で不変な多極子です。
- active multipole の名前は、`X` に軌道の名前を付けたものです（例：`Qdu`）。`multipole_rank` はそのランクを、`convert_mp_alias` はその別名を与えます（どちらも `multipie.util.util_response_tensor` にあります）。例えば、`Qdu` のランクは 2、別名は `Qu` です。

(ja-multipolar-harmonics)=
## 8. 点群多極子調和関数

点群多極子調和関数は、点群調和関数と、ランク `s = 1, 2, 3` で極性（`x = q`）または軸性（`x = g`）の内部自由度（例：変位やスピン）を組み合わせたものです。矢印と内部タイプを付けて書きます。例えば $\vec{\mathbb{G}}_{1}^{(1,0)}[q](T_{1g})$ は、O$_{\rm h}$ の電気トロイダル双極子 $\boldsymbol{Q} \times
\boldsymbol{r}$ です。`Group.multipolar_harmonics` は、それらをランク 3 まで、内部基底の係数と式とともに与えます。内部基底は、ランク `s` の立方晶系または六方晶系の軌道で、$Q_{x}$ や $G_{x}$ のような記号で書きます。
