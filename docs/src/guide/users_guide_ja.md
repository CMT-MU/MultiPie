---
orphan: true
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# ユーザーガイド

このガイドでは、MultiPie が提供する量と、その読み方・使い方を説明します。扱うのは、多極子のキーとその記号、点群調和関数、原子軌道、原子・クラスター・結合 SAMB（対称性適合多極子基底）、それらから作る強束縛模型、応答テンソルです。英語版は[こちら](../guide.md)です。

群のデータは、MultiPieGen が生成した群のデータベースから取ります。すべての群の表は[群の情報](https://cmt-mu.github.io/MultiPieGen/src/group_info_ja.html)にあります。コードの例では、`Group`（群のデータ）と `MaterialModel`（模型）を使います。

```{code-cell} ipython3
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
9. [簡易マニュアルとチュートリアル（PDF）](#ja-brief-manual)

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

```{code-cell} ipython3
Group.tag_multipole(("Q", 1, "T1u", -1, -1, 0, 0, "q"), 0, latex=True)
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

```{code-cell} ipython3
g = Group("Oh")
ex, U = g.harmonics[("Q", 1, "T1u", -1, -1, 0, 0, "q")]
ex
```

```{code-cell} ipython3
U
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

```{code-cell} ipython3
a = g.atomic_samb("lg", (1, 1))  # p-p block of Oh.
list(a.keys())
```

```{code-cell} ipython3
mat, ex = a[("M", 1, "T1g", -1, -1, 0, 0, "q")]
mat[0]  # x component, l_x / sqrt(2) in (px, py, pz).
```

O$_{\rm h}$ の p 軌道では、$\mathbb{Q}_{0}^{(a)}(A_{1g})$ は電荷（$\propto$ 単位行列）、$\mathbb{Q}_{2}^{(a)}(E_{g})$ と $\mathbb{Q}_{2}^{(a)}(T_{2g})$ は軌道の四極子、$\mathbb{M}_{1}^{(a)}(T_{1g})$ は軌道角運動量です。スピンを含めると（`lgs`）、スピンの磁気双極子は $\mathbb{M}_{1}^{(1,-1;a)}(T_{1g})$ で、p 軌道のスピン軌道相互作用 $\boldsymbol{l} \cdot \boldsymbol{s}$ は $\mathbb{Q}_{0}^{(1,1;a)}(A_{1g})$ に比例します。

(ja-atomic-multipole-operator)=
**演算子表現。** `Group.atomic_multipole_operator(X, l, m, s=0, k=0)` は、球面基底の成分 $m$ の原子多極子 $X_{lm}^{(s)}(k)$ を、H. Kusunose, R. Oiwa, and S. Hayami, J. Phys. Soc. Jpn. **89**, 104704 (2020) に従って演算子として与えます。スピンによらないもの（式 1-4）は、$\boldsymbol{l} = -i(\boldsymbol{r} \times \nabla)$ とし、$\nabla$ は $O_{lm}(\boldsymbol{r}) = \sqrt{4\pi/(2l+1)}\, r^l Y_{lm}(\hat{\boldsymbol{r}})$ にだけ作用するとして、

$$
\begin{aligned}
Q_{lm} &= O_{lm}, \\
M_{lm} &= \frac{1}{l+1}\,[(\nabla O_{lm})\cdot\boldsymbol{l} + \boldsymbol{l}\cdot(\nabla O_{lm})], \\
T_{lm} &= \frac{1}{(l+1)(l+2)}\,[(\nabla O_{lm})\cdot(\boldsymbol{r}\times\boldsymbol{l}) - (\boldsymbol{l}\times\boldsymbol{r})\cdot(\nabla O_{lm})], \\
G_{lm} &= \frac{2i}{(l+1)^2(l+2)}\,[(\nabla O_{lm})\cdot\boldsymbol{l}\,\boldsymbol{l}^2 - \boldsymbol{l}^2\,\boldsymbol{l}\cdot(\nabla O_{lm})]
\end{aligned}
$$

です。スピンに依存するもの（式 13）は、ランク $l + k$ の軌道の多極子 $X'$ とパウリ行列を結合したもので、

$$
X_{lm}^{(1)}(k) = i^{1+k}(-1)^{l+m}\sqrt{2l+1}\sum_{n}
\begin{pmatrix} l+k & l & 1\\ m-n & -m & n\end{pmatrix} X'_{l+k,m-n}\,\sigma_{1,n}
$$

です。$\sigma_{1,0} = \sigma_z$、$\sigma_{1,\pm1} = \mp(\sigma_x \pm i\sigma_y)/\sqrt2$ で、$X'$ は、`X = Q, M, T, G` に対して `k = 0` では `T`、`G`、`Q`、`M`、`k = ±1` では `M`、`Q`、`G`、`T` です（論文の表 I）。**軌道部分が常に極性の{ref}`点群多極子調和関数 <ja-multipolar-harmonics>`と違い、ここでの軌道部分 $X'$ は極性（`Q`、`T`）にも軸性（`G`、`M`。例えば軌道角運動量）にもなります。**

- 結果は、非可換な演算子、つまり位置 `x`、`y`、`z`、軌道角運動量 `l_x`、`l_y`、`l_z`（無次元）、パウリ行列 `sigma_x`、`sigma_y`、`sigma_z`（$\sigma_0 = 1$。スピンは $\boldsymbol{\sigma}/2$）の sympy の式です。
- これらは $[l_i, l_j] = i\epsilon_{ijk} l_k$、$[\sigma_i, \sigma_j] = 2i\epsilon_{ijk}\sigma_k$、$[l_i, x_j] = i\epsilon_{ijk} x_k$ に従い、$x_i$ どうしは交換し、$\sigma_i$ は $l_i$ とも $x_i$ とも交換します。
- 定義の演算子の並び順を保ちます。連続する位置の演算子の並びを $x, y, z$ の順に並べ替えることと、パウリ行列を右端に寄せることだけを行い、ほかの交換関係は使いません。$l < 2$ の $G_{lm}$ は 0 です。
- 演算子は、位置を通して、異なるランクの軌道の間（例えば p-d）も含むどの軌道にも作用します。その行列要素は、対称性への適合と規格化の前の `Group.atomic_samb` のもので、動径の因子は $\langle r^k \rangle = 1$ です。

```{code-cell} ipython3
Group.atomic_multipole_operator("Q", 0, 0, s=1, k=1)  # spin-orbit coupling, l.sigma / sqrt(3).
```

```{code-cell} ipython3
Group.atomic_multipole_operator("M", 1, 0)  # orbital angular momentum, l_z.
```

```{code-cell} ipython3
Group.atomic_multipole_operator("T", 1, 0)  # magnetic toroidal dipole, [(r x l) - (l x r)]_z / 6.
```

ランク 2 までの演算子の一覧は[原子多極子演算子](../atomic_multipole_table_ja.md)にあります。動く例は [docs/src/examples/atomic_multipole_operator.py](https://github.com/CMT-MU/MultiPie/tree/main/docs/src/examples/atomic_multipole_operator.py) です。

(ja-cluster-samb)=
## 5. クラスター SAMB

Wyckoff 位置のクラスター SAMB は、その $N$ 個の等価なサイト（またはボンド）上の関数で、群の既約表現の基底になります。ある量がクラスターのサイトやボンドにどう分布するかを表します。

- **サイト。** サイトの Wyckoff 位置（例：`2a`）の $N$ 個のサイト上の関数です。値は実数、タイプは `Q` で、完全な正規直交系をなします。
- **ボンド。** ボンドの Wyckoff 位置の $N$ 本のボンド上の関数です。ラベル `2b@1o` は、中心がサイトの Wyckoff 位置 `1o` にある、ボンドの Wyckoff 位置 `2b`（2 本のボンド）を表します。対称なボンドの関数（ボンドの 2 つの向きで値が同じ。実数の飛び移り）は、タイプが `Q` または `G` で、値は実数です。反対称なもの（2 つの向きで符号が逆。虚数の飛び移り）は、タイプが `T` または `M` で、値は純虚数です。

`Group.cluster_samb(wp)` と `Group.cluster_samb(wp, "bond")` は、各キーについて、サイトまたはボンド上の値 `(dim, N)` と調和関数の式を与えます。D$_{\rm 4h}$ では：

```{code-cell} ipython3
d = Group("D4h")
d.find_wyckoff_site("[0,0,1/4]")[0]  # two sites at (0, 0, ±z).
```

```{code-cell} ipython3
for key, (samb, ex) in d.cluster_samb("2a").items(): print(key[:3], samb, ex)
```

```{code-cell} ipython3
for key, (samb, ex) in d.cluster_samb("2b@1o", "bond").items(): print(key[:5], samb.tolist())
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

```{code-cell} ipython3
mm = MaterialModel()
mm.analyze({"model": "graphene", "group": 191, "cell": {"c": 4}, "site": {"C": ("[1/3,2/3,0]", "pz")},
            "bond": [("C", "C", [1, 2])], "SAMB_select": {"Gamma": "IR"}, "pdf": {"create": False}})
for z, (symbol, u, key, comp) in mm["combined_id"].items(): print(z, symbol, u.samb_type.wyckoff, u.neighbor)
```

- `mm["combined_id"]` は、各 $z_j$ について、その記号、サイトまたはボンドのクラスターとその近接の次数（`u.samb_type` と `u.neighbor`）、キー、成分を与えます。
- `mm["combined_samb"]` は、各結合 SAMB の定義を与えます。原子とクラスターの各ブロックについての、`(coefficient, atomic key, atomic component, cluster key, cluster component)` のリストです。
- `mm.get_combined_samb_matrix(fmt="value")` は、各 $Z_j$ の実空間での行列を、キー `(n1, n2, n3, m, n)` で与えます。キーは、格子ベクトル $\boldsymbol{R} = (n_1, n_2, n_3)$ と、`mm["full_matrix"]["ket"]` の bra と ket の添字です。例えば、グラフェンの最近接の飛び移り $Z_2$ の要素は $1/\sqrt{6}$ です（6 本のボンド、規格化済み）。
- 模型の PDF（`"pdf": {"create": True}` を付けた `mp_create`）には、すべての SAMB がその定義とともに載ります。

重み $z_j$ は、`ModelAnalyzer`（または `mp_analyze`）の制御ファイル（control）`{"samb": {"parameter": {"z1": ...}}}` で与えます。各 $Z_j$ は規格化されているので、$z_j$ はエネルギーの単位を持ちます。$z_j$ の番号付けは群のデータベースに依存しますが（[Ver.2 からの変更点](../changes_ja.md)を参照）、SAMB が張る空間は依存しません。完全な例は[物質の例](../materials_ja.md)を参照してください。

(ja-samb-order)=
### SAMB の順序

結合 SAMB $z_1, z_2, \ldots$ は、次の順に番号を付けます（3.2.0 から）。

1. 既約表現 $\Gamma$。群の指標表の順です（`SAMB_select["Gamma"]` で選んだもの）。
2. サイトまたはボンドのクラスター。サイトがボンドより先で、次に tail と head（名前の順）、近接、同じ近接の中の番号（`u.neighbor`、`u.n`）、軌道の bra-ket のブロック（bra と ket のランクと軌道）の順です。
3. 多極子。電気（$Q$、$G$）が磁気（$T$、$M$）より先で、次にそのキーの $s$、$k$、$l$、タイプ（$Q$ が $G$ より先、$T$ が $M$ より先）、$n$、$p$ の順です。
4. 既約表現の成分（成分ごとに 1 つの $z_j$）。

例えば、すべてのタイプ（`"SAMB_select": {"X": []}`）のグラフェンでは、既約表現は指標表の順の A1g, A2g, E2g, B1u, B2u, E1u です。

| $z_j$ | 既約表現 | クラスター | 多極子 |
| --- | --- | --- | --- |
| $z_1$, $z_2$, $z_3$ | A1g | サイト、第 1 近接、第 2 近接 | $Q$、$l = 0$ |
| $z_4$ | A2g | 第 2 近接 | $M$、$l = 1$ |
| $z_5$-$z_8$ | E2g | 第 1 近接、第 2 近接 | $Q$、$l = 2$（それぞれ 2 成分） |
| $z_9$, $z_{10}$ | E2g | 第 2 近接 | $T$、$l = 2$ |
| $z_{11}$-$z_{14}$ | B1u, B2u | サイト、第 1 近接、第 2 近接 | $Q$、$T$、$l = 3$ |
| $z_{15}$-$z_{20}$ | E1u | 第 1 近接、第 2 近接 | $T$、$Q$、$T$、$l = 1$ |

1 つのブロックの中では、例えばスピンのあるサイトの p-d のブロックで、$s$、$k$、$l$ をタイプより先に比べるので、$\mathbb{G}_{4}^{(1,-1)}$ は $\mathbb{Q}_{3}^{(1,-1)}$ のあと、$\mathbb{Q}_{1}^{(1,0)}$ の前に来ます。クラスター SAMB $y_j$ は、クラスター、$\Gamma$（指標表の順）、$l$、$k$、タイプ（$Q$、$G$、$T$、$M$）、$n$ の順に番号を付けます。

```{important}
3.2.0 より前は、既約表現は名前の順（例えば O$_{\rm h}$ では A1g, A1u, A2g, ...）で、多極子はタイプ（$Q$、$G$、$M$、$T$）、$s$、$k$、$l$ の順でした。SAMB は同じですが、3.2.0 以降で作り直した模型では $z_j$ の番号が変わります。保存した模型はその番号のまま使えます。`mp_convert` は、Ver.3 の以前の版の模型と、その制御ファイルの重みの番号を付け替えます（{ref}`3.2.0 より前の Ver.3 の模型の変換 <ja-converting-a-model-before-3-2>`を参照）。
```

(ja-response-tensors)=
## 7. 応答テンソル

ランク 0〜4 の応答テンソルは、2 つの物理量（例えば電流と電場）を結び付けます。157 の磁気点群それぞれについて、`Group.response_tensor_all()` は対称性で許されるテンソルを、キーが `(X, rank, opt)` の辞書で与えます（`Group.response_tensor_all(X)` ではタイプ `X` のものだけ）。

**タイプ。** テンソルのタイプ `X` は、空間反転と時間反転での偶奇で決まります。4 つのタイプは通常、i-polar、i-axial、c-polar、c-axial テンソルと呼ばれます。i（c）は、時間反転でテンソルが不変である（符号を変える）ことを表します。

| 名称 | `X` | 空間反転 | 時間反転 |
|:--|:--|:--|:--|
| i-polar テンソル | `Q` | 極性 | 偶 |
| i-axial テンソル | `G` | 軸性 | 偶 |
| c-polar テンソル | `T` | 極性 | 奇 |
| c-axial テンソル | `M` | 軸性 | 奇 |

応答は、そのテンソルに現れる多極子が担います（下の **Active multipole** を参照）。例えば次のとおりです。

| 応答 | テンソル `(X, rank, opt)` | 多極子 |
|:--|:--|:--|
| 電気伝導度（対称部分） | `("Q", 2, "s")` | 電気単極子、電気四極子 |
| ホール効果 | `("T", 2, "a")` | 磁気双極子（MD）。例：$T_{xy} = \mathbb{M}_{z}$ |
| 電気磁気効果 | `("M", 2, "s")`、`("M", 2, "a")` | 対称部分は磁気単極子（MM）と磁気四極子（MQ）、反対称部分は磁気トロイダル双極子（MTD） |
| エーデルシュタイン（電流磁気）効果 | `("G", 2, "s")`、`("G", 2, "a")` | 対称部分は電気トロイダル単極子（ETM）と電気トロイダル四極子（ETQ）、反対称部分は電気双極子（ED） |

**テンソルの部分（`opt`）。** ランク 2 以上のテンソルは、添字の対称性で部分に分けます。添字の組について `s` は対称、`a` は反対称です。下の表は `X = Q` の場合の関係です（ほかのタイプも同じ関係です）。`;` の前の添字が行列の行の組（または 3 つ組）、後ろの添字が列です。添字の組は

- `c` = `[x, y, z]`
- `v` = `[xx, yy, zz, yz, zx, xy]`（対称な組）
- `a` = `[yz, zx, xy]`（反対称な組）
- `u` = `[xx, yy, zz, yz, zx, xy, zy, xz, yx]`（すべての組）
- `t` = `[xxx, yyy, zzz, yyz, zzx, xxy, yzz, zxx, xyy, xyz]`（対称な 3 つ組）

で、「独立な数」は、対称性がないときのその部分の独立な成分の数です。

| ランク | `opt` | 成分間の関係 | 独立な数 | 行列（行 x 列） |
|:--|:--|:--|:--|:--|
| 0 | `""` | - | 1 | 1 x 1 |
| 1 | `""` | - | 3 | 1 x 3（`c`） |
| 2 | `s` | $X_{ab} = X_{ba}$ | 6 | 3 x 3（`c` x `c`） |
| 2 | `a` | $X_{ab} = -X_{ba}$ | 3 | 3 x 3（`c` x `c`） |
| 3 | `ss` | $X_{abc}$ は $a, b, c$ について完全対称 | 10 | 6 x 3（`v` x `c`） |
| 3 | `sa` | $X_{ab;c} = X_{ba;c} = -(X_{ac;b} + X_{bc;a})$ | 8 | 6 x 3（`v` x `c`） |
| 3 | `a` | $X_{ab;c} = -X_{ba;c}$ | 9 | 3 x 3（`a` x `c`） |
| 4 | `sss` | $X_{ab;cd} = X_{ba;cd} = X_{ab;dc} = X_{cd;ab}$ | 21 | 6 x 6（`v` x `v`） |
| 4 | `ssa` | $X_{ab;cd} = X_{ba;cd} = X_{ab;dc} = -X_{cd;ab}$ | 15 | 6 x 6（`v` x `v`） |
| 4 | `aas` | $X_{ab;cd} = -X_{ba;cd} = -X_{ab;dc} = X_{cd;ab}$ | 6 | 3 x 3（`a` x `a`） |
| 4 | `aaa` | $X_{ab;cd} = -X_{ba;cd} = -X_{ab;dc} = -X_{cd;ab}$ | 3 | 3 x 3（`a` x `a`） |
| 4 | `sa` | $X_{ab;cd} = X_{ba;cd} = -X_{ab;dc}$ | 18 | 6 x 3（`v` x `a`） |
| 4 | `as` | $X_{ab;cd} = -X_{ba;cd} = X_{ab;dc}$ | 18 | 3 x 6（`a` x `v`） |
| 4 | `s` | $X_{ab;c;d} = X_{ba;c;d}$（`sss` + `ssa` + `sa`） | 54 | 6 x 9（`v` x `u`） |
| 4 | `a` | $X_{ab;c;d} = -X_{ba;c;d}$（`aas` + `aaa` + `as`） | 27 | 3 x 9（`a` x `u`） |
| 4 | `t` | $X_{abc;d}$ は $a, b, c$ について完全対称 | 30 | 10 x 3（`t` x `c`） |

例えば、電気伝導度 $\sigma_{ab}$（$j_a = \sigma_{ab} E_b$）はランク 2 のテンソルです。その対称部分（`s`）が通常の伝導度で、反対称部分（`a`、時間反転で奇）がホール伝導度です。ホール伝導度は磁気双極子が担います。電気磁気テンソル $\alpha_{ab}$（$P_a = \alpha_{ab} H_b$）は軸性で時間反転に奇（`M`）です。その対称部分は磁気単極子と磁気四極子が、反対称部分は磁気トロイダル双極子が担います。エーデルシュタイン（電流磁気）効果のテンソル $\chi_{ab}$（$M_a = \chi_{ab} j_b$）は軸性で時間反転に偶（`G`）です。その対称部分はキラルな結晶のように電気トロイダル単極子と四極子が、反対称部分は極性結晶のように電気双極子が担います。

**なぜそのタイプと多極子になるのか。** 応答テンソルのタイプ `X` は、テンソルが結び付ける 2 つの量の偶奇から決まります。空間反転では、電場 $\boldsymbol{E}$、電流 $\boldsymbol{j}$、電気分極 $\boldsymbol{P}$ は極性ベクトル、磁場 $\boldsymbol{H}$ と磁化 $\boldsymbol{M}$ は軸性ベクトルです。時間反転では、$\boldsymbol{E}$ と $\boldsymbol{P}$ は偶、$\boldsymbol{j}$、$\boldsymbol{H}$、$\boldsymbol{M}$ は奇です。応答と外場を結び付けるテンソルは、その偶奇の積を持ちます。

| テンソル | 応答 = テンソル x 外場 | 空間反転 | 時間反転 | `X` |
|:--|:--|:--|:--|:--|
| 電気磁気 $\alpha_{ab}$ | $P_a = \alpha_{ab} H_b$ | 極性 x 軸性 = 軸性 | 偶 x 奇 = 奇 | `M` |
| エーデルシュタイン $\chi_{ab}$ | $M_a = \chi_{ab} j_b$ | 軸性 x 極性 = 軸性 | 奇 x 奇 = 偶 | `G` |
| 伝導度 $\sigma_{ab}$ | $j_a = \sigma_{ab} E_b$ | 極性 x 極性 = 極性 | 下を参照 | `Q`（`s`）、`T`（`a`） |

伝導度では、エネルギーを散逸する応答そのものは時間反転の対称性を持ちません。時間反転で奇な秩序変数 $\boldsymbol{X}$（例えば磁化）についてのオンサーガーの関係 $\sigma_{ab}(\boldsymbol{X}) = \sigma_{ba}(-\boldsymbol{X})$ により、伝導度は 2 つの部分に分かれます。対称部分は、緩和時間 $\tau$ に比例する散逸的な（オームの）伝導度で、$\boldsymbol{X}$ について偶です。散逸そのものと同じく時間反転で奇な因子 $\tau$ を除くと、物質の性質としては偶（`Q`）で、どの物質にもあります。反対称部分は、散逸のないホール伝導度で、$\boldsymbol{X}$ について奇（`T`）であり、物質の磁化か磁場によって時間反転が破れていることが必要です。この意味で `X` は、その応答のために物質が持つべき対称性、つまり応答を担う多極子のタイプです。

多極子は、ランク 2 のテンソル $X_{ab}$ を既約な部分に分けることで決まります。トレース（ランク 0）、反対称部分（ベクトル $\frac{1}{2}\epsilon_{abc} X_{bc}$、ランク 1）、対称でトレースのない部分（ランク 2）です。ランク 0 と 2 の部分はテンソルと同じタイプです。反対称部分のベクトルは、$\epsilon_{abc}$ が軸性なので空間反転の偶奇が逆になり、時間反転の偶奇は同じです。

| `X`（極性・軸性、時間反転） | ランク 0 | ランク 1（反対称） | ランク 2 |
|:--|:--|:--|:--|
| `Q`（極性、偶） | 電気単極子 $\mathbb{Q}_0$ | 電気トロイダル双極子 $\mathbb{G}_1$ | 電気四極子 $\mathbb{Q}_2$ |
| `G`（軸性、偶） | 電気トロイダル単極子 $\mathbb{G}_0$ | 電気双極子 $\mathbb{Q}_1$ | 電気トロイダル四極子 $\mathbb{G}_2$ |
| `T`（極性、奇） | 磁気トロイダル単極子 $\mathbb{T}_0$ | 磁気双極子 $\mathbb{M}_1$ | 磁気トロイダル四極子 $\mathbb{T}_2$ |
| `M`（軸性、奇） | 磁気単極子 $\mathbb{M}_0$ | 磁気トロイダル双極子 $\mathbb{T}_1$ | 磁気四極子 $\mathbb{M}_2$ |

このため、ホール効果（`T`、反対称）には磁気双極子（磁化または磁場）が、電気磁気効果には磁気単極子または磁気四極子（対称部分）か磁気トロイダル双極子（反対称部分）が、エーデルシュタイン効果には電気トロイダル単極子または四極子（対称部分。キラリティ）か電気双極子（反対称部分。極性）が必要です。より高いランクのテンソルも同じように分解します。

**Active multipole。** テンソルの各成分を多極子（既約なデカルト成分）に分解します。0 でない値をとれるのは、群の active multipole（時間反転を含むすべての操作で不変な多極子）だけです。これから、0 でない成分とその間の関係が決まります。`Group.active_multipole` がその一覧で、名前は `X` と軌道名をつなげたもの（例：`Qdu`）です。`multipole_rank` はその名前のランクを、`convert_mp_alias` は別名を与えます（どちらも `multipie.util.util_response_tensor`）。例えば `Qdu` のランクは 2、別名は `Qu` です。

**値。** 各キー `(X, rank, opt)` の値は、成分の行列（$Q_{xyz}$ のような記号。0 になる成分は `0`。行と列は上の表のとおり）と、0 でない各成分を active multipole で書いた式の辞書です。

```{code-cell} ipython3
m = Group("MPG:15.1.53")  # 4/mmm.
m.active_multipole
```

```{code-cell} ipython3
tensor, ex = m.response_tensor_all()[("Q", 2, "s")]  # e.g., the electric conductivity.
tensor
```

```{code-cell} ipython3
ex
```

```{code-cell} ipython3
m = Group("MPG:15.6.58")  # 4/mm'm', a ferromagnet along z.
tensor, ex = m.response_tensor_all()[("T", 2, "a")]  # e.g., the Hall effect.
tensor
```

```{code-cell} ipython3
ex
```

```{code-cell} ipython3
m = Group("MPG:20.4.74")  # -3'm', e.g., Cr2O3.
tensor, ex = m.response_tensor_all()[("M", 2, "s")]  # magneto-electric effect.
tensor
```

```{code-cell} ipython3
ex
```

```{code-cell} ipython3
m = Group("MPG:13.2.45")  # 4mm1', a polar crystal (e.g., a Rashba system).
tensor, ex = m.response_tensor_all()[("G", 2, "a")]  # Edelstein effect.
tensor
```

```{code-cell} ipython3
ex
```

```{code-cell} ipython3
m = Group("MPG:18.2.66")  # 321', a chiral crystal (e.g., Te).
tensor, ex = m.response_tensor_all()[("G", 2, "s")]  # Edelstein effect.
tensor
```

```{code-cell} ipython3
ex
```

4/mmm では、伝導度は $\mathrm{diag}(Q_{xx}, Q_{xx}, Q_{zz})$ で、active な単極子 $\mathbb{Q}_0$ と四極子 $\mathbb{Q}_u$ で与えられます。強磁性体 4/mm'm' では、0 でない反対称成分は $T_{xy} = -T_{yx}$ だけで、active な磁気双極子 $\mathbb{M}_z$、つまり $z$ 方向の磁化で与えられます。これがホール効果です。-3'm'（例えば Cr$_2$O$_3$）では、電気磁気テンソルは $\mathrm{diag}(M_{xx}, M_{xx}, M_{zz})$ で、active な磁気単極子 $\mathbb{M}_0$ と磁気四極子 $\mathbb{M}_u$ で与えられます。エーデルシュタイン効果は、極性の 4mm1' では反対称部分に現れ、active な電気双極子 $\mathbb{Q}_z$ で与えられます。キラルな 321' では対称部分に現れ、active な電気トロイダル単極子 $\mathbb{G}_0$ と四極子 $\mathbb{G}_u$ で与えられます。

(ja-multipolar-harmonics)=
## 8. 点群多極子調和関数

点群多極子調和関数は、点群調和関数と、各点の内部自由度（例えば変位、電気分極、スピン）を組み合わせたものです。

**内部自由度。** 内部自由度は、そのランク `s`（`s = 1`：ベクトル、`s = 2, 3`：より高いランクのテンソル）とタイプ `x` で決めます。タイプは極性（`x = q`。例：変位、電気分極）か軸性（`x = g`。例：スピン、回転）です。その基底 $e_\nu$ は、そのタイプのランク `s` の軌道の組で、立方晶系（O$_{\rm h}$）または六方晶系（D$_{\rm 6h}$）のものです。記号は $Q_\nu$（`x = q`）または $G_\nu$（`x = g`）で、例えば `s = 1` では $(Q_x, Q_y, Q_z)$ または $(G_x, G_y, G_z)$ です。基底は `Group.global_info()["harmonics"]["internal_basis"]` に、系（`"cubic"` または `"hexagonal"`）と `(s, x)` ごとに入っています。

**多極子球面調和関数。** 全ランク $l$ の多極子球面調和関数は、ランク $l + k$ の軌道部分と内部基底を結合したものです。

$$
O_{lm}^{(s)}(k) = i^{k+s} (-1)^{l+m} \sqrt{2l+1} \sum_{\nu}
\begin{pmatrix} l+k & l & s \\ m-\nu & -m & \nu \end{pmatrix} O_{l+k, m-\nu}\, e_\nu,
\qquad k = -s, \dots, s.
$$

```{important}
**軌道部分 $O_{l+k, m-\nu}(\boldsymbol{r}) = \sqrt{4\pi/(2l+2k+1)}\, r^{l+k} Y_{l+k, m-\nu}(\hat{\boldsymbol{r}})$ は、常に $\boldsymbol{r}$ の極性の関数です。** そのため、タイプ `X`（極性 `Q` か軸性 `G`）は、下の表のとおり、内部タイプ `x` と `s + k` の偶奇で決まります。これは MultiPie の点群多極子調和関数の定義で、H. Kusunose, R. Oiwa, and S. Hayami, J. Phys. Soc. Jpn. **89**, 104704 (2020) の原子多極子の定義とは異なります。原子多極子では、軌道部分が極性（`Q`、`T`）にも軸性（`G`、`M`。例えば軌道角運動量）にもなります。
```

**多極子のタイプ。** 内部タイプ `x` を固定すると、多極子のタイプ `X` は `s + k` の偶奇で決まります。

| 内部タイプ | `s + k` 偶 | `s + k` 奇 |
|:--|:--|:--|
| `x = q`（極性） | `Q` | `G` |
| `x = g`（軸性） | `G` | `Q` |

つまり、`X` の極性・軸性は、`s + k` が偶なら内部基底と同じ、奇なら逆です。点群の点群多極子調和関数は、点群調和関数 `(X, l, Gamma, n)` の $U$ を使った $\sum_m U_{m\gamma} O_{lm}^{(s)}(k)$ です。

**磁気的な内部自由度。** 内部基底 $Q_\nu$ と $G_\nu$ は時間反転で偶です。スピンや磁気モーメントのような、時間反転で奇な磁気的な内部自由度では、タイプを `Q` $\to$ `T`、`G` $\to$ `M` と読み替えてください。極性・軸性は変わらず、多極子が時間反転で奇になります。例えば、軸性ベクトル $\boldsymbol{G}$ をスピン $\boldsymbol{\sigma}$ とすると、$\boldsymbol{r} \cdot \boldsymbol{G}$（$\vec{\mathbb{G}}_{0}^{(1,1)}[g]$）は磁気単極子 $\boldsymbol{r} \cdot \boldsymbol{\sigma}$、$\boldsymbol{G} \times \boldsymbol{r}$（$\vec{\mathbb{Q}}_{1}^{(1,0)}[g]$）は磁気トロイダル双極子 $\boldsymbol{\sigma} \times \boldsymbol{r}$ です。

**キーと値。** キーは `PGMultipoleType(X, l, Gamma, n, p, s, k, x)` で、`s = 1, 2, 3`、`l <= 3`、47 の点群すべてについて与えられます。`Group.multipolar_harmonics` は、各キーについて、内部基底の係数（`dim` x 内部基底の関数の数。$x, y, z$ の多項式）と、内部基底の記号で書いた式を与えます。記号には矢印と上付きの $(s, k)$ が付き、内部タイプは $[q]$ または $[g]$ と書きます（例：$\vec{\mathbb{G}}_{1}^{(1,0)}[q](T_{1g})$）。`Group.tag_multipole(key, latex=True, vector=True)` は矢印の付いた記号を与えます。

```{code-cell} ipython3
coeff, ex = g.multipolar_harmonics[("G", 1, "T1g", -1, -1, 1, 0, "q")]  # Oh.
coeff  # rows: components of T1g, columns: Q_x, Q_y, Q_z.
```

```{code-cell} ipython3
ex
```

*例（O$_{\rm h}$、`s = 1`）。* 内部基底を極性ベクトル $\boldsymbol{Q} = (Q_x, Q_y, Q_z)$（`x = q`）または軸性ベクトル $\boldsymbol{G} = (G_x, G_y, G_z)$（`x = g`）とすると：

| キー | 記号 | 式 | 意味 |
|:--|:--|:--|:--|
| `("Q", 0, "A1g", -1, -1, 1, 1, "q")` | $\vec{\mathbb{Q}}_{0}^{(1,1)}[q](A_{1g})$ | $(xQ_x + yQ_y + zQ_z)/\sqrt{3}$ | $\boldsymbol{r} \cdot \boldsymbol{Q}$、電気単極子 |
| `("Q", 1, "T1u", -1, -1, 1, -1, "q")` | $\vec{\mathbb{Q}}_{1}^{(1,-1)}[q](T_{1u})$ | $(Q_x, Q_y, Q_z)$ | 極性ベクトルそのもの |
| `("G", 1, "T1g", -1, -1, 1, 0, "q")` | $\vec{\mathbb{G}}_{1}^{(1,0)}[q](T_{1g})$ | $(zQ_y - yQ_z, xQ_z - zQ_x, yQ_x - xQ_y)/\sqrt{2}$ | $\boldsymbol{Q} \times \boldsymbol{r}$、電気トロイダル双極子 |
| `("G", 0, "A1u", -1, -1, 1, 1, "g")` | $\vec{\mathbb{G}}_{0}^{(1,1)}[g](A_{1u})$ | $(xG_x + yG_y + zG_z)/\sqrt{3}$ | $\boldsymbol{r} \cdot \boldsymbol{G}$、電気トロイダル単極子 |

例えば、サイト $\boldsymbol{r}_i$ 上の電気双極子 $\boldsymbol{Q}(\boldsymbol{r}_i)$ が中心のまわりを回るような配置、$\boldsymbol{Q} \propto \hat{z} \times \boldsymbol{r}$ は、$z$ 方向の電気トロイダル双極子 $\vec{\mathbb{G}}_{1}^{(1,0)}[q](T_{1g})$ です。

*例：波数に依存するスピン軌道相互作用。* 軸性ベクトル $\boldsymbol{G}$ をスピン $\boldsymbol{\sigma}$ と、座標 $(x, y, z)$ を波数 $(k_x, k_y, k_z)$ と読み替えると、`x = g`、`s = 1` の点群多極子調和関数から、$\boldsymbol{k}$-$\boldsymbol{\sigma}$ の形のスピン軌道相互作用が読み取れます。$\boldsymbol{\sigma}$ も $\boldsymbol{k}$ も時間反転で奇なので、上の読み替え `Q` $\to$ `T`、`G` $\to$ `M` は打ち消し合い、タイプは表のとおり（`Q` または `G`。ハミルトニアンの項として時間反転で偶）です。スピン軌道相互作用は、その点群多極子調和関数が点群の恒等表現に属するとき、結晶で許されます。$\boldsymbol{k}$ について 1 次のものは次のとおりです。

| 点群 | キー | 式 | スピン軌道相互作用 |
|:--|:--|:--|:--|
| C$_{\rm 4v}$（極性） | `("Q", 1, "A1", -1, -1, 1, 0, "g")` | $(yG_x - xG_y)/\sqrt{2}$ | ラシュバ型、$k_y\sigma_x - k_x\sigma_y$ |
| D$_{\rm 2d}$ | `("G", 2, "A1", -1, -1, 1, -1, "g")` | $(xG_x - yG_y)/\sqrt{2}$ | ドレッセルハウス型、$k_x\sigma_x - k_y\sigma_y$ |
| O、D$_3$（カイラル） | `("G", 0, "A1", -1, -1, 1, 1, "g")` | $(xG_x + yG_y + zG_z)/\sqrt{3}$ | 放射状（ワイル型）、$\boldsymbol{k} \cdot \boldsymbol{\sigma}$ |

ラシュバ型の項は $z$ 方向の電気双極子 $\vec{\mathbb{Q}}_{1}^{(1,0)}[g]$ で、C$_{\rm 4v}$ の極性によって許されます。放射状の項は電気トロイダル単極子 $\vec{\mathbb{G}}_{0}^{(1,1)}[g]$ で、O や D$_3$ のカイラリティによって許されます（例えば Te。[物質の例](../materials_ja.md)を参照）。D$_3$ では、異方的な $\vec{\mathbb{G}}_{2}^{(1,-1)}[g](A_1)$、$(2zG_z - xG_x - yG_y)/\sqrt{6}$、つまり $2k_z\sigma_z - k_x\sigma_x - k_y\sigma_y$ も許されます。式は、例えば `Group("C4v").multipolar_harmonics[("Q", 1, "A1", -1, -1, 1, 0, "g")]` で得られます。

(ja-multipolar-harmonics-expression)=
**成分ごとの式。** `Group.multipolar_harmonics_expression(X, l, m, s=0, k=0, x="q", internal="cubic")` は、球面基底の成分 $m$（$-l \le m \le l$）の $O_{lm}^{(s)}(k)$ そのもの、つまり $x, y, z$ の多項式に内部基底の記号をかけた式を与えます。$l \le 3$ に限らず、内部ランクは `s = 0, 1, 2, 3`（`k = -s, ..., s`、`l + k >= 0`）で、内部基底は立方晶系（`"cubic"`）または六方晶系（`"hexagonal"`）のものです。

- `X`（`Q`、`G`、`T`、`M`）と `x` の両方を与えます。軌道部分は常に極性なので、偶奇の合わない `X` と `x` の組み合わせ（上の表にないもの）は 0 になります。
- `X = T, M` では、内部基底の記号を `Q` $\to$ `T`、`G` $\to$ `M` と読み替えます（磁気的な内部自由度）。
- `s = 0` では、内部基底の記号は $G_{0}$ と $M_{0}$（`x = g` の `X = G, M`）のときだけかけます。$Q_{0}$ と $T_{0}$ では、`Group.harmonics` と同じく多項式だけを与えます。
- 点群の点群多極子調和関数は、`Group.harmonics` の $U$（`T`/`M` では `Q`/`G` のもの）を使った $\sum_m U_{m\gamma} O_{lm}^{(s)}(k)$ です。

```{code-cell} ipython3
Group.multipolar_harmonics_expression("Q", 0, 0, s=1, k=1, x="q")  # r.Q / sqrt(3).
```

```{code-cell} ipython3
Group.multipolar_harmonics_expression("M", 1, 0, s=1, k=0, x="g")  # magnetic: G -> M, (y M_x - x M_y)/sqrt(2).
```

```{code-cell} ipython3
Group.multipolar_harmonics_expression("Q", 1, 0, s=1, k=0, x="q")  # the wrong parity: 0.
```

動く例は [docs/src/examples/multipolar_harmonics_expression.py](https://github.com/CMT-MU/MultiPie/tree/main/docs/src/examples/multipolar_harmonics_expression.py) です。

(ja-brief-manual)=

## 9. 簡易マニュアルとチュートリアル（PDF）

MultiPie を紹介する PDF の文書が二つあります。プログラムの簡易マニュアルと、MultiPie と QtDraw で多極子と群論を学ぶチュートリアル記事です。英語版は{ref}`英語のユーザーガイド <brief-manual>`にあります。

:::{note}
どちらも MultiPie Ver.2 向けに書かれています（簡易マニュアルは Version 2.3、チュートリアルは MultiPie 2.0.7 と QtDraw 3.1.0）。考え方と作業の流れは Ver.3 でも同じですが、名前や結果の一部が異なります。例えば、f 軌道 `f2`, `f1`, `f3x`, `f3y` の名前が変わり、`Group(tag, with_opt=True)` は廃止され、SAMB の番号や符号が異なることがあります。現在の版については [Ver.2 からの変更点](../changes_ja.md)と上の各節を見てください。
:::

### MultiPie 簡易マニュアル

[PDF](../brief_manual_ja.pdf)

MultiPie の概要をスライドでまとめたものです。群論的なデータ、SAMB（電気 `Q`、磁気 `M`、磁気トロイダル `T`、電気トロイダル `G` 多極子）、グラフェンの例、`mp_create` と `mp_analyze` による作業の流れ（入力ファイルと制御ファイルのすべてのキー）、`Group`、`MaterialModel`、`ModelAnalyzer` の各クラスの概要を扱います。

### QtDraw を用いた MultiPie のチュートリアル

[PDF](../tutorial_via_qtdraw_ja.pdf)

*MultiPie と QtDraw で学ぶ多極子表現と群論*（大槻純也、岡山大学）。QtDraw の画面を使い、利用者の立場から MultiPie を紹介します。基本的な使い方、軌道の結晶場分裂と既約表現の積の分解、p 軌道の原子多極子、URhSn（D$_{\rm 3h}$）の U サイトの電気四極子のクラスター多極子による分類、iPython からの利用を扱います。
