---
orphan: true
---

# Ver.2 からの変更点

MultiPie Ver.3 は、MultiPieGen が生成した新しい群のデータベースを使います。群そのもの（対称操作、Wyckoff 位置、指標）は変わりませんが、多極子のデータ、つまり調和関数とそれをもとに作る SAMB は、新しい方法で作られています。このページは、Ver.2 の利用者にとって何が変わるかをまとめたものです。英語版は[こちら](changes.md)です。

**要点：**

- モデルは Ver.3 で作り直してください。Ver.2 で作ったモデル（`.pkl`）は受け付けません。
- モデルの SAMB $Z_j$ が張る空間は Ver.2 と同じですが、番号、符号、値は違うことがあります。重み $z_j$ は、フィッティングなどでもう一度決めてください。
- 4 つの f 軌道の名前が変わりました：`f2`、`f1`、`f3x`、`f3y` は `fc`、`fs`、`fx`、`fy` になります。
- D$_{\rm 3h}$ など一部の群で、既約表現の並びが変わりました。
- `Group` の一部の関数が削除または変更されました（[Group](#ja-changes-group) を参照）。

## 群のデータベース

- データベース（`multipie/binary_data`）は別のプロジェクトである MultiPieGen が生成し、`tools/sync_binary_data.py` で MultiPie にコピーします。`multipie/binary_data/SOURCE.txt` に、元にした MultiPieGen のコミットが記録されています。Ver.2 でデータベースを作っていたツール（`others/tools`）は削除しました。
- データベースは約 20 MB です（Ver.2 は 81 MB）。点群のオプションのデータ（`PG/*_opt.pkl`）は、各点群のデータにまとめました。
- 群の表（PDF）も MultiPieGen が生成します。[Group Information](group_doc.md) からリンクしています。Ver.2 のドキュメントにあった群のページは削除しました。

## 規則

**既約表現の並び。** 既約表現は指標表の並びに従います。`g` が `u` より先、`'` の付く既約表現が `''` の付くものより先です。並びが変わるのは C$_{\rm 3h}$、D$_{\rm 3h}$ とその別の設定です。例えば D$_{\rm 3h}$ では次のようになります：

| Ver.2 | Ver.3 |
|:--|:--|
| $A_1', A_2', A_1'', A_2'', E', E''$ | $A_1', A_2', E', A_1'', A_2'', E''$ |

最初の既約表現は、これまでどおり恒等表現です。

**軌道名。** 軌道の名前は次のように変わりました。モデルの入力では古い名前は受け付けず、エラーで新しい名前を示します。

| Ver.2 | Ver.3 | 関数 | Wannier90 |
|:--|:--|:--|:--|
| `f2` | `fc` | $\frac{\sqrt{10}}{4}x(x^2-3y^2)$ | `fx(x2-3y2)` |
| `f1` | `fs` | $\frac{\sqrt{10}}{4}y(3x^2-y^2)$ | `fy(3x2-y2)` |
| `f3x` | `fx` | $\frac{\sqrt{6}}{4}x(4z^2-x^2-y^2)$ | `fxz2` |
| `f3y` | `fy` | $\frac{\sqrt{6}}{4}y(4z^2-x^2-y^2)$ | `fyz2` |

多極子の名前（例えば応答テンソルの active multipole）に現れる rank 4 の軌道の名前も変わりました：`ga`、`gb`、`gc` は`gc`、`gs`、`gc4` になります（`gc` は旧 `ga` を指すことに注意してください）。

**多重度。** キーの多重度 `n` は常に点群の調和関数 `(X, l, Gamma)` の多重度で、xyz の多項式を特定します。ほかのフィールドがすべて同じものは `p` で番号付けします。[ユーザーガイド](guide/users_guide_ja.md) を参照してください。

## 調和関数

- `Group.harmonics` の値は `(expression, U)` です。Ver.2 の 3 番目の項目（tesseral 調和関数による式）は削除しました。
- 式は展開した多項式です。例えば $\frac{\sqrt{10}}{4}y(3x^2 - y^2)$ ではなく$\frac{3\sqrt{10}}{4}x^2y - \frac{\sqrt{10}}{4}y^3$ となります。
- キーは 47 点群すべてで同じです。32 個の調和関数で符号が変わります。これらは三方晶系と六方晶系の 16 点群の、多重度 2 または 3 の rank 4 の $E$ の調和関数です（D$_{\rm 6h}$ の $E_{2g}$ に由来します）。式と $U$ の符号は一緒に変わります。
- multipolar harmonics（`Group.multipolar_harmonics`）は rank 3 まで（Ver.2 は rank 11 まで）で、内部基底は立方晶または六方晶の軌道、記号は $Q_{x}$ などです。

## 原子 SAMB

- キーと行列は Ver.2 と同じです。ただし D$_{\rm 6h}$ 系では符号が変わるものがあります（調和関数に由来します）。
- 複素の点群（`-c`）にも原子 SAMB があります。

## クラスター SAMB

- 各既約表現のサイトとボンドのクラスター SAMB が張る空間は Ver.2 と同じですが、別の代表点から作るため、値は違います。例えばサイトのクラスター SAMB では、約 69% が変わらず、14% が符号が変わり、17% が値が変わります。
- ボンドのクラスター SAMB は `Group.cluster_samb(wp, "bond")` で一組として与えます。対称（`Q`、`G`、実数）と反対称（`T`、`M`、純虚数）のボンド関数を合わせたものです。`"bond_s"`、`"bond_a"`、`"vector"` の種類は削除しました。
- 反対称のボンド関数のキーは、Ver.2 と `p`、`n`、または `(X, l, Gamma)` が違うことがあります。例えばキーの`p` はトップの群から引き継ぐので、$\mathbb{T}_{1}^{(b)}(E_{1u},a)$ のような記号になります。
- 複素の点群にはクラスター SAMB はありません（Ver.2 と同じ）。`Group.cluster_samb` はエラーを出します。

## 結合 SAMB とモデル

- モデルの結合 SAMB が張る行列の空間は Ver.2 と同じで、SAMB の数も同じですが、番号 $z_j$、符号、値は違うことがあります。そのため、与えた重みを持つモデルのバンド構造は、重みを変換しない限り違ったものになります。
- `MaterialModel.load`（したがって `mp_analyze`）は、Ver.2 で作ったモデルをエラーで拒否します。`mp_create` でモデルを作り直し、重みをもう一度決めてください。
- モデルのバンドは同じように対称性に従います。バージョンをまたいで成り立つテストについては[Checks by symmetry](materials.md) を参照してください。

(ja-changes-group)=
## Group

| Ver.2 | Ver.3 |
|:--|:--|
| `Group(tag, with_opt=True)`, `Group.opt` | 削除 |
| `Group.opt["representation_matrix"]` | `Group.representation_matrix` |
| `Group.opt["harmonics_multipole"]` | `Group.multipolar_harmonics`（rank <= 3） |
| `Group.opt["harmonics_irrep"]` | 削除（`Group.harmonics` のキーを使う） |
| `Group.harmonics` の値 `(expression, U, tesseral expression)` | `(expression, U)` |
| `Group.cluster_samb(wp, "bond_s" / "bond_a" / "vector")` | 削除。`Group.cluster_samb(wp, "bond")` |
| `Group.response_tensor_all(X)` がテンソルを計算する | データベースから読む（テンソルは同じ） |
| `wyckoff[...]["reference"]` | 別の値（別の代表点） |
| 点群の `symmetry_operation["mapping"]`、右剰余類 | 親の群をたどる左剰余類。空間群にもある |
| 複素の点群にはボンドがない | 複素の点群にも Wyckoff のボンドがある |
| - | `wyckoff["bond"][...]["symmetry"]`（ボンドの中心のサイト対称性） |

対称操作（とその並び）、Wyckoff のラベル・位置・mapping、指標表（既約表現の並びを除く）、MPG と MSG のデータは変わりません。

## 応答テンソル

- 応答テンソルは、多極子の名前と `gv` の符号を除いて Ver.2 と同じです。
- active multipole の名前は変更後の軌道名に従います。例えば `Qf1` は `Qfs` になります。その別名（`convert_mp_alias`）はデータベースの LaTeX の名前から作ります。例えば `Qfc` は $\mathbb{Q}_{3c}$ です。
- `multipole_rank(tag)` は、active multipole の名前からそのランクを返します。

## 関連パッケージ

- QtDraw 3.4 は MultiPie Ver.3（`multipie~=3.0.dev0`）で、QtDraw 3.3 は Ver.2（`multipie~=2.5`）で動きます。リリースされるまでは、MultiPie Ver.3 を入れたあとで、QtDraw 3.4 を開発ブランチから入れてください：`pip install "git+https://github.com/CMT-MU/QtDraw.git@develop-v3.4"`。
