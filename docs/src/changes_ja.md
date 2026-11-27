---
orphan: true
---

# Ver.2 からの変更点

MultiPie Ver.3 は、MultiPieGen が生成した新しい群のデータベースを使います。群そのもの（対称操作、Wyckoff 位置、指標）は変わりませんが、多極子のデータ、つまり調和関数とそれをもとに作る SAMB は、新しい方法で作られています。このページは、Ver.2 の利用者にとって何が変わるかをまとめたものです。英語版は[こちら](changes.md)です。

**要点：**

- 模型は Ver.3 で作り直してください。Ver.2 で作った模型（`.pkl`）は受け付けません。`mp_convert` は、Ver.2 の `.pkl` から模型を作り直し、その制御ファイルの重み $z_j$ を変換します（[Ver.2 の模型の変換](#ja-converting-a-model-of-ver-2)を参照）。
- 模型の SAMB $Z_j$ が張る空間は Ver.2 と同じですが、番号、符号、値は違うことがあります。重み $z_j$ は、フィッティングなどでもう一度決めてください。
- 4 つの f 軌道の名前が変わりました：`f2`、`f1`、`f3x`、`f3y` は `fc`、`fs`、`fx`、`fy` になります。
- D$_{\rm 3h}$ など一部の群で、既約表現の並びが変わりました。
- 模型のキー `max_neighbor` と `search_cell_range` は削除しました。模型の入力から消してください。
- `Group` の一部の関数が削除または変更されました（[Group](#ja-changes-group) を参照）。
- 3.2.0 から、結合 SAMB $z_j$ の番号の付け方が変わりました（既約表現は指標表の順。{ref}`SAMB の順序 <ja-samb-order>`を参照）。`mp_convert` は、Ver.3 の以前の版の模型と、その制御ファイルの重みの番号を付け替えます（{ref}`3.2.0 より前の Ver.3 の模型の変換 <ja-converting-a-model-before-3-2>`を参照）。

## 群のデータベース

- データベース（`multipie/binary_data`）は別のプロジェクトである MultiPieGen が生成し、`tools/sync_binary_data.py` で MultiPie にコピーします。`multipie/binary_data/SOURCE.txt` に、元にした MultiPieGen のコミットが記録されています。Ver.2 でデータベースを作っていたツール（`others/tools`）は削除しました。
- データベースは約 20 MB です（Ver.2 は 81 MB）。点群のオプションのデータ（`PG/*_opt.pkl`）は、各点群のデータにまとめました。ファイル名では群の ID の `:` を `_` に置き換えています（例：`PG:32` は `PG/PG_32.pkl`。Ver.2 は `PG/PG:32.pkl`）。Windows でも git から MultiPie を入れられるようにするためです。
- 群の表（PDF）も MultiPieGen が生成します。[群の情報](https://cmt-mu.github.io/MultiPieGen/src/group_info_ja.html) からリンクしています。Ver.2 のドキュメントにあった群のページは削除しました。

## 規則

**既約表現の並び。** 既約表現は指標表の並びに従います。`g` が `u` より先、`'` の付く既約表現が `''` の付くものより先です。並びが変わるのは C$_{\rm 3h}$、D$_{\rm 3h}$ とその別の設定です。例えば D$_{\rm 3h}$ では次のようになります：

| Ver.2 | Ver.3 |
|:--|:--|
| $A_1', A_2', A_1'', A_2'', E', E''$ | $A_1', A_2', E', A_1'', A_2'', E''$ |

最初の既約表現は、これまでどおり恒等表現です。

**軌道名。** 軌道の名前は次のように変わりました。模型の入力では古い名前は受け付けず、エラーで新しい名前を示します。

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
- multipolar harmonics（`Group.multipolar_harmonics`）は rank 3 まで（Ver.2 は rank 11 まで）で、内部基底は立方晶または六方晶の軌道、記号は $Q_{x}$ などです。任意の rank の成分 $m$ ごとの式は、`Group.multipolar_harmonics_expression()` で得られます（3.1.3 から。{ref}`成分ごとの式 <ja-multipolar-harmonics-expression>` を参照）。

## 原子 SAMB

- キーと行列は Ver.2 と同じです。ただし D$_{\rm 6h}$ 系では符号が変わるものがあります（調和関数に由来します）。
- 複素の点群（`-c`）にも原子 SAMB があります。

## クラスター SAMB

- 各既約表現のサイトとボンドのクラスター SAMB が張る空間は Ver.2 と同じですが、別の代表点から作るため、値は違います。例えばサイトのクラスター SAMB では、約 69% が変わらず、14% が符号が変わり、17% が値が変わります。
- ボンドのクラスター SAMB は `Group.cluster_samb(wp, "bond")` で一組として与えます。対称（`Q`、`G`、実数）と反対称（`T`、`M`、純虚数）のボンド関数を合わせたものです。`"bond_s"`、`"bond_a"`、`"vector"` の種類は削除しました。
- 反対称のボンド関数のキーは、Ver.2 と `p`、`n`、または `(X, l, Gamma)` が違うことがあります。例えばキーの`p` はトップの群から引き継ぐので、$\mathbb{T}_{1}^{(b)}(E_{1u},a)$ のような記号になります。
- 複素の点群にはクラスター SAMB はありません（Ver.2 と同じ）。`Group.cluster_samb` はエラーを出します。
- 3.3.0 から、SG 98、208、211、213、214 の 20 個の Wyckoff 位置のクラスター SAMB を、別の基準のサイトまたはボンドから平均して作ります。空間反転の剰余類の像のクラスター SAMB に、ラベルごとに写るようにするためです（{ref}`模型の剰余類 <ja-coset-model>`を参照）。これらの Wyckoff 位置を含む模型を作り直すと、これらのクラスター SAMB の値が変わります。

## 結合 SAMB と模型

- 模型の結合 SAMB が張る行列の空間は Ver.2 と同じで、SAMB の数も同じですが、番号 $z_j$、符号、値は違うことがあります。そのため、与えた重みを持つ模型のバンド構造は、重みを変換しない限り違ったものになります。
- `MaterialModel.load`（したがって `mp_analyze`）は、Ver.2 で作った模型をエラーで拒否します。`mp_create` で模型を作り直して重みをもう一度決めるか、`mp_convert` で模型と重みを変換してください。
- 模型の入力のキー `max_neighbor` と `search_cell_range` は削除しました。これらを含む入力はエラーになります。ボンドの近接は `"bond"` の指定だけで決まり、`max_neighbor`（Ver.2 の初期値は 10）による上限はありません。空間群では、指定した近接までのサイトがすべて入るまで、ボンドごとに探索範囲を広げます。Ver.2 では `search_cell_range` が狭すぎると、近接の殻が黙って抜け、近接の番号がずれていました。
- 模型のバンドは同じように対称性に従います。バージョンをまたいで成り立つテストについては{ref}`対称性によるチェック <ja-checks-by-symmetry>`を参照してください。
- 3.2.0 から、`info/<model>_var.py` は、各ボンドクラスターについてボンド 1 の head と tail のサイト（`"bond"`）を与え、原子の変数は軌道だけで `g_{head の軌道,tail の軌道}` と名付けます。例えば `g_{pz_C2,pz_C1}` は `g_{pz,pz}` になります（[はじめに](getting_started_ja.md)を参照）。
- 3.3.0 から、模型の入力キー `"coset"` で、親の群での剰余類の代表で写した SAMB と、親構造でのその固有値の重みを、`mm["coset"]` と PDF に出します（{ref}`模型の剰余類 <ja-coset-model>`を参照）。既定で計算し、`"coset": False` で省きます。
- 3.4.0 から、模型の入力キー `"create_inversion_model": True` で、最初の剰余類が空間反転であるカイラルな群の模型について、その像を作ります。サイトと代表ボンドは、`"coset"` で空間反転で写した SAMB が前提とするものです。`mm["coset"]` もこの像の代表を与えます（{ref}`模型の剰余類 <ja-coset-model>`を参照）。3.3.0 では、ラベル空間でのベクトルのノルムが 1/2 の SAMB（例えば E の s と p の混成）について、`mm["coset"]` の `"samb"` の係数が正しい値の 1/2 でした（回転の場合も同じ）。3.4.0 からは最小二乗で解きます。
- 3.5.0 から、群の任意の上位群（`Group.supergroups`、`embedding`、`coset_decomposition`、`parent_limit`、`map_site`。群のデータの直接の上位群から）と、模型の行列 $S_{kj} = \langle Z_k, g Z_j \rangle$ と $F_j = \sum_k (\delta_{kj} - S_{kj}) \langle Z_k \rangle$（`MaterialModel.transform_samb`、`image_overlap`、`parent_forbidden_part`（位数 $n > 2$ の $g$ では、$g$ の固有値への射影を使う）、`parent_shells`、`map_bond`）を与えます。`"coset"` は任意の上位群を受け付けます（{ref}`模型の剰余類 <ja-coset-model>`を参照）。重みは、親構造を持つすべての親について与えます。

(ja-converting-a-model-of-ver-2)=
## Ver.2 の模型の変換

`mp_convert`（または `convert_model`。3.1.0 から）は、模型の入力が手元にないときに、Ver.2 の模型 `<model>/<model>.pkl` を Ver.3 に変換します。模型の制御ファイル（`samb` モード）の重み $z_j$ も変換します。

```bash
mp_convert graphene graphene_ctrl.py
```

1. Ver.2 の `.pkl` から模型の入力を組み立て直し、`graphene_in.py` に書き出します。f 軌道の名前を変え、`max_neighbor` と `search_cell_range` を除きます。近接は、Ver.2 のボンドの距離で選びます。Ver.2 が近接の殻の一部を取りこぼしていたとき（`search_cell_range` が狭すぎたとき）は、近接の番号が Ver.3 と違い、Ver.3 の模型には Ver.2 が取りこぼしたボンドも入ります（重みは 0）。
2. Ver.3 で模型を作り、`graphene/graphene.pkl` とします。Ver.2 の模型は `graphene/graphene_v2.pkl` として残します。
3. ket が違うとき、または Ver.2 の SAMB が Ver.3 の SAMB の張る空間に入らないときは、変換を拒否し、Ver.2 の模型を元に戻します。
4. 重みは射影で変換します。Ver.2 の $H = \sum_j z_j Z_j$ から、Ver.3 の $z_j = {\rm Tr}[Z_j H]$ を求めます（SAMB が正規直交でないとき、例えば f 軌道の一部だけを持つサイトでは、最小二乗で求めます）。$H$ は丸め誤差の範囲で変わりません。変換した制御ファイルは `graphene_ctrl_v3.py` に書き出し、元のファイルは変えません。`samb/select` は変換しません。近接が変わったときは、その `"bond"` を確かめてください。

`mp_convert` は、変わった点（名前を変えた軌道、除いたキー、近接、SAMB の数、射影の残差）を表示します。`--no-pdf` と `--no-qtdraw` を付けると、PDF と QtDraw のファイルを作りません（既定では Ver.2 の模型と同じ）。変換済みの模型（Ver.3 の `<model>.pkl` と `<model>_v2.pkl`）では、制御ファイルだけを変換します。辞書で与えた重み（例えばフィッティングの結果）は、`convert_model("graphene", parameter={"z1": ...})` で変換し、`report["parameter"]` として返します。

(ja-converting-a-model-before-3-2)=
## 3.2.0 より前の Ver.3 の模型の変換

3.2.0 から、結合 SAMB は、指標表の既約表現の順に、1 つのブロックの中では電気（Q、G）/磁気（T、M）、$s$、$k$、$l$、タイプ、$n$、$p$ の順に番号を付けます（{ref}`SAMB の順序 <ja-samb-order>`を参照）。3.2.0 より前は、既約表現は名前の順で、多極子はタイプ（Q、G、M、T）、$s$、$k$、$l$ の順でした。SAMB そのものは同じです。

- 3.2.0 より前に作った模型（`.pkl`）は、その番号のままです。その制御ファイルで、これまでどおり解析できます。
- 同じ入力から 3.2.0 以降で作り直した模型は、番号が違います。古い制御ファイルの重み $z_j$ は、別の SAMB に掛かってしまいます。

`mp_convert` は、そのような模型と、その制御ファイルの SAMB の番号を付け替えます。

```bash
mp_convert graphene graphene_ctrl.py
```

1. `graphene_in.py` から、それがなければ古い模型から組み立て直した入力（`graphene_in.py` に書き出します）から、模型を `graphene/graphene.pkl` として作り直します。古い模型は `graphene/graphene_v3.1.pkl`（その版）として残します。
2. 古い模型と新しい模型の SAMB を、クラスターと bra-ket のブロック、キー、成分で一対一に対応させ、行列が一致することを確かめます。そうでないとき、例えば `graphene_in.py` を変えていたときは、変換を拒否し、古い模型を元に戻します。
3. 制御ファイルの `samb/parameter`、`output/callback_info`、`output/dispersion/z`、`output/z_expectation` の SAMB の番号を付け替え、変換した制御ファイルを `graphene_ctrl_v3.2.py` に書き出します。元のファイルは変えません。$H$ はまったく同じです。

`mp_convert` は、番号が変わった SAMB の数を表示します。変換済みの模型（3.2.0 以降の `<model>.pkl` と `<model>_v3.x.pkl`）では、制御ファイルだけを変換します。`convert_model` の `report["renumber"]` は、古い $z_j$ から新しい $z_j$ への対応を与えます。

(ja-changes-group)=
## Group

| Ver.2 | Ver.3 |
|:--|:--|
| `Group(tag, with_opt=True)`, `Group.opt` | 削除 |
| `Group.opt["representation_matrix"]` | `Group.representation_matrix` |
| `Group.opt["harmonics_multipole"]` | `Group.multipolar_harmonics`（rank <= 3）、`Group.multipolar_harmonics_expression()` |
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

- 応答テンソルは、多極子の名前と 4 階の `t` 部分を除いて Ver.2 と同じです。`t` 部分 $X_{abc;d}$ は $a, b, c$ について完全対称になり、独立な成分は 30 個です（Ver.2 では完全対称ではありませんでした）。
- 応答テンソルの多極子は、どれも同じ名前の調和関数（`Group.harmonics`）の正の倍数です。3.2.0 までは、応答テンソルの $\mathbb{Q}_{4v}$（`gv`）の符号が逆で、`t` 部分も完全対称ではありませんでした。
- active multipole の名前は変更後の軌道名に従います。例えば `Qf1` は `Qfs` になります。その別名（`convert_mp_alias`）はデータベースの LaTeX の名前から作ります。例えば `Qfc` は $\mathbb{Q}_{3c}$ です。
- `multipole_rank(tag)` は、active multipole の名前からそのランクを返します。
- `Group.response_tensor_all()` は、タイプ `X` を省略するとすべてのタイプのテンソルを返します（タイプはキー `(X, rank, opt)` の最初の要素です）。`response_tensor_all(X)` は Ver.2 と同じく `X` のものを返します。

## 関連パッケージ

- QtDraw 3.4 以降は MultiPie Ver.3 で動きます（2.5 以降の Ver.2 でも動きます）。QtDraw 3.3 以前は Ver.2 でのみ動きます。
- SymClosestWannier 2.1.3 以降は MultiPie Ver.3 で動きます。Wannier 関数に付ける f 軌道の名前は、入っている MultiPie に従います。対称化に使うモデルは Ver.3 で作り直してください。

(ja-brief-manual-v2)=
## Ver.2 の簡易マニュアル（PDF）

MultiPie の簡易マニュアルは Version 2.3 向けに書かれています。考え方と作業の流れは Ver.3 でも同じですが、このページに書いたように名前や結果の一部が異なります。例えば、f 軌道 `f2`, `f1`, `f3x`, `f3y` の名前が変わり、`Group(tag, with_opt=True)` は廃止され、SAMB の番号や符号が異なることがあります。英語版は{ref}`英語のページ <brief-manual-v2>`にあります。

[PDF](pdf/brief_manual_ja.pdf)

MultiPie の概要をスライドでまとめたものです。群論的なデータ、SAMB（電気 `Q`、磁気 `M`、磁気トロイダル `T`、電気トロイダル `G` 多極子）、グラフェンの例、`mp_create` と `mp_analyze` による作業の流れ（入力ファイルと制御ファイルのすべてのキー）、`Group`、`MaterialModel`、`ModelAnalyzer` の各クラスの概要を扱います。
