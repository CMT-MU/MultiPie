---
orphan: true
---

# 物質の例

英語版は[こちら](materials.md)です。

このページでは、MultiPie で 3 つの物質の強束縛模型を順を追って作ります。SrVO$_3$（立方晶ペロブスカイト）、単層 MoS$_2$（空間反転のない六方晶の層）、Te（スピン軌道相互作用のあるカイラル結晶）です。これらは模型の主な要素を一通り含みます。異なる軌道を持つ複数のサイト、異なるサイト間の遠方までのボンド、既約表現にプライムが付く点群（D$_{\rm 3h}$）、カイラルな群でのスピンを含む模型です。

同じ模型の入力はテスト（`tests/test_materials.py`）でも使われ、各節の最後に述べる性質を確かめています（[対称性によるチェック](#ja-checks-by-symmetry)を参照）。このページの出力は `docs/src/examples/materials/create_examples.py` で作っています。

以下のバンド分散に使う重み $z_j$ は、説明のために選んだものにすぎず、第一原理計算にフィットしたものではありません。現実的な模型を作るには、例えば `mp_analyze` の `symcw` モードで計算にフィットしてください（[後処理](post_processing_ja.md)を参照）。

## 手順の概要

各物質について、次のようにします。

1. 模型の入力ファイル `<name>_in.py` を書きます。空間群、セル、サイトとその軌道、近接ごとのボンドを与えます。`"SAMB_select": {"Gamma": "IR"}` は恒等表現の SAMB を選びます。これらがハミルトニアンに許されるものです。
2. 模型を作ります。

    ```bash
    mp_create SrVO3_in.py
    ```

    これで `SrVO3/SrVO3.pkl`（模型）と `SrVO3/SrVO3.pdf`（すべての SAMB とその定義をまとめたもの）が書き出されます。
3. 制御ファイル `<name>_ctrl.py` で、SAMB $Z_j$ の重み $z_j$ を選びます。ハミルトニアンは

    $$
    H = \sum_j z_j Z_j
    $$

    です。各 $Z_j$ は 1 に規格化されたエルミート行列なので、$z_j$ はエネルギーの単位を持ちます。
4. 模型を解析します。

    ```bash
    mp_analyze SrVO3_ctrl.py
    ```

    これで k パスに沿ったバンド分散が `SrVO3/output/` に書き出されます。

以下の SAMB の表で、"site/bond" は bra と ket のサイト、"Wyckoff" はサイトまたはボンドの Wyckoff 位置（`12b@12j` は、ボンドの中心が `12j` にあるボンドの Wyckoff 位置 `12b`）、"neighbor" はボンドの近接の次数（オンサイト項は `-`）、"bra-ket" は軌道のランクです。記号は結合 SAMB $\mathbb{X}_{l}^{(c)}(\Gamma)$ で、タイプ `X`、ランク `l`、既約表現 $\Gamma$ を表します。読み方は[ユーザーガイド](guide/users_guide_ja.md)で説明しています。各 $Z_j$ を原子 SAMB とクラスター SAMB で表した式は、模型の PDF にあります。

## SrVO$_3$：立方晶ペロブスカイト

SrVO$_3$ は立方晶ペロブスカイト構造（空間群 $Pm\bar{3}m$、No. 221、$a = 3.84$ Å）を持ちます。伝導バンドは V の $t_{2g}$ 軌道（$d_{yz}$、$d_{xz}$、$d_{xy}$）から来ており、O の $p$ 軌道と混成します。V のサイトを体心（Wyckoff `1b`）に、3 つの O のサイトを面心（`3c`）に置き、Sr は省きます。

```{literalinclude} ../../tests/data/materials/SrVO3_in.py
```

模型の SAMB は 10 個です（[SrVO3.pdf](examples/materials/SrVO3.pdf)）。

```{include} examples/materials/SrVO3_samb.md
```

各種類の SAMB の数は対称性から決まります。

- V のオンサイト（$z_3$）：$t_{2g}$ 軌道はサイト対称性 O$_{\rm h}$ の既約表現 $T_{2g}$ をなすので、エネルギーは 1 つです。
- O のオンサイト（$z_1$、$z_2$）：O のサイト対称性は D$_{\rm 4h}$ で、$p$ 軌道は V-O ボンドに沿った 1 つとそれに垂直な 2 つに分かれるので、エネルギーは 2 つです。
- V-O ホッピング（$z_8$）：$t_{2g}$ 軌道は $p$ と $\sigma$ 結合を作らないので、$pd\pi$ ホッピングだけが残ります。
- V-V ホッピング（$z_9$、$z_{10}$）：$t_{2g}$ 軌道間の $dd\pi$ と $dd\delta$ です。
- O-O ホッピング（$z_4$-$z_7$）：O-O ボンドの対称性が低いので、パラメータは 4 つです。

制御ファイルでは、O のオンサイトエネルギー、V-O ホッピング、O-O ホッピングの 1 つを与えます。

```{literalinclude} examples/materials/SrVO3_ctrl.py
```

```{image} figures/materials/SrVO3_dispersion.png
:width: 600px
```

3 本の $t_{2g}$ バンドは 9 本の O $p$ バンドの上にあります。$pd\pi$ ホッピングだけでは、$t_{2g}$ 軌道は自分の面内でしかホッピングしません（例えば $d_{yz}$ は $yz$ 面内だけ）。そのため、$t_{2g}$ バンドの 1 本は $\Gamma$-X に沿って平坦になります。$\Gamma$ と R では $t_{2g}$ バンドは 3 重に縮退し（$T_{2g}$）、$\Gamma$ では O $p$ バンドも 3 重に縮退します（$2T_{1u} + T_{2u}$）。

## MoS$_2$：2H-MoS$_2$ の単層

2H-MoS$_2$ の単層は空間群 $P\bar{6}m2$（No. 187、点群 D$_{\rm 3h}$）を持ち、空間反転はありません。Mo のサイトは `1a` にあり、その上下にある 2 つの S のサイトは `2h`、$(1/3, 2/3, \pm z)$ にあります。大きな $c$（20 Å）で層を隔てます。Mo の 5 つの $d$ 軌道と S の 3 つの $p$ 軌道を取ります。

```{literalinclude} ../../tests/data/materials/MoS2_in.py
```

模型の SAMB は 30 個です（[MoS2.pdf](examples/materials/MoS2.pdf)）。

```{include} examples/materials/MoS2_samb.md
```

D$_{\rm 3h}$ の既約表現は、$A_1'$（恒等表現）や $E''$ のようにプライムを付けて書きます。Ver.3 からは、既約表現は $A_1', A_2', E', A_1'', A_2'', E''$ の順に並びます（[Ver.2 からの変更点](changes_ja.md)を参照）。サイト対称性により、Mo の $d$ 軌道は $A_1'$（$d_{z^2}$）、$E'$（$d_{x^2-y^2}$、$d_{xy}$）、$E''$（$d_{xz}$、$d_{yz}$）に分かれ、3 つのオンサイトエネルギー $z_1$-$z_3$ を与えます。S の $p$ 軌道は $p_z$ と（$p_x$、$p_y$）に分かれ、$z_4$ と $z_5$ を与えます。

```{literalinclude} examples/materials/MoS2_ctrl.py
```

```{image} figures/materials/MoS2_dispersion.png
:width: 600px
```

$\Gamma$ では、バンドは縮退しない（$A_1'$、$A_2''$）か 2 重に縮退します（$E'$、$E''$）。11 本のバンドは、3 本の単独のバンドと 4 組の対になります。

## Te：スピン軌道相互作用のあるカイラル結晶

三方晶の Te はカイラルな空間群 $P3_121$（No. 152、点群 D$_3$）を持ちます。Te 原子は $c$ に沿ったらせん鎖をなし、結晶には鏡映も空間反転もありません。3 つの Te のサイトは `3a`、$(x, 0, 1/3)$（$x = 0.2636$）にあります。スピンを含めた $p$ 軌道（`"spinful": True`）と、鎖内のボンド（第 1 近接）および鎖間のボンド（第 2 近接）を取ります。

```{literalinclude} ../../tests/data/materials/Te_in.py
```

模型の SAMB は 66 個です（[Te.pdf](examples/materials/Te.pdf)）。

:::{admonition} Te の SAMB（66 個）
:class: dropdown

```{include} examples/materials/Te_samb.md
```
:::

2 つの特徴が現れます。

- **スピンに依存する SAMB。** スピンを含む模型では、原子 SAMB はスピンを含み、上付きの $(s, k)$（$s = 1$）を付けて書きます。例えば $z_6$、$\mathbb{Q}_{0}^{(1,1;c)}(A_1)$ はオンサイトのスピン軌道相互作用で、$\boldsymbol{l} \cdot \boldsymbol{s}$ に比例します。
- **電気トロイダル（`G`）SAMB。** D$_3$ は回転しか持たないので、恒等表現の軸性（`G`）SAMB が、極性（`Q`）のものと一緒にハミルトニアンに許されます。例えば $\mathbb{G}_{0}^{(c)}(A_1)$（$z_{19}$、$z_{20}$）はボンドの電気トロイダル単極子で、カイラリティを表します。

制御ファイルでは、スピン軌道相互作用と、鎖内および鎖間のホッピングを与え、k パスに沿った $S_z$ の期待値（`"local": ["Sz"]`）を計算して色で示します。

```{literalinclude} examples/materials/Te_ctrl.py
```

```{image} figures/materials/Te_dispersion.png
:width: 600px
```

空間反転がないので、スピン軌道相互作用は時間反転不変運動量（例えば $\Gamma$ と A）以外でバンドを分裂させます。時間反転不変運動量では、バンドは 2 重に縮退します（Kramers 縮退）。$\Gamma$-A（鎖の方向、$k_z$）に沿って、分裂したバンドの各対は互いに逆符号で $z$ 方向に偏極します（赤と青）。つまり、スピンは $\boldsymbol{k}$ に平行または反平行で、これがカイラル結晶の放射状のスピンテクスチャです。

(ja-checks-by-symmetry)=
## 対称性によるチェック

`tests/test_materials.py` は 3 つの模型を作り、すべての SAMB にランダムな重みを与えます。そして、SAMB の番号、符号、基底によらず、任意の完全な SAMB の組で成り立つ次の性質を確かめます。

- SAMB の数と、そのすべてが恒等表現の `Q` または `G` であること
- $H(\boldsymbol{k})$ がエルミートであること
- 空間群のすべての回転 $W$ について $E_n(\boldsymbol{k}) = E_n(\boldsymbol{k} W)$ であり、時間反転により $E_n(\boldsymbol{k}) = E_n(-\boldsymbol{k})$ であること
- $\Gamma$ での縮退が既約表現に従うこと：SrVO$_3$ は 4 組の 3 重縮退、MoS$_2$ は 3 つの非縮退と 4 組の 2 重縮退、Te は 9 組の Kramers 対
- 一般の $\boldsymbol{k}$ では縮退が残らないこと

これらの性質は、SAMB が同じハミルトニアンの空間を張る限り、群のデータベースが変わっても成り立ちます。
