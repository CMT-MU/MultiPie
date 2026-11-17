---
orphan: true
---

# 後処理ツール

英語版は[こちら](post_processing.md)です。

強束縛模型を（SAMB から、または Wannier90 の `seedname_hr.dat` から）作ったあと、`mp_analyze` はその固有値と固有ベクトルから物理量を計算します。
このページでは、制御ファイルの `output` で設定するツールを説明します。

| ツール | `output` のキー | 出力ファイル（`model_name/output/` 内、先頭に `model_name` が付く） |
|---|---|---|
| [状態密度](#ja-pp-density-of-states) | `dos`, `dos_e_num`, `dos_width`, `dos_method`, `dos_spin` | `_dos.txt`, `_dos.pdf`, `_dos.eps`, `plot_dos.gnu` |
| [占有数からのフェルミエネルギー](#ja-pp-fermi-energy) | `n_electron` | （すべての出力に使う） |
| [軌道の重み（ファットバンド）](#ja-pp-orbital-weight-fat-band) | `dispersion/weight` | `_weight.txt`, `_weight.pdf`, `plot_weight.gnu` |
| [固有値と固有ベクトル](#ja-pp-eigenvalues-and-eigenvectors) | `save_k_path`, `save_k_full` | `_eigen_k_path.npz`, `_eigen_k_full.npz` |
| [独自の計算のためのコールバック](#ja-pp-callbacks-for-custom-calculations) | `callback_k_path`, `callback_k_full`, `callback_info` | （コールバックが書き出す） |
| [局所演算子](#ja-pp-local-operators) | `dispersion/local`, `callback_info` | `_dispersion.txt`（`dispersion/local` を指定した場合）、またはコールバックに渡す |

いずれもデフォルトでは無効で、バンド分散（`_dispersion.txt`, `.pdf`, `.eps`, `plot_band.gnu`）はこれまでどおり計算されます。

## 例：グラフェン

ファイルは [docs/src/examples](https://github.com/CMT-MU/MultiPie/tree/main/docs/src/examples) にあります。模型の入力 `graphene_in.py`（[はじめに](getting_started_ja.md)を参照）、制御ファイル `graphene_post_ctrl.py`、コールバック `graphene_callback.py` です。

```{literalinclude} examples/graphene_post_ctrl.py
```

3 つのファイルを 1 つのディレクトリに置き、次を実行します。

```bash
$ mp_create graphene_in.py
$ mp_analyze -v graphene_post_ctrl.py
```

すると `graphene/output/` には、バンド分散に加えて次のファイルができます。

- `graphene_dos.txt`, `graphene_dos.pdf`, `graphene_dos.eps`, `plot_dos.gnu`：DOS。
- `graphene_weight.txt`, `graphene_weight.pdf`, `plot_weight.gnu`：`pz@C(1)` の軌道の重み。
- `graphene_eigen_k_path.npz`：k パスに沿った固有値と固有ベクトル。
- `my_band.txt`, `my_occupation.txt`：コールバックが書き出したもの。

`"n_electron": 1`（2 本のバンドのハーフフィリング）ではフェルミエネルギーは 0 で、`-v` は `fermi energy = 0.0 for n_electron = 1 (gaussian).` と表示します。

(ja-pp-density-of-states)=
## 状態密度

`output/dos` が `False` でないとき、基本単位胞あたりの状態密度（DOS）を、`grid` で与えた k 点のグリッド上で計算します（点群の模型では Γ 点のみ）。

- `"auto"`（または `True`）：サイトのクラスターと軌道ごとの部分 DOS で、スピンについては和を取ります。例えば `s@A` は `(s,u)@A(1)`, `(s,d)@A(1)`, `(s,u)@A(2)`, ... の和です。`"dos_spin": True` とすると、スピンを含む基底ではスピンを分解します。例：`(s,u)@A` と `(s,d)@A`。
- `[[ket name]]`：ケット名のリストごとに、それらのケットの部分 DOS の和を計算します。例：`[["pz@C(1)"], ["pz@C(1)", "pz@C(2)"]]`。ラベルは、ケット名を `+` でつないだものです。ケット名は、`info/model_name_hr.dat`（wannier モードでは `info/model_name_info_wannier.py`）にある基底の名前です。

`model_name_dos.txt` の列は、フェルミエネルギーから測ったエネルギー（バンド分散と同じ）、全 DOS、全 DOS の積分 N(E)、部分 DOS です。先頭のコメント行に各列の内容が書かれています。例えばグラフェンでは次のとおりです。

```text
# DOS per primitive cell (summed over spin).
# shifted by fermi energy = 0.0
# number of k points = 14400 (grid = [120, 120, 1])
# Gaussian width = 0.02449489742783178
# column 1: E
# column 2: total
# column 3: N(E), integrated total DOS
# column 4: pz@C = pz@C(1) + pz@C(2)
```

DOS は、全 DOS の積分がバンドの数になるように規格化しています。エネルギーの範囲はすべての固有値を含み、点の数は `dos_e_num`（デフォルトは 1000）です。

計算方法は `dos_method` で選びます。

- `"gaussian"`（デフォルト）：デルタ関数を、標準偏差 `dos_width` のガウス関数で広げます（`None` ならバンド幅の 1% で、エネルギーの刻み幅より小さくはしません）。刻み幅より小さい幅を与えるとエラーになります。その場合は `dos_e_num` を増やしてください（`dos_e_num` は 12 以上でなければなりません）。
- `"tetrahedron"`：[libtetrabz](https://mitsuaki1987.github.io/libtetrabz/) による最適化四面体法です。libtetrabz は別途インストールします（[インストール](install_ja.md)を参照）。N(E) を四面体法で計算し、DOS はエネルギーの各刻みでのその差分 $(N(E+dE/2)-N(E-dE/2))/dE$ とします。これにより、平坦なバンドでも規格化が保たれます。広がりのない鋭い DOS が得られますが（`dos_width` は使いません）、点群の模型には使えません。`"gaussian"` より遅く、メモリも多く使います。大きな模型では、`grid` と `dos_e_num` をほどほどの値にしてください。

`"gaussian"`（左）と `"tetrahedron"`（右）によるグラフェンの DOS です。四面体法では、van Hove 特異点とバンド端が鋭く分解されます。

```{image} figures/post_processing/graphene_dos_gaussian.png
:width: 48%
```
```{image} figures/post_processing/graphene_dos_tetrahedron.png
:width: 48%
```

(ja-pp-fermi-energy)=
## フェルミエネルギー

フェルミエネルギーは、samb モードでは 0、wannier モードと symcw モードでは `seedname.win` のものです。`output/n_electron` を与えると、基本単位胞あたり `n_electron` 個の状態が占有されるように、全 k グリッド上の固有値からフェルミエネルギーを決めます。このフェルミエネルギーは、すべての出力（バンド分散、DOS、軌道の重み、保存する固有値と固有ベクトル、コールバックに渡す値）に使います。

各バンドには 1 つの状態が入るので、`n_electron` は、スピンを含む模型では電子数、スピンを含まない模型ではその半分です。例えば、ハーフフィリングのグラフェンでは `"n_electron": 1` です。
`"dos_method": "gaussian"` では、状態をエネルギーの低い順に占有させ（T = 0）、占有された最高の状態と非占有の最低の状態の中間にフェルミエネルギーを置きます。`"tetrahedron"` では、libtetrabz で計算します。絶縁体（ギャップがあり `n_electron` が整数）では、どちらの方法でもフェルミエネルギーはギャップの中央になります。

`output/dispersion/power` に正の数 $p$ を与えると、エネルギーをそのべき $E_k^p$ にします。例えば、ボソン励起では $\omega_k = E_k^{1/2}$ です。フェルミエネルギーは 0 になり、`output/n_electron` とは併用できません（wannier モードでも、`seedname.win` のフェルミエネルギーは使いません）。バンド分散と軌道の重みは $E_k^p$ として書き出し、プロットします。DOS は $E_k^p$ のもので、`"gaussian"` と `"tetrahedron"` のどちらでも、$\delta(E - E_k)$ の代わりに $\delta(E - E_k^p)$ を使います。$E_k^p$ が実数にならないところ、つまり $E_k$ が負で $p$ が整数でないところでは、負のエネルギーが負のままになるように $-|E_k|^p$ を使います。$p$ が整数なら $E_k^p$ をそのまま使います。例えば $p$ が偶数なら正になります。`.npz` ファイルの固有値とコールバックに渡す固有値は、べきにしません。

(ja-pp-orbital-weight-fat-band)=
## 軌道の重み（ファットバンド）

`output/dispersion/weight` を `output/dos` と同じように（`"auto"` または `[[ket name]]` で）与えると、各固有ベクトルでのケットの重み、つまりケット $m$ についての $|U_k[m,n]|^2$ の和を、k パスに沿って計算します。
`model_name_weight.txt` の列は、k、フェルミエネルギーから測ったエネルギー（バンド分散と同じ）、重みです。`_weight.pdf` と `plot_weight.gnu` では、重みごとに 1 ページずつ、バンドを重みの色でプロットします。
`"auto"` では常にスピンについて和を取ります（`dos_spin` は使いません）。スピンを分解した重みには `[[ket name]]` を与えてください。例：`[["(pz,u)@C(1)", "(pz,u)@C(2)"]]`。

グラフェンでは、2 つの副格子の間の対称性により、`pz@C(1)` の重みはどのバンドでも 1/2 です。

```{image} figures/post_processing/graphene_weight.png
:width: 60%
```

## 複数の軌道を持つ例

複数の軌道を持つ模型では、部分 DOS と軌道の重みからバンドの性格がわかります。次の例は、A に s 軌道と p 軌道、B に s 軌道を持つスピンを含む模型（CsCl 型、空間群 221、恒等既約表現の SAMB に任意のパラメータを与えたもの）で、`"dos": "auto"`, `"dos_method": "tetrahedron"`, `"n_electron": 4`, `"weight": "auto"` で計算しました。

バンド分散と DOS（鋭いピークは B の s の狭いバンド）：

```{image} figures/post_processing/sp_dispersion.png
:width: 48%
```
```{image} figures/post_processing/sp_dos.png
:width: 48%
```

`s@A` と `s@B` の軌道の重み（`_weight.pdf` の 5 ページのうち 2 ページ）：

```{image} figures/post_processing/sp_weight_s_A.png
:width: 48%
```
```{image} figures/post_processing/sp_weight_s_B.png
:width: 48%
```

(ja-pp-eigenvalues-and-eigenvectors)=
## 固有値と固有ベクトル

`"save_k_path": True` と `"save_k_full": True` の一方または両方を指定すると、k パスに沿った固有値と固有ベクトル、全 k グリッド上の固有値と固有ベクトルを、それぞれ `model_name_eigen_k_path.npz` と `model_name_eigen_k_full.npz` に保存します。中身は次のとおりです。

- `Ek[k, n]`, `Uk[k, :, n]`, `k_grid`, `ket`, `fermi_energy`：後述のコールバックの説明のとおり。
- `k_linear`（k パスに沿った場合）または `grid`（全 k グリッドの場合）。

これらは `numpy.load` で読めます。例えば、Jupyter で固有ベクトルを解析するには次のようにします。

```python
import numpy as np

z = np.load("graphene/output/graphene_eigen_k_path.npz")
Ek, Uk, ket = z["Ek"], z["Uk"], list(z["ket"])
weight = np.abs(Uk[:, ket.index("pz@C(1)"), :]) ** 2  # weight[k, n] of pz@C(1).
```

コールバックと同様に、`save_k_path` にはバンド分散が必要です（`output/dispersion/k_path` が `None` でなく、模型が空間群のものであること）。また、`save_k_full` では全グリッド上のすべての固有ベクトルをメモリに保持します。

(ja-pp-callbacks-for-custom-calculations)=
## 独自の計算のためのコールバック

MultiPie が提供しない量は、H(k) の固有値と固有ベクトルを受け取るユーザー定義の関数で計算できます。
`output/callback_k_path` と `output/callback_k_full` の一方または両方に、トップディレクトリ（制御ファイルがある場所）からの相対パスで Python ファイルを与え、そのファイルで同じ名前の関数を定義します。1 つのファイルで両方を定義でき、ファイルは 1 回だけ読み込まれます。

例の `graphene_callback.py` は、k パスに沿ったバンドのエネルギーと、$T = 0$ での各ケットの占有数 $n_m = (1/N_k) \sum_{k,n} \theta(E_F - E_{kn}) |U_k[m,n]|^2$ を書き出します。

```{literalinclude} examples/graphene_callback.py
```

ハーフフィリングのグラフェンでは、`my_occupation.txt` は次のようになります。

```text
pz@C(1) 0.500000
pz@C(2) 0.500000
```

コールバックの引数は次のとおりです。

- `Ek[k, n]`：昇順の固有値（フェルミエネルギーでずらさず、`output/dispersion/power` のべきにもしません）。
- `Uk[k, :, n]`：`Ek[k, n]` の固有ベクトル。つまり各 k について `H(k) = Uk @ diag(Ek) @ Uk^†` で、`info["ket"]` の順、`output/fourier/tb_gauge` で与えたゲージです。
- `info`（dict）：
  - `k_path`（str）、`k_point`（dict、ラベルから分数座標へ）：バンド分散の k パスと k 点（バンド分散を計算しない場合は `None`）。
  - `k_grid`（ndarray）：`Ek` と `Uk` の k 点（基本単位胞での分数座標）、`[number of k, 3]`。
  - `ket`（list）：ケット名。
  - `fermi_energy`（float）、`name`（str）：フェルミエネルギーと模型の名前（模型のない wannier モードでは seedname）。
  - `_topdir`（str）、`_verbose`（bool）：トップディレクトリ（絶対パス）と verbose オプション。
  - `callback_k_path` のみ：`k_linear`（ndarray）と `k_dis_pos`（dict）。バンド分散のプロットに使う、k パスに沿った位置とラベルの位置です。
  - `callback_k_full` のみ：`grid`（list）。`[N1, N2, N3]` です。
  - `callback_k_full` のみ：`output/callback_info` で指定した量（後述）。

これらの関数は、ほかの出力のあと、解析の最後に呼ばれます。`callback_k_path` にはバンド分散が必要です。`callback_k_full` では、全グリッド上のすべての固有ベクトル、つまり `N1 x N2 x N3 x (number of kets)^2` 個の複素数をメモリに保持するので、大きな模型では粗い `grid` が必要になることがあります。

`output/callback_info` を設定すると、全 k グリッド上の各固有状態での局所演算子の期待値を、`info[keyword]`（`[k, n]`）として渡します。キーワードは `output/dispersion/local` と同じで、スピン `Sx`, `Sy`, `Sz`、軌道角運動量 `Lx`, `Ly`, `Lz`、四極子 `Qu`, `Qv`, `Qyz`, `Qxz`, `Qxy` です（定義は[局所演算子](#ja-pp-local-operators)を参照）。例えば

```python
"output": {
    "callback_k_full": "my_callback.py",
    "callback_info": ["Sx", "Sy", "Sz"],
},
```

とすると、`callback_k_full(Ek, Uk, info)` で `info["Sx"]`、`info["Sy"]`、`info["Sz"]` が得られます。`output/callback_info` には `output/callback_k_full` と `lg` または `lgs` 基底の模型が必要で、スピンの演算子にはスピンを含む模型が必要です。

```{note}
縮退したバンドの中では、固有ベクトル、したがって各バンドの値は任意です。例えば、空間反転対称性と時間反転対称性の両方を持つスピンを含む模型では、バンドはすべての k で 2 重に縮退しています。意味があるのは、縮退したバンドについての和だけです。`output/dispersion/local` についても同様です。

演算子は模型の軸に沿って定義されます。Wannier の入力では、`seedname.nnkp` の射影の軸（z 軸、x 軸、スピンの量子化軸）は考慮しません。
```

```{note}
入力ファイルと違って、コールバックのファイルは Python のコードとして実行されます。信頼できるコールバックのファイルだけを使ってください。
```

(ja-pp-local-operators)=
## 局所演算子

`output/dispersion/local` と `output/callback_info` の局所演算子は、模型の原子基底、つまり方域調和関数の軌道 `s`, `px`, `py`, `pz`, `dv`, `dxy`, `dxz`, `dyz`, `du`, `fc`, `fs`, `fbz`, `f3`, `fx`, `fy`, `faz`（スピンを含む模型ではスピン `u`, `d` 付き）での行列 $O_{mm'}$ です。$k$ での固有状態 $n$ の値は、期待値 $\langle O \rangle_{kn} = \sum_{mm'} U_k[m,n]^* \, O_{mm'} \, U_k[m',n]$ です。

- **スピン** `Sx`, `Sy`, `Sz`：$\boldsymbol{S} = \boldsymbol{\sigma}/2$（$\hbar$ 単位）。$\boldsymbol{\sigma}$ は (`u`, `d`) の順のパウリ行列です。
- **軌道角運動量** `Lx`, `Ly`, `Lz`：ランク $l$ の軌道に対する通常の角運動量演算子（$\hbar$ 単位）で、$[L_x, L_y] = iL_z$ と $\boldsymbol{L}^2 = l(l+1)$ を満たします。
- **四極子** `Qu`, `Qv`, `Qyz`, `Qxz`, `Qxy`：角度関数の行列要素 $\langle l m | Q(\hat{\boldsymbol{r}}) | l m' \rangle$ です。

  | キーワード | $Q(\hat{\boldsymbol{r}})$ |
  |---|---|
  | `Qu` | $(3z^2 - r^2)/(2r^2)$ |
  | `Qv` | $\sqrt{3}(x^2 - y^2)/(2r^2)$ |
  | `Qyz` | $\sqrt{3}\,yz/r^2$ |
  | `Qxz` | $\sqrt{3}\,zx/r^2$ |
  | `Qxy` | $\sqrt{3}\,xy/r^2$ |

  例えば、$\langle p_z | Q_u | p_z \rangle = 2/5$、$\langle d_u | Q_u | d_u \rangle = 2/7$、$\langle f_{az} | Q_u | f_{az} \rangle = 4/15$ です。

演算子は、各サイト（原子と副格子）の各軌道ランク $l$ の中で作用します。異なるサイトの間や、p と d のような異なるランクの間の行列要素は 0 です。$\boldsymbol{L}$ と $Q$ はスピンについて対角で、$\boldsymbol{S}$ は軌道について対角です。s 軌道（$l = 0$）では、$\boldsymbol{L}$ と $Q$ は 0 です。局所演算子は `lg` または `lgs` 基底の模型で使え、`jml` 基底では使えません。

## gnuplot によるプロット

PDF と EPS のプロットは matplotlib で書き出します。gnuplot スクリプトも同じプロットを書き出し、プロットの調整に便利です。

```bash
$ cd graphene/output
$ gnuplot plot_dos.gnu      # graphene_dos.pdf and .eps
$ gnuplot plot_weight.gnu   # graphene_weight.pdf (one page for each weight)
```

PDF と EPS のファイルは上書きされます。プロットの範囲とフェルミ準位は、スクリプトの先頭で変数として与えています。例えば `plot_dos.gnu` では次のとおりです。

```text
Ef = 0.0  # fermi energy (data are already shifted by Ef).
Ef_line = 0.0  # position of fermi level in plot.
xmin = -1.347219358530748
xmax = 1.347219358530748
ymin = 0.0
ymax = 1.9823887755144438
```

また、`plot_weight.gnu` には、プロットする重みのリスト `labels` があります。
