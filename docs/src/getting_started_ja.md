---
orphan: true
---

# はじめに

英語版は[こちら](getting_started.md)です。

## チュートリアル

チュートリアルとして、**グラフェン** の場合にモデルを生成する手順を説明します。

1. [examples](https://github.com/CMT-MU/MultiPie/tree/main/docs/src/examples) ディレクトリから `graphene_in.py` と `graphene_ctrl.py` をダウンロードします。
2. 上の二つのファイルがある場所で、次を実行します：

    ```bash
    $ mp_create -v graphene_in.py
    ```

    `graphene` ディレクトリの下に、次のファイルが作られます：

     - `graphene.pkl`：モデルの情報のファイル（バイナリ）
     - `graphene.tex`、`graphene.pdf`：モデルの情報のまとめ
     - `graphene.qtdw`：モデルの構造の QtDraw ファイル

    `.tex`/`.pdf` ファイルは、`"pdf": {"create": True}`（デフォルト）で、かつ[インストール](install_ja.md)に挙げた TeX の構成要素（pLaTeX、`ptex2pdf`、LaTeX パッケージ）が使えるときだけ作られます：

     - `latex` が見つからない場合は、どちらのファイルも書き出しません（`-v` を付けると、その旨を表示します）。
     - そうでなければ、まず `.tex` ファイルを書き出します。そのあと `ptex2pdf` や LaTeX パッケージがないとわかった場合や、コンパイルが失敗したり終わらなかったりした場合は、警告を表示し、以前の実行で作られたものや不完全なものである `.pdf` ファイルを削除します（例えば権限がないなどの理由で削除できない場合を除きます）。
     - `.pkl` ファイルとほかの出力は、どの場合にも書き出します。

    `.qtdw` ファイルは、`"qtdraw": {"create": True}`（デフォルト）で、かつ QtDraw がインストールされているときだけ作られます（`-v` を付けると、QtDraw が見つからない場合にその旨を表示します）。

3. モデルを解析する（例えば分散を描く）には、次を実行します：

    ```bash
    $ mp_analyze -v graphene_ctrl.py
    ```

    `graphene` ディレクトリの下に、次のファイルが作られます：

     - `info/graphene_matrix.py`、`info/graphene_hr.dat`：選んだ SAMB の行列の完全な情報
     - `output/`：物理量のさまざまな出力
     - `samb/`：SAMB の QtDraw ファイル。`samb_figure` が `True` で、QtDraw がインストールされていて、モデルが `"qtdraw": {"create": True}` で作られたときに書き出します。`samb_figure` が `True` でもほかの条件のどれかが満たされない場合は、警告を表示し、ディレクトリは作りません。

4. `graphene.pkl` を対話的に扱うには、IPython のインターフェースを使います。詳しくは [analyze_model.ipynb](examples/analyze_model.ipynb) を参照してください。

## モデルの入力ファイル

グラフェンの場合、モデルの入力ファイルは次のとおりです。

```{literalinclude} examples/graphene_in.py
```

ほかの設定はデフォルト値で与えられ、それは次のとおりです：

```{literalinclude} examples/default_model.py
```

グラフェンの場合のように、どの設定でも上書きできます。

:::{note}
MultiPie Ver.3 では、4 つの f 軌道の名前を変えました。古い名前はもう受け付けないので、モデルの入力ファイルで置き換えてください。

| Ver.2 | Ver.3 | Wannier90 |
|:---|:---|:---|
| `f2` | `fc` | `fx(x2-3y2)` |
| `f1` | `fs` | `fy(3x2-y2)` |
| `f3x` | `fx` | `fxz2` |
| `f3y` | `fy` | `fyz2` |
:::

## 制御ファイル

グラフェンの場合、制御ファイルは次のとおりです。

```{literalinclude} examples/graphene_ctrl.py
```

制御ファイルのデフォルト値は次のとおりです：

```{literalinclude} examples/default_control.py
```

典型的な使い方では、まずすべての既約表現の SAMB を作り、それから制御ファイルで `samb/select` や `samb/parameter` を指定して、恒等表現に加えて対称性を破る項などの必要な SAMB を選びます。

### 後処理ツール

状態密度、フィリングからのフェルミエネルギー、軌道の重み（ファットバンド）、固有値と固有ベクトルの出力、独自の計算のためのコールバックは、制御ファイルの `output` で設定します。詳しい説明と例は[後処理ツール](post_processing_ja.md)を参照してください。

## 入力ファイルの形式

モデルの入力ファイルと制御ファイルは、辞書を含む Python ファイルです（例：`graphene_in = {...}`）。

- これらは Python のコードとしては実行されません。各辞書は `ast.literal_eval` で読むので、リテラル（文字列、数値、リスト、タプル、辞書、`True`/`False`/`None`）だけが使えます。`2*2` のような算術、変数、`import` は使えません。
- 分数や記号的な値は、文字列で与えます（例：`"[1/3,2/3,0]"`）。数式を表す文字列は SymPy で解析します。
- `mp_create`（または `create_model`）では、モデルの入力ファイル中のすべての辞書 `name = {...}` をそれぞれ別のモデルとして扱い、`mp_analyze`（または `analyze_model`）では、制御ファイル中のすべての辞書をそれぞれ別の制御として扱います。`MaterialModel.analyze(filename)` と `ModelAnalyzer.analyze(filename)` は、1 つのファイルに 1 つの辞書を想定しています。
- 未知のキー（例えば `"spinfull"` のようなタイプミス）はエラーになります：`ValueError: unknown key(s) in model 'x': 'spinfull' (did you mean 'spinful'?).` 受け付けるキーは、上のデフォルトのモデルとデフォルトの制御にあるキーと、次のものです。
  - `cell`：`a`、`b`、`c`、`alpha`、`beta`、`gamma`
  - `SAMB_select`、`atomic_select`、`site_select`、`bond_select`：`X`、`l`、`Gamma`、`s`
  - `samb/select`：`site`、`bond`、`X`、`l`、`Gamma`、`s`

  次の辞書のキーはデフォルトと照合しませんが、それぞれの規則があります：

  - `site` のサイト名：`;` を含まない空でない文字列で、`bond` で使うサイトでは `_` も含まないもの（`;` と `_` はサイトとボンドのクラスターの名前に使います）。名前は PDF（LaTeX）やファイル名にも現れるので、すべてのサイトで英数字（例：`"Fe1"`）を推奨します。モデル名とサイト名に含まれる LaTeX の特殊文字（例：`#` や `%`）は、PDF ではエスケープされます。
  - `samb/parameter` の SAMB 名：生成した SAMB の名前（例：`"z1"`）。`samb/parameter` には、同じ形の辞書を 1 つ含む Python ファイルの名前を与えることもできます。そのときは `model_name/info/` からの相対パスで、拡張子を付けて与えます（例：`model_name/info/my_z.py` なら `"my_z.py"`。`"info/my_z.py"`、`"my_z"`、絶対パスではありません）。パラメータが空でないとき、`mp_analyze` はそれを `model_name/info/model_name_z.py` に書き出します。このファイルは実行のたびに上書きされるので、編集する場合は別の名前にコピーしてください。
  - `output/dispersion/k_point` の k 点のラベル：`k_path` で使うラベルで、`-`、`|`、空白を含まないもの（これらは `k_path` で点を区切るのに使います）。
- `wannier/read_KS` は削除しました：`"read_KS": True` はエラーになり、古い制御ファイルにある `"read_KS": False` は警告を出して無視します（制御ファイルから削除してください）。

```{note}
入力ファイル中の数式を表す文字列は、SymPy の `parse_expr` で解析します。これは内部で `eval` を使います。また、モデルのファイル `model_name.pkl` は Python の `pickle` で保存しています。
入力ファイルと `.pkl` ファイルは、信頼できる入手元のものだけを使ってください。
```

## コマンドラインツール

- `mp_create` と `mp_analyze` は、エラーが起きると終了ステータス 1 で、使い方の誤り（例えば入力ファイルがない、未知のオプション）では終了ステータス 2 で終了します。エラーは、入力ファイルと、該当する場合はモデルとサイトとともに、1 つのメッセージで表示します。`-v` を付けると、代わりに完全なトレースバックを表示します。
- 複数の入力ファイルを与えた場合は、どのモデルを作成・解析するよりも前に、すべてのファイルを読み込みます。
- `mp_create -i` と `mp_analyze -i` は、デフォルトのモデルと制御を表示します。これは入力ファイルのテンプレートとして使えます。
- `python -m multipie.scripts.mp_create` と `python -m multipie.scripts.mp_analyze` も同じように動きます。
- 経過時間や警告などのメッセージは、`multipie` ロガーを通して標準出力に書き出します。`create_model` と `analyze_model`（したがって `mp_create` と `mp_analyze`）は、このロガーにまだハンドラがなければ設定します。アプリケーションのログの設定（例えば Jupyter のルートロガー）は変更しません。メッセージを別の場所に集めるには、これらを呼ぶ前に `multipie` ロガーにハンドラを加えてください。

## 並列計算

原子多極子の行列（`multipie.util.util_atomic_multipole` の `create_atomic_multipole_matrix`）と原子 SAMB（`Group.create_atomic_samb_L`）は、joblib で並列に計算します。
これらはデータベースの生成に使います。
`mp_create`（`MaterialModel.analyze`）では、大きなモデル（例えば d 軌道を持つスピンありのモデル）について、原子 SAMB と結合 SAMB の直交化も並列に行います。小さなモデルは逐次に処理します。プロセスの起動のほうが計算より時間がかかるためです。
出力はプロセス数によりません。
プロセス数は環境変数 `MULTIPIE_N_JOBS` で与えます：

- 設定なし、または空：すべてのコア（`-1`、デフォルト）
- `1`：逐次
- `n > 1`：`n` プロセス、`n < 0`：コア数 + 1 + n プロセス（最低 1）。joblib と同じです。
- `0` や整数でない値はエラーです。

```bash
$ MULTIPIE_N_JOBS=4 python my_script.py   # a script calling the functions above
```

Windows の PowerShell では、スクリプトを実行する前に `$env:MULTIPIE_N_JOBS = "4"` とします。この設定は、`Remove-Item env:MULTIPIE_N_JOBS` で削除するまでセッションの間残ります。

## 出力ファイル

出力ファイルの説明は次のとおりです（先頭に model_name が付きます）：

- **model_name**：すべてのファイルは `model_name` の下に作られます。
  - **.pdf**：モデルの情報（PDF）。
  - **.pkl**：モデルのデータ（バイナリ）。
  - **.qtdw**：モデルの構造（QtDraw）。
  - **.tex**：モデルの情報のソース。
  - **samb**/
    - **_atomic_samb.qtdw**：原子 SAMB。
    - **_A-B-..._def.qtdw**：サイトのクラスターの定義（すべてのサイトを 1 つのファイルに）。
    - **_A-B-....qtdw**：サイトのクラスター SAMB（すべてのサイトを 1 つのファイルに）。
    - **_A;B_nnn_m_def.qtdw**：A;B（nnn 近接、m）のボンドのクラスターの定義（nnn：近接の次数、3 桁。例：001）。
    - **_A;B_nnn_m.qtdw**：A;B（nnn 近接、m）のボンドのクラスター SAMB。
  - **info/**
    - **_hr.dat**：H[R] の行列データ。
    - **_info_output.py**：出力の情報。
    - **_info_samb.py**：モデルと SAMB の情報。
    - **_info_wannier.py**：wannier の情報。
    - **_info.py**：全体の情報。
    - **_k.py**：運動量表示の SAMB。
    - **_matrix.py**：SAMB の行列データ。
    - **_var.py**：zj と、ボンド 1 での原子のパラメータとの関係。
    - **_z.py**：zj パラメータ。
  - **output/**
    - **_dispersion.eps**：分散（EPS）。
    - **_dispersion.pdf**：分散（PDF）。
    - **_dispersion.txt**：分散のデータ。
    - **_dos.eps**：DOS（EPS）。
    - **_dos.pdf**：DOS（PDF）。
    - **_dos.txt**：DOS のデータ（`output/dos` が `False` でないとき）。
    - **_eigen_k_path.npz**、**_eigen_k_full.npz**：固有値と固有ベクトル（`output/save_k_path` または `output/save_k_full` が `True` のとき）。
    - **_weight.pdf**：軌道の重みで色付けした分散。
    - **_weight.txt**：軌道の重み付きの分散（`output/dispersion/weight` が `False` でないとき）。
    - **plot_dos.gnu**、**plot_weight.gnu**：DOS と軌道の重みを gnuplot で描くためのコマンドスクリプト。
    - **plot_band.gnu**：gnuplot で描くためのコマンドスクリプト。
