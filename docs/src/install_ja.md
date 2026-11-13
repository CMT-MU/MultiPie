---
orphan: true
---

# MultiPie のインストール

英語版は[こちら](install.md)です。

**MultiPie** は、Python 3.12 以降で pip を使って PyPI からインストールできます。

## 必要なもの

- Python ≥ 3.12
- [オプション] LaTeX と PDF のファイルを生成するには、[TeX Live](https://www.tug.org/texlive/) の環境を推奨します。
  - PDF ファイルは `ptex2pdf` を通して pLaTeX で生成するので、日本語の TeX 環境（pLaTeX と `ptex2pdf`）が必要です。
  - 次の LaTeX パッケージを使います：`color`、`amsmath`、`breqn`、`graphicx`、`amssymb`、`bm`、`enumitem`、`braket`、`multirow`、`booktabs`、`longtable`、`caption`。
    これらの一部は、TeX Live の "basic" のような小さなインストールには含まれていません。パッケージがあるかどうかは、例えば `kpsewhich braket.sty` で確認できます。
- [オプション] QtDraw（`.qtdw`）ファイルを生成するには、[QtDraw](https://github.com/CMT-MU/QtDraw) が必要です。

## インストール

1. Python の環境を準備します

- **macOS**

  - [Homebrew](https://brew.sh) をインストールします（まだ入っていない場合）。
  - Python をインストールします：
    ```bash
    brew install python@3.13
    ```
  - Python（と、必要なら LaTeX）のパスを `.bashrc` または `.zshrc` に加えます：
    ```bash
    export PATH=/opt/homebrew/opt/python@3.13/libexec/bin:$PATH
    ```
  - 仮想環境（例：`~/.venv`）を使う場合は、その bin ディレクトリを PATH の先頭に加えます。
  - シェルを再起動します。

- **Linux**

  - ディストリビューションのパッケージマネージャを使うか、ソースから Python をインストールします。
  - Python（と、必要なら LaTeX）のパスをシェルの設定ファイル（例：`.bashrc` または `.zshrc`）に加えます：
    ```bash
    export PATH=/path/to/python/bin:$PATH
    ```
  - 仮想環境を使う場合は、その bin ディレクトリがグローバルな Python のパスより前に来るようにします。
  - シェルを再起動します。

- **Windows**

  - 次のページの説明に従って、PowerShell と Python をインストールします：
  https://microsoft.com/PowerShell
  - インストールの際に、Python のパスがシステムの PATH に加えられるようにします。

2. MultiPie をインストールします。必要な依存パッケージはすべて自動でインストールされます：
    ```bash
    pip install multipie
    ```
3. 可視化ツール QtDraw をインストールします：
    ```bash
    pip install qtdraw
    playwright install chromium  # on Linux, also run `playwright install-deps chromium`.
    ```
4. （オプション）最適化四面体法で DOS を計算するには（制御ファイルで `"dos_method": "tetrahedron"`）、[libtetrabz](https://mitsuaki1987.github.io/libtetrabz/) をインストールします。C コンパイラが必要です：
    ```bash
    pip install "multipie[tetrahedron]"  # or pip install libtetrabz
    ```

- **Linux**（WSL2 上の Ubuntu 22.04.4 LTS）
  - `.bashrc` に次の行を加えます：
    ```bash
    export QT_QPA_PLATFORM=xcb
    ```
  - 必要なシステムライブラリをインストールします：
    ```bash
    sudo apt update
    sudo apt upgrade  # optional, but recommended
    sudo apt install libxcb-cursor0
    ```

## シェルコマンド
MultiPie は、次のコマンドラインユーティリティを提供します：
- `mp_create [option] input1 input2 ...`

  モデルの入力ファイルからモデルを作成します。
- `mp_analyze [option] control1 control2 ...`

  制御ファイルによってモデルを解析します。

例については[はじめに](getting_started_ja.md)のガイドを参照してください。

## ソースコード
- PyPI: https://pypi.org/project/multipie/
- GitHub: https://github.com/CMT-MU/MultiPie
