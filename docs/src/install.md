# Installing MultiPie

**MultiPie** can be installed from PyPI using pip on Python 3.11 or later.

## Requirements

- Python ≥ 3.11
- [Optional] A [TeX Live](https://www.tug.org/texlive/) environment is recommended for generating LaTeX and PDF files.

## Installation

1. Set up a Python environment

- **macOS**

  - Install [Homebrew](https://brew.sh) (if not already installed):
  - Install Python:
    ```bash
    brew install python@3.13
    ```
  - Add Python (and LaTeX, if needed) to your `.bashrc` or `.zshrc`:
    ```bash
    export PATH=/opt/homebrew/opt/python@3.13/libexec/bin:$PATH
    ```
  - If you use a virtual environment (e.g. `~/.venv`), prepend its bin directory to PATH.
  - Restart your shell.

- **Linux**

  - Install Python using your distribution’s package manager or from source.
  - Add Python (and LaTeX, if needed) to your shell configuration file (e.g. `.bashrc` or `.zshrc`):
    ```bash
    export PATH=/path/to/python/bin:$PATH
    ```
  - If you use a virtual environment, ensure its bin directory appears before the global Python path.
  - Restart your shell.

- **Windows**

  - Install PowerShell and Python by following the instructions at:
  https://microsoft.com/PowerShell
  - Ensure that Python is added to your system PATH during installation.

2. Install MultiPie
All required dependencies will be installed automatically.
    ```bash
    pip install multipie
    ```
3. Install QtDraw
Install the visualization tool QtDraw:
    ```bash
    pip install qtdraw
    playwright install chromium  # for Linux use `playwright install-deps chromium` instead.
    ```

- **Linux** (Ubuntu 22.04.4 LTS on WSL2)
  - Add the following line to your `.bashrc`:
    ```bash
    export QT_QPA_PLATFORM=xcb
    ```
  - Install required system libraries:
    ```bash
    sudo apt update
    sudo apt upgrade  # optional, but recommended
    sudo apt install libxcb-cursor0
    ```

## Shell Commands
MultiPie provides the following command-line utilities:
- `mp_create_samb file1, file2, ...`

  create the SAMB for the specified model files.
- `mp_create_samb_matrix file1, file2, ...`

  create the SAMB matrix for the specified select_parameter files.
- `mp_create_samb_qtdraw file1, file2, ...`

  generate QtDraw files of the SAMB for the specified model files.

See the [Getting Started](getting_started.md) guide for examples.

## Source Code
- PyPI: https://pypi.org/project/multipie/
- GitHub: https://github.com/CMT-MU/MultiPie
