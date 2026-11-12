# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

MultiPie is a Python library (Python >= 3.12) that provides a database of crystallographic point groups and space groups, both ordinary and magnetic (PG/SG/MPG/MSG). It builds Symmetry-Adapted Multipole Bases (SAMB) and tight-binding models from that database. Structures and SAMB figures are written as QtDraw files (`.qtdw`), and PDF output goes through LaTeX (TeXLive, optional).

## Branch workflow

- `main` is the default branch and is used for releases. Every PR branches from `develop` and targets `develop`.
- After merging a PR into `develop`, delete the feature branch, then run the test suite on the updated `develop`.
- The maintainer opens and merges the `develop` → `main` PR by hand. Never delete `develop`.
- A GitHub release whose tag is `v` + `multipie.__version__` (e.g. `v2.4.1`) publishes the package to PyPI (`.github/workflows/release.yml`). The version is set in `multipie/core/multipie_info.py`.

## Commands

```bash
pip install -e ".[dev,tetrahedron]"        # dev install (tetrahedron = optional libtetrabz for DOS)
pytest -n auto                              # all tests in parallel (testpaths = tests)
pytest tests/test_model.py::test_model_files  # a single test
MPLBACKEND=Agg pytest -n auto               # headless, as in CI
black .                                     # formatter; line length 130 (pyproject.toml)
make -C docs build                          # Jupyter Book docs -> docs/_build/html
mp_create model_in.py                       # CLI: create model(s); `mp_create -i` prints the input template
mp_analyze model_ctrl.py                    # CLI: analyze model(s) from control file(s)
```

CI (`.github/workflows/test.yml`) installs the package non-editably on ubuntu and macOS with Python 3.12 and 3.13. It runs the tests from outside the checkout, so they must import the installed package rather than the source tree. As a result, anything the package needs at runtime has to be covered by `package-data` in `pyproject.toml`.

## Architecture

The public API is lazily imported in `multipie/__init__.py` (`Group`, `MaterialModel`, `ModelAnalyzer`, `create_model`, `analyze_model`). The same file defines the namedtuple types used throughout the package (`PGMultipoleType`, `RepSiteType`, `CellBondType`, `SAMBType`, ...). Logic lives in `multipie/core/`, and most of the heavy lifting is in the helper modules under `multipie/util/`.

The pipeline has three layers:

1. **`Group`** (`core/group.py`) is a read-only view of precomputed group data. A group is identified by a tag such as `PG:1-47`, `SG:1-230` (or a bare integer), an MPG `PG.no.ID`, an MSG BNS id, or a Schoenflies symbol. Data loads from pickles in `multipie/binary_data/{PG,SG,MPG,MSG}/` and `info.pkl` (about 81 MB) through `BinaryManager` (`util/util_binary.py`). Loads are cached per process (`Group.clear_cache()`). NumPy arrays are shared and read-only, while dicts and lists are copied for each `Group`. `Group` provides symmetry operations, Wyckoff positions, harmonics, atomic, cluster and combined SAMB, and response tensors.
2. **`MaterialModel`** (`core/material_model.py`, also a `BinaryManager`) turns a *model input* dict into a model. `analyze()` finds Wyckoff sites and bonds, builds the cell sites and bonds, and builds the atomic, cluster and combined SAMB. `save()`, `save_view()` and `save_pdf()` then write `<topdir>/<model>/<model>.pkl`, a `.qtdw` file and a PDF. The input keys and their defaults are in `core/default_model.py`, which is also what `mp_create -i` prints.
3. **`ModelAnalyzer`** (`core/model_analyzer.py`) reads a *control* dict (defaults in `core/default_control.py`) and runs in one of three modes:
   - `samb`: builds a tight-binding Hamiltonian from SAMB weights `z_j` on a saved model `.pkl`.
   - `wannier`: reads Wannier90 output (`.win`, `.nnkp`, `_hr.dat`, `_wsvec.dat`) through `util/util_wannier.py`.
   - `symcw`: combines the two modes above by decomposing the Wannier Hamiltonian onto SAMB.

   It then computes physical quantities: band dispersion along a seekpath k-path, DOS (Gaussian or tetrahedron), the Fermi energy from filling, fat bands, local operators, k-multipoles, and user callbacks (`callback_k_path` / `callback_k_full`). Results go to `<topdir>/<model>/output/`, including gnuplot scripts.

`core/cmd.py` (`create_model`, `analyze_model`) wraps layers 2 and 3. The click CLIs in `scripts/` are thin front ends to it.

Input and control files are Python files holding exactly one dict literal. They are read with `util.util.read_dict` and `read_dict_file`, never executed as code. Output dict files are formatted with black if it is installed (`do_black`). Missing keys are filled in from the defaults with `deep_update`, and unknown keys are rejected by `check_dict_keys`.

## Data, examples, and tests

- `docs/src/examples/` holds the worked example models (`*_in.py`, `*_ctrl.py`, and their generated outputs). `tests/test_model.py` uses the graphene example as a regression fixture.
- **Do not modify `multipie/binary_data/`.** That means no edits, no regeneration, and no deletion. It holds the group database, and only the maintainers change it. If a task seems to need a change there, stop and ask. The binaries and the group PDFs are generated by MultiPieGen, a separate repository, and are not built in this repository. `docs/src/group_doc.md` links to the Group Information of MultiPieGen; its URL is defined once at the bottom of that file.
- There is no `others/` in the repository. `others/uc/` and `others/doc_src/` are gitignored for local material only.
- Every source file starts with the MPL-2.0 SPDX header. Earlier versions, up to 2.3.4, were released under MIT.
