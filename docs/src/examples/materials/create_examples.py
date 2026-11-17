# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create the examples of materials for the documentation (docs/src/materials.md).

    python docs/src/examples/materials/create_examples.py

For each material, the model input in tests/data/materials (also used by tests/test_materials.py) and the control in
this directory give:
    - <name>.pdf: the summary of the model with all the SAMB (linked from materials.md).
    - <name>_samb.md: the table of the SAMB (included in materials.md).
    - ../../figures/materials/<name>_dispersion.png: the band dispersion.
"""

import os
import shutil
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from multipie.core.cmd import create_model
from multipie.core.model_analyzer import ModelAnalyzer
from multipie.util.util import read_dict
from multipie.util.util_model_analyzer import k_label

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir, os.pardir, os.pardir))
MATERIAL_DIR = os.path.join(TOP, "tests", "data", "materials")
FIGURE_DIR = os.path.join(HERE, os.pardir, os.pardir, "figures", "materials")
MATERIALS = ["SrVO3", "MoS2", "Te"]


# ==================================================
def samb_table(mm):
    """
    Markdown table of the combined SAMB of a model.
    """
    rows = ["| $z_j$ | site/bond | Wyckoff | neighbor | bra-ket | symbol |", "|:--|:--|:--|:--|:--|:--|"]
    for z, (symbol, u, tag, comp) in mm["combined_id"].items():
        st = u.samb_type
        bk = st.bk_info
        kind = st.head if st.head == st.tail and "@" not in st.wyckoff else f"{st.head}-{st.tail}"
        neighbor = "-" if u.neighbor == 0 else str(u.neighbor)
        braket = f"{'spdf'[bk.bh_rank]}-{'spdf'[bk.kt_rank]}"
        rows.append(f"| {z} | {kind} | `{st.wyckoff}` | {neighbor} | {braket} | ${symbol}$ |")
    return "\n".join(rows) + "\n"


# ==================================================
def plot_dispersion(name, ma, ylim=None, color=None):
    """
    Band dispersion in PNG (with the weight of the first local operator as color, if any).
    """
    data = np.loadtxt(os.path.join(name, "output", f"{name}_dispersion.txt"))
    k_dis_pos = ma._k_path_info["k_dis_pos"]
    n_band = len(ma.get_ket())
    k, E = data[:, 0].reshape(n_band, -1), data[:, 1].reshape(n_band, -1)

    fig, ax = plt.subplots(figsize=(6, 3.6))
    if color is not None and data.shape[1] > 2:
        c = data[:, 2].reshape(n_band, -1)
        sc = None
        for i in range(n_band):
            sc = ax.scatter(k[i], E[i], c=c[i], cmap="coolwarm", vmin=-1, vmax=1, s=2, lw=0)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
        cbar.set_label(color)
    else:
        for i in range(n_band):
            ax.plot(k[i], E[i], color="tab:blue", lw=1)
    for x in k_dis_pos:
        ax.axvline(x, color="black", lw=0.5)
    ax.set_xticks(list(k_dis_pos))
    ax.set_xticklabels([k_label(v, "latex") for v in k_dis_pos.values()])
    ax.set_xlim(0, k.max())
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_ylabel("Energy")
    fig.tight_layout()
    os.makedirs(FIGURE_DIR, exist_ok=True)
    fig.savefig(os.path.join(FIGURE_DIR, f"{name}_dispersion.png"), dpi=150)
    plt.close(fig)


# ==================================================
def main():
    for name in MATERIALS:
        print(f"* {name}")
        model = read_dict(os.path.join(MATERIAL_DIR, f"{name}_in.py"))
        model["qtdraw"] = {"create": False}
        control = read_dict(os.path.join(HERE, f"{name}_ctrl.py"))
        local = control.get("output", {}).get("dispersion", {}).get("local", [])

        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as workdir:  # the model and the analysis (only the results are kept).
            try:
                os.chdir(workdir)
                create_model(model, topdir=workdir)
                ma = ModelAnalyzer(workdir)
                ma.analyze(control)
                plot_dispersion(name, ma, color=local[0] if local else None)
            finally:
                os.chdir(cwd)
            shutil.copy2(os.path.join(workdir, name, f"{name}.pdf"), os.path.join(HERE, f"{name}.pdf"))

        with open(os.path.join(HERE, f"{name}_samb.md"), "w", encoding="utf-8") as f:
            f.write(samb_table(ma.model))
    return 0


if __name__ == "__main__":
    sys.exit(main())
