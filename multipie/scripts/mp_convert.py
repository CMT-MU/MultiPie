# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Convert a model (.pkl) of MultiPie Ver.2 into Ver.3.
"""

import click
from multipie.core.cmd import convert_model


# ==================================================
def report_text(report):
    """
    Report of the conversion as text.

    Args:
        report (dict): report of convert_model.

    Returns:
        - (str) -- report.
    """
    lst = [f"model '{report['model']}':"]
    if report["input"] is not None:
        lst.append(f"  model input of Ver.3: {report['input']}")
        lst.append(f"  model of Ver.2 kept as: {report['model_v2']}")
    for site, ren in report["renamed_orbital"].items():
        lst.append(f"  renamed orbitals of {site}: " + ", ".join(f"{a} -> {b}" for a, b in ren))
    for key, value in report["removed_key"].items():
        lst.append(f"  removed key: {key} = {value}")
    for (tail, head), (old, new) in report["neighbor"].items():
        lst.append(f"  neighbors of {tail};{head} (chosen by the distances of the bonds of Ver.2): {old} -> {new}")
    n2, n3 = report["SAMB"]
    lst.append(f"  SAMB: {n2} (Ver.2) -> {n3} (Ver.3), residual of the projection {report['span_residual']:.2g}")
    for c, r in report["control"].items():
        p2, p3 = r["parameter"]
        lst.append(f"  control '{c}' -> '{r['output']}': {p2} -> {p3} weights, residual of H {r['residual']:.2g}")
    return "\n".join(lst)


# ================================================== mp_convert
@click.command()
@click.option("-v", "--verbose", is_flag=True, help="verbose on (and show full traceback for an error).")
@click.option("--pdf/--no-pdf", default=None, help="create PDF ? [default: as in the model of Ver.2]")
@click.option("--qtdraw/--no-qtdraw", default=None, help="create QtDraw file ? [default: as in the model of Ver.2]")
@click.argument("model")
@click.argument("controls", nargs=-1)
def cmd(model, controls, verbose, pdf, qtdraw):
    """
    Convert MODEL (MODEL/MODEL.pkl) of MultiPie Ver.2 into Ver.3, and the SAMB weights of CONTROLS (mp_analyze).

    \b
    - MODEL_in.py: model input rebuilt for Ver.3.
    - MODEL/MODEL.pkl: model created with Ver.3 (the model of Ver.2 is kept as MODEL/MODEL_v2.pkl).
    - CONTROL_v3.py: control with the weights z_j converted to the SAMB of Ver.3.
    """
    try:
        report = convert_model(model, list(controls), pdf=pdf, qtdraw=qtdraw, verbose=verbose)
    except Exception as e:
        if verbose:
            raise
        msg = f"{model}: {type(e).__name__}: {e}"
        for note in getattr(e, "__notes__", []):
            msg += f"\n  {note}"
        msg += "\n  (use -v to show the full traceback.)"
        raise click.ClickException(msg) from e
    click.echo(report_text(report))


# ==================================================
if __name__ == "__main__":
    cmd()
