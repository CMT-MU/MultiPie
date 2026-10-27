# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Analyze model by control file.
"""

import click
from multipie.core.multipie_info import __top_dir__
from multipie.core.cmd import analyze_model
from multipie.scripts.mp_create import extract_dict, run_files

DEFAULT_CONTROL = __top_dir__ + "/multipie/core/default_control.py"


# ================================================== mp_analyze
@click.command()
@click.option("-v", "--verbose", is_flag=True, help="verbose on (and show full traceback for an error).")
@click.option("-i", "--input", "input", is_flag=True, help="show input format, and exit.")
@click.argument("controls", nargs=-1)
def cmd(controls, verbose, input):
    """
    Analyze model by control files (CONTROLS w or w/o '.py').
    """
    if input:
        input_str = "default_control = " + extract_dict(DEFAULT_CONTROL, "default_control")
        click.echo(input_str)
        return
    if len(controls) < 1:
        raise click.UsageError("Missing argument 'CONTROLS'.")

    # analyze model.
    run_files(analyze_model, controls, verbose)


# ==================================================
if __name__ == "__main__":
    cmd()
