"""
Analyze model by control file.
"""

import click
import sys
from multipie.core.multipie_info import __top_dir__
from multipie.core.cmd import analyze_model
from multipie.scripts.mp_create import extract_dict

DEFAULT_CONTROL = __top_dir__ + "/multipie/core/default_control.py"


# ================================================== mp_analyze
@click.command()
@click.option("-v", "--verbose", is_flag=True, help="verbose on.")
@click.option("-i", "--input", "input", is_flag=True, help="show input format, and exit.")
@click.argument("controls", nargs=-1)
def cmd(controls, verbose, input):
    """
    Analyze model by control files (CONTROLS w or w/o '.py').
    """
    if input:
        input_str = "default_control = " + extract_dict(DEFAULT_CONTROL, "default_control")
        click.echo(input_str)
        exit()
    if len(controls) < 1:
        click.echo("Usage: mp_analyze [OPTIONS] [CONTROLs]...")
        click.echo("Try 'mp_analyze --help' for help.\n")
        click.echo("Error: Missing argument 'CONTROLS'.")
        sys.exit(1)

    # analyze model.
    analyze_model(controls, verbose=verbose)
