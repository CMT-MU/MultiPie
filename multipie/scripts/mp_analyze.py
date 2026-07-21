"""
Analyze model by control file.
"""

import click
import sys
from multipie.core.multipie_info import __top_dir__
from multipie.core.cmd import analyze_model
from multipie.scripts.mp_create import extract_dict

DEFAULT_CONTROL = __top_dir__ + "/multipie/core/default_control.py"


# ==================================================
def parse_grid(ctx, param, value):
    if value is None:
        return None
    try:
        parts = [int(v) for v in value.replace(" ", "").split(",")]
    except ValueError:
        raise click.BadParameter("grid values must be integers, e.g. '-g 1,2,3'.")
    if not (1 <= len(parts) <= 3):
        raise click.BadParameter("grid must have 1 to 3 values, e.g. '-g 1,2,3'.")
    return tuple(parts)


# ================================================== mp_analyze
@click.command()
@click.option("-v", "--verbose", is_flag=True, help="verbose on.")
@click.option("-i", "--input", "input", is_flag=True, help="show input format, and exit.")
@click.option("-g", "--grid", callback=parse_grid, default=None, help="grid points, N1[,N2[,N3]].")
@click.argument("controls", nargs=-1)
def cmd(controls, verbose, input, grid):
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
    if grid is not None:
        analyze_model(controls, *grid, verbose=verbose)
    else:
        analyze_model(controls, verbose=verbose)
