"""
Create SAMB from model input.
"""

import click
from multipie.core.cmd import create_samb_matrix


# ================================================== mp_create_samb_matrix
@click.command()
@click.option("-v", "--verbose", is_flag=True, help="verbose off.")
@click.argument("models", nargs=-1)
def cmd(models, verbose):
    """
    Create SAMB matrix (and hr) from select input files (MODELS w or w/o '.py').
    """
    if len(models) < 1:
        click.echo("Usage: mp_create_samb_matrix [OPTIONS] [MODELS]...")
        click.echo("Try 'mp_create_samb_matrix --help' for help.\n")
        click.echo("Error: Missing argument 'MODELS'.")
        exit()

    # create all models.
    create_samb_matrix(models, verbose=verbose)
