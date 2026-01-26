"""
Create SAMB qtdraw from model binary.
"""

import click
from multipie.core.cmd import create_samb_qtdraw


# ================================================== mp_create_samb_qtdraw
@click.command()
@click.option("-v", "--verbose", is_flag=True, help="verbose off.")
@click.argument("models", nargs=-1)
def cmd(models, verbose):
    """
    Create SAMB QtDraw files from model names (MODELS).
    """
    if len(models) < 1 or create_samb_qtdraw(models, verbose=verbose):
        click.echo("Usage: mp_create_samb_qtdraw [OPTIONS] [MODELS]...")
        click.echo("Try 'mp_create_samb_qtdraw --help' for help.\n")
        click.echo("Error: Missing argument 'MODELS'.")
        exit()
