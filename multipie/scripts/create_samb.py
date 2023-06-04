"""
create SAMB from model input.
"""
import os
import click
from multipie.model.create_model import create_model
from multipie.model.material_model import input_str


# ================================================== create_samb
@click.command()
@click.option("-p", "--parallel", is_flag=True, help="parallel off.")
@click.option("-v", "--verbose", is_flag=True, help="verbose off.")
@click.option("-l", "--latex", is_flag=True, help="LaTeX and PDF off.")
@click.option("-f", "--formatter", is_flag=True, help="formatter off.")
@click.option("-m", "--mode", type=click.Choice(["arrow", "debug"]), default=None, help="view mode.")
@click.option("-d", "--qtdraw", is_flag=True, help="QtDraw off.")
@click.option("-o", "--output", default=".", help="output top-directory name.")
@click.option("-i", "--input", is_flag=True, help="show input format, and exit.")
@click.argument("models", nargs=-1)
def cmd(models, parallel, verbose, latex, formatter, mode, qtdraw, output, input):
    """
    create SAMB from input files.

        MODELS : file names (w or w/o `.py`).
    """
    if input:
        click.echo(input_str)
        exit()
    if len(models) < 1:
        exit()

    rel = os.path.relpath(".", "./" + output)

    models = [i[:-3] if i[-3:] == ".py" else i for i in models]
    models = [os.path.join(rel, i + ".py") for i in models]
    create_model(
        models,
        topdir=output,
        symbolic=True,
        parallel=not parallel,
        verbose=not verbose,
        pdf=not latex,
        formatter=not formatter,
        view_mode=mode,
        qtdraw=not qtdraw,
    )
