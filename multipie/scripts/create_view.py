"""
create QtDraw file from model input.
"""
import click
from multipie.model.create_view import create_view
from multipie.model.material_model import input_str


# ================================================== create_view
@click.command()
@click.option("-v", "--verbose", is_flag=True, help="verbose off.")
@click.option("-f", "--formatter", is_flag=True, help="formatter off.")
@click.option("-m", "--mode", type=click.Choice(["arrow", "debug"]), default=None, help="view mode.")
@click.option("-i", "--input", is_flag=True, help="show input format, and exit.")
@click.argument("models", nargs=-1)
def cmd(models, verbose, formatter, mode, input):
    """
    create QtDraw files from input files.

        MODELS : base file names without `.py`.
    """
    if input:
        click.echo(input_str)
        exit()
    if len(models) < 1:
        exit()

    models = [i + ".py" for i in models]
    create_view(models, verbose=not verbose, formatter=not formatter, view_mode=mode)
