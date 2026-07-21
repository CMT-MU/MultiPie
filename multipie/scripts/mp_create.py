"""
Create model by model input.
"""

import click
from multipie.core.multipie_info import __top_dir__
from multipie.core.cmd import create_model

DEFAULT_MODEL = __top_dir__ + "/multipie/core/default_model.py"


# ==================================================
def extract_dict(filepath, key):
    with open(filepath, encoding="utf-8") as f:
        text = f.read()

    start = text.find(key)
    if start == -1:
        return None

    start_brace = text.find("{", start)
    if start_brace == -1:
        return None

    count = 0
    for i in range(start_brace, len(text)):
        if text[i] == "{":
            count += 1
        elif text[i] == "}":
            count -= 1
            if count == 0:
                return text[start_brace : i + 1]

    return None


# ================================================== mp_create
@click.command()
@click.option("-v", "--verbose", is_flag=True, help="verbose on.")
@click.option("-i", "--input", is_flag=True, help="show input format, and exit.")
@click.argument("models", nargs=-1)
def cmd(models, verbose, input):
    """
    Create models by input files (MODELS w or w/o '.py').
    """
    if input:
        input_str = "default_model = " + extract_dict(DEFAULT_MODEL, "default_model")
        click.echo(input_str)
        exit()
    if len(models) < 1:
        click.echo("Usage: mp_create [OPTIONS] [MODELS]...")
        click.echo("Try 'mp_create --help' for help.\n")
        click.echo("Error: Missing argument 'MODELS'.")
        exit()

    # create all models.
    create_model(models, verbose=verbose)
