"""
Create model by model input.
"""

import click
from multipie import __version__
from multipie.core.multipie_info import __top_dir__
from multipie.core.cmd import create_model
from multipie.util.util import read_dict_file

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


# ==================================================
def run_files(func, files, verbose):
    """
    Read all files, and run func for each file. An error is reported in one message.

    Args:
        func (function): create_model or analyze_model.
        files ([str]): input files.
        verbose (bool): verbose ? (show full traceback for an error).

    Raises:
        click.ClickException: if an error occurs (exit status 1).

    Note:
        - all files are read before func is called, so that no output is written when one of the files cannot be read.
    """

    def run(file, func, *args):
        try:
            return func(*args)
        except Exception as e:
            if verbose:
                raise
            msg = f"{file}: {type(e).__name__}: {e}"
            for note in getattr(e, "__notes__", []):
                msg += f"\n  {note}"
            msg += "\n  (use -v to show the full traceback.)"
            raise click.ClickException(msg) from e

    data = [(file, run(file, read_dict_file, file, None, verbose)) for file in files]
    for file, dic in data:
        run(file, func, list(dic.values()), None, verbose)


# ================================================== mp_create
@click.command()
@click.option("-v", "--verbose", is_flag=True, help="verbose on (and show full traceback for an error).")
@click.option("-i", "--input", is_flag=True, help="show input format, and exit.")
@click.argument("models", nargs=-1)
def cmd(models, verbose, input):
    """
    Create models by input files (MODELS w or w/o '.py').
    """
    if input:
        input_str = "default_model = " + extract_dict(DEFAULT_MODEL, "default_model")
        input_str = input_str.replace("__version__", repr(__version__))
        click.echo(input_str)
        return
    if len(models) < 1:
        raise click.UsageError("Missing argument 'MODELS'.")

    # create all models.
    run_files(create_model, models, verbose)


# ==================================================
if __name__ == "__main__":
    cmd()
