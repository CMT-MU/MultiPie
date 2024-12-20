"""
create energy plot.
"""

import os
import click
from multipie.scripts.plot_utility import plot_model


# ================================================== plot energy
@click.command()
@click.option("-v", "--verbose", is_flag=True, help="verbose off.")
@click.argument("params", nargs=-1)
def cmd(params, verbose):
    """
    create energy plot.

        PARAMS : parameter file names (w/ or w/o `.py`).
    """
    if len(params) < 1:
        exit()

    rel = os.path.relpath(".", "./")

    params = [i[:-3] if i[-3:] == ".py" else i for i in params]
    params = [os.path.join(rel, i + ".py") for i in params]

    plot_model(params, verbose=verbose)
