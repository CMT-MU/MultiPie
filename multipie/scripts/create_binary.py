"""
create binary files.
"""
import click
import multipie as mp


# ==================================================
def create_bin():
    mp.create_binary(clean=True, verbose=True)


# ==================================================
def read_bin():
    mp.get_binary(load=True, verbose=True)


# ================================================== create_binary
@click.command()
def cmd():
    """
    create binary files.
    """
    create_bin()
    read_bin()


# ================================================== create
if __name__ == "__main__":
    cmd()
