# Example usage of full matrix

The example use of the output file, `_matrix.py` is given here.

1. Prepare parameter file. Example parameter files for C3v molecule and graphene are given as follows.
    - C3v molecule
    ```{literalinclude} ../examples/p_C3v.py
    ```
    - graphene
    ```{literalinclude} ../examples/p_graphene.py
    ```

2. At the folder where the parameter file exists, do the following to create energy plot.
See for more detail, try `create_plot --help` command.
    ```
    $ create_plot p_C3v p_graphene
    ```
1. The `.png` files are created in the current folder.
