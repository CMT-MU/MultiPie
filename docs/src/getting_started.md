# Getting Started

## Tutorial

As a tutorial, we describe the procedure for generating model in the case of **graphene**.

1. Download `graphene_in.py` and `graphene_ctrl.py` from [examples](https://github.com/CMT-MU/MultiPie/tree/main/docs/src/examples) directory.
2. At the same place where the above two files are, run the following:

    ```bash
    $ mp_create graphene_in.py
    ```

    It creates four files under `graphene` directory:

     - `graphene.pkl` : Model information file (binary)
     - `graphene.tex`, `graphene.pdf` : Summary of the model information
     - `graphene.qtdw` : QtDraw file of the model structure

3. To analyze the model, e.g., draw dispersion, run the following:

    ```bash
    $ mp_analyze graphene_ctrl.py
    ```

    It creates the following files under `graphene` directory:

     - `graphene_matrix.py`, `graphene_hr.dat` : Full matrix information for selected SAMBs
     - `output/` : Various output for physical quantities
     - `samb/` : QtDraw files for SAMBs

4. In order to handle `graphene.pkl` interactively, use IPython interface. See in detail [analyze_model.ipynb](examples/analyze_model.ipynb)

## Model input file

In case of graphene, the model input file is

```{literalinclude} examples/graphene_in.py
```

The other setting are provided as default values, which are given as follows:

```{literalinclude} examples/default_model.py
```

You can overwrite whatever you want as in the case of graphene.

## Control file

In case of graphene, the control file is

```{literalinclude} examples/graphene_ctrl.py
```

The default values of control file are provided as follows:

```{literalinclude} examples/default_control.py
```

The typical use of generating SAMBs, you first create the SAMBs for all irreps., and then choose the necessary SAMBs, such as the symmetry-breaking terms in addition to the identity irreps., by specifying `samb/select` and/or `samb/parameter` in the control file.
