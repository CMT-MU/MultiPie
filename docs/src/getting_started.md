# Getting Started

## Tutorial

As a tutorial, we describe the procedure for generating model in case of graphene.

1. Download the contents of the [examples](https://github.com/CMT-MU/MultiPie/tree/main/docs/src/examples) directory.
2. Inside the `examples` directory, run the following:

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

     - `graphene_matrix.py`, `graphene_hr.dat` : Full matrix information for given selected SAMBs
     - `output/` : Various output for physical quantities
     - `samb/` : QtDraw files for SAMBs

4. In order to handle `graphene.pkl` interactively, use IPython interface. See in detail [analyze_model.ipynb](examples/analyze_model.ipynb)

## Model input file

The default values of model input file is the following:

```{literalinclude} examples/default_model.py
```

You can overwrite whatever you want. In case of graphene, the model input file is

```{literalinclude} examples/graphene_in.py
```

## Control file

The default values of control file is the following:

```{literalinclude} examples/default_control.py
```

You can overwrite whatever you want. In case of graphene, the control file is

```{literalinclude} examples/graphene_ctrl.py
```

The typical usage is the following: you first create the SAMBs for all irreps., and then select the desired SAMBs in order to create matrix elements, such as the symmetry-breaking terms in addition to the identity irreps.
