# Getting Started

## Tutorial

As a tutorial, we describe the procedure for generating models and SAMBs.

1. Download the contents of the [examples](https://github.com/CMT-MU/MultiPie/tree/main/docs/examples) directory.
2. Inside the `examples` directory, run the `create_model.py` script:

    ```bash
    python -m create_model.py
    ```

    ```{literalinclude} examples/create_model.py
    ```

    This will create the `C3v` and `graphene` directories. Each directory will contain six files and a samb subdirectory.
    For `C3v`, the generated files are:

     - `C3v.pkl` : Model information file (binary)
     - `C3v.tex`, `C3v.pdf` : Summary of the model information
     - `C3v.qtdw` : QtDraw file of the model structure
     - `C3v_matrix.py` : Matrix element data of the SAMB
     - `C3v_hr.dat` : Numerical matrix element data with SAMB weight parameters applied
     - `samb` : QtDraw files for the SAMB

    The same set of files is generated for `graphene`.

3. For details on how to load and analyze the model information binary file, `.pkl`, see [read_model.ipynb](examples/read_model.ipynb)

## SAMB construction (detail)

### model input file

The default values of model input file is the following:

```{literalinclude} examples/default_model.py
```

You can overwrite whatever you want.

### model information files

Using the following input file, you can construct the model information files by running:

```bash
mp_create_samb C3v
```

This command generates `C3v.pkl`, `C3v.tex`, `C3v.pdf`, and `C3v.qtdw`.
Example input file:

```{literalinclude} examples/C3v.py
```

### selection and parameter file

The constructed SAMBs can be selected by selection file:

```{literalinclude} examples/selection_ex.py
```

The typical usage is the following: you first create the SAMBs for all irreps., and then select the desired SAMBs in order to create matrix elements with necessary bases, such as the symmetry-breaking terms in addition to the identity irreps.

### matrix element files

To generate full matrix file, `_matrix.py` by selection (and `_hr.dat` for given parameters) file, run:

```bash
mp_create_samb_matrix C3v_param
```

This command generates two files, `C3v_matrix.py`, and `C3v_hr.dat`.
Example input file:

```{literalinclude} examples/C3v_param.py
```

### QtDraw SAMB files

To generate QtDraw files for the constructed SAMB, run:

```bash
mp_create_samb_qtdraw C3v graphene
```

The generated SAMB QtDraw files are saved in the `samb` directory.
