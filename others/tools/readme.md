# Creation of group database and PDF

The `tools` folder is to create various binary database and PDF for crystallographic point and space groups.

## Create core data

1. Run `create_data_info.py` and `create_data_group.py` in `create_binary`.
2. In `pre`, run `create_data_group.py` and `create_wyckoff_bond_data.py`. Then, the output files are created in `output`.
3. Move the created files in `output` to `data`.

## Create binary data

1. Run `create_data.py` or the files in `create_binary` one by one.
2. Comment out `from .core.group import Group` in `__init__.py`, and Run `create_final_data.py`. Then, uncomment it again.
3. The binary files are created in `tools/binary_data` as temporary ones.
4. Run `create_final_data.py`, and the final binary files in `multipie/binary_data`.

In the above procedures, files in `data` are used.
The binary creation codes refer utility functions in `utils`.

## Create PDF

1. Run `create_pdf.py` or the files in `create_pdf` one by one.

## Binary Data

The temporary binary data is created in `tools/binary_data/`, and the final ones in `multipie/binary_data`.

The temporary binary files are the following:

- `info`: information for all crystallographic and magnetic point and space group.
- `group`: symmetry operation, character, Wyckoff position for all crystallographic and magnetic point and space group.
- `harmonics`: harmonics for all point group.
- `harmonics_multipole`: multipolar harmonics for all point group.
- `harmonics_spherical`: spherical harmonics.
- `harmonics_root_cluster`: root-cluster harmonics.
- `atomic_multipole`: matrix element of atomic multipole.
- `atomic_multipole_group`: matrix element of atomic multipole for all point group.
- `cluster_samb`: cluster SAMB for all point and space group.
- `representation_matrix`: representation matix for all point group.
- `root_cluster`: information for root cluster and site mapping.
- `symmetry_operation_matrix`: symmetry operation matrix for point-group multiopole.

## PDF

PDF is created in `docs/pdf/`.

- `PG`: data for each point group.
- `SG`: data for each space group.
- `MPG`: data for each magnetic point group.
- `MSG`: data for each magnetic space group.
- `info`: overall information.
- `misc`: detailed overall information.

Technical documents are in `docs/technical_note`.

## Test & Check

- The test and check codes are in `tests` and `check`.
- To check temporary binary files, use codes in `check_binary`.
