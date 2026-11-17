# Models of MultiPie Ver.2

Models (`<model>/<model>.pkl`) created with MultiPie Ver.2, for the tests of `convert_model` (`tests/test_convert.py`).
They were created with the tags `v2.5.2` and `v2.3.4` of this repository, with `"pdf": {"create": False}` and
`"qtdraw": {"create": False}`:

| model | version | input |
|---|---|---|
| `graphene` | 2.5.2 | `{"model": "graphene", "group": 191, "cell": {"a": 2.435, "c": 10}, "site": {"C": ("[1/3,2/3,0]", "pz")}, "bond": [("C", "C", [1, 2, 3])], "SAMB_select": {"Gamma": "IR"}}` |
| `ch4` | 2.5.2 | `{"model": "ch4", "group": "Td", "site": {"C": ("[0,0,0]", ["s", "p"]), "H": ("[1/3,1/3,1/3]", "s")}, "bond": [("C", "H", 1), ("H", "H", 1)], "spinful": True}` |
| `pt_jml` | 2.5.2 | `{"model": "pt_jml", "group": "D4h", "site": {"Pt": ("[0,0,0]", ["(1/2,p)", "(3/2,p)"]), "O": ("[1,0,0]", "(1/2,s)")}, "bond": [("Pt", "O", 1)]}` |
| `ce_f` | 2.5.2 | `{"model": "ce_f", "group": "Oh", "site": {"Ce": ("[0,0,0]", "f")}, "bond": []}` |
| `ce_f_part` | 2.5.2 | `{"model": "ce_f_part", "group": "Oh", "site": {"Ce": ("[0,0,0]", ["f1", "f2", "f3x", "f3y"])}, "bond": []}` |
| `cubic_s` | 2.3.4 | `{"model": "cubic_s", "group": 221, "cell": {"a": 3.0}, "site": {"A": ("[0,0,0]", "s")}, "bond": [("A", "A", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])], "max_neighbor": 10, "SAMB_select": {"Gamma": "IR"}}` |

`cubic_s` was created with 2.3.4, which searched the bonds only in the cells within `search_cell_range`
(`(-2, 3, -2, 3, -2, 3)` by default). It missed a part of the shells beyond the 7th neighbor, so that its 9th neighbor is
the 11th neighbor of Ver.3 (and its 10th neighbor is not found).

`ce_f_part` has only a part of the f orbitals, which are not a complete set of irreps of the site symmetry, so that its
SAMB are not orthonormal (some of them vanish, and the others are linearly dependent).

`graphene_ctrl.py` and `cubic_s_ctrl.py` are controls (`samb` mode) with weights of the SAMB of Ver.2.
