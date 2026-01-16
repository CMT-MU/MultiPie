import os

header = """
Group information for all point groups.

- pg_info.
    { PG_id: field }
    PG_id: PG:{ID} (ID:1-37; complex:38-47).
- pg_id_set.
    { all/crystal/complex/irrep: [PG_id]}
    irrep: 1d, 2d, 3d.
- pg_tag.
    { PG_id: tag }
- pg_id.
    { tag: PG_id }

- field.
    - tag: Schoenflies in text.
    - international: international short symbol in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting.
    - PG: = id.
    - SG: associated SG, 1st SG in the same PGs.
    - MPG: associated MPG, 1st MPG with type II in the same PGs.
    - MSG: associated MSG, PG -> MPG -> MSG.
    - lattice: "0".
    - hexagonal_g: trigonal or hexagonal ?
    - SO: symmetry operations.
    - SO_gen: generator of SOs.
"""
import re
from multipie import PGInfoType
from data_pg import _pg_info, _sg_complex_to_real
from output.data_group_sg import sg_id_set
from output.data_group_mpg import mpg_id_set
from output.data_group_msg import msg_id_set


def replace_bar(s):
    return re.sub(r"-(\d)", r"\\bar{\1}", s)


def create_pg():
    pg_info = {}
    for PG_id, (schoenflies, international, crystal, setting, SO, generator) in _pg_info.items():
        s = schoenflies.split("-")
        S_latex = s[0][0] + r"_{\rm " + s[0][1:] + "}"
        if len(s) > 1:
            S_latex += f"({s[1]})"
        tag = schoenflies
        hex_g = crystal in ["trigonal", "hexagonal"]
        PG = f"PG:{_sg_complex_to_real.get(PG_id, PG_id)}"
        SG = sg_id_set["PG"][PG][0]
        MPG_lst = mpg_id_set["PG"][PG]
        MPG = next(x for x in MPG_lst if x in mpg_id_set["type"]["II"])
        MSG = msg_id_set["MPG"][MPG][0]

        pg_info[f"PG:{PG_id}"] = PGInfoType(
            tag, replace_bar(international), S_latex, crystal, setting, f"PG:{PG_id}", SG, MPG, MSG, "0", hex_g, SO, generator
        )

    pg_tag = {id_s: v.tag for id_s, v in pg_info.items()}
    pg_id = {v.tag: id_s for id_s, v in pg_info.items()}

    # point group ID.
    pg_id_set = {
        # point group (real).
        "all": [f"PG:{no}" for no in range(1, 38)],
        "crystal": {},
        # point group with complex character.
        "complex": [f"PG:{no}" for no in range(38, 48)],
        "irrep": {
            # point group wiht 1d irrep.
            "1d": [f"PG:{no}" for no in [1, 2, 3, 4, 5, 6, 7, 8, 38, 39, 40, 41, 42, 43, 44, 45]],
            # point group with 1d and 2d irreps.
            "2d": [
                f"PG:{no}"
                for no in [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 33, 34, 35, 36, 37]
            ],
            # point group with 1d, 2d, and 3d irreps.
            "3d": [f"PG:{no}" for no in [28, 29, 30, 31, 32, 46, 47]],
        },
    }
    for id_s, v in pg_info.items():
        crystal = v.crystal
        pg_id_set["crystal"][crystal] = pg_id_set["crystal"].get(crystal, []) + [id_s]

    cur_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(cur_dir, "output/data_group_pg.py")
    with open(out_file, "w", encoding="utf-8") as f:
        print('"""' + header + '"""', file=f)
        print("from multipie import PGInfoType", file=f)
        print("pg_info =", {k: tuple(v) for k, v in pg_info.items()}, file=f)
        print("pg_id_set =", pg_id_set, file=f)
        print("pg_tag =", pg_tag, file=f)
        print("pg_id =", pg_id, file=f)
        print("pg_info = {k: PGInfoType(*v) for k, v in pg_info.items()}", file=f)
