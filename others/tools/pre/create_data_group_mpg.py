import os

header = """
Group information for all magnetic point groups.

- mpg_info.
    { MPG_id: field }
    MPG_id: MPG:{PG}.{no}.{ID} (PG:1-37, ID:1-155).
- mpg_id_set
    { all/crystal/PG/type: [MPG_id] }.
- mpg_tag.
    { MPG_id: tag }
- mpg_id.
    { tag: MPG_id }

- field.
    - tag: MPG_id (PG.no.ID).
    - international: international short symbol in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting.
    - PG: associated PG, unique.
    - SG: associated SG, MPG -> MSG -> SG.
    - MPG: = id.
    - MSG: associated MSG, 1st MSG in the same MPGs.
    - lattice: "0".
    - hexagonal_g: trigonal or hexagonal ?
    - type: type I, II, III.
    - SO: symmetry operations.
"""
import re
from multipie import MPGInfoType
from data_mpg import _mpg_info
from data_pg import _pg_info
from output.data_group_msg import msg_id_set


def replace_bar(s):
    return re.sub(r"-(\d)", r"\\bar{\1}", s)


def create_mpg():
    mpg_info = {}
    for MPG_id, (international, T_type, setting, SO) in _mpg_info.items():
        PG, no, ID = MPG_id.split(".")
        PG, no, ID = int(PG), int(no), int(ID)
        tag = MPG_id
        schoenflies = _pg_info[PG][0]
        crystal = _pg_info[PG][2]
        hex_g = crystal in ["trigonal", "hexagonal"]
        MSG = msg_id_set["MPG"][f"MPG:{MPG_id}"][0]
        SG = MSG.split(":")[1].split(".")[0]

        mpg_info[f"MPG:{MPG_id}"] = MPGInfoType(
            tag,
            replace_bar(international),
            schoenflies,
            crystal,
            setting,
            f"PG:{PG}",
            f"SG:{SG}",
            f"MPG:{MPG_id}",
            MSG,
            "0",
            hex_g,
            T_type,
            SO,
        )

    mpg_tag = {id_s: v.tag for id_s, v in mpg_info.items()}
    mpg_id = {v.tag: id_s for id_s, v in mpg_info.items()}

    # magnetic point group ID.
    mpg_id_set = {}
    mpg_id_set["all"] = list(mpg_info.keys())
    mpg_id_set["crystal"] = {}
    mpg_id_set["PG"] = {}
    mpg_id_set["type"] = {}
    for id_s, v in mpg_info.items():
        crystal = v.crystal
        pg = v.PG
        tp = v.type
        mpg_id_set["crystal"][crystal] = mpg_id_set["crystal"].get(crystal, []) + [id_s]
        mpg_id_set["PG"][pg] = mpg_id_set["PG"].get(pg, []) + [id_s]
        mpg_id_set["type"][tp] = mpg_id_set["type"].get(tp, []) + [id_s]

    cur_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(cur_dir, "output/data_group_mpg.py")
    with open(out_file, "w", encoding="utf-8") as f:
        print('"""' + header + '"""', file=f)
        print("from multipie import MPGInfoType", file=f)
        print("mpg_info =", {k: tuple(v) for k, v in mpg_info.items()}, file=f)
        print("mpg_id_set =", mpg_id_set, file=f)
        print("mpg_tag =", mpg_tag, file=f)
        print("mpg_id =", mpg_id, file=f)
        print("mpg_info = {k: MPGInfoType(*v) for k, v in mpg_info.items()}", file=f)
