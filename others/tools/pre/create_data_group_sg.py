import os

header = """
Group information for all space groups.

- sg_info.
    { SG_id: field }
    SG_id: SG:{ID} (ID:1-230).
- sg_id_set.
    { all/crystal/PG: [SG_id]}
- sg_tag.
    { SG_id: tag }
- sg_id.
    { tag: SG_id }

- field.
    - tag: Schoenflies symbol in text.
    - international: international short symbol in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting comment.
    - PG: associated PG, unique.
    - SG: = id.
    - MPG: associated MPG, SG -> MSG -> MPG.
    - MSG: associated MSG, 1st MSG with type II in the same SGs.
    - lattice: A, B, C, P, I, F, R.
    - hexagonal_g: trigonal or hexagonal ?
    - SO: symmetry operations.
    - SO_gen: generator of SOs.
"""
import re
from multipie import SGInfoType
from data_sg import _sg_info
from data_pg import _pg_info
from data_msg import _msg_info
from output.data_group_msg import msg_id_set


def replace_bar(s):
    return re.sub(r"-(\d)", r"\\bar{\1}", s)


def create_sg():
    sg_info = {}
    for SG_id, (schoenflies, international, setting, PG_id, SO, generator) in _sg_info.items():
        tag = schoenflies
        s = schoenflies.split("^")
        S_latex = s[0][0] + r"_{\rm " + s[0][1:] + "}" + "^{" + s[1] + "}"
        crystal = _pg_info[PG_id][2]
        lattice = international[0]
        hex_g = crystal in ["trigonal", "hexagonal"]
        MSG_lst = msg_id_set["SG"][f"SG:{SG_id}"]
        MSG = next(x for x in MSG_lst if x in msg_id_set["type"]["II"])
        MPG = _msg_info[MSG.split(":")[1]][3]

        sg_info[f"SG:{SG_id}"] = SGInfoType(
            tag,
            replace_bar(international),
            S_latex,
            crystal,
            setting,
            f"PG:{PG_id}",
            f"SG:{SG_id}",
            f"MPG:{MPG}",
            MSG,
            lattice,
            hex_g,
            SO,
            generator,
        )

    sg_tag = {id_s: v.tag for id_s, v in sg_info.items()}
    sg_id = {v.tag: id_s for id_s, v in sg_info.items()}

    # space group ID.
    sg_id_set = {}
    sg_id_set["all"] = list(sg_info.keys())
    sg_id_set["crystal"] = {}
    PG_dict = {}
    for id_s, v in sg_info.items():
        crystal = v.crystal
        pg = v.PG
        sg_id_set["crystal"][crystal] = sg_id_set["crystal"].get(crystal, []) + [id_s]
        PG_dict[pg] = PG_dict.get(pg, []) + [id_s]
    sg_id_set["PG"] = PG_dict

    cur_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(cur_dir, "output/data_group_sg.py")
    with open(out_file, "w", encoding="utf-8") as f:
        print('"""' + header + '"""', file=f)
        print("from multipie import SGInfoType", file=f)
        print("sg_info =", {k: tuple(v) for k, v in sg_info.items()}, file=f)
        print("sg_id_set =", sg_id_set, file=f)
        print("sg_tag =", sg_tag, file=f)
        print("sg_id =", sg_id, file=f)
        print("sg_info = {k: SGInfoType(*v) for k, v in sg_info.items()}", file=f)
