import os

header = """
Group information for all magnetic space groups.

- msg_info.
    { MSG_id: field }
    MSG_id: MSG:{SG}.{no} (SG:1-230).
- msg_id_set
    { all/crystal/PG/SG/MPG/type: [MPG_id] }.
- msg_tag.
    { MSG_id: tag }
- msg_id.
    { tag: MSG_id }

- field.
    - tag: BNS_id (SG.no).
    - BNS: Belov-Neronova-Smirnova (international) notation in LaTeX.
    - OG: Opechowski-Guccione notation in LaTeX.
    - schoenflies: Schoenflies symbol in LaTeX.
    - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
    - setting: setting.
    - PG: associated PG, unique.
    - SG: associated SG, unique.
    - MPG: associated MPG, unique.
    - MSG: = id.
    - lattice: A, B, C, P, I, F, R.
    - hexagonal_g: trigonal or hexagonal ?
    - type: type I, II, III, IV.
    - SO: symmetry operations.
"""
import re
from multipie import MSGInfoType
from data_msg import _msg_info
from data_pg import _pg_info
from data_sg import _sg_info


def replace_bar(s):
    return re.sub(r"-(\d)", r"\\bar{\1}", s)


def create_msg():
    msg_info = {}
    for MSG_id, (BNS, T_type, OG, MPG_id, setting, SO) in _msg_info.items():
        SG, no = MSG_id.split(".")
        SG, no = int(SG), int(no)
        PG = int(MPG_id.split(".")[0])

        tag = MSG_id
        schoenflies = _sg_info[SG][0]
        crystal = _pg_info[PG][2]
        lattice = BNS[0]
        hex_g = crystal in ["trigonal", "hexagonal"]

        msg_info[f"MSG:{MSG_id}"] = MSGInfoType(
            tag,
            replace_bar(BNS),
            replace_bar(OG),
            schoenflies,
            crystal,
            setting,
            f"PG:{PG}",
            f"SG:{SG}",
            f"MPG:{MPG_id}",
            f"MSG:{MSG_id}",
            lattice,
            hex_g,
            T_type,
            SO,
        )
    msg_tag = {id_s: v.tag for id_s, v in msg_info.items()}
    msg_id = {v.tag: id_s for id_s, v in msg_info.items()}

    # space group ID.
    msg_id_set = {}
    msg_id_set["all"] = list(msg_info.keys())
    msg_id_set["crystal"] = {}
    msg_id_set["PG"] = {}
    msg_id_set["SG"] = {}
    msg_id_set["MPG"] = {}
    msg_id_set["type"] = {}
    for id_s, v in msg_info.items():
        crystal = v.crystal
        pg = v.PG
        sg = v.SG
        mpg = v.MPG
        tp = v.type
        msg_id_set["crystal"][crystal] = msg_id_set["crystal"].get(crystal, []) + [id_s]
        msg_id_set["PG"][pg] = msg_id_set["PG"].get(pg, []) + [id_s]
        msg_id_set["SG"][sg] = msg_id_set["SG"].get(sg, []) + [id_s]
        msg_id_set["MPG"][mpg] = msg_id_set["MPG"].get(mpg, []) + [id_s]
        msg_id_set["type"][tp] = msg_id_set["type"].get(tp, []) + [id_s]

    cur_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(cur_dir, "output/data_group_msg.py")
    with open(out_file, "w", encoding="utf-8") as f:
        print('"""' + header + '"""', file=f)
        print("from multipie import MSGInfoType", file=f)
        print("msg_info =", {k: tuple(v) for k, v in msg_info.items()}, file=f)
        print("msg_id_set =", msg_id_set, file=f)
        print("msg_tag =", msg_tag, file=f)
        print("msg_id =", msg_id, file=f)
        print("msg_info = {k: MSGInfoType(*v) for k, v in msg_info.items()}", file=f)
