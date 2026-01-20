"""
Create group name list.
Created "output/data_group_name_list.py" should be moved to "multipie/data".
"""

import os
from multipie import Group


# ==================================================
header = """
Group name list for all magnetic space groups.

id : (dir_name, md_name, no, tag, name, PG_id, SG_id, MPG_id, MSG_id).
"""


# ==================================================
def create_group_name_list():
    g_info = Group.global_info()
    dic = {"PG": {}, "SG": {}, "MPG": {}, "MSG": {}}

    for id_s in g_info["id_set"]["PG"]["all"] + g_info["id_set"]["PG"]["complex"]:
        g = Group(id_s)
        tag, no = g.info.tag, int(id_s.split(":")[1])
        name = f"${g.latex()}$ (${g.info.international}$)"
        d = f"{no:03d}-{tag}"
        md = f"No_{no}.md"

        dic["PG"][id_s] = (d, md, no, tag, name, g.info.PG, g.info.SG, g.info.MPG, g.info.MSG)

    for id_s in g_info["id_set"]["SG"]["all"]:
        g = Group(id_s)
        tag, no = g.info.tag.replace("^", "_"), int(id_s.split(":")[1])
        name = f"${g.latex()}$ (${g.info.international}$)"
        d = f"{no:03d}-{tag}"
        md = f"No_{no}.md"

        dic["SG"][id_s] = (d, md, no, tag, name, g.info.PG, g.info.SG, g.info.MPG, g.info.MSG)

    for id_s in g_info["id_set"]["MPG"]["all"]:
        g = Group(id_s)
        tag, no = g.info.tag, id_s.split(":")[1]
        name = f"${g.latex()}$"
        d = f"{no}"
        md = f"No_{no}.md"

        dic["MPG"][id_s] = (d, md, no, tag, name, g.info.PG, g.info.SG, g.info.MPG, g.info.MSG)

    for id_s in g_info["id_set"]["MSG"]["all"]:
        g = Group(id_s)
        tag, no = g.info.tag, id_s.split(":")[1]
        name = f"${g.info.BNS}$"
        d = f"{no}"
        md = f"No_{no}.md"

        dic["MSG"][id_s] = (d, md, no, tag, name, g.info.PG, g.info.SG, g.info.MPG, g.info.MSG)

    return dic


# ==================================================
dic = create_group_name_list()
cur_dir = os.path.dirname(os.path.abspath(__file__))
out_file = os.path.join(cur_dir, "output/data_group_name_list.py")
with open(out_file, "w", encoding="utf-8") as f:
    print('"""' + header + '"""', file=f)
    print("group_name_list =", dic, file=f)
