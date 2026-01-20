"""
Create markdown for all information.
"""

import os
import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
from multipie.core.group import Group
from multipie.util.util import timer
from multipie.util.util_binary import BinaryManager
from others.tools.data.data_group_name_list import group_name_list


# ==================================================
def dict_to_markdown(d, indent=0):
    md = ""
    for key, value in d.items():
        mk = "**" if indent == 0 else ""
        prefix = "  " * indent + f"- {mk}{key}{mk}"
        if isinstance(value, dict):
            md += f"{prefix}\n" + dict_to_markdown(value, indent + 1)
        else:
            md += f"{prefix} : {value}\n"
    return md


# ==================================================
def list_to_table(data, cols=4, header=True, tag="tag"):
    """Return markdown table from list."""
    if header:
        data = [tag] + [""] * (cols - 1) + data
        header_row = "| " + " | ".join(data[:cols]) + " |"
        separator = "|" + ":---|" * cols
        rows = []
        for i in range(cols, len(data), cols):
            row = "| " + " | ".join(data[i : i + cols]) + " |"
            rows.append(row)
        return "\n".join([header_row, separator] + rows)
    else:
        # append rows only (no header)
        while len(data) % cols != 0:
            data.append("")
        rows = []
        for i in range(0, len(data), cols):
            row = "| " + " | ".join(data[i : i + cols]) + " |"
            rows.append(row)
        return "\n".join(rows)


# ==================================================
@timer
def create_markdown():
    info = BinaryManager("info")
    dir_pdf = "../../others/pdf/"
    dir_pdf2 = "../../../others/pdf/"
    dir_md = os.path.join(__top_dir__, "docs/src")

    os.makedirs(os.path.join(dir_md, "PG"), exist_ok=True)
    os.makedirs(os.path.join(dir_md, "SG"), exist_ok=True)
    os.makedirs(os.path.join(dir_md, "MPG"), exist_ok=True)
    os.makedirs(os.path.join(dir_md, "MSG"), exist_ok=True)

    # ---------- Overview ----------
    overview_dict = {
        "Group information": {
            "Point group": f"\n[PDF]({dir_pdf}info/PG.pdf)",
            "Space group": f"\n[PDF]({dir_pdf}info/SG.pdf)",
            "Magnetic point group": f"\n[PDF]({dir_pdf}info/MPG.pdf)",
            "Magnetic space group": f"\n[PDF]({dir_pdf}info/MSG.pdf)",
            "Compatibility relation": f"\n[PDF]({dir_pdf}info/compatibility_relation.pdf)",
            "Active multipole": f"\n[PDF]({dir_pdf}info/active_multipole.pdf)",
        },
        "Atomic multipole": {
            "JML basis": f"[PDF]({dir_pdf}misc/atomic_multipole_jml.pdf)",
            "LM basis": f"[PDF]({dir_pdf}misc/atomic_multipole_lm.pdf)",
        },
        "Harmonics (in terms of tesseral harmonics)": {
            r"$D_{\rm 6h}$": {
                "polar": f"[PDF]({dir_pdf}misc/harmonics_D6h_polar.pdf)",
                "axial": f"[PDF]({dir_pdf}misc/harmonics_D6h_axial.pdf)",
            },
            r"$O_{\rm h}$": {
                "polar": f"[PDF]({dir_pdf}misc/harmonics_Oh_polar.pdf)",
                "axial": f"[PDF]({dir_pdf}misc/harmonics_Oh_axial.pdf)",
            },
        },
    }

    overview_md = dict_to_markdown(overview_dict)

    # ---------- Helper ----------
    def write_group_markdown(group_dict, dir_out, title, md):
        md_content = title + "\n" + dict_to_markdown(group_dict) + "\n"
        out_path = os.path.join(dir_md, dir_out, md)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(md_content)

    # ---------- Containers ----------
    point_tbl, complex_tbl, space_tbl, mpg_tbl, msg_tbl = [], [], [], [], []

    # ---------- Point group (real) ----------
    for id_s, (d, md, no, tag, name, pg, sg, mpg, msg) in group_name_list["PG"].items():
        if no > 37:
            continue

        point_tbl.append(f"{no}:  [{name}](PG/{md})")

        xl = ["Q", "G"]
        mno = group_name_list["MPG"][mpg][2]  # no.
        title = f"# PG #{no}:  {name}\n"
        dd = f"{dir_pdf2}PG/{d}"

        sg_md = group_name_list["SG"][sg][1]
        mpg_md = group_name_list["MPG"][mpg][1]
        sg_name = group_name_list["SG"][sg][4]
        mpg_name = group_name_list["MPG"][mpg][4]

        pg_dict = {
            "Related group": {
                "SG": f"[{sg_name}](../SG/{sg_md})",
                "MPG": f"[{mpg_name}](../MPG/{mpg_md})",
                "MSG (SG)": f"[{sg_name}](../MSG/{sg_md})",
            },
            "Symmetry operation": f"[PDF]({dd}/symmetry_operation.pdf)",
            "Character table": f"[PDF]({dd}/character_table.pdf)",
            "Wyckoff": {
                "site": f"[PDF]({dd}/wyckoff_site.pdf)",
                "bond": f"[PDF]({dd}/wyckoff_bond.pdf)",
            },
            "Atomic multipole": {
                "JML basis": f"[PDF]({dd}/atomic_multipole_jml.pdf)",
                r"L$\Gamma\sigma$ basis": f"[PDF]({dd}/atomic_multipole_lgs.pdf)",
                r"L$\Gamma$ basis": f"[PDF]({dd}/atomic_multipole_lg.pdf)",
            },
            "Harmonics": {
                "polar": f"[PDF]({dd}/harmonics_polar.pdf)",
                "axial": f"[PDF]({dd}/harmonics_axial.pdf)",
            },
            "Multipolar Harmonics (internal)": {
                "dipolar internal polar(Q)/axial(G) varialble": {
                    "polar (Q)": f"[PDF]({dd}/harmonics_s1_polar_q.pdf)",
                    "axial (Q)": f"[PDF]({dd}/harmonics_s1_axial_q.pdf)",
                    "polar (G)": f"[PDF]({dd}/harmonics_s1_polar_g.pdf)",
                    "axial (G)": f"[PDF]({dd}/harmonics_s1_axial_g.pdf)",
                },
                "quadrupolar internal polar(Q)/axial(G) variable": {
                    "polar (Q)": f"[PDF]({dd}/harmonics_s2_polar_q.pdf)",
                    "axial (Q)": f"[PDF]({dd}/harmonics_s2_axial_q.pdf)",
                    "polar (G)": f"[PDF]({dd}/harmonics_s2_polar_g.pdf)",
                    "axial (G)": f"[PDF]({dd}/harmonics_s2_axial_g.pdf)",
                },
                "octupolar internal polar(Q)/axial(G) variable": {
                    "polar (Q)": f"[PDF]({dd}/harmonics_s3_polar_q.pdf)",
                    "axial (Q)": f"[PDF]({dd}/harmonics_s3_axial_q.pdf)",
                    "polar (G)": f"[PDF]({dd}/harmonics_s3_polar_g.pdf)",
                    "axial (G)": f"[PDF]({dd}/harmonics_s3_axial_g.pdf)",
                },
            },
            "Response tensor": {X: f"[PDF]({dir_pdf2}MPG/{mno}/response_tensor_" + X + ".pdf)" for X in xl},
        }

        write_group_markdown(pg_dict, "PG", title, md)

    # ---------- Point group (complex) ----------
    for id_s, (d, md, no, tag, name, pg, sg, mpg, msg) in group_name_list["PG"].items():
        if no < 38:
            continue

        complex_tbl.append(f"{no}:  [{name}](PG/{md})")

        title = f"# PG #{no}:  {name}\n"
        dd = f"{dir_pdf2}PG/{d}"

        sg_md = group_name_list["SG"][sg][1]
        mpg_md = group_name_list["MPG"][mpg][1]
        sg_name = group_name_list["SG"][sg][4]
        mpg_name = group_name_list["MPG"][mpg][4]

        cg_dict = {
            "Related group": {
                "SG": f"[{sg_name}](../SG/{sg_md})",
                "MPG": f"[{mpg_name}](../MPG/{mpg_md})",
                "MSG (SG)": f"[{sg_name}](../MSG/{sg_md})",
            },
            "Symmetry operation": f"[PDF]({dd}/symmetry_operation.pdf)",
            "Character table": f"[PDF]({dd}/character_table.pdf)",
            "Harmonics": {
                "polar": f"[PDF]({dd}/harmonics_polar.pdf)",
                "axial": f"[PDF]({dd}/harmonics_axial.pdf)",
            },
            "Representation matrix": f"[PDF]({dd}/representation_matrix.pdf)",
        }

        write_group_markdown(cg_dict, "PG", title, md)

    # ---------- Space group ----------
    for id_s, (d, md, no, tag, name, pg, sg, mpg, msg) in group_name_list["SG"].items():
        space_tbl.append(f"{no}:  [{name}](SG/{md})")

        title = f"# SG #{no}:  {name}\n"
        dd = f"{dir_pdf2}SG/{d}"

        pg_md = group_name_list["PG"][pg][1]
        mpg_md = group_name_list["MPG"][mpg][1]
        pg_name = group_name_list["PG"][pg][4]
        mpg_name = group_name_list["MPG"][mpg][4]

        sg_dict = {
            "Related group": {
                "PG": f"[{pg_name}](../PG/{pg_md})",
                "MPG": f"[{mpg_name}](../MPG/{mpg_md})",
                "MSG (SG)": f"[{name}](../MSG/{md})",
            },
            "Symmetry operation": f"[PDF]({dd}/symmetry_operation.pdf)",
            "Wyckoff": {
                "site": f"[PDF]({dd}/wyckoff_site.pdf)",
                "bond": f"[PDF]({dd}/wyckoff_bond.pdf)",
            },
        }

        write_group_markdown(sg_dict, "SG", title, md)

    # ---------- Magnetic point group ----------
    for id_s, (d, md, no, tag, name, pg, sg, mpg, msg) in group_name_list["MPG"].items():
        mpg_tbl.append(f"{no}:  [{name}](MPG/{md})")

        xl = ["Q", "G"] if name.endswith("1'$") else ["Q", "G", "T", "M"]
        title = f"# MPG #{no}:  {name}\n"
        dd = f"{dir_pdf2}MPG/{d}"

        pg_md = group_name_list["PG"][pg][1]
        sg_md = group_name_list["SG"][sg][1]
        pg_name = group_name_list["PG"][pg][4]
        sg_name = group_name_list["SG"][sg][4]

        cg_dict = {
            "Related group": {
                "PG": f"[{pg_name}](../PG/{pg_md})",
                "SG": f"[{sg_name}](../SG/{sg_md})",
                "MSG (SG)": f"[{sg_name}](../MSG/{sg_md})",
            },
            "Symmetry operation": f"[PDF]({dd}/symmetry_operation.pdf)",
            "Wyckoff": {"site": f"[PDF]({dd}/wyckoff_site.pdf)"},
            "Response tensor": {X: f"[PDF]({dd}/response_tensor_" + X + ".pdf)" for X in xl},
        }

        write_group_markdown(cg_dict, "MPG", title, md)

    # ---------- Magnetic space group ----------
    for SG, lst in info["id_set"]["MSG"]["SG"].items():
        d, md, no, tag, name, pg, sg, mpg, msg = group_name_list["SG"][SG]

        msg_tbl.append(f"{no}:  [{name}](MSG/{md})")

        title = f"# MSG (SG) #{no}:  {name}\n"

        cg_dict = {}
        for id_s in lst:
            msg_info = group_name_list["MSG"][id_s]
            no_c = msg_info[2]  # no.
            tag_c = f"#{no_c}:  {msg_info[4]}"  # name.
            dd = f"{dir_pdf2}MSG/{no_c}"

            pg_md = group_name_list["PG"][pg][1]
            sg_md = group_name_list["SG"][sg][1]
            mpg_md = group_name_list["MPG"][mpg][1]
            pg_name = group_name_list["PG"][pg][4]
            sg_name = group_name_list["SG"][sg][4]
            mpg_name = group_name_list["MPG"][mpg][4]

            cg_dict[f"{tag_c}"] = {
                "Related group": {
                    "PG": f"[{pg_name}](../PG/{pg_md})",
                    "SG": f"[{sg_name}](../SG/{sg_md})",
                    "MPG": f"[{mpg_name}](../MPG/{mpg_md})",
                },
                "Symmetry operation": f"[PDF]({dd}/symmetry_operation.pdf)",
                "Wyckoff": {"site": f"[PDF]({dd}/wyckoff_site.pdf)"},
            }

        write_group_markdown(cg_dict, "MSG", title, md)

    # ---------- Merge real + complex ----------
    point_tbl_md = list_to_table(point_tbl, header=True)
    point_tbl_md += "\n" + list_to_table(complex_tbl, header=False)

    space_tbl_md = list_to_table(space_tbl)
    mpg_tbl_md = list_to_table(mpg_tbl)
    msg_tbl_md = list_to_table(msg_tbl, tag="tag (SG)")

    # ---------- Main markdown ----------
    main_md = [
        "# Group information\n",
        "## Overview\n",
        overview_md,
        "## Point group\n",
        point_tbl_md,
        "## Space group\n",
        space_tbl_md,
        "## Magnetic point group\n",
        mpg_tbl_md,
        "## Magnetic space group (by SG)\n",
        msg_tbl_md,
    ]

    with open(os.path.join(dir_md, "group_doc.md"), "w", encoding="utf-8") as f:
        for line in main_md:
            print(line, file=f)


# ================================================== main
if __name__ == "__main__":
    create_markdown()
