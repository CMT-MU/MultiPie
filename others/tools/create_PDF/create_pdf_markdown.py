"""
Create markdown for all information.
"""

import os

from multipie.core.multipie_info import __top_dir__
from multipie.core.group import Group
from multipie.util.util import timer
from multipie.util.util_binary import BinaryManager


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
    def write_group_markdown(group_dict, dir_out, g, no, with_inter=False):
        if with_inter:
            title = f"# #{no}: ${g.info.international}$\n"
        else:
            title = f"# #{no}: ${g.latex()}$\n"
        md_content = title + "\n" + dict_to_markdown(group_dict) + "\n"
        out_path = os.path.join(dir_md, dir_out, f"No_{no}.md")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(md_content)

    # ---------- Containers ----------
    point_tbl, complex_tbl, space_tbl, mpg_tbl, msg_tbl = [], [], [], [], []

    # ---------- Space group ----------
    for id_s in info["id_set"]["SG"]["all"]:
        g = Group(id_s)
        tag, no = g.info.tag.replace("^", "_"), int(id_s.split(":")[1])
        space_tbl.append(f"{no}: [$\\,\\,{g.latex()}$ (${g.info.international}$)](SG/No_{no}.md)")

        sg_dict = {
            "Symmetry operation": f"[PDF]({dir_pdf2}SG/{no:03d}-{tag}/symmetry_operation.pdf)",
            "Wyckoff": {
                "site": f"[PDF]({dir_pdf2}SG/{no:03d}-{tag}/wyckoff_site.pdf)",
                "bond": f"[PDF]({dir_pdf2}SG/{no:03d}-{tag}/wyckoff_bond.pdf)",
            },
        }
        write_group_markdown(sg_dict, "SG", g, no)

    # ---------- Point group (real) ----------
    for id_s in info["id_set"]["PG"]["all"]:
        g = Group(id_s)
        tag, no = g.info.tag, int(id_s.split(":")[1])
        point_tbl.append(f"{no}: [$\\,\\,{g.latex()}$ (${g.info.international}$)](PG/No_{no}.md)")

        xl = ["Q", "G"]
        mno = g.info.MPG.split(":")[1]

        pg_dict = {
            "Symmetry operation": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/symmetry_operation.pdf)",
            "Character table": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/character_table.pdf)",
            "Wyckoff": {
                "site": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/wyckoff_site.pdf)",
                "bond": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/wyckoff_bond.pdf)",
            },
            "Atomic multipole": {
                "JML basis": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/atomic_multipole_jml.pdf)",
                r"L$\Gamma\sigma$ basis": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/atomic_multipole_lgs.pdf)",
                r"L$\Gamma$ basis": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/atomic_multipole_lg.pdf)",
            },
            "Harmonics": {
                "polar": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_polar.pdf)",
                "axial": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_axial.pdf)",
            },
            "Multipolar Harmonics (internal)": {
                "dipolar internal polar(Q)/axial(G) varialble": {
                    "polar (Q)": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_s1_polar_q.pdf)",
                    "axial (Q)": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_s1_axial_q.pdf)",
                    "polar (G)": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_s1_polar_g.pdf)",
                    "axial (G)": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_s1_axial_g.pdf)",
                },
                "quadrupolar internal polar(Q)/axial(G) variable": {
                    "polar (Q)": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_s2_polar_q.pdf)",
                    "axial (Q)": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_s2_axial_q.pdf)",
                    "polar (G)": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_s2_polar_g.pdf)",
                    "axial (G)": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_s2_axial_g.pdf)",
                },
                "octupolar internal polar(Q)/axial(G) variable": {
                    "polar (Q)": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_s3_polar_q.pdf)",
                    "axial (Q)": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_s3_axial_q.pdf)",
                    "polar (G)": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_s3_polar_g.pdf)",
                    "axial (G)": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_s3_axial_g.pdf)",
                },
            },
            "Response tensor": {X: f"[PDF]({dir_pdf2}MPG/{mno}/response_tensor_" + X + ".pdf)" for X in xl},
        }
        write_group_markdown(pg_dict, "PG", g, no)

    # ---------- Point group (complex) ----------
    for id_s in info["id_set"]["PG"]["complex"]:
        g = Group(id_s)
        tag, no = g.info.tag, int(id_s.split(":")[1])
        complex_tbl.append(f"{no} :[$\\,\\,{g.latex()}$ (${g.info.international}$)](PG/No_{no}.md)")

        cg_dict = {
            "Symmetry operation": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/symmetry_operation.pdf)",
            "Character table": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/character_table.pdf)",
            "Harmonics": {
                "polar": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_polar.pdf)",
                "axial": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/harmonics_axial.pdf)",
            },
            "Representation matrix": f"[PDF]({dir_pdf2}PG/{no:03d}-{tag}/representation_matrix.pdf)",
        }
        write_group_markdown(cg_dict, "PG", g, no)

    # ---------- Magnetic point group ----------
    for id_s in info["id_set"]["MPG"]["all"]:
        g = Group(id_s)
        tag, no = g.info.tag, id_s.split(":")[1]
        mpg_tbl.append(f"{no} :[$\\,\\,{g.latex()}$](MPG/No_{no}.md)")
        xl = ["Q", "G"] if g.info.type == "II" else ["Q", "G", "T", "M"]

        cg_dict = {
            "Symmetry operation": f"[PDF]({dir_pdf2}MPG/{no}/symmetry_operation.pdf)",
            "Wyckoff": {"site": f"[PDF]({dir_pdf2}MPG/{no}/wyckoff_site.pdf)"},
            "Response tensor": {X: f"[PDF]({dir_pdf2}MPG/{no}/response_tensor_" + X + ".pdf)" for X in xl},
        }
        write_group_markdown(cg_dict, "MPG", g, no)

    # ---------- Magnetic space group ----------
    for SG, lst in info["id_set"]["MSG"]["SG"].items():
        gs = Group(SG)
        tag, no = gs.info.tag, SG.split(":")[1]
        msg_tbl.append(f"{no} :[$\\,\\,{gs.info.international}$](MSG/No_{no}.md)")
        cg_dict = {}
        for id_s in lst:
            g = Group(id_s)
            no_c = id_s.split(":")[1]
            tag_c = "#" + no_c + ": $" + g.info.BNS + "$"
            cg_dict[f"{tag_c}"] = {
                "Symmetry operation": f"[PDF]({dir_pdf2}MSG/{no_c}/symmetry_operation.pdf)",
                "Wyckoff": {"site": f"[PDF]({dir_pdf2}MSG/{no_c}/wyckoff_site.pdf)"},
            }
        write_group_markdown(cg_dict, "MSG", gs, no, True)

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
