"""
create virtual-cluster basis.
"""
import os
import sympy as sp
from gcoreutils.nsarray import NSArray
from gcoreutils.string_util import remove_space

from multipie.data.data_harmonics import _data_harmonics_polar, _data_harmonics_axial
from util.output_util import write_data

create_dir = __file__[: __file__.rfind("/")]


# ==================================================
header_polar = """
data of point-group harmonics (polar).
{ tag_harmonics: data }
    data = ( tesseral_harmonics ex., cartesian ex., unitary matrix in descending order, [l,l-1,...,-l] )
Ea and Eb are promoted to E 0:(Ea+Eb)/sqrt(2), 1:(Ea-Eb)/sqrt(2)i.
"""

header_axial = """
data of point-group harmonics (axial).
{ tag_harmonics: data }
    data = ( tesseral_harmonics ex., cartesian ex., unitary matrix in descending order, [l,l-1,...,-l] )
Ea and Eb are promoted to E 0:(Ea+Eb)/sqrt(2), 1:(Ea-Eb)/sqrt(2)i.
"""


# ==================================================
def create_real_harmonics(output_dir=None):
    if output_dir is not None:
        os.chdir(output_dir)

    ofile = "data_harmonics_real.py"

    print("=== create data_harmonics_real ===")

    cpg = ["C4", "S4", "C4h", "C3", "C3i", "C6", "C3h", "C6h", "T", "Th"]

    data = []
    for harm in [_data_harmonics_polar, _data_harmonics_axial]:
        dic = {}
        for pg in cpg:
            dic[pg] = {}
            for k, v in harm[pg].items():
                if "E" in k and "a" in k:
                    ea_def, ea_ex, ea_u = v
                    eb_def, eb_ex, eb_u = harm[pg][k.replace("a", "b")]
                    x_def = remove_space(str(((NSArray(ea_def) + NSArray(eb_def)) / sp.sqrt(2)).simplify()))
                    y_def = remove_space(str(((NSArray(ea_def) - NSArray(eb_def)) / sp.sqrt(2) / sp.I).simplify()))
                    x_ex = remove_space(str(((NSArray(ea_ex) + NSArray(eb_ex)) / sp.sqrt(2)).simplify()))
                    y_ex = remove_space(str(((NSArray(ea_ex) - NSArray(eb_ex)) / sp.sqrt(2) / sp.I).simplify()))
                    x_u = remove_space(str(((NSArray(ea_u) + NSArray(eb_u)) / sp.sqrt(2)).simplify()))
                    y_u = remove_space(str(((NSArray(ea_u) - NSArray(eb_u)) / sp.sqrt(2) / sp.I).simplify()))
                    dic[pg][k.replace("a", "")[:-1] + "0)"] = (x_def, x_ex, x_u)
                    dic[pg][k.replace("a", "")[:-1] + "1)"] = (y_def, y_ex, y_u)
                elif "b" not in k:
                    dic[pg][k] = v
        data.append(dic)

    write_data(ofile, data[0], header_polar, "_data_harmonics_polar_real", mode="w")
    write_data(ofile, data[1], header_axial, "_data_harmonics_axial_real", mode="a")


# ================================================== main
create_real_harmonics(output_dir=create_dir)
