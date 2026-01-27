# graphene_param.py
"""
selection_parameter file for graphene.
if parameter is omitted, no hr.dat file is created.
"""
graphene_sel_par = {
    # model pkl name (w/o .pkl).
    "model": "graphene",
    # selection for SAMB.
    "select": {
        "bond": [[1, 2]],
    },
    # weight for SAMBs, non-specified parameters are zero. (float or sympy const.)
    "parameter": {
        "z1": 1,
        "z2": 1,
    },
}
