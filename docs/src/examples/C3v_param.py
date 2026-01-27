# C3v_param.py
"""
selection_parameter file for C3v.
if parameter is omitted, no hr.dat file is created.
"""
C3v_sel_par = {
    # model pkl name (w/o .pkl).
    "model": "C3v",
    # selection for SAMB.
    "select": {
        "site": "A",
        "bond": [[1]],
    },
    # weight for SAMBs, non-specified parameters are zero. (float or sympy const.)
    "parameter": {
        "z1": 1,
        "z6": 1,
    },
}
