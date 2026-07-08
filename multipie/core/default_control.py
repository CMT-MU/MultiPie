"""
Default control for ModelAnalyzer.
"""

from multipie import __version__

# ==================================================
# default control.
default_control = {
    "samb": {
        "model": None,  # model name for .pkl.
        "select": {  # SAMB select condition, where S=site name, R=orbital rank, N=neighbor bond.
            # "site": [("A", [1, 2])],  # S or (S, [R]) in list.
            # "bond": [[0, 1, 2]],  # S1;S2, [N], (S1;S2, [N]), (S1;S2, R;R2), (S1;S2, [N]) or (S1;S2, R1;R2, [N]) in list.
            # "X": [],  # SAMB type Q/G/T/M in list. empty []=all.
            # "l": [],  # SAMB rank, 0, 1, ..., 11 in list. empty []=all.
            # "Gamma": ["A1g"],  # SAMB irreps. in list. "IR"=identity irrep., empty []=all.
            # "s": [],  # SAMB internal rank, 0,1 in list. empty []=all
        },
        "parameter": {  # SAMB with finite weight (float or sympy const.).
            "z1": 1.0,
        },
        "samb_figure": False,  # save SAMB QtDraw files ?
    },
    "wannier": {  # Closest Wannier (CW) or SymWannier setting.
        "cw": None,  # CW or SymWannier file.
    },
    "output": {  # physical quantity setting.
        "dir": "output",  #  output directory.
        "fourier": {"atom_phase": True},  # phase option to use atomic position.
        "dispersion": {  # dispersion info.
            "k_point": {"Γ": "[0,0,0]", "X": "[1,0,0]"},  # k-point definition (primitive, fractional).
            "k_path": "Γ",  # symmetry line, separated by "-". disconnected points by "|".
        },
        "dos": False,  # compute DOS ?
    },
}
