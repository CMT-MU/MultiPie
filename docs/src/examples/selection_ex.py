# selection and parameter example.
selection_paramter = {
    # model pkl name (w/o .pkl).
    "model": "model_pkl_name",
    # selection for SAMB.
    "select": {  # S=site, R=orbital rank.
        "site": [("A", [1, 2])],  # S or (S, [R]).
        "bond": [
            [0, 1, 2]
        ],  # S1;S2, [neighbor], (S1;S2, [neighbor]), (S1;S2, R;R2), (S1;S2, [neighbor]) or (S1;S2, R1;R2, [neighbor]).
        "X": ["Q", "G"],  # combined SAMB type Q/G/T/M.
        "l": [0, 1, 2],  # combined SAMB rank.
        "Gamma": "A1g",  # combined SAMB irreps.
        "s": [0, 1],  # combined SAMB s=0, 1.
    },
    # weight for SAMBs, non-specified parameters are zero. (float or sympy const.)
    "parameter": {
        "z1": 1,
        "z2": 1,
    },
}
