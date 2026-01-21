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
        "X": [],  # combined SAMB type Q/G/T/M, []=all.
        "l": [],  # combined SAMB rank, []=all.
        "Gamma": "A1g",  # combined SAMB irreps., e.g., "A1g", ["A1g","A1u"], "IR"=identity irrep., []=all.
        "s": [],  # combined SAMB s=0, 1, []=all.
    },
    # weight for SAMBs, non-specified parameters are zero. (float or sympy const.)
    "parameter": {
        "z1": 1,
        "z2": 1,
    },
}
