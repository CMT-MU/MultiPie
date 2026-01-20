# selection and parameter example.
selection_paramter = {
    # model pkl name (w/o .pkl).
    "model": "model_pkl_name",
    #
    #
    # select.
    "SAMB_select": {  # select combined SAMB.
        "X": ["Q", "G"],  # type, "Q/G/T/M", []=all.
        "l": [],  # rank, []=all.
        "Gamma": [],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
        "s": [],  # spin, 0/1, []=all.
    },
    "atomic_select": {  # select atomic SAMB.
        "X": [],  # type, "Q/G/T/M", []=all.
        "l": [],  # rank, []=all.
        "Gamma": [],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
        "s": [],  # spin, 0/1, []=all.
    },
    "site_select": {  # select site-cluster SAMB.
        "l": [],  # rank, []=all.
        "Gamma": [],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
        "site": [("A", [1, 2])],  # Site or (Site, [Rank]), []=all.
    },
    "bond_select": {  # select bond-cluster SAMB.
        "X": [],  # type, "Q/T/M", []=all.
        "l": [],  # rank, []=all.
        "Gamma": [],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
        "bond": [
            [0, 1, 2]
        ],  # S1;S2, [neighbor], (S1;S2, [neighbor]), (S1;S2, R1;R2), (S1;S2, [neighbor]) or (S1;S2, R1;R2, [neighbor]), []=all.
    },
    #
    #
    # weight for SAMBs, non-specified parameters are zero. (float or sympy const.)
    "parameter": {
        "z1": 1,
        "z2": 1,
    },
}
