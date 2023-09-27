"""
=== SAMB (* only for crystal with fourier_transform==True) ===
- info
    - atomic : { "M_#" : ["amp_#"] }
    - site_cluster : { "S_#" : ["smp_#"] }
    - bond_cluster : { "B_#" : ["bmp_#"] }
    - uniform : { "S_#"/"B_#" : ["ump_#"] }
    - structure* : { "B_#" : ["kmp_#"] }
    - Z : { ("M_#", "S_#"/"B_#") : ["z_#"] }
    - version : MultiPie version
    - harmonics : { head : [TagMultipole] }

- data
    - atomic : { "amp_#" : ( TagMultipole, shape, [(i, j, matrix element)] ) }
    - site_cluster : { "smp_#" : ( TagMultipole, [vector component] ) }
    - bond_cluster : { "bmp_#" : ( TagMultipole, [vector component] ) }
    - uniform : { "ump_#" : ( TagMultipole, shape, [(i, j, matrix element)] ) }
    - structure* : { "kmp_#" : (TagMultipole, "structure factor") }
    - Z : {"z_#" : ( TagMultipole, [(coeff, "amp_#", "smp_#"/"bmp_#/ump_#")] ) }
    - Zk* : {"z_#" : ( TagMultipole, [(coeff, "amp_#", "ump_#", "kmp_#")] ) }
"""
GaAs = {
    "info": {
        "atomic": {"M_001": ["amp_001", "amp_002", "amp_003", "amp_004"]},
        "site_cluster": {"S_001": ["smp_001"], "S_002": ["smp_002"]},
        "bond_cluster": {"B_001": ["bmp_003", "bmp_004", "bmp_005", "bmp_006"]},
        "Z": {
            ("A1", "M_001", "S_001"): ["z_001"],
            ("A1", "M_001", "S_002"): ["z_002"],
            ("A1", "M_001", "B_001"): ["z_003", "z_004"],
        },
        "version": "1.1.14",
        "harmonics": {
            "Q": ["Qh(0,A1,,)", "Qh(1,T2,,0)", "Qh(1,T2,,1)", "Qh(1,T2,,2)", "Qh(2,T2,,0)", "Qh(2,T2,,1)", "Qh(2,T2,,2)"],
            "G": [],
        },
    },
    "data": {
        "atomic": {
            "amp_001": ("Qa(0,A1,,)", (3, 3), [(0, 0, "sqrt(3)/3"), (1, 1, "sqrt(3)/3"), (2, 2, "sqrt(3)/3")]),
            "amp_002": ("Qa(2,T2,,0)", (3, 3), [(1, 2, "sqrt(2)/2"), (2, 1, "sqrt(2)/2")]),
            "amp_003": ("Qa(2,T2,,1)", (3, 3), [(0, 2, "sqrt(2)/2"), (2, 0, "sqrt(2)/2")]),
            "amp_004": ("Qa(2,T2,,2)", (3, 3), [(0, 1, "sqrt(2)/2"), (1, 0, "sqrt(2)/2")]),
        },
        "site_cluster": {"smp_001": ("Qs(0,A1,,)", "[1]"), "smp_002": ("Qs(0,A1,,)", "[1]")},
        "bond_cluster": {
            "bmp_003": ("Qb(0,A1,,)", "[1/2, 1/2, 1/2, 1/2]"),
            "bmp_004": ("Qb(1,T2,,0)", "[1/2, -1/2, 1/2, -1/2]"),
            "bmp_005": ("Qb(1,T2,,1)", "[1/2, -1/2, -1/2, 1/2]"),
            "bmp_006": ("Qb(1,T2,,2)", "[1/2, 1/2, -1/2, -1/2]"),
        },
        "Z": {
            "z_001": ("Q(0,A1,,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(0,A1,,)", [("1", "amp_001", "smp_002")]),
            "z_003": ("Q(0,A1,,)", [("1", "amp_001", "bmp_003")]),
            "z_004": (
                "Q(3,A1,,)",
                [("sqrt(3)/3", "amp_002", "bmp_004"), ("sqrt(3)/3", "amp_003", "bmp_005"), ("sqrt(3)/3", "amp_004", "bmp_006")],
            ),
        },
    },
}
