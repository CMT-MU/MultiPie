"""
=== SAMB (* only for crystal) ===
- info
    - atomic : { "M_#" : ["amp_#"] }
    - site_cluster : { "S_#" : ["smp_#"] }
    - bond_cluster : { "B_#" : ["bmp_#"] }
    - uniform : { "S_#"/"B_#" : ["ump_#"] }
    - structure* : { "B_#" : ["kmp_#"] }
    - Z : { ("M_#", "S_#"/"B_#") : ["z_#"] }
    - version : MultiPie version
    - harmonics : { head : { "harm_tag" } }

- data
    - atomic : { "amp_#" : ( TagMultipole, shape, [(i, j, matrix element)] ) }
    - site_cluster : { "smp_#" : ( TagMultipole, [vector component] ) }
    - bond_cluster : { "bmp_#" : ( TagMultipole, [vector component] ) }
    - uniform : { "ump_#" : ( TagMultipole, shape, [(i, j, matrix element)] ) }
    - structure* : { "kmp_#" : (TagMultipole, "formfactor") }
    - Z : {"z_#" : ( TagMultipole, [(coeff, "amp_#", "smp_#"/"bmp_#/ump_#")] ) }
    - Zk* : {"z_#" : ( TagMultipole, [(coeff, "amp_#", "ump_#", "kmp_#")] ) }
"""
GaAs = {
    "info": {
        "atomic": {"M_001": ["amp_001", "amp_002", "amp_003", "amp_004"]},
        "site_cluster": {"S_001": ["smp_001"], "S_002": ["smp_002"]},
        "bond_cluster": {"B_001": ["bmp_003", "bmp_004", "bmp_005", "bmp_006"]},
        "uniform": {"S_001": ["ump_001"], "S_002": ["ump_002"], "B_001": ["ump_003", "ump_004"]},
        "Z": {
            ("A1", "M_001", "S_001"): ["z_001"],
            ("A1", "M_001", "S_002"): ["z_002"],
            ("A1", "M_001", "B_001"): ["z_003", "z_004"],
        },
        "version": "1.1.1",
        "structure": {"B_001": ["kmp_001", "kmp_002", "kmp_003", "kmp_004", "kmp_005", "kmp_006", "kmp_007", "kmp_008"]},
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
        "uniform": {
            "ump_001": ("Qs(0,A1,,)", (2, 2), [(0, 0, "1")]),
            "ump_002": ("Qs(0,A1,,)", (2, 2), [(1, 1, "1")]),
            "ump_003": ("Qu(0,A1,,)", (2, 2), [(0, 1, "sqrt(2)/2"), (1, 0, "sqrt(2)/2")]),
            "ump_004": ("Tu(0,A1,,)", (2, 2), [(0, 1, "sqrt(2)*I/2"), (1, 0, "-sqrt(2)*I/2")]),
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
        "structure": {
            "kmp_001": ("Qk(0,A1,,)", "sqrt(2)*c001/2 + sqrt(2)*c002/2 + sqrt(2)*c003/2 + sqrt(2)*c004/2"),
            "kmp_002": ("Qk(1,T2,,0)", "sqrt(2)*c001/2 - sqrt(2)*c002/2 + sqrt(2)*c003/2 - sqrt(2)*c004/2"),
            "kmp_003": ("Qk(1,T2,,1)", "sqrt(2)*c001/2 - sqrt(2)*c002/2 - sqrt(2)*c003/2 + sqrt(2)*c004/2"),
            "kmp_004": ("Qk(1,T2,,2)", "sqrt(2)*c001/2 + sqrt(2)*c002/2 - sqrt(2)*c003/2 - sqrt(2)*c004/2"),
            "kmp_005": ("Tk(0,A1,,)", "sqrt(2)*s001/2 + sqrt(2)*s002/2 + sqrt(2)*s003/2 + sqrt(2)*s004/2"),
            "kmp_006": ("Tk(1,T2,,0)", "sqrt(2)*s001/2 - sqrt(2)*s002/2 + sqrt(2)*s003/2 - sqrt(2)*s004/2"),
            "kmp_007": ("Tk(1,T2,,1)", "sqrt(2)*s001/2 - sqrt(2)*s002/2 - sqrt(2)*s003/2 + sqrt(2)*s004/2"),
            "kmp_008": ("Tk(1,T2,,2)", "sqrt(2)*s001/2 + sqrt(2)*s002/2 - sqrt(2)*s003/2 - sqrt(2)*s004/2"),
        },
        "Zk": {
            "z_001": ("Q(0,A1,,)", [("1", "amp_001", "ump_001")]),
            "z_002": ("Q(0,A1,,)", [("1", "amp_001", "ump_002")]),
            "z_003": (
                "Q(0,A1,,)",
                [("sqrt(2)/2", "amp_001", "ump_003", "kmp_001"), ("-sqrt(2)/2", "amp_001", "ump_004", "kmp_005")],
            ),
            "z_004": (
                "Q(3,A1,,)",
                [
                    ("sqrt(6)/6", "amp_002", "ump_003", "kmp_002"),
                    ("-sqrt(6)/6", "amp_002", "ump_004", "kmp_006"),
                    ("sqrt(6)/6", "amp_003", "ump_003", "kmp_003"),
                    ("-sqrt(6)/6", "amp_003", "ump_004", "kmp_007"),
                    ("sqrt(6)/6", "amp_004", "ump_003", "kmp_004"),
                    ("-sqrt(6)/6", "amp_004", "ump_004", "kmp_008"),
                ],
            ),
        },
    },
}
