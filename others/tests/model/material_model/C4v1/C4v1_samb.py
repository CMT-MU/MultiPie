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
C4v1 = {
    "info": {
        "atomic": {"M_001": ["amp_001", "amp_002", "amp_003"]},
        "site_cluster": {"S_001": ["smp_001"]},
        "bond_cluster": {"B_001": ["bmp_002", "bmp_003", "bmp_004"], "B_002": ["bmp_005"]},
        "uniform": {"S_001": ["ump_001"]},
        "Z": {
            ("A1", "M_001", "S_001"): ["z_001"],
            ("A1", "M_001", "B_001"): ["z_002", "z_003"],
            ("A1", "M_001", "B_002"): ["z_004"],
        },
        "version": "1.1.10",
        "structure": {"B_001": ["kmp_001", "kmp_002", "kmp_003"], "B_002": ["kmp_004"]},
        "harmonics": {"Q": ["Qh(0,A1,,)", "Qh(1,E,,0)", "Qh(1,E,,1)"], "G": ["Gh(1,E,,0)", "Gh(1,E,,1)"]},
    },
    "data": {
        "atomic": {
            "amp_001": ("Qa(0,A1,,)", (2, 2), [(0, 0, "sqrt(2)/2"), (1, 1, "sqrt(2)/2")]),
            "amp_002": ("Ma(1,E,,0|1,-1)", (2, 2), [(0, 1, "sqrt(2)/2"), (1, 0, "sqrt(2)/2")]),
            "amp_003": ("Ma(1,E,,1|1,-1)", (2, 2), [(0, 1, "-sqrt(2)*I/2"), (1, 0, "sqrt(2)*I/2")]),
        },
        "site_cluster": {"smp_001": ("Qs(0,A1,,)", "[1]")},
        "bond_cluster": {
            "bmp_002": ("Qb(0,A1,,)", "[sqrt(2)/2, sqrt(2)/2]"),
            "bmp_003": ("Tb(1,E,,0)", "[0, I]"),
            "bmp_004": ("Tb(1,E,,1)", "[I, 0]"),
            "bmp_005": ("Qb(0,A1,,)", "[1]"),
        },
        "uniform": {"ump_001": ("Qs(0,A1,,)", (1, 1), [(0, 0, "1")])},
        "Z": {
            "z_001": ("Q(0,A1,,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(0,A1,,)", [("1", "amp_001", "bmp_002")]),
            "z_003": ("Q(1,A1,,|1,-1)", [("sqrt(2)/2", "amp_002", "bmp_004"), ("-sqrt(2)/2", "amp_003", "bmp_003")]),
            "z_004": ("Q(0,A1,,)", [("1", "amp_001", "bmp_005")]),
        },
        "structure": {
            "kmp_001": ("Qk(0,A1,,)", "c001 + c002"),
            "kmp_002": ("Tk(1,E,,0)", "sqrt(2)*s002"),
            "kmp_003": ("Tk(1,E,,1)", "sqrt(2)*s001"),
            "kmp_004": ("Qk(0,A1,,)", "sqrt(2)*c003"),
        },
        "Zk": {
            "z_001": ("Q(0,A1,,)", [("1", "amp_001", "ump_001")]),
            "z_002": ("Q(0,A1,,)", [("1", "amp_001", "ump_001", "kmp_001")]),
            "z_003": (
                "Q(1,A1,,|1,-1)",
                [("sqrt(2)/2", "amp_002", "ump_001", "kmp_003"), ("-sqrt(2)/2", "amp_003", "ump_001", "kmp_002")],
            ),
            "z_004": ("Q(0,A1,,)", [("1", "amp_001", "ump_001", "kmp_004")]),
        },
    },
}
