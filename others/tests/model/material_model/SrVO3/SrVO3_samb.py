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
SrVO3 = {
    "info": {
        "atomic": {"M_001": ["amp_001", "amp_002", "amp_003", "amp_004", "amp_005", "amp_006"]},
        "site_cluster": {"S_001": ["smp_001"]},
        "bond_cluster": {
            "B_001": ["bmp_002", "bmp_003", "bmp_004"],
            "B_002": ["bmp_005", "bmp_006", "bmp_007", "bmp_008", "bmp_009", "bmp_010"],
        },
        "uniform": {"S_001": ["ump_001"]},
        "Z": {
            ("A1g", "M_001", "S_001"): ["z_001"],
            ("A1g", "M_001", "B_001"): ["z_002", "z_003"],
            ("A1g", "M_001", "B_002"): ["z_004", "z_005", "z_006"],
        },
        "version": "1.1.2",
        "structure": {
            "B_001": ["kmp_001", "kmp_002", "kmp_003"],
            "B_002": ["kmp_004", "kmp_005", "kmp_006", "kmp_007", "kmp_008", "kmp_009"],
        },
        "harmonics": {
            "Q": ["Qh(0,A1g,,)", "Qh(2,Eg,,0)", "Qh(2,Eg,,1)", "Qh(2,T2g,,0)", "Qh(2,T2g,,1)", "Qh(2,T2g,,2)"],
            "G": [],
        },
    },
    "data": {
        "atomic": {
            "amp_001": ("Qa(0,A1g,,)", (3, 3), [(0, 0, "sqrt(3)/3"), (1, 1, "sqrt(3)/3"), (2, 2, "sqrt(3)/3")]),
            "amp_002": ("Qa(2,Eg,,0)", (3, 3), [(0, 0, "sqrt(6)/6"), (1, 1, "sqrt(6)/6"), (2, 2, "-sqrt(6)/3")]),
            "amp_003": ("Qa(2,Eg,,1)", (3, 3), [(0, 0, "-sqrt(2)/2"), (1, 1, "sqrt(2)/2")]),
            "amp_004": ("Qa(2,T2g,,0)", (3, 3), [(1, 2, "sqrt(2)/2"), (2, 1, "sqrt(2)/2")]),
            "amp_005": ("Qa(2,T2g,,1)", (3, 3), [(0, 2, "sqrt(2)/2"), (2, 0, "sqrt(2)/2")]),
            "amp_006": ("Qa(2,T2g,,2)", (3, 3), [(0, 1, "sqrt(2)/2"), (1, 0, "sqrt(2)/2")]),
        },
        "site_cluster": {"smp_001": ("Qs(0,A1g,,)", "[1]")},
        "bond_cluster": {
            "bmp_002": ("Qb(0,A1g,,)", "[sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]"),
            "bmp_003": ("Qb(2,Eg,,0)", "[-sqrt(6)/3, sqrt(6)/6, sqrt(6)/6]"),
            "bmp_004": ("Qb(2,Eg,,1)", "[0, -sqrt(2)/2, sqrt(2)/2]"),
            "bmp_005": ("Qb(0,A1g,,)", "[sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6]"),
            "bmp_006": ("Qb(2,Eg,,0)", "[-sqrt(3)/6, -sqrt(3)/6, -sqrt(3)/6, sqrt(3)/3, -sqrt(3)/6, sqrt(3)/3]"),
            "bmp_007": ("Qb(2,Eg,,1)", "[1/2, 1/2, -1/2, 0, -1/2, 0]"),
            "bmp_008": ("Qb(2,T2g,,0)", "[sqrt(2)/2, -sqrt(2)/2, 0, 0, 0, 0]"),
            "bmp_009": ("Qb(2,T2g,,1)", "[0, 0, -sqrt(2)/2, 0, sqrt(2)/2, 0]"),
            "bmp_010": ("Qb(2,T2g,,2)", "[0, 0, 0, -sqrt(2)/2, 0, sqrt(2)/2]"),
        },
        "uniform": {"ump_001": ("Qs(0,A1g,,)", (1, 1), [(0, 0, "1")])},
        "Z": {
            "z_001": ("Q(0,A1g,,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(0,A1g,,)", [("1", "amp_001", "bmp_002")]),
            "z_003": ("Q(0,A1g,,)", [("sqrt(2)/2", "amp_002", "bmp_003"), ("sqrt(2)/2", "amp_003", "bmp_004")]),
            "z_004": ("Q(0,A1g,,)", [("1", "amp_001", "bmp_005")]),
            "z_005": (
                "Q(0,A1g,,)",
                [
                    ("sqrt(5)/5", "amp_002", "bmp_006"),
                    ("sqrt(5)/5", "amp_003", "bmp_007"),
                    ("sqrt(5)/5", "amp_004", "bmp_008"),
                    ("sqrt(5)/5", "amp_005", "bmp_009"),
                    ("sqrt(5)/5", "amp_006", "bmp_010"),
                ],
            ),
            "z_006": (
                "Q(4,A1g,,)",
                [
                    ("sqrt(30)/10", "amp_002", "bmp_006"),
                    ("sqrt(30)/10", "amp_003", "bmp_007"),
                    ("-sqrt(30)/15", "amp_004", "bmp_008"),
                    ("-sqrt(30)/15", "amp_005", "bmp_009"),
                    ("-sqrt(30)/15", "amp_006", "bmp_010"),
                ],
            ),
        },
        "structure": {
            "kmp_001": ("Qk(0,A1g,,)", "sqrt(6)*c001/3 + sqrt(6)*c002/3 + sqrt(6)*c003/3"),
            "kmp_002": ("Qk(2,Eg,,0)", "-2*sqrt(3)*c001/3 + sqrt(3)*c002/3 + sqrt(3)*c003/3"),
            "kmp_003": ("Qk(2,Eg,,1)", "-c002 + c003"),
            "kmp_004": (
                "Qk(0,A1g,,)",
                "sqrt(3)*c004/3 + sqrt(3)*c005/3 + sqrt(3)*c006/3 + sqrt(3)*c007/3 + sqrt(3)*c008/3 + sqrt(3)*c009/3",
            ),
            "kmp_005": (
                "Qk(2,Eg,,0)",
                "-sqrt(6)*c004/6 - sqrt(6)*c005/6 - sqrt(6)*c006/6 + sqrt(6)*c007/3 - sqrt(6)*c008/6 + sqrt(6)*c009/3",
            ),
            "kmp_006": ("Qk(2,Eg,,1)", "sqrt(2)*c004/2 + sqrt(2)*c005/2 - sqrt(2)*c006/2 - sqrt(2)*c008/2"),
            "kmp_007": ("Qk(2,T2g,,0)", "c004 - c005"),
            "kmp_008": ("Qk(2,T2g,,1)", "-c006 + c008"),
            "kmp_009": ("Qk(2,T2g,,2)", "-c007 + c009"),
        },
        "Zk": {
            "z_001": ("Q(0,A1g,,)", [("1", "amp_001", "ump_001")]),
            "z_002": ("Q(0,A1g,,)", [("1", "amp_001", "ump_001", "kmp_001")]),
            "z_003": (
                "Q(0,A1g,,)",
                [("sqrt(2)/2", "amp_002", "ump_001", "kmp_002"), ("sqrt(2)/2", "amp_003", "ump_001", "kmp_003")],
            ),
            "z_004": ("Q(0,A1g,,)", [("1", "amp_001", "ump_001", "kmp_004")]),
            "z_005": (
                "Q(0,A1g,,)",
                [
                    ("sqrt(5)/5", "amp_002", "ump_001", "kmp_005"),
                    ("sqrt(5)/5", "amp_003", "ump_001", "kmp_006"),
                    ("sqrt(5)/5", "amp_004", "ump_001", "kmp_007"),
                    ("sqrt(5)/5", "amp_005", "ump_001", "kmp_008"),
                    ("sqrt(5)/5", "amp_006", "ump_001", "kmp_009"),
                ],
            ),
            "z_006": (
                "Q(4,A1g,,)",
                [
                    ("sqrt(30)/10", "amp_002", "ump_001", "kmp_005"),
                    ("sqrt(30)/10", "amp_003", "ump_001", "kmp_006"),
                    ("-sqrt(30)/15", "amp_004", "ump_001", "kmp_007"),
                    ("-sqrt(30)/15", "amp_005", "ump_001", "kmp_008"),
                    ("-sqrt(30)/15", "amp_006", "ump_001", "kmp_009"),
                ],
            ),
        },
    },
}
