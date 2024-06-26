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
C3v1 = {
    "info": {
        "atomic": {
            "M_001": [
                "amp_001",
                "amp_002",
                "amp_003",
                "amp_004",
                "amp_005",
                "amp_006",
                "amp_007",
                "amp_008",
                "amp_009",
                "amp_010",
                "amp_011",
                "amp_012",
                "amp_013",
            ]
        },
        "site_cluster": {"S_001": ["smp_001"], "S_002": ["smp_002"]},
        "bond_cluster": {"B_001": ["bmp_003", "bmp_004", "bmp_005", "bmp_006", "bmp_007", "bmp_008"]},
        "Z": {
            ("A1", "M_001", "S_001"): ["z_001", "z_002"],
            ("A1", "M_001", "S_002"): ["z_003", "z_004"],
            ("A1", "M_001", "B_001"): ["z_005", "z_006", "z_007", "z_008", "z_009", "z_010", "z_011", "z_012"],
        },
        "version": "1.1.14",
        "harmonics": {
            "Q": ["Qh(0,A1,,)", "Qh(1,E,,0)", "Qh(1,E,,1)", "Qh(2,E,1,0)", "Qh(2,E,1,1)", "Qh(2,E,2,0)", "Qh(2,E,2,1)"],
            "G": ["Gh(1,E,,0)", "Gh(1,E,,1)", "Gh(3,A1,,)", "Gh(3,E,2,0)", "Gh(3,E,2,1)"],
        },
    },
    "data": {
        "atomic": {
            "amp_001": ("Qa(0,A1,,)", (4, 4), [(0, 0, "1/2"), (1, 1, "1/2"), (2, 2, "1/2"), (3, 3, "1/2")]),
            "amp_002": ("Qa(0,A1,,|1,1)", (4, 4), [(0, 2, "-I/2"), (1, 3, "I/2"), (2, 0, "I/2"), (3, 1, "-I/2")]),
            "amp_003": ("Qa(2,E,2,0)", (4, 4), [(0, 2, "-1/2"), (1, 3, "-1/2"), (2, 0, "-1/2"), (3, 1, "-1/2")]),
            "amp_004": ("Qa(2,E,2,1)", (4, 4), [(0, 0, "-1/2"), (1, 1, "-1/2"), (2, 2, "1/2"), (3, 3, "1/2")]),
            "amp_005": ("Qa(2,E,1,0|1,-1)", (4, 4), [(0, 3, "-I/2"), (1, 2, "-I/2"), (2, 1, "I/2"), (3, 0, "I/2")]),
            "amp_006": ("Qa(2,E,1,1|1,-1)", (4, 4), [(0, 3, "-1/2"), (1, 2, "1/2"), (2, 1, "1/2"), (3, 0, "-1/2")]),
            "amp_007": (
                "Ma(1,E,,0|1,1)",
                (4, 4),
                [
                    (0, 1, "-sqrt(19)*I/19"),
                    (0, 3, "-3*sqrt(19)/38"),
                    (1, 0, "sqrt(19)*I/19"),
                    (1, 2, "-3*sqrt(19)/38"),
                    (2, 1, "-3*sqrt(19)/38"),
                    (2, 3, "2*sqrt(19)*I/19"),
                    (3, 0, "-3*sqrt(19)/38"),
                    (3, 2, "-2*sqrt(19)*I/19"),
                ],
            ),
            "amp_008": (
                "Ma(1,E,,1|1,1)",
                (4, 4),
                [
                    (0, 1, "2*sqrt(19)/19"),
                    (0, 3, "-3*sqrt(19)*I/38"),
                    (1, 0, "2*sqrt(19)/19"),
                    (1, 2, "3*sqrt(19)*I/38"),
                    (2, 1, "-3*sqrt(19)*I/38"),
                    (2, 3, "-sqrt(19)/19"),
                    (3, 0, "3*sqrt(19)*I/38"),
                    (3, 2, "-sqrt(19)/19"),
                ],
            ),
            "amp_009": (
                "Ma(1,E,,0|1,-1)",
                (4, 4),
                [
                    (0, 1, "7*sqrt(38)*I/76"),
                    (0, 3, "sqrt(38)/76"),
                    (1, 0, "-7*sqrt(38)*I/76"),
                    (1, 2, "sqrt(38)/76"),
                    (2, 1, "sqrt(38)/76"),
                    (2, 3, "5*sqrt(38)*I/76"),
                    (3, 0, "sqrt(38)/76"),
                    (3, 2, "-5*sqrt(38)*I/76"),
                ],
            ),
            "amp_010": (
                "Ma(1,E,,1|1,-1)",
                (4, 4),
                [
                    (0, 1, "5*sqrt(38)/76"),
                    (0, 3, "sqrt(38)*I/76"),
                    (1, 0, "5*sqrt(38)/76"),
                    (1, 2, "-sqrt(38)*I/76"),
                    (2, 1, "sqrt(38)*I/76"),
                    (2, 3, "7*sqrt(38)/76"),
                    (3, 0, "-sqrt(38)*I/76"),
                    (3, 2, "7*sqrt(38)/76"),
                ],
            ),
            "amp_011": ("Ma(3,E,2,0|1,-1)", (4, 4), [(0, 0, "1/2"), (1, 1, "-1/2"), (2, 2, "-1/2"), (3, 3, "1/2")]),
            "amp_012": ("Ma(3,E,2,1|1,-1)", (4, 4), [(0, 2, "-1/2"), (1, 3, "1/2"), (2, 0, "-1/2"), (3, 1, "1/2")]),
            "amp_013": (
                "Ma(3,A1,,|1,-1)",
                (4, 4),
                [
                    (0, 1, "sqrt(2)/4"),
                    (0, 3, "sqrt(2)*I/4"),
                    (1, 0, "sqrt(2)/4"),
                    (1, 2, "-sqrt(2)*I/4"),
                    (2, 1, "sqrt(2)*I/4"),
                    (2, 3, "-sqrt(2)/4"),
                    (3, 0, "-sqrt(2)*I/4"),
                    (3, 2, "-sqrt(2)/4"),
                ],
            ),
        },
        "site_cluster": {"smp_001": ("Qs(0,A1,,)", "[1]"), "smp_002": ("Qs(0,A1,,)", "[1]")},
        "bond_cluster": {
            "bmp_003": ("Qb(0,A1,,)", "[sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]"),
            "bmp_004": ("Qb(1,E,,0)", "[0, sqrt(2)/2, -sqrt(2)/2]"),
            "bmp_005": ("Qb(1,E,,1)", "[-sqrt(6)/3, sqrt(6)/6, sqrt(6)/6]"),
            "bmp_006": ("Tb(0,A1,,)", "[sqrt(3)*I/3, sqrt(3)*I/3, sqrt(3)*I/3]"),
            "bmp_007": ("Tb(1,E,,0)", "[0, sqrt(2)*I/2, -sqrt(2)*I/2]"),
            "bmp_008": ("Tb(1,E,,1)", "[-sqrt(6)*I/3, sqrt(6)*I/6, sqrt(6)*I/6]"),
        },
        "Z": {
            "z_001": ("Q(0,A1,,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(0,A1,,|1,1)", [("1", "amp_002", "smp_001")]),
            "z_003": ("Q(0,A1,,)", [("1", "amp_001", "smp_002")]),
            "z_004": ("Q(0,A1,,|1,1)", [("1", "amp_002", "smp_002")]),
            "z_005": ("Q(0,A1,,)", [("1", "amp_001", "bmp_003")]),
            "z_006": ("Q(0,A1,,|1,1)", [("1", "amp_002", "bmp_003")]),
            "z_007": ("Q(3,A1,2,)", [("-sqrt(2)/2", "amp_003", "bmp_004"), ("-sqrt(2)/2", "amp_004", "bmp_005")]),
            "z_008": ("Q(1,A1,,|1,-1)", [("sqrt(2)/2", "amp_005", "bmp_004"), ("sqrt(2)/2", "amp_006", "bmp_005")]),
            "z_009": ("Q(1,A1,,|1,1)", [("sqrt(2)/2", "amp_007", "bmp_007"), ("sqrt(2)/2", "amp_008", "bmp_008")]),
            "z_010": ("G(3,A1,,|1,-1)", [("1", "amp_013", "bmp_006")]),
            "z_011": ("Q(3,A1,2,|1,-1)", [("sqrt(2)/2", "amp_011", "bmp_007"), ("sqrt(2)/2", "amp_012", "bmp_008")]),
            "z_012": ("Q(1,A1,,|1,-1)", [("sqrt(2)/2", "amp_009", "bmp_007"), ("sqrt(2)/2", "amp_010", "bmp_008")]),
        },
    },
}
