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
MoS2 = {
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
            ],
            "M_002": ["amp_014", "amp_015", "amp_016", "amp_017", "amp_018", "amp_019", "amp_020", "amp_021", "amp_022"],
            "M_003": [
                "amp_023",
                "amp_024",
                "amp_025",
                "amp_026",
                "amp_027",
                "amp_028",
                "amp_029",
                "amp_030",
                "amp_031",
                "amp_032",
                "amp_033",
                "amp_034",
                "amp_035",
            ],
        },
        "site_cluster": {"S_001": ["smp_001"], "S_002": ["smp_002"]},
        "bond_cluster": {
            "B_001": ["bmp_003", "bmp_004", "bmp_005", "bmp_006", "bmp_007", "bmp_008"],
            "B_002": ["bmp_009", "bmp_010", "bmp_011", "bmp_012", "bmp_013", "bmp_014"],
            "B_003": ["bmp_015", "bmp_016", "bmp_017", "bmp_018", "bmp_019", "bmp_020", "bmp_021", "bmp_022"],
        },
        "Z": {
            ("A1'", "M_001", "S_001"): ["z_001", "z_002", "z_003"],
            ("A1'", "M_002", "S_002"): ["z_004", "z_005"],
            ("A1'", "M_001", "B_001"): ["z_006", "z_007", "z_008", "z_009", "z_010", "z_011", "z_012", "z_013", "z_014"],
            ("A1'", "M_003", "B_002"): ["z_015", "z_016", "z_017", "z_018", "z_019", "z_020", "z_021", "z_022"],
            ("A1'", "M_002", "B_003"): ["z_023", "z_024", "z_025", "z_026", "z_027", "z_028"],
        },
        "version": "1.1.14",
        "harmonics": {
            "Q": [
                "Qh(0,A1',,)",
                "Qh(1,A2'',,)",
                "Qh(1,E',,0)",
                "Qh(1,E',,1)",
                "Qh(2,A1',,)",
                "Qh(2,E'',,0)",
                "Qh(2,E'',,1)",
                "Qh(2,E',,0)",
                "Qh(2,E',,1)",
                "Qh(3,A1',,)",
                "Qh(3,A2'',,)",
                "Qh(3,A2',,)",
                "Qh(3,E'',,0)",
                "Qh(3,E'',,1)",
                "Qh(3,E',,0)",
                "Qh(3,E',,1)",
                "Qh(4,A1',,)",
                "Qh(4,E',1,0)",
                "Qh(4,E',1,1)",
                "Qh(4,E',2,0)",
                "Qh(4,E',2,1)",
            ],
            "G": [
                "Gh(1,A2',,)",
                "Gh(1,E'',,0)",
                "Gh(1,E'',,1)",
                "Gh(2,E'',,0)",
                "Gh(2,E'',,1)",
                "Gh(2,E',,0)",
                "Gh(2,E',,1)",
                "Gh(3,A2',,)",
                "Gh(3,E',,0)",
                "Gh(3,E',,1)",
            ],
        },
    },
    "data": {
        "atomic": {
            "amp_001": (
                "Qa(0,A1',,)",
                (5, 5),
                [(0, 0, "sqrt(5)/5"), (1, 1, "sqrt(5)/5"), (2, 2, "sqrt(5)/5"), (3, 3, "sqrt(5)/5"), (4, 4, "sqrt(5)/5")],
            ),
            "amp_002": (
                "Qa(2,A1',,)",
                (5, 5),
                [
                    (0, 0, "sqrt(14)/7"),
                    (1, 1, "-sqrt(14)/7"),
                    (2, 2, "sqrt(14)/14"),
                    (3, 3, "sqrt(14)/14"),
                    (4, 4, "-sqrt(14)/7"),
                ],
            ),
            "amp_003": (
                "Qa(4,A1',,)",
                (5, 5),
                [
                    (0, 0, "3*sqrt(70)/35"),
                    (1, 1, "sqrt(70)/70"),
                    (2, 2, "-2*sqrt(70)/35"),
                    (3, 3, "-2*sqrt(70)/35"),
                    (4, 4, "sqrt(70)/70"),
                ],
            ),
            "amp_004": (
                "Qa(2,E',,0)",
                (5, 5),
                [(0, 4, "sqrt(14)/7"), (2, 3, "-sqrt(42)/14"), (3, 2, "-sqrt(42)/14"), (4, 0, "sqrt(14)/7")],
            ),
            "amp_005": (
                "Qa(2,E',,1)",
                (5, 5),
                [(0, 1, "sqrt(14)/7"), (1, 0, "sqrt(14)/7"), (2, 2, "sqrt(42)/14"), (3, 3, "-sqrt(42)/14")],
            ),
            "amp_006": ("Qa(4,E',1,0)", (5, 5), [(1, 4, "sqrt(2)/2"), (4, 1, "sqrt(2)/2")]),
            "amp_007": ("Qa(4,E',1,1)", (5, 5), [(1, 1, "-sqrt(2)/2"), (4, 4, "sqrt(2)/2")]),
            "amp_008": (
                "Qa(4,E',2,0)",
                (5, 5),
                [(0, 4, "-sqrt(42)/14"), (2, 3, "-sqrt(14)/7"), (3, 2, "-sqrt(14)/7"), (4, 0, "-sqrt(42)/14")],
            ),
            "amp_009": (
                "Qa(4,E',2,1)",
                (5, 5),
                [(0, 1, "-sqrt(42)/14"), (1, 0, "-sqrt(42)/14"), (2, 2, "sqrt(14)/7"), (3, 3, "-sqrt(14)/7")],
            ),
            "amp_010": (
                "Ma(1,A2',,)",
                (5, 5),
                [(1, 4, "-sqrt(10)*I/5"), (2, 3, "sqrt(10)*I/10"), (3, 2, "-sqrt(10)*I/10"), (4, 1, "sqrt(10)*I/5")],
            ),
            "amp_011": (
                "Ma(3,A2',,)",
                (5, 5),
                [(1, 4, "sqrt(10)*I/10"), (2, 3, "sqrt(10)*I/5"), (3, 2, "-sqrt(10)*I/5"), (4, 1, "-sqrt(10)*I/10")],
            ),
            "amp_012": ("Ma(3,E',,0)", (5, 5), [(0, 4, "sqrt(2)*I/2"), (4, 0, "-sqrt(2)*I/2")]),
            "amp_013": ("Ma(3,E',,1)", (5, 5), [(0, 1, "sqrt(2)*I/2"), (1, 0, "-sqrt(2)*I/2")]),
            "amp_014": ("Qa(0,A1',,)", (3, 3), [(0, 0, "sqrt(3)/3"), (1, 1, "sqrt(3)/3"), (2, 2, "sqrt(3)/3")]),
            "amp_015": ("Qa(2,A1',,)", (3, 3), [(0, 0, "-sqrt(6)/6"), (1, 1, "-sqrt(6)/6"), (2, 2, "sqrt(6)/3")]),
            "amp_016": ("Qa(2,E',,0)", (3, 3), [(0, 1, "-sqrt(2)/2"), (1, 0, "-sqrt(2)/2")]),
            "amp_017": ("Qa(2,E',,1)", (3, 3), [(0, 0, "-sqrt(2)/2"), (1, 1, "sqrt(2)/2")]),
            "amp_018": ("Qa(2,E'',,0)", (3, 3), [(0, 2, "sqrt(2)/2"), (2, 0, "sqrt(2)/2")]),
            "amp_019": ("Qa(2,E'',,1)", (3, 3), [(1, 2, "sqrt(2)/2"), (2, 1, "sqrt(2)/2")]),
            "amp_020": ("Ma(1,A2',,)", (3, 3), [(0, 1, "-sqrt(2)*I/2"), (1, 0, "sqrt(2)*I/2")]),
            "amp_021": ("Ma(1,E'',,0)", (3, 3), [(0, 2, "-sqrt(2)*I/2"), (2, 0, "sqrt(2)*I/2")]),
            "amp_022": ("Ma(1,E'',,1)", (3, 3), [(1, 2, "-sqrt(2)*I/2"), (2, 1, "sqrt(2)*I/2")]),
            "amp_023": ("Qa(1,A2'',,)", (3, 5), [(0, 3, "sqrt(30)/10"), (1, 2, "sqrt(30)/10"), (2, 0, "sqrt(10)/5")]),
            "amp_024": ("Qa(3,A2'',,)", (3, 5), [(0, 3, "-sqrt(5)/5"), (1, 2, "-sqrt(5)/5"), (2, 0, "sqrt(15)/5")]),
            "amp_025": (
                "Qa(1,E',,0)",
                (3, 5),
                [(0, 0, "-sqrt(10)/10"), (0, 1, "sqrt(30)/10"), (1, 4, "sqrt(30)/10"), (2, 3, "sqrt(30)/10")],
            ),
            "amp_026": (
                "Qa(1,E',,1)",
                (3, 5),
                [(0, 4, "sqrt(30)/10"), (1, 0, "-sqrt(10)/10"), (1, 1, "-sqrt(30)/10"), (2, 2, "sqrt(30)/10")],
            ),
            "amp_027": (
                "Qa(3,E',,0)",
                (3, 5),
                [(0, 0, "sqrt(10)/5"), (0, 1, "-sqrt(30)/30"), (1, 4, "-sqrt(30)/30"), (2, 3, "2*sqrt(30)/15")],
            ),
            "amp_028": (
                "Qa(3,E',,1)",
                (3, 5),
                [(0, 4, "-sqrt(30)/30"), (1, 0, "sqrt(10)/5"), (1, 1, "sqrt(30)/30"), (2, 2, "2*sqrt(30)/15")],
            ),
            "amp_029": (
                "Ga(2,E',,0)",
                (3, 5),
                [(0, 0, "sqrt(2)/2"), (0, 1, "sqrt(6)/6"), (1, 4, "sqrt(6)/6"), (2, 3, "-sqrt(6)/6")],
            ),
            "amp_030": (
                "Ga(2,E',,1)",
                (3, 5),
                [(0, 4, "sqrt(6)/6"), (1, 0, "sqrt(2)/2"), (1, 1, "-sqrt(6)/6"), (2, 2, "-sqrt(6)/6")],
            ),
            "amp_031": ("Qa(3,A1',,)", (3, 5), [(0, 4, "sqrt(2)/2"), (1, 1, "sqrt(2)/2")]),
            "amp_032": ("Qa(3,E'',,0)", (3, 5), [(0, 2, "-sqrt(3)/3"), (1, 3, "-sqrt(3)/3"), (2, 4, "-sqrt(3)/3")]),
            "amp_033": ("Qa(3,E'',,1)", (3, 5), [(0, 3, "-sqrt(3)/3"), (1, 2, "sqrt(3)/3"), (2, 1, "-sqrt(3)/3")]),
            "amp_034": ("Ga(2,E'',,0)", (3, 5), [(0, 2, "sqrt(6)/6"), (1, 3, "sqrt(6)/6"), (2, 4, "-sqrt(6)/3")]),
            "amp_035": ("Ga(2,E'',,1)", (3, 5), [(0, 3, "sqrt(6)/6"), (1, 2, "-sqrt(6)/6"), (2, 1, "-sqrt(6)/3")]),
        },
        "site_cluster": {"smp_001": ("Qs(0,A1',,)", "[1]"), "smp_002": ("Qs(0,A1',,)", "[sqrt(2)/2, sqrt(2)/2]")},
        "bond_cluster": {
            "bmp_003": ("Qb(0,A1',,)", "[sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]"),
            "bmp_004": ("Qb(1,E',,0)", "[-sqrt(2)/2, sqrt(2)/2, 0]"),
            "bmp_005": ("Qb(1,E',,1)", "[-sqrt(6)/6, -sqrt(6)/6, sqrt(6)/3]"),
            "bmp_006": ("Tb(1,E',,0)", "[sqrt(6)*I/6, -sqrt(6)*I/6, -sqrt(6)*I/3]"),
            "bmp_007": ("Tb(1,E',,1)", "[-sqrt(2)*I/2, -sqrt(2)*I/2, 0]"),
            "bmp_008": ("Tb(3,A2',,)", "[sqrt(3)*I/3, -sqrt(3)*I/3, sqrt(3)*I/3]"),
            "bmp_009": ("Qb(0,A1',,)", "[sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6]"),
            "bmp_010": ("Qb(1,A2'',,)", "[sqrt(6)/6, -sqrt(6)/6, -sqrt(6)/6, -sqrt(6)/6, sqrt(6)/6, sqrt(6)/6]"),
            "bmp_011": ("Qb(1,E',,0)", "[-1/2, 1/2, -1/2, 0, 1/2, 0]"),
            "bmp_012": ("Qb(1,E',,1)", "[-sqrt(3)/6, -sqrt(3)/6, -sqrt(3)/6, sqrt(3)/3, -sqrt(3)/6, sqrt(3)/3]"),
            "bmp_013": ("Qb(2,E'',,0)", "[-1/2, -1/2, 1/2, 0, 1/2, 0]"),
            "bmp_014": ("Qb(2,E'',,1)", "[-sqrt(3)/6, sqrt(3)/6, sqrt(3)/6, -sqrt(3)/3, -sqrt(3)/6, sqrt(3)/3]"),
            "bmp_015": ("Qb(0,A1',,)", "[sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6]"),
            "bmp_016": ("Qb(1,E',,0)", "[0, 0, -1/2, 1/2, 1/2, -1/2]"),
            "bmp_017": ("Qb(1,E',,1)", "[-sqrt(3)/3, -sqrt(3)/3, sqrt(3)/6, sqrt(3)/6, sqrt(3)/6, sqrt(3)/6]"),
            "bmp_018": ("Qb(2,E'',,0)", "[0, 0, 1/2, -1/2, 1/2, -1/2]"),
            "bmp_019": ("Qb(2,E'',,1)", "[-sqrt(3)/3, sqrt(3)/3, -sqrt(3)/6, -sqrt(3)/6, sqrt(3)/6, sqrt(3)/6]"),
            "bmp_020": ("Tb(2,E'',,0)", "[-sqrt(3)*I/3, sqrt(3)*I/3, sqrt(3)*I/6, -sqrt(3)*I/6, sqrt(3)*I/6, -sqrt(3)*I/6]"),
            "bmp_021": ("Tb(2,E'',,1)", "[0, 0, I/2, I/2, -I/2, -I/2]"),
            "bmp_022": ("Tb(3,A2',,)", "[sqrt(6)*I/6, sqrt(6)*I/6, -sqrt(6)*I/6, sqrt(6)*I/6, sqrt(6)*I/6, -sqrt(6)*I/6]"),
        },
        "Z": {
            "z_001": ("Q(0,A1',,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(2,A1',,)", [("1", "amp_002", "smp_001")]),
            "z_003": ("Q(4,A1',,)", [("1", "amp_003", "smp_001")]),
            "z_004": ("Q(0,A1',,)", [("1", "amp_014", "smp_002")]),
            "z_005": ("Q(2,A1',,)", [("1", "amp_015", "smp_002")]),
            "z_006": ("Q(0,A1',,)", [("1", "amp_001", "bmp_003")]),
            "z_007": ("Q(2,A1',,)", [("1", "amp_002", "bmp_003")]),
            "z_008": ("Q(3,A1',,)", [("-sqrt(2)/2", "amp_004", "bmp_004"), ("-sqrt(2)/2", "amp_005", "bmp_005")]),
            "z_009": ("Q(4,A1',,)", [("1", "amp_003", "bmp_003")]),
            "z_010": (
                "Q(3,A1',,)",
                [
                    ("sqrt(406)/29", "amp_006", "bmp_004"),
                    ("sqrt(406)/29", "amp_007", "bmp_005"),
                    ("sqrt(58)/58", "amp_008", "bmp_004"),
                    ("sqrt(58)/58", "amp_009", "bmp_005"),
                ],
            ),
            "z_011": (
                "G(4,A1',,)",
                [
                    ("-sqrt(58)/58", "amp_006", "bmp_004"),
                    ("-sqrt(58)/58", "amp_007", "bmp_005"),
                    ("sqrt(406)/29", "amp_008", "bmp_004"),
                    ("sqrt(406)/29", "amp_009", "bmp_005"),
                ],
            ),
            "z_012": ("Q(3,A1',,)", [("1", "amp_010", "bmp_008")]),
            "z_013": ("Q(3,A1',,)", [("sqrt(2)/2", "amp_012", "bmp_006"), ("sqrt(2)/2", "amp_013", "bmp_007")]),
            "z_014": ("Q(3,A1',,)", [("-1", "amp_011", "bmp_008")]),
            "z_015": (
                "Q(0,A1',,)",
                [("sqrt(3)/3", "amp_023", "bmp_010"), ("sqrt(3)/3", "amp_025", "bmp_011"), ("sqrt(3)/3", "amp_026", "bmp_012")],
            ),
            "z_016": (
                "Q(2,A1',,)",
                [("sqrt(6)/3", "amp_023", "bmp_010"), ("-sqrt(6)/6", "amp_025", "bmp_011"), ("-sqrt(6)/6", "amp_026", "bmp_012")],
            ),
            "z_017": ("Q(3,A1',,)", [("1", "amp_031", "bmp_009")]),
            "z_018": (
                "Q(2,A1',,)",
                [
                    ("sqrt(21)/7", "amp_024", "bmp_010"),
                    ("sqrt(14)/7", "amp_027", "bmp_011"),
                    ("sqrt(14)/7", "amp_028", "bmp_012"),
                ],
            ),
            "z_019": (
                "Q(4,A1',,)",
                [
                    ("2*sqrt(7)/7", "amp_024", "bmp_010"),
                    ("-sqrt(42)/14", "amp_027", "bmp_011"),
                    ("-sqrt(42)/14", "amp_028", "bmp_012"),
                ],
            ),
            "z_020": ("Q(3,A1',,)", [("-sqrt(2)/2", "amp_032", "bmp_013"), ("-sqrt(2)/2", "amp_033", "bmp_014")]),
            "z_021": ("Q(2,A1',,)", [("sqrt(2)/2", "amp_029", "bmp_011"), ("sqrt(2)/2", "amp_030", "bmp_012")]),
            "z_022": ("Q(3,A1',,)", [("-sqrt(2)/2", "amp_034", "bmp_013"), ("-sqrt(2)/2", "amp_035", "bmp_014")]),
            "z_023": ("Q(0,A1',,)", [("1", "amp_014", "bmp_015")]),
            "z_024": ("Q(2,A1',,)", [("1", "amp_015", "bmp_015")]),
            "z_025": ("Q(3,A1',,)", [("-sqrt(2)/2", "amp_016", "bmp_016"), ("-sqrt(2)/2", "amp_017", "bmp_017")]),
            "z_026": ("Q(0,A1',,)", [("sqrt(2)/2", "amp_018", "bmp_018"), ("sqrt(2)/2", "amp_019", "bmp_019")]),
            "z_027": ("Q(2,A1',,)", [("sqrt(2)/2", "amp_021", "bmp_020"), ("sqrt(2)/2", "amp_022", "bmp_021")]),
            "z_028": ("Q(3,A1',,)", [("1", "amp_020", "bmp_022")]),
        },
    },
}