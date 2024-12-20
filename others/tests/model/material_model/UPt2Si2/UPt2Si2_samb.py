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
UPt2Si2 = {
    "info": {
        "atomic": {
            "M_001": ["amp_001", "amp_002", "amp_003", "amp_004", "amp_005", "amp_006"],
            "M_002": ["amp_007", "amp_008", "amp_009", "amp_010"],
            "M_003": ["amp_011", "amp_012"],
            "M_004": [
                "amp_013",
                "amp_014",
                "amp_015",
                "amp_016",
                "amp_017",
                "amp_018",
                "amp_019",
                "amp_020",
                "amp_021",
                "amp_022",
                "amp_023",
                "amp_024",
            ],
        },
        "site_cluster": {
            "S_001": ["smp_001"],
            "S_002": ["smp_002"],
            "S_003": ["smp_003"],
            "S_004": ["smp_004"],
            "S_005": ["smp_005"],
        },
        "bond_cluster": {
            "B_001": ["bmp_006", "bmp_007", "bmp_008", "bmp_009"],
            "B_002": ["bmp_010", "bmp_011", "bmp_012", "bmp_013"],
        },
        "Z": {
            ("A1g", "M_001", "S_001"): ["z_001", "z_002", "z_003", "z_004", "z_005", "z_006"],
            ("A1g", "M_002", "S_002"): ["z_007", "z_008", "z_009", "z_010"],
            ("A1g", "M_002", "S_003"): ["z_011", "z_012", "z_013", "z_014"],
            ("A1g", "M_003", "S_004"): ["z_015", "z_016"],
            ("A1g", "M_003", "S_005"): ["z_017", "z_018"],
            ("A1g", "M_004", "B_001"): ["z_019", "z_020", "z_021", "z_022", "z_023", "z_024", "z_025", "z_026"],
            ("A1g", "M_004", "B_002"): ["z_027", "z_028", "z_029", "z_030", "z_031", "z_032", "z_033", "z_034"],
        },
        "version": "1.1.14",
        "harmonics": {
            "Q": [
                "Qh(0,A1g,,)",
                "Qh(1,A2u,,)",
                "Qh(1,Eu,,0)",
                "Qh(1,Eu,,1)",
                "Qh(2,A1g,,)",
                "Qh(3,A2u,,)",
                "Qh(3,B2u,,)",
                "Qh(3,Eu,1,0)",
                "Qh(3,Eu,1,1)",
                "Qh(3,Eu,2,0)",
                "Qh(3,Eu,2,1)",
                "Qh(4,A1g,1,)",
                "Qh(4,A1g,2,)",
                "Qh(6,A1g,1,)",
                "Qh(6,A1g,2,)",
            ],
            "G": ["Gh(2,B2u,,)", "Gh(2,Eu,,0)", "Gh(2,Eu,,1)"],
        },
    },
    "data": {
        "atomic": {
            "amp_001": (
                "Qa(0,A1g,,)",
                (7, 7),
                [
                    (0, 0, "sqrt(7)/7"),
                    (1, 1, "sqrt(7)/7"),
                    (2, 2, "sqrt(7)/7"),
                    (3, 3, "sqrt(7)/7"),
                    (4, 4, "sqrt(7)/7"),
                    (5, 5, "sqrt(7)/7"),
                    (6, 6, "sqrt(7)/7"),
                ],
            ),
            "amp_002": (
                "Qa(2,A1g,,)",
                (7, 7),
                [
                    (1, 1, "-sqrt(21)/21"),
                    (1, 4, "sqrt(35)/14"),
                    (2, 2, "-sqrt(21)/21"),
                    (2, 5, "-sqrt(35)/14"),
                    (3, 3, "2*sqrt(21)/21"),
                    (4, 1, "sqrt(35)/14"),
                    (5, 2, "-sqrt(35)/14"),
                ],
            ),
            "amp_003": (
                "Qa(4,A1g,1,)",
                (7, 7),
                [
                    (0, 0, "-sqrt(66)/11"),
                    (1, 1, "sqrt(66)/22"),
                    (2, 2, "sqrt(66)/22"),
                    (3, 3, "sqrt(66)/22"),
                    (4, 4, "-sqrt(66)/66"),
                    (5, 5, "-sqrt(66)/66"),
                    (6, 6, "-sqrt(66)/66"),
                ],
            ),
            "amp_004": (
                "Qa(4,A1g,2,)",
                (7, 7),
                [
                    (1, 1, "-sqrt(2310)/308"),
                    (1, 4, "-3*sqrt(154)/308"),
                    (2, 2, "-sqrt(2310)/308"),
                    (2, 5, "3*sqrt(154)/308"),
                    (3, 3, "sqrt(2310)/154"),
                    (4, 1, "-3*sqrt(154)/308"),
                    (4, 4, "sqrt(2310)/132"),
                    (5, 2, "3*sqrt(154)/308"),
                    (5, 5, "sqrt(2310)/132"),
                    (6, 6, "-sqrt(2310)/66"),
                ],
            ),
            "amp_005": (
                "Qa(6,A1g,1,)",
                (7, 7),
                [
                    (0, 0, "2*sqrt(462)/77"),
                    (1, 1, "5*sqrt(462)/462"),
                    (2, 2, "5*sqrt(462)/462"),
                    (3, 3, "5*sqrt(462)/462"),
                    (4, 4, "-3*sqrt(462)/154"),
                    (5, 5, "-3*sqrt(462)/154"),
                    (6, 6, "-3*sqrt(462)/154"),
                ],
            ),
            "amp_006": (
                "Qa(6,A1g,2,)",
                (7, 7),
                [
                    (1, 1, "-5*sqrt(66)/132"),
                    (1, 4, "-sqrt(110)/44"),
                    (2, 2, "-5*sqrt(66)/132"),
                    (2, 5, "sqrt(110)/44"),
                    (3, 3, "5*sqrt(66)/66"),
                    (4, 1, "-sqrt(110)/44"),
                    (4, 4, "-sqrt(66)/44"),
                    (5, 2, "sqrt(110)/44"),
                    (5, 5, "-sqrt(66)/44"),
                    (6, 6, "sqrt(66)/22"),
                ],
            ),
            "amp_007": (
                "Qa(0,A1g,,)",
                (5, 5),
                [(0, 0, "sqrt(5)/5"), (1, 1, "sqrt(5)/5"), (2, 2, "sqrt(5)/5"), (3, 3, "sqrt(5)/5"), (4, 4, "sqrt(5)/5")],
            ),
            "amp_008": (
                "Qa(2,A1g,,)",
                (5, 5),
                [
                    (0, 0, "sqrt(14)/7"),
                    (1, 1, "-sqrt(14)/7"),
                    (2, 2, "sqrt(14)/14"),
                    (3, 3, "sqrt(14)/14"),
                    (4, 4, "-sqrt(14)/7"),
                ],
            ),
            "amp_009": (
                "Qa(4,A1g,1,)",
                (5, 5),
                [
                    (0, 0, "sqrt(30)/10"),
                    (1, 1, "sqrt(30)/10"),
                    (2, 2, "-sqrt(30)/15"),
                    (3, 3, "-sqrt(30)/15"),
                    (4, 4, "-sqrt(30)/15"),
                ],
            ),
            "amp_010": (
                "Qa(4,A1g,2,)",
                (5, 5),
                [
                    (0, 0, "sqrt(42)/14"),
                    (1, 1, "-sqrt(42)/14"),
                    (2, 2, "-sqrt(42)/21"),
                    (3, 3, "-sqrt(42)/21"),
                    (4, 4, "2*sqrt(42)/21"),
                ],
            ),
            "amp_011": ("Qa(0,A1g,,)", (3, 3), [(0, 0, "sqrt(3)/3"), (1, 1, "sqrt(3)/3"), (2, 2, "sqrt(3)/3")]),
            "amp_012": ("Qa(2,A1g,,)", (3, 3), [(0, 0, "-sqrt(6)/6"), (1, 1, "-sqrt(6)/6"), (2, 2, "sqrt(6)/3")]),
            "amp_013": ("Qa(1,A2u,,)", (3, 5), [(0, 3, "sqrt(30)/10"), (1, 2, "sqrt(30)/10"), (2, 0, "sqrt(10)/5")]),
            "amp_014": ("Qa(3,A2u,,)", (3, 5), [(0, 3, "-sqrt(5)/5"), (1, 2, "-sqrt(5)/5"), (2, 0, "sqrt(15)/5")]),
            "amp_015": (
                "Qa(1,Eu,,0)",
                (3, 5),
                [(0, 0, "-sqrt(10)/10"), (0, 1, "sqrt(30)/10"), (1, 4, "sqrt(30)/10"), (2, 3, "sqrt(30)/10")],
            ),
            "amp_016": (
                "Qa(1,Eu,,1)",
                (3, 5),
                [(0, 4, "sqrt(30)/10"), (1, 0, "-sqrt(10)/10"), (1, 1, "-sqrt(30)/10"), (2, 2, "sqrt(30)/10")],
            ),
            "amp_017": (
                "Qa(3,Eu,1,0)",
                (3, 5),
                [(0, 0, "-sqrt(15)/10"), (0, 1, "3*sqrt(5)/10"), (1, 4, "-sqrt(5)/5"), (2, 3, "-sqrt(5)/5")],
            ),
            "amp_018": (
                "Qa(3,Eu,1,1)",
                (3, 5),
                [(0, 4, "-sqrt(5)/5"), (1, 0, "-sqrt(15)/10"), (1, 1, "-3*sqrt(5)/10"), (2, 2, "-sqrt(5)/5")],
            ),
            "amp_019": (
                "Qa(3,Eu,2,0)",
                (3, 5),
                [(0, 0, "-1/2"), (0, 1, "-sqrt(3)/6"), (1, 4, "sqrt(3)/3"), (2, 3, "-sqrt(3)/3")],
            ),
            "amp_020": ("Qa(3,Eu,2,1)", (3, 5), [(0, 4, "sqrt(3)/3"), (1, 0, "-1/2"), (1, 1, "sqrt(3)/6"), (2, 2, "-sqrt(3)/3")]),
            "amp_021": (
                "Ga(2,Eu,,0)",
                (3, 5),
                [(0, 0, "-sqrt(2)/2"), (0, 1, "-sqrt(6)/6"), (1, 4, "-sqrt(6)/6"), (2, 3, "sqrt(6)/6")],
            ),
            "amp_022": (
                "Ga(2,Eu,,1)",
                (3, 5),
                [(0, 4, "sqrt(6)/6"), (1, 0, "sqrt(2)/2"), (1, 1, "-sqrt(6)/6"), (2, 2, "-sqrt(6)/6")],
            ),
            "amp_023": ("Qa(3,B2u,,)", (3, 5), [(0, 3, "sqrt(3)/3"), (1, 2, "-sqrt(3)/3"), (2, 1, "sqrt(3)/3")]),
            "amp_024": ("Ga(2,B2u,,)", (3, 5), [(0, 3, "-sqrt(6)/6"), (1, 2, "sqrt(6)/6"), (2, 1, "sqrt(6)/3")]),
        },
        "site_cluster": {
            "smp_001": ("Qs(0,A1g,,)", "[sqrt(2)/2, sqrt(2)/2]"),
            "smp_002": ("Qs(0,A1g,,)", "[sqrt(2)/2, sqrt(2)/2]"),
            "smp_003": ("Qs(0,A1g,,)", "[sqrt(2)/2, sqrt(2)/2]"),
            "smp_004": ("Qs(0,A1g,,)", "[sqrt(2)/2, sqrt(2)/2]"),
            "smp_005": ("Qs(0,A1g,,)", "[sqrt(2)/2, sqrt(2)/2]"),
        },
        "bond_cluster": {
            "bmp_006": (
                "Qb(1,A2u,,)",
                "[sqrt(2)/4, sqrt(2)/4, -sqrt(2)/4, -sqrt(2)/4, -sqrt(2)/4, -sqrt(2)/4, sqrt(2)/4, sqrt(2)/4]",
            ),
            "bmp_007": ("Qb(1,Eu,,0)", "[1/2, -1/2, 1/2, -1/2, 0, 0, 0, 0]"),
            "bmp_008": ("Qb(1,Eu,,1)", "[0, 0, 0, 0, 1/2, -1/2, 1/2, -1/2]"),
            "bmp_009": (
                "Qb(3,B2u,,)",
                "[sqrt(2)/4, sqrt(2)/4, -sqrt(2)/4, -sqrt(2)/4, sqrt(2)/4, sqrt(2)/4, -sqrt(2)/4, -sqrt(2)/4]",
            ),
            "bmp_010": (
                "Qb(1,A2u,,)",
                "[sqrt(2)/4, sqrt(2)/4, -sqrt(2)/4, -sqrt(2)/4, -sqrt(2)/4, -sqrt(2)/4, sqrt(2)/4, sqrt(2)/4]",
            ),
            "bmp_011": ("Qb(1,Eu,,0)", "[1/2, -1/2, 1/2, -1/2, 0, 0, 0, 0]"),
            "bmp_012": ("Qb(1,Eu,,1)", "[0, 0, 0, 0, 1/2, -1/2, 1/2, -1/2]"),
            "bmp_013": (
                "Qb(3,B2u,,)",
                "[sqrt(2)/4, sqrt(2)/4, -sqrt(2)/4, -sqrt(2)/4, sqrt(2)/4, sqrt(2)/4, -sqrt(2)/4, -sqrt(2)/4]",
            ),
        },
        "Z": {
            "z_001": ("Q(0,A1g,,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(2,A1g,,)", [("1", "amp_002", "smp_001")]),
            "z_003": ("Q(4,A1g,1,)", [("1", "amp_003", "smp_001")]),
            "z_004": ("Q(4,A1g,2,)", [("1", "amp_004", "smp_001")]),
            "z_005": ("Q(6,A1g,1,)", [("1", "amp_005", "smp_001")]),
            "z_006": ("Q(6,A1g,2,)", [("1", "amp_006", "smp_001")]),
            "z_007": ("Q(0,A1g,,)", [("1", "amp_007", "smp_002")]),
            "z_008": ("Q(2,A1g,,)", [("1", "amp_008", "smp_002")]),
            "z_009": ("Q(4,A1g,1,)", [("1", "amp_009", "smp_002")]),
            "z_010": ("Q(4,A1g,2,)", [("1", "amp_010", "smp_002")]),
            "z_011": ("Q(0,A1g,,)", [("1", "amp_007", "smp_003")]),
            "z_012": ("Q(2,A1g,,)", [("1", "amp_008", "smp_003")]),
            "z_013": ("Q(4,A1g,1,)", [("1", "amp_009", "smp_003")]),
            "z_014": ("Q(4,A1g,2,)", [("1", "amp_010", "smp_003")]),
            "z_015": ("Q(0,A1g,,)", [("1", "amp_011", "smp_004")]),
            "z_016": ("Q(2,A1g,,)", [("1", "amp_012", "smp_004")]),
            "z_017": ("Q(0,A1g,,)", [("1", "amp_011", "smp_005")]),
            "z_018": ("Q(2,A1g,,)", [("1", "amp_012", "smp_005")]),
            "z_019": (
                "Q(0,A1g,,)",
                [("sqrt(3)/3", "amp_013", "bmp_006"), ("sqrt(3)/3", "amp_015", "bmp_007"), ("sqrt(3)/3", "amp_016", "bmp_008")],
            ),
            "z_020": (
                "Q(2,A1g,,)",
                [("sqrt(6)/3", "amp_013", "bmp_006"), ("-sqrt(6)/6", "amp_015", "bmp_007"), ("-sqrt(6)/6", "amp_016", "bmp_008")],
            ),
            "z_021": (
                "Q(2,A1g,,)",
                [
                    ("sqrt(21)/7", "amp_014", "bmp_006"),
                    ("-sqrt(21)/14", "amp_017", "bmp_007"),
                    ("-sqrt(21)/14", "amp_018", "bmp_008"),
                    ("-sqrt(35)/14", "amp_019", "bmp_007"),
                    ("-sqrt(35)/14", "amp_020", "bmp_008"),
                ],
            ),
            "z_022": (
                "Q(4,A1g,1,)",
                [("sqrt(3)/3", "amp_014", "bmp_006"), ("sqrt(3)/3", "amp_017", "bmp_007"), ("sqrt(3)/3", "amp_018", "bmp_008")],
            ),
            "z_023": (
                "Q(4,A1g,2,)",
                [
                    ("sqrt(105)/21", "amp_014", "bmp_006"),
                    ("-sqrt(105)/42", "amp_017", "bmp_007"),
                    ("-sqrt(105)/42", "amp_018", "bmp_008"),
                    ("3*sqrt(7)/14", "amp_019", "bmp_007"),
                    ("3*sqrt(7)/14", "amp_020", "bmp_008"),
                ],
            ),
            "z_024": ("Q(0,A1g,,)", [("1", "amp_023", "bmp_009")]),
            "z_025": ("Q(2,A1g,,)", [("-sqrt(2)/2", "amp_021", "bmp_007"), ("sqrt(2)/2", "amp_022", "bmp_008")]),
            "z_026": ("Q(2,A1g,,)", [("-1", "amp_024", "bmp_009")]),
            "z_027": (
                "Q(0,A1g,,)",
                [("sqrt(3)/3", "amp_013", "bmp_010"), ("sqrt(3)/3", "amp_015", "bmp_011"), ("sqrt(3)/3", "amp_016", "bmp_012")],
            ),
            "z_028": (
                "Q(2,A1g,,)",
                [("sqrt(6)/3", "amp_013", "bmp_010"), ("-sqrt(6)/6", "amp_015", "bmp_011"), ("-sqrt(6)/6", "amp_016", "bmp_012")],
            ),
            "z_029": (
                "Q(2,A1g,,)",
                [
                    ("sqrt(21)/7", "amp_014", "bmp_010"),
                    ("-sqrt(21)/14", "amp_017", "bmp_011"),
                    ("-sqrt(21)/14", "amp_018", "bmp_012"),
                    ("-sqrt(35)/14", "amp_019", "bmp_011"),
                    ("-sqrt(35)/14", "amp_020", "bmp_012"),
                ],
            ),
            "z_030": (
                "Q(4,A1g,1,)",
                [("sqrt(3)/3", "amp_014", "bmp_010"), ("sqrt(3)/3", "amp_017", "bmp_011"), ("sqrt(3)/3", "amp_018", "bmp_012")],
            ),
            "z_031": (
                "Q(4,A1g,2,)",
                [
                    ("sqrt(105)/21", "amp_014", "bmp_010"),
                    ("-sqrt(105)/42", "amp_017", "bmp_011"),
                    ("-sqrt(105)/42", "amp_018", "bmp_012"),
                    ("3*sqrt(7)/14", "amp_019", "bmp_011"),
                    ("3*sqrt(7)/14", "amp_020", "bmp_012"),
                ],
            ),
            "z_032": ("Q(0,A1g,,)", [("1", "amp_023", "bmp_013")]),
            "z_033": ("Q(2,A1g,,)", [("-sqrt(2)/2", "amp_021", "bmp_011"), ("sqrt(2)/2", "amp_022", "bmp_012")]),
            "z_034": ("Q(2,A1g,,)", [("-1", "amp_024", "bmp_013")]),
        },
    },
}
