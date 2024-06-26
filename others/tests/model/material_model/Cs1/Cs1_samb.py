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
Cs1 = {
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
                "amp_014",
                "amp_015",
                "amp_016",
            ]
        },
        "site_cluster": {"S_001": ["smp_001"]},
        "bond_cluster": {
            "B_001": ["bmp_002", "bmp_003"],
            "B_002": ["bmp_004", "bmp_005"],
            "B_003": ["bmp_006", "bmp_007"],
            "B_004": ["bmp_008", "bmp_009"],
            "B_005": ["bmp_010", "bmp_011", "bmp_012", "bmp_013"],
            "B_006": ["bmp_014", "bmp_015", "bmp_016", "bmp_017"],
            "B_007": ["bmp_018", "bmp_019"],
        },
        "Z": {
            ("A'", "M_001", "S_001"): ["z_001", "z_002", "z_003", "z_004"],
            ("A'", "M_001", "B_001"): ["z_005", "z_006", "z_007", "z_008", "z_009", "z_010", "z_011", "z_012", "z_013", "z_014"],
            ("A'", "M_001", "B_002"): ["z_015", "z_016", "z_017", "z_018", "z_019", "z_020", "z_021", "z_022"],
            ("A'", "M_001", "B_003"): ["z_023", "z_024", "z_025", "z_026", "z_027", "z_028", "z_029", "z_030"],
            ("A'", "M_001", "B_004"): ["z_031", "z_032", "z_033", "z_034", "z_035", "z_036", "z_037", "z_038"],
            ("A'", "M_001", "B_005"): [
                "z_039",
                "z_040",
                "z_041",
                "z_042",
                "z_043",
                "z_044",
                "z_045",
                "z_046",
                "z_047",
                "z_048",
                "z_049",
                "z_050",
                "z_051",
                "z_052",
                "z_053",
                "z_054",
            ],
            ("A'", "M_001", "B_006"): [
                "z_055",
                "z_056",
                "z_057",
                "z_058",
                "z_059",
                "z_060",
                "z_061",
                "z_062",
                "z_063",
                "z_064",
                "z_065",
                "z_066",
                "z_067",
                "z_068",
                "z_069",
                "z_070",
            ],
            ("A'", "M_001", "B_007"): ["z_071", "z_072", "z_073", "z_074", "z_075", "z_076", "z_077", "z_078"],
        },
        "version": "1.1.14",
        "harmonics": {
            "Q": ["Qh(0,A',,)", "Qh(1,A'',,)", "Qh(2,A'',1,)", "Qh(2,A'',2,)", "Qh(2,A',2,)", "Qh(2,A',3,)"],
            "G": ["Gh(1,A'',1,)", "Gh(1,A'',2,)", "Gh(1,A',,)", "Gh(3,A'',1,)", "Gh(3,A'',4,)", "Gh(3,A',1,)", "Gh(3,A',2,)"],
        },
    },
    "data": {
        "atomic": {
            "amp_001": ("Qa(0,A',,)", (4, 4), [(0, 0, "1/2"), (1, 1, "1/2"), (2, 2, "1/2"), (3, 3, "1/2")]),
            "amp_002": ("Qa(2,A',2,)", (4, 4), [(0, 0, "1/2"), (1, 1, "1/2"), (2, 2, "-1/2"), (3, 3, "-1/2")]),
            "amp_003": ("Qa(0,A',,|1,1)", (4, 4), [(0, 2, "-I/2"), (1, 3, "I/2"), (2, 0, "I/2"), (3, 1, "-I/2")]),
            "amp_004": ("Qa(2,A',3,|1,-1)", (4, 4), [(0, 3, "-I/2"), (1, 2, "-I/2"), (2, 1, "I/2"), (3, 0, "I/2")]),
            "amp_005": ("Qa(2,A'',2,)", (4, 4), [(0, 2, "1/2"), (1, 3, "1/2"), (2, 0, "1/2"), (3, 1, "1/2")]),
            "amp_006": ("Qa(2,A'',1,|1,-1)", (4, 4), [(0, 3, "-1/2"), (1, 2, "1/2"), (2, 1, "1/2"), (3, 0, "-1/2")]),
            "amp_007": (
                "Ma(1,A',,|1,1)",
                (4, 4),
                [
                    (0, 1, "sqrt(19)*I/19"),
                    (0, 3, "3*sqrt(19)/38"),
                    (1, 0, "-sqrt(19)*I/19"),
                    (1, 2, "3*sqrt(19)/38"),
                    (2, 1, "3*sqrt(19)/38"),
                    (2, 3, "-2*sqrt(19)*I/19"),
                    (3, 0, "3*sqrt(19)/38"),
                    (3, 2, "2*sqrt(19)*I/19"),
                ],
            ),
            "amp_008": (
                "Ma(1,A',,|1,-1)",
                (4, 4),
                [
                    (0, 1, "-7*sqrt(38)*I/76"),
                    (0, 3, "-sqrt(38)/76"),
                    (1, 0, "7*sqrt(38)*I/76"),
                    (1, 2, "-sqrt(38)/76"),
                    (2, 1, "-sqrt(38)/76"),
                    (2, 3, "-5*sqrt(38)*I/76"),
                    (3, 0, "-sqrt(38)/76"),
                    (3, 2, "5*sqrt(38)*I/76"),
                ],
            ),
            "amp_009": ("Ma(3,A',1,|1,-1)", (4, 4), [(0, 2, "1/2"), (1, 3, "-1/2"), (2, 0, "1/2"), (3, 1, "-1/2")]),
            "amp_010": (
                "Ma(3,A',2,|1,-1)",
                (4, 4),
                [
                    (0, 1, "sqrt(2)*I/4"),
                    (0, 3, "-sqrt(2)/4"),
                    (1, 0, "-sqrt(2)*I/4"),
                    (1, 2, "-sqrt(2)/4"),
                    (2, 1, "-sqrt(2)/4"),
                    (2, 3, "-sqrt(2)*I/4"),
                    (3, 0, "-sqrt(2)/4"),
                    (3, 2, "sqrt(2)*I/4"),
                ],
            ),
            "amp_011": (
                "Ma(1,A'',1,|1,1)",
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
            "amp_012": ("Ma(1,A'',2,|1,1)", (4, 4), [(0, 0, "-1/2"), (1, 1, "1/2"), (2, 2, "-1/2"), (3, 3, "1/2")]),
            "amp_013": (
                "Ma(1,A'',1,|1,-1)",
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
            "amp_014": (
                "Ma(3,A'',1,|1,-1)",
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
            "amp_015": ("Ma(3,A'',4,|1,-1)", (4, 4), [(0, 0, "1/2"), (1, 1, "-1/2"), (2, 2, "-1/2"), (3, 3, "1/2")]),
            "amp_016": ("Ma(1,A'',2,)", (4, 4), [(0, 2, "-I/2"), (1, 3, "-I/2"), (2, 0, "I/2"), (3, 1, "I/2")]),
        },
        "site_cluster": {"smp_001": ("Qs(0,A',,)", "[1]")},
        "bond_cluster": {
            "bmp_002": ("Qb(0,A',,)", "[1]"),
            "bmp_003": ("Tb(1,A'',,)", "[I]"),
            "bmp_004": ("Qb(0,A',,)", "[1]"),
            "bmp_005": ("Tb(0,A',,)", "[I]"),
            "bmp_006": ("Qb(0,A',,)", "[1]"),
            "bmp_007": ("Tb(0,A',,)", "[I]"),
            "bmp_008": ("Qb(0,A',,)", "[1]"),
            "bmp_009": ("Tb(0,A',,)", "[I]"),
            "bmp_010": ("Qb(0,A',,)", "[sqrt(2)/2, sqrt(2)/2]"),
            "bmp_011": ("Qb(1,A'',,)", "[sqrt(2)/2, -sqrt(2)/2]"),
            "bmp_012": ("Tb(0,A',,)", "[sqrt(2)*I/2, sqrt(2)*I/2]"),
            "bmp_013": ("Tb(1,A'',,)", "[sqrt(2)*I/2, -sqrt(2)*I/2]"),
            "bmp_014": ("Qb(0,A',,)", "[sqrt(2)/2, sqrt(2)/2]"),
            "bmp_015": ("Qb(1,A'',,)", "[sqrt(2)/2, -sqrt(2)/2]"),
            "bmp_016": ("Tb(0,A',,)", "[sqrt(2)*I/2, sqrt(2)*I/2]"),
            "bmp_017": ("Tb(1,A'',,)", "[sqrt(2)*I/2, -sqrt(2)*I/2]"),
            "bmp_018": ("Qb(0,A',,)", "[1]"),
            "bmp_019": ("Tb(0,A',,)", "[I]"),
        },
        "Z": {
            "z_001": ("Q(0,A',,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(2,A',2,)", [("1", "amp_002", "smp_001")]),
            "z_003": ("Q(0,A',,|1,1)", [("1", "amp_003", "smp_001")]),
            "z_004": ("Q(2,A',3,|1,-1)", [("1", "amp_004", "smp_001")]),
            "z_005": ("Q(0,A',,)", [("1", "amp_001", "bmp_002")]),
            "z_006": ("Q(2,A',2,)", [("1", "amp_002", "bmp_002")]),
            "z_007": ("Q(0,A',,|1,1)", [("1", "amp_003", "bmp_002")]),
            "z_008": ("Q(2,A',3,|1,-1)", [("1", "amp_004", "bmp_002")]),
            "z_009": ("Q(1,A',1,|1,1)", [("-1", "amp_012", "bmp_003")]),
            "z_010": ("Q(1,A',2,|1,1)", [("1", "amp_011", "bmp_003")]),
            "z_011": ("Q(1,A',2,|1,-1)", [("1", "amp_013", "bmp_003")]),
            "z_012": ("G(2,A',1,|1,-1)", [("-1", "amp_015", "bmp_003")]),
            "z_013": ("G(2,A',2,|1,-1)", [("-1", "amp_014", "bmp_003")]),
            "z_014": ("Q(1,A',1,)", [("-1", "amp_016", "bmp_003")]),
            "z_015": ("Q(0,A',,)", [("1", "amp_001", "bmp_004")]),
            "z_016": ("Q(2,A',2,)", [("1", "amp_002", "bmp_004")]),
            "z_017": ("Q(0,A',,|1,1)", [("1", "amp_003", "bmp_004")]),
            "z_018": ("Q(2,A',3,|1,-1)", [("1", "amp_004", "bmp_004")]),
            "z_019": ("G(1,A',,|1,1)", [("1", "amp_007", "bmp_005")]),
            "z_020": ("G(1,A',,|1,-1)", [("1", "amp_008", "bmp_005")]),
            "z_021": ("G(3,A',1,|1,-1)", [("1", "amp_009", "bmp_005")]),
            "z_022": ("G(3,A',2,|1,-1)", [("1", "amp_010", "bmp_005")]),
            "z_023": ("Q(0,A',,)", [("1", "amp_001", "bmp_006")]),
            "z_024": ("Q(2,A',2,)", [("1", "amp_002", "bmp_006")]),
            "z_025": ("Q(0,A',,|1,1)", [("1", "amp_003", "bmp_006")]),
            "z_026": ("Q(2,A',3,|1,-1)", [("1", "amp_004", "bmp_006")]),
            "z_027": ("G(1,A',,|1,1)", [("1", "amp_007", "bmp_007")]),
            "z_028": ("G(1,A',,|1,-1)", [("1", "amp_008", "bmp_007")]),
            "z_029": ("G(3,A',1,|1,-1)", [("1", "amp_009", "bmp_007")]),
            "z_030": ("G(3,A',2,|1,-1)", [("1", "amp_010", "bmp_007")]),
            "z_031": ("Q(0,A',,)", [("1", "amp_001", "bmp_008")]),
            "z_032": ("Q(2,A',2,)", [("1", "amp_002", "bmp_008")]),
            "z_033": ("Q(0,A',,|1,1)", [("1", "amp_003", "bmp_008")]),
            "z_034": ("Q(2,A',3,|1,-1)", [("1", "amp_004", "bmp_008")]),
            "z_035": ("G(1,A',,|1,1)", [("1", "amp_007", "bmp_009")]),
            "z_036": ("G(1,A',,|1,-1)", [("1", "amp_008", "bmp_009")]),
            "z_037": ("G(3,A',1,|1,-1)", [("1", "amp_009", "bmp_009")]),
            "z_038": ("G(3,A',2,|1,-1)", [("1", "amp_010", "bmp_009")]),
            "z_039": ("Q(0,A',,)", [("1", "amp_001", "bmp_010")]),
            "z_040": ("Q(2,A',2,)", [("1", "amp_002", "bmp_010")]),
            "z_041": ("Q(1,A',1,)", [("1", "amp_005", "bmp_011")]),
            "z_042": ("Q(0,A',,|1,1)", [("1", "amp_003", "bmp_010")]),
            "z_043": ("Q(2,A',3,|1,-1)", [("1", "amp_004", "bmp_010")]),
            "z_044": ("Q(1,A',2,|1,-1)", [("1", "amp_006", "bmp_011")]),
            "z_045": ("G(1,A',,|1,1)", [("1", "amp_007", "bmp_012")]),
            "z_046": ("Q(1,A',1,|1,1)", [("-1", "amp_012", "bmp_013")]),
            "z_047": ("Q(1,A',2,|1,1)", [("1", "amp_011", "bmp_013")]),
            "z_048": ("G(1,A',,|1,-1)", [("1", "amp_008", "bmp_012")]),
            "z_049": ("Q(1,A',2,|1,-1)", [("1", "amp_013", "bmp_013")]),
            "z_050": ("G(3,A',1,|1,-1)", [("1", "amp_009", "bmp_012")]),
            "z_051": ("G(3,A',2,|1,-1)", [("1", "amp_010", "bmp_012")]),
            "z_052": ("G(2,A',1,|1,-1)", [("-1", "amp_015", "bmp_013")]),
            "z_053": ("G(2,A',2,|1,-1)", [("-1", "amp_014", "bmp_013")]),
            "z_054": ("Q(1,A',1,)", [("-1", "amp_016", "bmp_013")]),
            "z_055": ("Q(0,A',,)", [("1", "amp_001", "bmp_014")]),
            "z_056": ("Q(2,A',2,)", [("1", "amp_002", "bmp_014")]),
            "z_057": ("Q(1,A',1,)", [("1", "amp_005", "bmp_015")]),
            "z_058": ("Q(0,A',,|1,1)", [("1", "amp_003", "bmp_014")]),
            "z_059": ("Q(2,A',3,|1,-1)", [("1", "amp_004", "bmp_014")]),
            "z_060": ("Q(1,A',2,|1,-1)", [("1", "amp_006", "bmp_015")]),
            "z_061": ("G(1,A',,|1,1)", [("1", "amp_007", "bmp_016")]),
            "z_062": ("Q(1,A',1,|1,1)", [("-1", "amp_012", "bmp_017")]),
            "z_063": ("Q(1,A',2,|1,1)", [("1", "amp_011", "bmp_017")]),
            "z_064": ("G(1,A',,|1,-1)", [("1", "amp_008", "bmp_016")]),
            "z_065": ("Q(1,A',2,|1,-1)", [("1", "amp_013", "bmp_017")]),
            "z_066": ("G(3,A',1,|1,-1)", [("1", "amp_009", "bmp_016")]),
            "z_067": ("G(3,A',2,|1,-1)", [("1", "amp_010", "bmp_016")]),
            "z_068": ("G(2,A',1,|1,-1)", [("-1", "amp_015", "bmp_017")]),
            "z_069": ("G(2,A',2,|1,-1)", [("-1", "amp_014", "bmp_017")]),
            "z_070": ("Q(1,A',1,)", [("-1", "amp_016", "bmp_017")]),
            "z_071": ("Q(0,A',,)", [("1", "amp_001", "bmp_018")]),
            "z_072": ("Q(2,A',2,)", [("1", "amp_002", "bmp_018")]),
            "z_073": ("Q(0,A',,|1,1)", [("1", "amp_003", "bmp_018")]),
            "z_074": ("Q(2,A',3,|1,-1)", [("1", "amp_004", "bmp_018")]),
            "z_075": ("G(1,A',,|1,1)", [("1", "amp_007", "bmp_019")]),
            "z_076": ("G(1,A',,|1,-1)", [("1", "amp_008", "bmp_019")]),
            "z_077": ("G(3,A',1,|1,-1)", [("1", "amp_009", "bmp_019")]),
            "z_078": ("G(3,A',2,|1,-1)", [("1", "amp_010", "bmp_019")]),
        },
    },
}
