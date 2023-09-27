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
grapheneAB = {
    "info": {
        "atomic": {"M_001": ["amp_001"], "M_002": ["amp_002", "amp_003", "amp_004", "amp_005"], "M_003": ["amp_006", "amp_007"]},
        "site_cluster": {"S_001": ["smp_001"], "S_002": ["smp_002"]},
        "bond_cluster": {
            "B_001": ["bmp_003", "bmp_004"],
            "B_002": ["bmp_005"],
            "B_003": ["bmp_006", "bmp_007", "bmp_008", "bmp_009"],
        },
        "Z": {
            ("A1'", "M_001", "S_001"): ["z_001"],
            ("A1'", "M_002", "S_002"): ["z_002"],
            ("A1'", "M_003", "B_001"): ["z_003"],
            ("A1'", "M_001", "B_002"): ["z_004"],
            ("A1'", "M_002", "B_003"): ["z_005", "z_006", "z_007"],
        },
        "version": "1.1.14",
        "harmonics": {
            "Q": ["Qh(0,A1',,)", "Qh(1,E',,0)", "Qh(1,E',,1)", "Qh(2,E',,0)", "Qh(2,E',,1)", "Qh(3,A2',,)"],
            "G": ["Gh(1,A2',,)"],
        },
    },
    "data": {
        "atomic": {
            "amp_001": ("Qa(0,A1',,)", (1, 1), [(0, 0, "1")]),
            "amp_002": ("Qa(0,A1',,)", (2, 2), [(0, 0, "sqrt(2)/2"), (1, 1, "sqrt(2)/2")]),
            "amp_003": ("Qa(2,E',,0)", (2, 2), [(0, 1, "-sqrt(2)/2"), (1, 0, "-sqrt(2)/2")]),
            "amp_004": ("Qa(2,E',,1)", (2, 2), [(0, 0, "-sqrt(2)/2"), (1, 1, "sqrt(2)/2")]),
            "amp_005": ("Ma(1,A2',,)", (2, 2), [(0, 1, "-sqrt(2)*I/2"), (1, 0, "sqrt(2)*I/2")]),
            "amp_006": ("Qa(1,E',,0)", (2, 1), [(0, 0, "1")]),
            "amp_007": ("Qa(1,E',,1)", (2, 1), [(1, 0, "1")]),
        },
        "site_cluster": {"smp_001": ("Qs(0,A1',,)", "[1]"), "smp_002": ("Qs(0,A1',,)", "[1]")},
        "bond_cluster": {
            "bmp_003": ("Qb(1,E',,0)", "[0, -sqrt(2)/2, sqrt(2)/2]"),
            "bmp_004": ("Qb(1,E',,1)", "[-sqrt(6)/3, sqrt(6)/6, sqrt(6)/6]"),
            "bmp_005": ("Qb(0,A1',,)", "[sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]"),
            "bmp_006": ("Qb(0,A1',,)", "[sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]"),
            "bmp_007": ("Qb(1,E',,0)", "[0, -sqrt(2)/2, sqrt(2)/2]"),
            "bmp_008": ("Qb(1,E',,1)", "[-sqrt(6)/3, sqrt(6)/6, sqrt(6)/6]"),
            "bmp_009": ("Tb(3,A2',,)", "[sqrt(3)*I/3, -sqrt(3)*I/3, sqrt(3)*I/3]"),
        },
        "Z": {
            "z_001": ("Q(0,A1',,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(0,A1',,)", [("1", "amp_002", "smp_002")]),
            "z_003": ("Q(0,A1',,)", [("sqrt(2)/2", "amp_006", "bmp_003"), ("sqrt(2)/2", "amp_007", "bmp_004")]),
            "z_004": ("Q(0,A1',,)", [("1", "amp_001", "bmp_005")]),
            "z_005": ("Q(0,A1',,)", [("1", "amp_002", "bmp_006")]),
            "z_006": ("Q(3,A1',,)", [("-sqrt(2)/2", "amp_003", "bmp_007"), ("-sqrt(2)/2", "amp_004", "bmp_008")]),
            "z_007": ("Q(3,A1',,)", [("1", "amp_005", "bmp_009")]),
        },
    },
}
