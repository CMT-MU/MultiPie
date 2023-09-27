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
Te = {
    "info": {
        "atomic": {"M_001": ["amp_001", "amp_002", "amp_003", "amp_004", "amp_005", "amp_006"]},
        "site_cluster": {"S_001": ["smp_001", "smp_002", "smp_003"]},
        "bond_cluster": {"B_001": ["bmp_004", "bmp_005", "bmp_006"]},
        "Z": {
            ("A1", "M_001", "S_001"): ["z_001", "z_002", "z_003", "z_004"],
            ("A1", "M_001", "B_001"): ["z_005", "z_006", "z_007", "z_008"],
        },
        "version": "1.1.15",
        "harmonics": {
            "Q": [
                "Qh(0,A1,,)",
                "Qh(1,E,,0)",
                "Qh(1,E,,1)",
                "Qh(2,A1,,)",
                "Qh(2,E,1,0)",
                "Qh(2,E,1,1)",
                "Qh(2,E,2,0)",
                "Qh(2,E,2,1)",
            ],
            "G": [],
        },
    },
    "data": {
        "atomic": {
            "amp_001": ("Qa(0,A1,,)", (3, 3), [(0, 0, "sqrt(3)/3"), (1, 1, "sqrt(3)/3"), (2, 2, "sqrt(3)/3")]),
            "amp_002": ("Qa(2,A1,,)", (3, 3), [(0, 0, "-sqrt(6)/6"), (1, 1, "-sqrt(6)/6"), (2, 2, "sqrt(6)/3")]),
            "amp_003": ("Qa(2,E,1,0)", (3, 3), [(1, 2, "sqrt(2)/2"), (2, 1, "sqrt(2)/2")]),
            "amp_004": ("Qa(2,E,1,1)", (3, 3), [(0, 2, "-sqrt(2)/2"), (2, 0, "-sqrt(2)/2")]),
            "amp_005": ("Qa(2,E,2,0)", (3, 3), [(0, 0, "sqrt(2)/2"), (1, 1, "-sqrt(2)/2")]),
            "amp_006": ("Qa(2,E,2,1)", (3, 3), [(0, 1, "-sqrt(2)/2"), (1, 0, "-sqrt(2)/2")]),
        },
        "site_cluster": {
            "smp_001": ("Qs(0,A1,,)", "[sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]"),
            "smp_002": ("Qs(1,E,,0)", "[sqrt(6)/3, -sqrt(6)/6, -sqrt(6)/6]"),
            "smp_003": ("Qs(1,E,,1)", "[0, -sqrt(2)/2, sqrt(2)/2]"),
        },
        "bond_cluster": {
            "bmp_004": ("Qb(0,A1,,)", "[sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]"),
            "bmp_005": ("Qb(1,E,,0)", "[sqrt(6)/6, sqrt(6)/6, -sqrt(6)/3]"),
            "bmp_006": ("Qb(1,E,,1)", "[-sqrt(2)/2, sqrt(2)/2, 0]"),
        },
        "Z": {
            "z_001": ("Q(0,A1,,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(2,A1,,)", [("1", "amp_002", "smp_001")]),
            "z_003": ("G(2,A1,,)", [("-sqrt(2)/2", "amp_003", "smp_002"), ("-sqrt(2)/2", "amp_004", "smp_003")]),
            "z_004": ("Q(3,A1,,)", [("sqrt(2)/2", "amp_005", "smp_002"), ("sqrt(2)/2", "amp_006", "smp_003")]),
            "z_005": ("Q(0,A1,,)", [("1", "amp_001", "bmp_004")]),
            "z_006": ("Q(2,A1,,)", [("1", "amp_002", "bmp_004")]),
            "z_007": ("G(2,A1,,)", [("-sqrt(2)/2", "amp_003", "bmp_005"), ("-sqrt(2)/2", "amp_004", "bmp_006")]),
            "z_008": ("Q(3,A1,,)", [("sqrt(2)/2", "amp_005", "bmp_005"), ("sqrt(2)/2", "amp_006", "bmp_006")]),
        },
    },
}
