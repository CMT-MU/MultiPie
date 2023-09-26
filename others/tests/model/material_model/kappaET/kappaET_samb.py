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
kappaET = {
    "info": {
        "atomic": {"M_001": ["amp_001", "amp_002", "amp_003", "amp_004"]},
        "site_cluster": {"S_001": ["smp_001"]},
        "bond_cluster": {
            "B_001": ["bmp_002", "bmp_003", "bmp_004"],
            "B_002": ["bmp_005", "bmp_006", "bmp_007", "bmp_008"],
            "B_003": ["bmp_009", "bmp_010", "bmp_011", "bmp_012"],
        },
        "Z": {
            ("A1", "M_001", "S_001"): ["z_001"],
            ("A1", "M_001", "B_001"): ["z_002", "z_003", "z_004"],
            ("A1", "M_001", "B_002"): ["z_005", "z_006", "z_007", "z_008"],
            ("A1", "M_001", "B_003"): ["z_009", "z_010", "z_011", "z_012"],
        },
        "version": "1.1.14",
        "harmonics": {
            "Q": ["Qh(0,A1,,)", "Qh(1,B1,,)", "Qh(1,B2,,)", "Qh(2,A2,,)"],
            "G": ["Gh(1,A2,,)", "Gh(1,B1,,)", "Gh(1,B2,,)"],
        },
    },
    "data": {
        "atomic": {
            "amp_001": ("Qa(0,A1,,)", (2, 2), [(0, 0, "sqrt(2)/2"), (1, 1, "sqrt(2)/2")]),
            "amp_002": ("Ma(1,A2,,|1,-1)", (2, 2), [(0, 0, "sqrt(2)/2"), (1, 1, "-sqrt(2)/2")]),
            "amp_003": ("Ma(1,B1,,|1,-1)", (2, 2), [(0, 1, "-sqrt(2)*I/2"), (1, 0, "sqrt(2)*I/2")]),
            "amp_004": ("Ma(1,B2,,|1,-1)", (2, 2), [(0, 1, "sqrt(2)/2"), (1, 0, "sqrt(2)/2")]),
        },
        "site_cluster": {"smp_001": ("Qs(0,A1,,)", "[1/2, 1/2, 1/2, 1/2]")},
        "bond_cluster": {
            "bmp_002": ("Qb(0,A1,,)", "[sqrt(2)/2, sqrt(2)/2]"),
            "bmp_003": ("Tb(1,B1,,)", "[sqrt(2)*I/2, sqrt(2)*I/2]"),
            "bmp_004": ("Tb(1,B2,,)", "[sqrt(2)*I/2, -sqrt(2)*I/2]"),
            "bmp_005": ("Qb(0,A1,,)", "[1/2, 1/2, 1/2, 1/2]"),
            "bmp_006": ("Tb(1,B1,,)", "[I/2, -I/2, -I/2, I/2]"),
            "bmp_007": ("Tb(1,B2,,)", "[I/2, -I/2, I/2, -I/2]"),
            "bmp_008": ("Tb(2,A2,,)", "[I/2, I/2, I/2, I/2]"),
            "bmp_009": ("Qb(0,A1,,)", "[1/2, 1/2, 1/2, 1/2]"),
            "bmp_010": ("Tb(1,B1,,)", "[I/2, -I/2, -I/2, I/2]"),
            "bmp_011": ("Tb(1,B2,,)", "[I/2, -I/2, I/2, -I/2]"),
            "bmp_012": ("Tb(2,A2,,)", "[I/2, I/2, I/2, I/2]"),
        },
        "Z": {
            "z_001": ("Q(0,A1,,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(0,A1,,)", [("1", "amp_001", "bmp_002")]),
            "z_003": ("Q(1,A1,,|1,-1)", [("-sqrt(2)/2", "amp_003", "bmp_003"), ("sqrt(2)/2", "amp_004", "bmp_004")]),
            "z_004": ("G(2,A1,,|1,-1)", [("sqrt(2)/2", "amp_003", "bmp_003"), ("sqrt(2)/2", "amp_004", "bmp_004")]),
            "z_005": ("Q(0,A1,,)", [("1", "amp_001", "bmp_005")]),
            "z_006": ("Q(1,A1,,|1,-1)", [("-sqrt(2)/2", "amp_003", "bmp_006"), ("sqrt(2)/2", "amp_004", "bmp_007")]),
            "z_007": ("G(2,A1,,|1,-1)", [("sqrt(2)/2", "amp_003", "bmp_006"), ("sqrt(2)/2", "amp_004", "bmp_007")]),
            "z_008": ("Q(2,A1,2,|1,-1)", [("-1", "amp_002", "bmp_008")]),
            "z_009": ("Q(0,A1,,)", [("1", "amp_001", "bmp_009")]),
            "z_010": ("Q(1,A1,,|1,-1)", [("-sqrt(2)/2", "amp_003", "bmp_010"), ("sqrt(2)/2", "amp_004", "bmp_011")]),
            "z_011": ("G(2,A1,,|1,-1)", [("sqrt(2)/2", "amp_003", "bmp_010"), ("sqrt(2)/2", "amp_004", "bmp_011")]),
            "z_012": ("Q(2,A1,2,|1,-1)", [("-1", "amp_002", "bmp_012")]),
        },
    },
}
