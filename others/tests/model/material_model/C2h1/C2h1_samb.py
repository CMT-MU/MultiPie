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
C2h1 = {
    "info": {
        "atomic": {"M_001": ["amp_001", "amp_002", "amp_003", "amp_004"]},
        "site_cluster": {"S_001": ["smp_001"], "S_002": ["smp_002"]},
        "bond_cluster": {
            "B_001": ["bmp_003"],
            "B_002": ["bmp_004"],
            "B_003": ["bmp_005"],
            "B_004": ["bmp_006"],
            "B_005": ["bmp_007", "bmp_008", "bmp_009"],
        },
        "Z": {
            ("Ag", "M_001", "S_001"): ["z_001"],
            ("Ag", "M_001", "S_002"): ["z_002"],
            ("Ag", "M_001", "B_001"): ["z_003"],
            ("Ag", "M_001", "B_002"): ["z_004"],
            ("Ag", "M_001", "B_003"): ["z_005"],
            ("Ag", "M_001", "B_004"): ["z_006"],
            ("Ag", "M_001", "B_005"): ["z_007", "z_008", "z_009", "z_010"],
        },
        "version": "1.1.14",
        "harmonics": {"Q": ["Qh(0,Ag,,)", "Qh(2,Bg,2,)"], "G": ["Gh(1,Ag,,)", "Gh(1,Bg,1,)", "Gh(1,Bg,2,)"]},
    },
    "data": {
        "atomic": {
            "amp_001": ("Qa(0,Ag,,)", (2, 2), [(0, 0, "sqrt(2)/2"), (1, 1, "sqrt(2)/2")]),
            "amp_002": ("Ma(1,Ag,,|1,-1)", (2, 2), [(0, 1, "-sqrt(2)*I/2"), (1, 0, "sqrt(2)*I/2")]),
            "amp_003": ("Ma(1,Bg,1,|1,-1)", (2, 2), [(0, 1, "sqrt(2)/2"), (1, 0, "sqrt(2)/2")]),
            "amp_004": ("Ma(1,Bg,2,|1,-1)", (2, 2), [(0, 0, "sqrt(2)/2"), (1, 1, "-sqrt(2)/2")]),
        },
        "site_cluster": {"smp_001": ("Qs(0,Ag,,)", "[1]"), "smp_002": ("Qs(0,Ag,,)", "[1]")},
        "bond_cluster": {
            "bmp_003": ("Qb(0,Ag,,)", "[1]"),
            "bmp_004": ("Qb(0,Ag,,)", "[1]"),
            "bmp_005": ("Qb(0,Ag,,)", "[1]"),
            "bmp_006": ("Qb(0,Ag,,)", "[1]"),
            "bmp_007": ("Qb(0,Ag,,)", "[1/2, 1/2, 1/2, 1/2]"),
            "bmp_008": ("Tb(0,Ag,,)", "[I/2, I/2, I/2, I/2]"),
            "bmp_009": ("Tb(2,Bg,2,)", "[I/2, -I/2, I/2, -I/2]"),
        },
        "Z": {
            "z_001": ("Q(0,Ag,,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(0,Ag,,)", [("1", "amp_001", "smp_002")]),
            "z_003": ("Q(0,Ag,,)", [("1", "amp_001", "bmp_003")]),
            "z_004": ("Q(0,Ag,,)", [("1", "amp_001", "bmp_004")]),
            "z_005": ("Q(0,Ag,,)", [("1", "amp_001", "bmp_005")]),
            "z_006": ("Q(0,Ag,,)", [("1", "amp_001", "bmp_006")]),
            "z_007": ("Q(0,Ag,,)", [("1", "amp_001", "bmp_007")]),
            "z_008": ("G(1,Ag,,|1,-1)", [("1", "amp_002", "bmp_008")]),
            "z_009": ("G(1,Ag,,|1,-1)", [("1", "amp_003", "bmp_009")]),
            "z_010": ("Q(2,Ag,2,|1,-1)", [("-1", "amp_004", "bmp_009")]),
        },
    },
}
