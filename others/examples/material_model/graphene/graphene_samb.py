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
graphene = {
    "info": {
        "atomic": {"M_001": ["amp_001", "amp_002", "amp_003", "amp_004"]},
        "site_cluster": {"S_001": ["smp_001", "smp_002"]},
        "bond_cluster": {"B_001": ["bmp_003", "bmp_004", "bmp_005", "bmp_006", "bmp_007", "bmp_008"]},
        "Z": {
            ("A1g", "M_001", "S_001"): ["z_001"],
            ("A1g", "M_001", "B_001"): ["z_002"],
            ("E2g", "M_001", "B_001"): ["z_003", "z_004"],
            ("A1u", "M_001", "B_001"): ["z_005"],
            ("A2u", "M_001", "B_001"): ["z_006"],
            ("B1u", "M_001", "S_001"): ["z_007"],
            ("B2u", "M_001", "B_001"): ["z_008"],
            ("E1u", "M_001", "B_001"): ["z_009", "z_010"],
            ("E2u", "M_001", "B_001"): ["z_011", "z_012", "z_013", "z_014"],
        },
        "version": "1.1.15",
        "harmonics": {
            "Q": ["Qh(0,A1g,,)", "Qh(1,E1u,,0)", "Qh(1,E1u,,1)", "Qh(2,E2g,,0)", "Qh(2,E2g,,1)", "Qh(3,B1u,,)"],
            "G": ["Gh(1,A2g,,)", "Gh(1,E1g,,0)", "Gh(1,E1g,,1)"],
        },
    },
    "data": {
        "atomic": {
            "amp_001": ("Qa(0,A1g,,)", (2, 2), [(0, 0, "sqrt(2)/2"), (1, 1, "sqrt(2)/2")]),
            "amp_002": ("Ma(1,A2g,,|1,1)", (2, 2), [(0, 0, "sqrt(2)/2"), (1, 1, "-sqrt(2)/2")]),
            "amp_003": ("Ma(1,E1g,,0|1,1)", (2, 2), [(0, 1, "-sqrt(2)*I/2"), (1, 0, "sqrt(2)*I/2")]),
            "amp_004": ("Ma(1,E1g,,1|1,1)", (2, 2), [(0, 1, "-sqrt(2)/2"), (1, 0, "-sqrt(2)/2")]),
        },
        "site_cluster": {
            "smp_001": ("Qs(0,A1g,,)", "[sqrt(2)/2, sqrt(2)/2]"),
            "smp_002": ("Qs(3,B1u,,)", "[sqrt(2)/2, -sqrt(2)/2]"),
        },
        "bond_cluster": {
            "bmp_003": ("Qb(0,A1g,,)", "[sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]"),
            "bmp_004": ("Qb(2,E2g,,0)", "[sqrt(6)/3, -sqrt(6)/6, -sqrt(6)/6]"),
            "bmp_005": ("Qb(2,E2g,,1)", "[0, -sqrt(2)/2, sqrt(2)/2]"),
            "bmp_006": ("Tb(1,E1u,,0)", "[0, sqrt(2)*I/2, -sqrt(2)*I/2]"),
            "bmp_007": ("Tb(1,E1u,,1)", "[sqrt(6)*I/3, -sqrt(6)*I/6, -sqrt(6)*I/6]"),
            "bmp_008": ("Tb(3,B1u,,)", "[sqrt(3)*I/3, sqrt(3)*I/3, sqrt(3)*I/3]"),
        },
        "Z": {
            "z_001": ("Q(0,A1g,,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(0,A1g,,)", [("1", "amp_001", "bmp_003")]),
            "z_003": ("Q(2,E2g,,0)", [("1", "amp_001", "bmp_004")]),
            "z_004": ("Q(2,E2g,,1)", [("1", "amp_001", "bmp_005")]),
            "z_005": ("G(0,A1u,,|1,1)", [("-sqrt(2)/2", "amp_003", "bmp_007"), ("sqrt(2)/2", "amp_004", "bmp_006")]),
            "z_006": ("Q(1,A2u,,|1,1)", [("sqrt(2)/2", "amp_003", "bmp_006"), ("sqrt(2)/2", "amp_004", "bmp_007")]),
            "z_007": ("Q(3,B1u,,)", [("1", "amp_001", "smp_002")]),
            "z_008": ("Q(3,B2u,,|1,1)", [("-1", "amp_002", "bmp_008")]),
            "z_009": ("Q(1,E1u,,0|1,1)", [("-1", "amp_002", "bmp_007")]),
            "z_010": ("Q(1,E1u,,1|1,1)", [("1", "amp_002", "bmp_006")]),
            "z_011": ("G(2,E2u,,0|1,1)", [("-sqrt(2)/2", "amp_003", "bmp_006"), ("sqrt(2)/2", "amp_004", "bmp_007")]),
            "z_012": ("G(2,E2u,,1|1,1)", [("sqrt(2)/2", "amp_003", "bmp_007"), ("sqrt(2)/2", "amp_004", "bmp_006")]),
            "z_013": ("G(2,E2u,,0|1,1)", [("1", "amp_004", "bmp_008")]),
            "z_014": ("G(2,E2u,,1|1,1)", [("-1", "amp_003", "bmp_008")]),
        },
    },
}
