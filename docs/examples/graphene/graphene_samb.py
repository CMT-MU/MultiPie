"""
=== SAMB (* only for crystal with fourier transform) ===
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
        "atomic": {"M_001": ["amp_001"]},
        "site_cluster": {"S_001": ["smp_001"]},
        "bond_cluster": {
            "B_001": ["bmp_002"],
            "B_002": ["bmp_003"],
            "B_003": ["bmp_004"],
            "B_004": ["bmp_005"],
            "B_005": ["bmp_006"],
            "B_006": ["bmp_007"],
        },
        "Z": {
            ("A1g", "M_001", "S_001"): ["z_001"],
            ("A1g", "M_001", "B_001"): ["z_002"],
            ("A1g", "M_001", "B_002"): ["z_003"],
            ("A1g", "M_001", "B_003"): ["z_004"],
            ("A1g", "M_001", "B_004"): ["z_005"],
            ("A1g", "M_001", "B_005"): ["z_006"],
            ("A1g", "M_001", "B_006"): ["z_007"],
        },
        "version": "1.1.15",
        "harmonics": {"Q": ["Qh(0,A1g,,)"], "G": []},
    },
    "data": {
        "atomic": {"amp_001": ("Qa(0,A1g,,)", (1, 1), [(0, 0, "1")])},
        "site_cluster": {"smp_001": ("Qs(0,A1g,,)", "[sqrt(2)/2, sqrt(2)/2]")},
        "bond_cluster": {
            "bmp_002": ("Qb(0,A1g,,)", "[sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]"),
            "bmp_003": ("Qb(0,A1g,,)", "[sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6]"),
            "bmp_004": ("Qb(0,A1g,,)", "[sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]"),
            "bmp_005": ("Qb(0,A1g,,)", "[sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6]"),
            "bmp_006": ("Qb(0,A1g,,)", "[sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6]"),
            "bmp_007": ("Qb(0,A1g,,)", "[sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6]"),
        },
        "Z": {
            "z_001": ("Q(0,A1g,,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(0,A1g,,)", [("1", "amp_001", "bmp_002")]),
            "z_003": ("Q(0,A1g,,)", [("1", "amp_001", "bmp_003")]),
            "z_004": ("Q(0,A1g,,)", [("1", "amp_001", "bmp_004")]),
            "z_005": ("Q(0,A1g,,)", [("1", "amp_001", "bmp_005")]),
            "z_006": ("Q(0,A1g,,)", [("1", "amp_001", "bmp_006")]),
            "z_007": ("Q(0,A1g,,)", [("1", "amp_001", "bmp_007")]),
        },
    },
}
