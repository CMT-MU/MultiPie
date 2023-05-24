"""
=== SAMB (* only for crystal) ===
- info
    - atomic : { "M_#" : ["amp_#"] }
    - site_cluster : { "S_#" : ["smp_#"] }
    - bond_cluster : { "B_#" : ["bmp_#"] }
    - uniform : { "S_#"/"B_#" : ["ump_#"] }
    - structure* : { "B_#" : ["kmp_#"] }
    - Z : { ("M_#", "S_#"/"B_#") : ["z_#"] }
    - version : MultiPie version
    - harmonics : { head : { "harm_tag" } }

- data
    - atomic : { "amp_#" : ( TagMultipole, shape, [(i, j, matrix element)] ) }
    - site_cluster : { "smp_#" : ( TagMultipole, [vector component] ) }
    - bond_cluster : { "bmp_#" : ( TagMultipole, [vector component] ) }
    - uniform : { "ump_#" : ( TagMultipole, shape, [(i, j, matrix element)] ) }
    - structure* : { "kmp_#" : (TagMultipole, "formfactor") }
    - Z : {"z_#" : ( TagMultipole, [(coeff, "amp_#", "smp_#"/"bmp_#/ump_#")] ) }
    - Zk* : {"z_#" : ( TagMultipole, [(coeff, "amp_#", "ump_#", "kmp_#")] ) }
"""
graphene = {
    "info": {
        "atomic": {"M_001": ["amp_001"]},
        "site_cluster": {"S_001": ["smp_001"]},
        "bond_cluster": {"B_001": ["bmp_002"], "B_002": ["bmp_003"]},
        "uniform": {"S_001": ["ump_001"], "B_001": ["ump_002", "ump_003"]},
        "Z": {("A1g", "M_001", "S_001"): ["z_001"], ("A1g", "M_001", "B_001"): ["z_002"], ("A1g", "M_001", "B_002"): ["z_003"]},
        "version": "1.1.1",
        "structure": {"B_001": ["kmp_001", "kmp_002"], "B_002": ["kmp_003"]},
        "harmonics": {"Q": ["Qh(0,A1g,,)", "Qh(3,B1u,,)"], "G": []},
    },
    "data": {
        "atomic": {"amp_001": ("Qa(0,A1g,,)", (1, 1), [(0, 0, "1")])},
        "site_cluster": {"smp_001": ("Qs(0,A1g,,)", "[sqrt(2)/2, sqrt(2)/2]")},
        "bond_cluster": {
            "bmp_002": ("Qb(0,A1g,,)", "[sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]"),
            "bmp_003": ("Qb(0,A1g,,)", "[sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6, sqrt(6)/6]"),
        },
        "uniform": {
            "ump_001": ("Qs(0,A1g,,)", (2, 2), [(0, 0, "sqrt(2)/2"), (1, 1, "sqrt(2)/2")]),
            "ump_002": ("Qu(0,A1g,,)", (2, 2), [(0, 1, "sqrt(2)/2"), (1, 0, "sqrt(2)/2")]),
            "ump_003": ("Tu(3,B1u,,)", (2, 2), [(0, 1, "sqrt(2)*I/2"), (1, 0, "-sqrt(2)*I/2")]),
        },
        "Z": {
            "z_001": ("Q(0,A1g,,)", [("1", "amp_001", "smp_001")]),
            "z_002": ("Q(0,A1g,,)", [("1", "amp_001", "bmp_002")]),
            "z_003": ("Q(0,A1g,,)", [("1", "amp_001", "bmp_003")]),
        },
        "structure": {
            "kmp_001": ("Qk(0,A1g,,)", "sqrt(6)*c001/3 + sqrt(6)*c002/3 + sqrt(6)*c003/3"),
            "kmp_002": ("Tk(3,B1u,,)", "sqrt(6)*s001/3 + sqrt(6)*s002/3 + sqrt(6)*s003/3"),
            "kmp_003": ("Qk(0,A1g,,)", "sqrt(6)*c004/3 + sqrt(6)*c006/3 + sqrt(6)*c007/3"),
        },
        "Zk": {
            "z_001": ("Q(0,A1g,,)", [("1", "amp_001", "ump_001")]),
            "z_002": (
                "Q(0,A1g,,)",
                [("sqrt(2)/2", "amp_001", "ump_002", "kmp_001"), ("-sqrt(2)/2", "amp_001", "ump_003", "kmp_002")],
            ),
            "z_003": ("Q(0,A1g,,)", [("1", "amp_001", "ump_001", "kmp_003")]),
        },
    },
}
