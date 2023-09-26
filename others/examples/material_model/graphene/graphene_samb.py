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
        "atomic": {"M_001": ["amp_001"]},
        "site_cluster": {"S_001": ["smp_001"]},
        "bond_cluster": {"B_001": ["bmp_002"]},
        "Z": {("A1g", "M_001", "S_001"): ["z_001"], ("A1g", "M_001", "B_001"): ["z_002"]},
        "version": "1.1.14",
        "harmonics": {"Q": ["Qh(0,A1g,,)"], "G": []},
    },
    "data": {
        "atomic": {"amp_001": ("Qa(0,A1g,,)", (2, 2), [(0, 0, "sqrt(2)/2"), (1, 1, "sqrt(2)/2")])},
        "site_cluster": {"smp_001": ("Qs(0,A1g,,)", "[sqrt(2)/2, sqrt(2)/2]")},
        "bond_cluster": {"bmp_002": ("Qb(0,A1g,,)", "[sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]")},
        "Z": {"z_001": ("Q(0,A1g,,)", [("1", "amp_001", "smp_001")]), "z_002": ("Q(0,A1g,,)", [("1", "amp_001", "bmp_002")])},
    },
}
