"""
example of molecule and crystal models.
"""


# ================================================== cluster models.
C3h = {
    "model": "C3h",
    "group": "C3h",
    #
    "site": {"H1": ("[0,0,0]", "s"), "O": ("[1/3,0,0]", "s p"), "H2": ("[1/2,1/6,0]", "s")},
    "bond": [("H1", "O", 1), ("O", "H2", 1)],
    #
    "spinful": True,
    #
    "option": {"view": [0, 0, 1]},
}

C3v = {
    "model": "C3v",
    "group": "C3v-1",
    "cell": {"c": 10},
    #
    "site": {"A": ("[-1/6,-1/6,0]", "s"), "B": ("[-2/3,0,0]", ["px", "py", "pz"])},
    "bond": [("A", "A", 1), ("A", "B", 1)],
    #
    "spinful": False,
    #
    "option": {"view": [0, 0, 1]},
}

CH4 = {
    "model": "CH4",
    "group": "Td",
    #
    "site": {"C": ("[0,0,0]", "s p"), "H": ("[1/3,1/3,1/3]", "s")},
    "bond": [("C", "H", 1)],
    #
    "spinful": True,
    "generate": {"irrep": ["A1", "A2"]},
}

D3 = {
    "model": "D3",
    "group": "D3-1",
    #
    "site": {"A": ("[1,0,1]", "s p")},
    "bond": [("A", "A", [1, 2])],
    #
    "spinful": True,
}

Th = {
    "model": "Th",
    "group": "Th",
    #
    "site": {"A": ("[1,1,1]", "s p")},
    "bond": [("A", "A", 1)],
    #
    "spinful": True,
    #
}

cluster_models = [C3h, C3v, CH4, D3, Th]


# ================================================== primitive models.
C2h1 = {
    "model": "C2h1",
    "group": "C2h^1",
    "cell": {"b": 1.2},
    #
    "site": {"A": ("[0,0,0]", "s"), "B": ("[1/2,1/2,0]", "s")},
    "bond": [("A", "A", 1), ("B", "B", 1), ("A", "B", 1)],
    #
    "spinful": True,
    #
    "option": {"view": [0, 0, 1]},
}

C3v1 = {
    "model": "C3v1",
    "group": "C3v^1",
    #
    "site": {"A": ("[1/3,2/3,0]", ["px", "py"]), "B": ("[2/3,1/3,0]", ["px", "py"])},
    "bond": [("A", "B", 1)],
    #
    "spinful": True,
    #
    "option": {"view": [0, 0, 1]},
}

C4v1 = {
    "model": "C4v1",
    "group": "C4v^1",
    #
    "site": {"A": ("[0,0,0]", "s")},
    "bond": [("A", "A", 1)],
    #
    "spinful": True,
}

CeCoSi = {
    "model": "CeCoSi",
    "group": 129,
    "cell": {"a": 4.057, "c": 6.987},
    #
    "site": {"Ce": ("[1/4,1/4,0.678]", "p"), "Co": ("[1/4,3/4,0]", "p"), "Si": ("[1/4,1/4,0.178]", "p")},
    "bond": [("Ce", "Co", [1]), ("Ce", "Si", [1]), ("Co", "Si", [1])],
    #
    "spinful": True,
}

Cs1 = {
    "model": "Cs1",
    "group": 6,
    #
    "site": {"A": ("[0,0,0]", ["px", "py"])},
    "bond": [("A", "A", [1, 2])],
    #
    "spinful": True,
    #
    "option": {"view": [0, 0, 1]},
}

D2h1 = {
    "model": "D2h1",
    "group": 47,
    "cell": {"a": 1.0, "b": 1.2, "c": 1.5},
    #
    "site": {"A": ("[0,0,0]", "s p")},
    "bond": [("A", "A", [1, 2, 3, 4])],
    #
    "spinful": True,
}

D4h1 = {
    "model": "D4h1",
    "group": 123,
    "cell": {"a": 1.0, "c": 1.5},
    #
    "site": {"A": ("[0,0,0]", "s p")},
    "bond": [("A", "A", [1, 2, 3, 4])],
    #
    "spinful": True,
}

grapheneAB = {
    "model": "grapheneAB",
    "group": 187,
    "cell": {"a": 2.435, "b": 2.435, "c": 10},
    #
    "site": {"A": ("[1/3,2/3,0]", "s"), "B": ("[2/3,1/3,0]", ["px", "py"])},
    "bond": [("A", "B", 1), ("A", "A", 1), ("B", "B", 1)],
    #
    "spinful": False,
    #
    "option": {"view": [0, 0, 1]},
}

graphene = {
    "model": "graphene",
    "group": 191,
    "cell": {"a": 2.435, "b": 2.435, "c": 10},
    #
    "site": {"A": ("[1/3,2/3,0]", "pz")},
    "bond": [("A", "A", [1, 2])],
    #
    "spinful": False,
    #
    "k_point": {"Γ": "[0, 0, 0]", "M": "[1/2, 0, 0]", "K": "[1/3, 1/3, 0]", "K'": "[-1/3, -1/3, 0]"},
    "k_path": "Γ-M-K-Γ-K'",
    #
    # "option": {"view": [0, 0, 1]},
    "option": {"view": [0, 0, 1]},
}


kappaET = {
    "model": "kappaET",
    "group": 32,
    "cell": {"b": 1.2},
    #
    "site": {"A": ("[9/10,1/20,0]", "s")},
    "bond": [("A", "A", [1, 2, 3])],
    #
    "spinful": True,
    #
    "option": {"view": [0, 0, 1]},
}

Mn3Sn = {
    "model": "Mn3Sn",
    "group": 194,
    #
    "site": {"Mn": ("[0.8388,0.6776,1/4]", "d"), "Sn": ("[1/3,2/3,1/4]", "d")},
    "bond": [("Mn", "Mn", 1), ("Sn", "Sn", 1), ("Mn", "Sn", 1)],
    #
    "spinful": False,
    #
    "option": {"view": [0, 0, 1]},
}

MoS2 = {
    "model": "MoS2",
    "group": 187,
    "cell": {"a": 3.1661, "b": 3.1661, "c": 20},
    #
    "site": {"Mo": ("[0,0,0]", ["du", "dv", "dxy", "dyz", "dzx"]), "S": ("[2/3, 1/3, 0.12425]", ["px", "py", "pz"])},
    "bond": [("Mo", "Mo", [1]), ("Mo", "S", [1]), ("S", "S", [1])],
    #
    "spinful": False,
    #
    "k_point": {"Γ": "[0, 0, 0]", "M": "[1/2, 0, 0]", "K": "[1/3, 1/3, 0]", "K'": "[-1/3, -1/3, 0]"},
    "k_path": "Γ-M-K-Γ-K'",
}

O1 = {
    "model": "O1",
    "group": 207,
    "cell": {"c": 1},
    #
    "site": {"A": ("[0,0,0]", "s p")},
    "bond": [("A", "A", [1, 2])],
    #
    "spinful": True,
}

Oh1 = {
    "model": "Oh1",
    "group": 221,
    "cell": {"c": 1},
    #
    "site": {"A": ("[0,0,0]", "s p")},
    "bond": [("A", "A", [1, 2])],
    #
    "spinful": True,
}

SnTe = {
    "model": "SnTe",
    "group": 31,
    "cell": {"a": 4.559, "b": 6, "c": 4.570},
    #
    "site": {"Sn": ("[1/2,6/10,6/10]", ["px", "py", "pz"]), "Te": ("[1/2,7/20,4/10]", ["px", "py", "pz"])},
    "bond": [("Sn", "Te", [1, 2]), ("Sn", "Sn", 1), ("Te", "Te", 1)],
    #
    "spinful": True,
    #
    "option": {"view": [0, 0, 1]},
}

SrVO3 = {
    "model": "SrVO3",
    "group": 221,
    "cell": {"a": 3.8409, "b": 3.8409, "c": 3.8409},
    #
    "site": {"V": ("[0.0, 0.0, 0.0]", ["dyz", "dzx", "dxy"])},
    "bond": [("V", "V", [1, 2])],
    #
    "spinful": False,
    #
    "k_point": {
        "R": "[1/2,1/2,1/2]",
        "Γ": "[0,0,0]",
        "X": "[1/2,0,0]",
        "M": "[1/2,1/2,0]",
    },
    "k_path": "R-Γ-X-R-M-Γ",
}

Te = {
    "model": "Te",
    "group": 152,
    "cell": {"a": 4.4580000, "b": 4.4580000, "c": 5.925},
    #
    "site": {"A": ("[0.274,0,1/3]", ["px", "py", "pz"])},
    "bond": [("A", "A", [1])],
    #
    "spinful": False,
    #
    "generate": {"model_type": "phonon"},
    #
    "k_path": "A-Γ-H-A-L-H-K-Γ-M-K",  # high-symmetry line
    "k_point": {  # k points
        "Γ": "[0, 0, 0]",
        "A": "[0, 0, 1/2]",
        "M": "[1/2, 0, 0]",
        "K": "[1/3, 1/3, 0]",
        "H": "[1/3, 1/3, 1/2]",
        "L": "[1/2, 0, 1/2]",
    },
}

UPt2Si2 = {
    "model": "UPt2Si2",
    "group": 129,
    "cell": {"a": 4.19720, "c": 9.69060},
    #
    "site": {
        "U": ("[1/4,1/4,0.7484]", "f"),
        "Pt1": ("[3/4,1/4,0]", "d"),
        "Pt2": ("[1/4,1/4,0.3785]", "d"),
        "Si1": ("[3/4,1/4,1/2]", "p"),
        "Si2": ("[1/4,1/4,0.1330]", "p"),
    },
    "bond": [("Pt1", "Si2", 1), ("Pt2", "Si1", 1)],
    #
    "spinful": False,
    #
    "option": {"view": [0, 0, 1]},
}


primitive_models = [
    C2h1,
    C3v1,
    C4v1,
    CeCoSi,
    Cs1,
    D2h1,
    D4h1,
    grapheneAB,
    graphene,
    kappaET,
    Mn3Sn,
    MoS2,
    O1,
    Oh1,
    SnTe,
    SrVO3,
    Te,
    UPt2Si2,
]

# ================================================== conventional models.
BCT = {
    "model": "BCT",
    "group": 139,
    "cell": {"c": 2.32},
    #
    "site": {"A": ("[0,0,0]", ["s", "px", "py"])},
    "bond": [("A", "A", [1, 2, 7])],
    #
    "spinful": True,
}

GaAs = {
    "model": "GaAs",
    "group": "Td^2",
    #
    "site": {"Ga": ("[0,0,0]", ["px", "py", "pz"]), "As": ("[1/4,1/4,1/4]", ["px", "py", "pz"])},
    "bond": [("Ga", "As", [1])],
    #
    "spinful": False,
}

C3v5 = {
    "model": "C3v5",
    "group": "C3v^5",
    #
    "site": {"A": ("[1/2,1/2,0]", "p"), "B": ("[1/2,1/3,0]", "p")},
    "bond": [("A", "B", 1)],
    #
    "spinful": True,
}

conventional_models = [BCT, GaAs, C3v5]


# ================================================== complex models.
kagome = {
    "model": "kagome",
    "group": 147,
    #
    "site": {"A": ("[1/2,0,0]", "s p")},
    "bond": [("A", "A", 1)],
    #
    "spinful": True,
    #
    "option": {"view": [0, 0, 1]},
}

Th1 = {
    "model": "Th1",
    "group": 200,
    "cell": {"c": 1},
    #
    "site": {"A": ("[1/2,0,0]", "s p")},
    "bond": [("A", "A", [1, 2])],
    #
    "spinful": True,
}

complex_models = [kagome, Th1]

# ================================================== all models.
models = cluster_models + primitive_models + conventional_models + complex_models
