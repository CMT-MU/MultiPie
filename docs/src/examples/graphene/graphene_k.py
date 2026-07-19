"""
Selected SAMB matrix in momentum representation.
- dimension (int): matrix size.
- ket_site (list): ket info., [ket_name].
- index (dict): ket index, dict[(site,sublattice,rank), (top_index,size)].
- cluster_vector (dict): cluster vector, dict[site/bond name, dict[kb, expression] ].
- k_multipole (dict): momentum multipole in terms of k.b_n, dict[wyckoff, dict[idx, (k_multipole, symmetry)] ].
- k_matrix (dict): momentum matrix, dict[tag, (site/bond_name, wyckoff, dict[(m,n), value]) ].
"""

graphene_k = {
    "dimension": 2,
    "ket_site": ["pz@C(1)", "pz@C(2)"],
    "index": {("C", 1, 1): (0, 1), ("C", 2, 1): (1, 1)},
    "cluster_vector": {"C": {}, "C;C_001_1": {"kb_1": "0.33333333*k_1+0.66666666*k_2", "kb_2": "-0.66666666*k_1-0.33333333*k_2", "kb_3": "0.33333333*k_1-0.33333333*k_2"}, "C;C_002_1": {"kb_1": "1.0*k_1", "kb_2": "1.0*k_2", "kb_3": "-1.0*k_1-1.0*k_2"}},
    "k_multipole": {
        "2c": {("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[[sqrt(2)/2,sqrt(2)/2]]", "[1]"), ("Q", 3, "B1u", -1, -1, 0, 0, "q"): ("[[sqrt(2)/2,-sqrt(2)/2]]", "[sqrt(10)*y*(3*x**2-y**2)/4]")},
        "3a@3f": {
            ("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[sqrt(6)*cos(kb_1)/3+sqrt(6)*cos(kb_2)/3+sqrt(6)*cos(kb_3)/3]", "[1]"),
            ("T", 3, "B1u", -1, -1, 0, 0, "q"): ("[-sqrt(6)*sin(kb_1)/3-sqrt(6)*sin(kb_2)/3-sqrt(6)*sin(kb_3)/3]", "[sqrt(10)*y*(3*x**2-y**2)/4]"),
            ("T", 1, "E1u", -1, -1, 0, 0, "q"): ("[sin(kb_2)-sin(kb_3),-2*sqrt(3)*sin(kb_1)/3+sqrt(3)*sin(kb_2)/3+sqrt(3)*sin(kb_3)/3]", "[x,y]"),
            ("Q", 2, "E2g", -1, -1, 0, 0, "q"): ("[2*sqrt(3)*cos(kb_1)/3-sqrt(3)*cos(kb_2)/3-sqrt(3)*cos(kb_3)/3,cos(kb_2)-cos(kb_3)]", "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]"),
        },
        "6b@6l": {
            ("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[sqrt(3)*cos(kb_1)/3+sqrt(3)*cos(kb_2)/3+sqrt(3)*cos(kb_3)/3+sqrt(3)*cos(kb_4)/3+sqrt(3)*cos(kb_5)/3+sqrt(3)*cos(kb_6)/3]", "[1]"),
            ("M", 1, "A2g", -1, -1, 0, 0, "q"): ("[-sqrt(3)*sin(kb_1)/3-sqrt(3)*sin(kb_2)/3-sqrt(3)*sin(kb_3)/3-sqrt(3)*sin(kb_4)/3-sqrt(3)*sin(kb_5)/3-sqrt(3)*sin(kb_6)/3]", "[z]"),
            ("Q", 3, "B1u", -1, -1, 0, 0, "q"): ("[sqrt(3)*cos(kb_1)/3+sqrt(3)*cos(kb_2)/3+sqrt(3)*cos(kb_3)/3-sqrt(3)*cos(kb_4)/3-sqrt(3)*cos(kb_5)/3-sqrt(3)*cos(kb_6)/3]", "[sqrt(10)*y*(3*x**2-y**2)/4]"),
            ("T", 3, "B2u", -1, -1, 0, 0, "q"): ("[-sqrt(3)*sin(kb_1)/3-sqrt(3)*sin(kb_2)/3-sqrt(3)*sin(kb_3)/3+sqrt(3)*sin(kb_4)/3+sqrt(3)*sin(kb_5)/3+sqrt(3)*sin(kb_6)/3]", "[sqrt(10)*x*(x**2-3*y**2)/4]"),
            ("Q", 1, "E1u", -1, -1, 0, 0, "q"): ("[-sqrt(2)*cos(kb_2)/2+sqrt(2)*cos(kb_3)/2+sqrt(2)*cos(kb_5)/2-sqrt(2)*cos(kb_6)/2,sqrt(6)*cos(kb_1)/3-sqrt(6)*cos(kb_2)/6-sqrt(6)*cos(kb_3)/6-sqrt(6)*cos(kb_4)/3+sqrt(6)*cos(kb_5)/6+sqrt(6)*cos(kb_6)/6]", "[x,y]"),
            ("T", 1, "E1u", -1, -1, 0, 0, "q"): ("[-sqrt(6)*sin(kb_1)/3+sqrt(6)*sin(kb_2)/6+sqrt(6)*sin(kb_3)/6+sqrt(6)*sin(kb_4)/3-sqrt(6)*sin(kb_5)/6-sqrt(6)*sin(kb_6)/6,-sqrt(2)*sin(kb_2)/2+sqrt(2)*sin(kb_3)/2+sqrt(2)*sin(kb_5)/2-sqrt(2)*sin(kb_6)/2]", "[x,y]"),
            ("Q", 2, "E2g", -1, -1, 0, 0, "q"): ("[sqrt(6)*cos(kb_1)/3-sqrt(6)*cos(kb_2)/6-sqrt(6)*cos(kb_3)/6+sqrt(6)*cos(kb_4)/3-sqrt(6)*cos(kb_5)/6-sqrt(6)*cos(kb_6)/6,sqrt(2)*cos(kb_2)/2-sqrt(2)*cos(kb_3)/2+sqrt(2)*cos(kb_5)/2-sqrt(2)*cos(kb_6)/2]", "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]"),
            ("T", 2, "E2g", -1, -1, 0, 0, "q"): ("[-sqrt(2)*sin(kb_2)/2+sqrt(2)*sin(kb_3)/2-sqrt(2)*sin(kb_5)/2+sqrt(2)*sin(kb_6)/2,sqrt(6)*sin(kb_1)/3-sqrt(6)*sin(kb_2)/6-sqrt(6)*sin(kb_3)/6+sqrt(6)*sin(kb_4)/3-sqrt(6)*sin(kb_5)/6-sqrt(6)*sin(kb_6)/6]", "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]"),
        },
    },
    "k_matrix": {
        "z1": ("C", "2c", {(0, 0): "sqrt(2)/2", (1, 1): "sqrt(2)/2"}),
        "z2": ("C;C_001_1", "3a@3f", {(1, 0): "sqrt(6)*exp(I*kb_1)/6+sqrt(6)*exp(I*kb_2)/6+sqrt(6)*exp(I*kb_3)/6", (0, 1): "sqrt(6)*exp(-I*kb_3)/6+sqrt(6)*exp(-I*kb_2)/6+sqrt(6)*exp(-I*kb_1)/6"}),
        "z3": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "sqrt(3)*exp(I*kb_1)/6+sqrt(3)*exp(I*kb_2)/6+sqrt(3)*exp(I*kb_3)/6+sqrt(3)*exp(-I*kb_3)/6+sqrt(3)*exp(-I*kb_2)/6+sqrt(3)*exp(-I*kb_1)/6",
                (1, 1): "sqrt(3)*exp(I*kb_4)/6+sqrt(3)*exp(I*kb_5)/6+sqrt(3)*exp(I*kb_6)/6+sqrt(3)*exp(-I*kb_6)/6+sqrt(3)*exp(-I*kb_5)/6+sqrt(3)*exp(-I*kb_4)/6",
            },
        ),
        "z4": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "sqrt(3)*I*exp(I*kb_1)/6+sqrt(3)*I*exp(I*kb_2)/6+sqrt(3)*I*exp(I*kb_3)/6-sqrt(3)*I*exp(-I*kb_3)/6-sqrt(3)*I*exp(-I*kb_2)/6-sqrt(3)*I*exp(-I*kb_1)/6",
                (1, 1): "sqrt(3)*I*exp(I*kb_4)/6+sqrt(3)*I*exp(I*kb_5)/6+sqrt(3)*I*exp(I*kb_6)/6-sqrt(3)*I*exp(-I*kb_6)/6-sqrt(3)*I*exp(-I*kb_5)/6-sqrt(3)*I*exp(-I*kb_4)/6",
            },
        ),
        "z5": ("C", "2c", {(0, 0): "sqrt(2)/2", (1, 1): "-sqrt(2)/2"}),
        "z6": ("C;C_001_1", "3a@3f", {(1, 0): "sqrt(6)*I*exp(I*kb_1)/6+sqrt(6)*I*exp(I*kb_2)/6+sqrt(6)*I*exp(I*kb_3)/6", (0, 1): "-sqrt(6)*I*exp(-I*kb_3)/6-sqrt(6)*I*exp(-I*kb_2)/6-sqrt(6)*I*exp(-I*kb_1)/6"}),
        "z7": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "sqrt(3)*exp(I*kb_1)/6+sqrt(3)*exp(I*kb_2)/6+sqrt(3)*exp(I*kb_3)/6+sqrt(3)*exp(-I*kb_3)/6+sqrt(3)*exp(-I*kb_2)/6+sqrt(3)*exp(-I*kb_1)/6",
                (1, 1): "-sqrt(3)*exp(I*kb_4)/6-sqrt(3)*exp(I*kb_5)/6-sqrt(3)*exp(I*kb_6)/6-sqrt(3)*exp(-I*kb_6)/6-sqrt(3)*exp(-I*kb_5)/6-sqrt(3)*exp(-I*kb_4)/6",
            },
        ),
        "z8": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "sqrt(3)*I*exp(I*kb_1)/6+sqrt(3)*I*exp(I*kb_2)/6+sqrt(3)*I*exp(I*kb_3)/6-sqrt(3)*I*exp(-I*kb_3)/6-sqrt(3)*I*exp(-I*kb_2)/6-sqrt(3)*I*exp(-I*kb_1)/6",
                (1, 1): "-sqrt(3)*I*exp(I*kb_4)/6-sqrt(3)*I*exp(I*kb_5)/6-sqrt(3)*I*exp(I*kb_6)/6+sqrt(3)*I*exp(-I*kb_6)/6+sqrt(3)*I*exp(-I*kb_5)/6+sqrt(3)*I*exp(-I*kb_4)/6",
            },
        ),
        "z9": ("C;C_001_1", "3a@3f", {(1, 0): "-I*exp(I*kb_2)/2+I*exp(I*kb_3)/2", (0, 1): "-I*exp(-I*kb_3)/2+I*exp(-I*kb_2)/2"}),
        "z10": ("C;C_001_1", "3a@3f", {(1, 0): "sqrt(3)*I*exp(I*kb_1)/3-sqrt(3)*I*exp(I*kb_2)/6-sqrt(3)*I*exp(I*kb_3)/6", (0, 1): "sqrt(3)*I*exp(-I*kb_3)/6+sqrt(3)*I*exp(-I*kb_2)/6-sqrt(3)*I*exp(-I*kb_1)/3"}),
        "z11": ("C;C_002_1", "6b@6l", {(0, 0): "-sqrt(2)*exp(I*kb_2)/4+sqrt(2)*exp(I*kb_3)/4+sqrt(2)*exp(-I*kb_3)/4-sqrt(2)*exp(-I*kb_2)/4", (1, 1): "sqrt(2)*exp(I*kb_5)/4-sqrt(2)*exp(I*kb_6)/4-sqrt(2)*exp(-I*kb_6)/4+sqrt(2)*exp(-I*kb_5)/4"}),
        "z12": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "sqrt(6)*exp(I*kb_1)/6-sqrt(6)*exp(I*kb_2)/12-sqrt(6)*exp(I*kb_3)/12-sqrt(6)*exp(-I*kb_3)/12-sqrt(6)*exp(-I*kb_2)/12+sqrt(6)*exp(-I*kb_1)/6",
                (1, 1): "-sqrt(6)*exp(I*kb_4)/6+sqrt(6)*exp(I*kb_5)/12+sqrt(6)*exp(I*kb_6)/12+sqrt(6)*exp(-I*kb_6)/12+sqrt(6)*exp(-I*kb_5)/12-sqrt(6)*exp(-I*kb_4)/6",
            },
        ),
        "z13": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "sqrt(6)*I*exp(I*kb_1)/6-sqrt(6)*I*exp(I*kb_2)/12-sqrt(6)*I*exp(I*kb_3)/12+sqrt(6)*I*exp(-I*kb_3)/12+sqrt(6)*I*exp(-I*kb_2)/12-sqrt(6)*I*exp(-I*kb_1)/6",
                (1, 1): "-sqrt(6)*I*exp(I*kb_4)/6+sqrt(6)*I*exp(I*kb_5)/12+sqrt(6)*I*exp(I*kb_6)/12-sqrt(6)*I*exp(-I*kb_6)/12-sqrt(6)*I*exp(-I*kb_5)/12+sqrt(6)*I*exp(-I*kb_4)/6",
            },
        ),
        "z14": ("C;C_002_1", "6b@6l", {(0, 0): "sqrt(2)*I*exp(I*kb_2)/4-sqrt(2)*I*exp(I*kb_3)/4+sqrt(2)*I*exp(-I*kb_3)/4-sqrt(2)*I*exp(-I*kb_2)/4", (1, 1): "-sqrt(2)*I*exp(I*kb_5)/4+sqrt(2)*I*exp(I*kb_6)/4-sqrt(2)*I*exp(-I*kb_6)/4+sqrt(2)*I*exp(-I*kb_5)/4"}),
        "z15": ("C;C_001_1", "3a@3f", {(1, 0): "sqrt(3)*exp(I*kb_1)/3-sqrt(3)*exp(I*kb_2)/6-sqrt(3)*exp(I*kb_3)/6", (0, 1): "-sqrt(3)*exp(-I*kb_3)/6-sqrt(3)*exp(-I*kb_2)/6+sqrt(3)*exp(-I*kb_1)/3"}),
        "z16": ("C;C_001_1", "3a@3f", {(1, 0): "exp(I*kb_2)/2-exp(I*kb_3)/2", (0, 1): "-exp(-I*kb_3)/2+exp(-I*kb_2)/2"}),
        "z17": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "sqrt(6)*exp(I*kb_1)/6-sqrt(6)*exp(I*kb_2)/12-sqrt(6)*exp(I*kb_3)/12-sqrt(6)*exp(-I*kb_3)/12-sqrt(6)*exp(-I*kb_2)/12+sqrt(6)*exp(-I*kb_1)/6",
                (1, 1): "sqrt(6)*exp(I*kb_4)/6-sqrt(6)*exp(I*kb_5)/12-sqrt(6)*exp(I*kb_6)/12-sqrt(6)*exp(-I*kb_6)/12-sqrt(6)*exp(-I*kb_5)/12+sqrt(6)*exp(-I*kb_4)/6",
            },
        ),
        "z18": ("C;C_002_1", "6b@6l", {(0, 0): "sqrt(2)*exp(I*kb_2)/4-sqrt(2)*exp(I*kb_3)/4-sqrt(2)*exp(-I*kb_3)/4+sqrt(2)*exp(-I*kb_2)/4", (1, 1): "sqrt(2)*exp(I*kb_5)/4-sqrt(2)*exp(I*kb_6)/4-sqrt(2)*exp(-I*kb_6)/4+sqrt(2)*exp(-I*kb_5)/4"}),
        "z19": ("C;C_002_1", "6b@6l", {(0, 0): "sqrt(2)*I*exp(I*kb_2)/4-sqrt(2)*I*exp(I*kb_3)/4+sqrt(2)*I*exp(-I*kb_3)/4-sqrt(2)*I*exp(-I*kb_2)/4", (1, 1): "sqrt(2)*I*exp(I*kb_5)/4-sqrt(2)*I*exp(I*kb_6)/4+sqrt(2)*I*exp(-I*kb_6)/4-sqrt(2)*I*exp(-I*kb_5)/4"}),
        "z20": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "-sqrt(6)*I*exp(I*kb_1)/6+sqrt(6)*I*exp(I*kb_2)/12+sqrt(6)*I*exp(I*kb_3)/12-sqrt(6)*I*exp(-I*kb_3)/12-sqrt(6)*I*exp(-I*kb_2)/12+sqrt(6)*I*exp(-I*kb_1)/6",
                (1, 1): "-sqrt(6)*I*exp(I*kb_4)/6+sqrt(6)*I*exp(I*kb_5)/12+sqrt(6)*I*exp(I*kb_6)/12-sqrt(6)*I*exp(-I*kb_6)/12-sqrt(6)*I*exp(-I*kb_5)/12+sqrt(6)*I*exp(-I*kb_4)/6",
            },
        ),
    },
}
