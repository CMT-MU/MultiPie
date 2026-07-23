"""
Selected SAMB matrix in momentum representation.
- dimension (int): matrix size.
- ket_site (list): ket info., [ket_name].
- index (dict): ket index, dict[(site,sublattice,rank), (top_index,size)].
- cluster_vector (dict): cluster vector, dict[site/bond name, dict[kb, expression] ].
- k_multipole (dict): momentum multipole in terms of p_n=k.b_n, dict[wyckoff, dict[idx, (k_multipole, symmetry)] ].
- k_matrix (dict): momentum matrix, dict[tag, (site/bond_name, wyckoff, dict[(m,n), value]) ].
"""

graphene_k = {
    "dimension": 2,
    "ket_site": ["pz@C(1)", "pz@C(2)"],
    "index": {("C", 1, 1): (0, 1), ("C", 2, 1): (1, 1)},
    "cluster_vector": {"C": {}, "C;C_001_1": {"p_1": "0.33333333*k_1+0.66666666*k_2", "p_2": "-0.66666666*k_1-0.33333333*k_2", "p_3": "0.33333333*k_1-0.33333333*k_2"}, "C;C_002_1": {"p_1": "1.0*k_1", "p_2": "1.0*k_2", "p_3": "-1.0*k_1-1.0*k_2"}},
    "k_multipole": {
        "2c": {("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[[sqrt(2)/2,sqrt(2)/2]]", "[1]"), ("Q", 3, "B1u", -1, -1, 0, 0, "q"): ("[[sqrt(2)/2,-sqrt(2)/2]]", "[sqrt(10)*y*(3*x**2-y**2)/4]")},
        "3a@3f": {
            ("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[0.816496580927726*cos(p_1)+0.816496580927726*cos(p_2)+0.816496580927726*cos(p_3)]", "[1]"),
            ("T", 3, "B1u", -1, -1, 0, 0, "q"): ("[0.816496580927726*I*sin(p_1)+0.816496580927726*I*sin(p_2)+0.816496580927726*I*sin(p_3)]", "[sqrt(10)*y*(3*x**2-y**2)/4]"),
            ("T", 1, "E1u", -1, -1, 0, 0, "q"): ("[-1.0*I*sin(p_2)+1.0*I*sin(p_3),1.15470053837925*I*sin(p_1)-0.577350269189626*I*sin(p_2)-0.577350269189626*I*sin(p_3)]", "[x,y]"),
            ("Q", 2, "E2g", -1, -1, 0, 0, "q"): ("[1.15470053837925*cos(p_1)-0.577350269189626*cos(p_2)-0.577350269189626*cos(p_3),1.0*cos(p_2)-1.0*cos(p_3)]", "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]"),
        },
        "6b@6l": {
            ("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[0.577350269189626*cos(p_1)+0.577350269189626*cos(p_2)+0.577350269189626*cos(p_3)+0.577350269189626*cos(p_4)+0.577350269189626*cos(p_5)+0.577350269189626*cos(p_6)]", "[1]"),
            ("M", 1, "A2g", -1, -1, 0, 0, "q"): ("[0.577350269189626*I*sin(p_1)+0.577350269189626*I*sin(p_2)+0.577350269189626*I*sin(p_3)+0.577350269189626*I*sin(p_4)+0.577350269189626*I*sin(p_5)+0.577350269189626*I*sin(p_6)]", "[z]"),
            ("Q", 3, "B1u", -1, -1, 0, 0, "q"): ("[0.577350269189626*cos(p_1)+0.577350269189626*cos(p_2)+0.577350269189626*cos(p_3)-0.577350269189626*cos(p_4)-0.577350269189626*cos(p_5)-0.577350269189626*cos(p_6)]", "[sqrt(10)*y*(3*x**2-y**2)/4]"),
            ("T", 3, "B2u", -1, -1, 0, 0, "q"): ("[0.577350269189626*I*sin(p_1)+0.577350269189626*I*sin(p_2)+0.577350269189626*I*sin(p_3)-0.577350269189626*I*sin(p_4)-0.577350269189626*I*sin(p_5)-0.577350269189626*I*sin(p_6)]", "[sqrt(10)*x*(x**2-3*y**2)/4]"),
            ("Q", 1, "E1u", -1, -1, 0, 0, "q"): (
                "[-0.707106781186548*cos(p_2)+0.707106781186548*cos(p_3)+0.707106781186548*cos(p_5)-0.707106781186548*cos(p_6),0.816496580927726*cos(p_1)-0.408248290463863*cos(p_2)-0.408248290463863*cos(p_3)-0.816496580927726*cos(p_4)+0.408248290463863*cos(p_5)+0.408248290463863*cos(p_6)]",
                "[x,y]",
            ),
            ("T", 1, "E1u", -1, -1, 0, 0, "q"): (
                "[0.816496580927726*I*sin(p_1)-0.408248290463863*I*sin(p_2)-0.408248290463863*I*sin(p_3)-0.816496580927726*I*sin(p_4)+0.408248290463863*I*sin(p_5)+0.408248290463863*I*sin(p_6),0.707106781186548*I*sin(p_2)-0.707106781186548*I*sin(p_3)-0.707106781186548*I*sin(p_5)+0.707106781186548*I*sin(p_6)]",
                "[x,y]",
            ),
            ("Q", 2, "E2g", -1, -1, 0, 0, "q"): (
                "[0.816496580927726*cos(p_1)-0.408248290463863*cos(p_2)-0.408248290463863*cos(p_3)+0.816496580927726*cos(p_4)-0.408248290463863*cos(p_5)-0.408248290463863*cos(p_6),0.707106781186548*cos(p_2)-0.707106781186548*cos(p_3)+0.707106781186548*cos(p_5)-0.707106781186548*cos(p_6)]",
                "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]",
            ),
            ("T", 2, "E2g", -1, -1, 0, 0, "q"): (
                "[0.707106781186548*I*sin(p_2)-0.707106781186548*I*sin(p_3)+0.707106781186548*I*sin(p_5)-0.707106781186548*I*sin(p_6),-0.816496580927726*I*sin(p_1)+0.408248290463863*I*sin(p_2)+0.408248290463863*I*sin(p_3)-0.816496580927726*I*sin(p_4)+0.408248290463863*I*sin(p_5)+0.408248290463863*I*sin(p_6)]",
                "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]",
            ),
        },
    },
    "k_matrix": {
        "z1": ("C", "2c", {(0, 0): "sqrt(2)/2", (1, 1): "sqrt(2)/2"}),
        "z2": ("C;C_001_1", "3a@3f", {(1, 0): "sqrt(6)*exp(I*p_1)/6+sqrt(6)*exp(I*p_2)/6+sqrt(6)*exp(I*p_3)/6", (0, 1): "sqrt(6)*exp(-I*p_3)/6+sqrt(6)*exp(-I*p_2)/6+sqrt(6)*exp(-I*p_1)/6"}),
        "z3": (
            "C;C_002_1",
            "6b@6l",
            {(0, 0): "sqrt(3)*exp(I*p_1)/6+sqrt(3)*exp(I*p_2)/6+sqrt(3)*exp(I*p_3)/6+sqrt(3)*exp(-I*p_3)/6+sqrt(3)*exp(-I*p_2)/6+sqrt(3)*exp(-I*p_1)/6", (1, 1): "sqrt(3)*exp(I*p_4)/6+sqrt(3)*exp(I*p_5)/6+sqrt(3)*exp(I*p_6)/6+sqrt(3)*exp(-I*p_6)/6+sqrt(3)*exp(-I*p_5)/6+sqrt(3)*exp(-I*p_4)/6"},
        ),
        "z4": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "sqrt(3)*I*exp(I*p_1)/6+sqrt(3)*I*exp(I*p_2)/6+sqrt(3)*I*exp(I*p_3)/6-sqrt(3)*I*exp(-I*p_3)/6-sqrt(3)*I*exp(-I*p_2)/6-sqrt(3)*I*exp(-I*p_1)/6",
                (1, 1): "sqrt(3)*I*exp(I*p_4)/6+sqrt(3)*I*exp(I*p_5)/6+sqrt(3)*I*exp(I*p_6)/6-sqrt(3)*I*exp(-I*p_6)/6-sqrt(3)*I*exp(-I*p_5)/6-sqrt(3)*I*exp(-I*p_4)/6",
            },
        ),
        "z5": ("C", "2c", {(0, 0): "sqrt(2)/2", (1, 1): "-sqrt(2)/2"}),
        "z6": ("C;C_001_1", "3a@3f", {(1, 0): "sqrt(6)*I*exp(I*p_1)/6+sqrt(6)*I*exp(I*p_2)/6+sqrt(6)*I*exp(I*p_3)/6", (0, 1): "-sqrt(6)*I*exp(-I*p_3)/6-sqrt(6)*I*exp(-I*p_2)/6-sqrt(6)*I*exp(-I*p_1)/6"}),
        "z7": (
            "C;C_002_1",
            "6b@6l",
            {(0, 0): "sqrt(3)*exp(I*p_1)/6+sqrt(3)*exp(I*p_2)/6+sqrt(3)*exp(I*p_3)/6+sqrt(3)*exp(-I*p_3)/6+sqrt(3)*exp(-I*p_2)/6+sqrt(3)*exp(-I*p_1)/6", (1, 1): "-sqrt(3)*exp(I*p_4)/6-sqrt(3)*exp(I*p_5)/6-sqrt(3)*exp(I*p_6)/6-sqrt(3)*exp(-I*p_6)/6-sqrt(3)*exp(-I*p_5)/6-sqrt(3)*exp(-I*p_4)/6"},
        ),
        "z8": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "sqrt(3)*I*exp(I*p_1)/6+sqrt(3)*I*exp(I*p_2)/6+sqrt(3)*I*exp(I*p_3)/6-sqrt(3)*I*exp(-I*p_3)/6-sqrt(3)*I*exp(-I*p_2)/6-sqrt(3)*I*exp(-I*p_1)/6",
                (1, 1): "-sqrt(3)*I*exp(I*p_4)/6-sqrt(3)*I*exp(I*p_5)/6-sqrt(3)*I*exp(I*p_6)/6+sqrt(3)*I*exp(-I*p_6)/6+sqrt(3)*I*exp(-I*p_5)/6+sqrt(3)*I*exp(-I*p_4)/6",
            },
        ),
        "z9": ("C;C_001_1", "3a@3f", {(1, 0): "-I*exp(I*p_2)/2+I*exp(I*p_3)/2", (0, 1): "-I*exp(-I*p_3)/2+I*exp(-I*p_2)/2"}),
        "z10": ("C;C_001_1", "3a@3f", {(1, 0): "sqrt(3)*I*exp(I*p_1)/3-sqrt(3)*I*exp(I*p_2)/6-sqrt(3)*I*exp(I*p_3)/6", (0, 1): "sqrt(3)*I*exp(-I*p_3)/6+sqrt(3)*I*exp(-I*p_2)/6-sqrt(3)*I*exp(-I*p_1)/3"}),
        "z11": ("C;C_002_1", "6b@6l", {(0, 0): "-sqrt(2)*exp(I*p_2)/4+sqrt(2)*exp(I*p_3)/4+sqrt(2)*exp(-I*p_3)/4-sqrt(2)*exp(-I*p_2)/4", (1, 1): "sqrt(2)*exp(I*p_5)/4-sqrt(2)*exp(I*p_6)/4-sqrt(2)*exp(-I*p_6)/4+sqrt(2)*exp(-I*p_5)/4"}),
        "z12": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "sqrt(6)*exp(I*p_1)/6-sqrt(6)*exp(I*p_2)/12-sqrt(6)*exp(I*p_3)/12-sqrt(6)*exp(-I*p_3)/12-sqrt(6)*exp(-I*p_2)/12+sqrt(6)*exp(-I*p_1)/6",
                (1, 1): "-sqrt(6)*exp(I*p_4)/6+sqrt(6)*exp(I*p_5)/12+sqrt(6)*exp(I*p_6)/12+sqrt(6)*exp(-I*p_6)/12+sqrt(6)*exp(-I*p_5)/12-sqrt(6)*exp(-I*p_4)/6",
            },
        ),
        "z13": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "sqrt(6)*I*exp(I*p_1)/6-sqrt(6)*I*exp(I*p_2)/12-sqrt(6)*I*exp(I*p_3)/12+sqrt(6)*I*exp(-I*p_3)/12+sqrt(6)*I*exp(-I*p_2)/12-sqrt(6)*I*exp(-I*p_1)/6",
                (1, 1): "-sqrt(6)*I*exp(I*p_4)/6+sqrt(6)*I*exp(I*p_5)/12+sqrt(6)*I*exp(I*p_6)/12-sqrt(6)*I*exp(-I*p_6)/12-sqrt(6)*I*exp(-I*p_5)/12+sqrt(6)*I*exp(-I*p_4)/6",
            },
        ),
        "z14": ("C;C_002_1", "6b@6l", {(0, 0): "sqrt(2)*I*exp(I*p_2)/4-sqrt(2)*I*exp(I*p_3)/4+sqrt(2)*I*exp(-I*p_3)/4-sqrt(2)*I*exp(-I*p_2)/4", (1, 1): "-sqrt(2)*I*exp(I*p_5)/4+sqrt(2)*I*exp(I*p_6)/4-sqrt(2)*I*exp(-I*p_6)/4+sqrt(2)*I*exp(-I*p_5)/4"}),
        "z15": ("C;C_001_1", "3a@3f", {(1, 0): "sqrt(3)*exp(I*p_1)/3-sqrt(3)*exp(I*p_2)/6-sqrt(3)*exp(I*p_3)/6", (0, 1): "-sqrt(3)*exp(-I*p_3)/6-sqrt(3)*exp(-I*p_2)/6+sqrt(3)*exp(-I*p_1)/3"}),
        "z16": ("C;C_001_1", "3a@3f", {(1, 0): "exp(I*p_2)/2-exp(I*p_3)/2", (0, 1): "-exp(-I*p_3)/2+exp(-I*p_2)/2"}),
        "z17": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "sqrt(6)*exp(I*p_1)/6-sqrt(6)*exp(I*p_2)/12-sqrt(6)*exp(I*p_3)/12-sqrt(6)*exp(-I*p_3)/12-sqrt(6)*exp(-I*p_2)/12+sqrt(6)*exp(-I*p_1)/6",
                (1, 1): "sqrt(6)*exp(I*p_4)/6-sqrt(6)*exp(I*p_5)/12-sqrt(6)*exp(I*p_6)/12-sqrt(6)*exp(-I*p_6)/12-sqrt(6)*exp(-I*p_5)/12+sqrt(6)*exp(-I*p_4)/6",
            },
        ),
        "z18": ("C;C_002_1", "6b@6l", {(0, 0): "sqrt(2)*exp(I*p_2)/4-sqrt(2)*exp(I*p_3)/4-sqrt(2)*exp(-I*p_3)/4+sqrt(2)*exp(-I*p_2)/4", (1, 1): "sqrt(2)*exp(I*p_5)/4-sqrt(2)*exp(I*p_6)/4-sqrt(2)*exp(-I*p_6)/4+sqrt(2)*exp(-I*p_5)/4"}),
        "z19": ("C;C_002_1", "6b@6l", {(0, 0): "sqrt(2)*I*exp(I*p_2)/4-sqrt(2)*I*exp(I*p_3)/4+sqrt(2)*I*exp(-I*p_3)/4-sqrt(2)*I*exp(-I*p_2)/4", (1, 1): "sqrt(2)*I*exp(I*p_5)/4-sqrt(2)*I*exp(I*p_6)/4+sqrt(2)*I*exp(-I*p_6)/4-sqrt(2)*I*exp(-I*p_5)/4"}),
        "z20": (
            "C;C_002_1",
            "6b@6l",
            {
                (0, 0): "-sqrt(6)*I*exp(I*p_1)/6+sqrt(6)*I*exp(I*p_2)/12+sqrt(6)*I*exp(I*p_3)/12-sqrt(6)*I*exp(-I*p_3)/12-sqrt(6)*I*exp(-I*p_2)/12+sqrt(6)*I*exp(-I*p_1)/6",
                (1, 1): "-sqrt(6)*I*exp(I*p_4)/6+sqrt(6)*I*exp(I*p_5)/12+sqrt(6)*I*exp(I*p_6)/12-sqrt(6)*I*exp(-I*p_6)/12-sqrt(6)*I*exp(-I*p_5)/12+sqrt(6)*I*exp(-I*p_4)/6",
            },
        ),
    },
}
