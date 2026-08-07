"""
Selected SAMB matrix in momentum representation.
- model (str): model name.
- source (str): source binary.
- created (str): binary created date.
- dimension (int): matrix size.
- ket (list): ket name of full matrix, [name].
- index (dict): ket index, dict[(site,sublattice,rank), (top_index,size)].
- cluster_vector (dict): cluster vector, dict[site/bond name, dict[kb, expression] ].
- k_multipole (dict): momentum multipole in terms of p_n=k.b_n, dict[wyckoff, dict[idx, (k_multipole, symmetry)] ].
- k_matrix (dict): momentum matrix, dict[tag, (site/bond_name, wyckoff, dict[(m,n), value]) ].
"""

graphene_k = {
    "model": "graphene",
    "source": "graphene.pkl",
    "created": "2026-08-07 20:24:47",
    "dimension": 2,
    "ket": ["pz@A(1)", "pz@A(2)"],
    "index": {("A", 1, 1): (0, 1), ("A", 2, 1): (1, 1)},
    "cluster_vector": {
        "A": {},
        "A;A_001_1": {"p_1": "0.33333333*k_1+0.66666666*k_2", "p_2": "-0.66666666*k_1-0.33333333*k_2", "p_3": "0.33333333*k_1-0.33333333*k_2"},
        "A;A_002_1": {"p_1": "1.0*k_1", "p_2": "1.0*k_2", "p_3": "-1.0*k_1-1.0*k_2", "p_4": "-1.0*k_1", "p_5": "-1.0*k_2", "p_6": "1.0*k_1+1.0*k_2"},
        "A;A_003_1": {"p_1": "-0.66666667*k_1-1.33333334*k_2", "p_2": "1.33333334*k_1+0.66666667*k_2", "p_3": "-0.66666667*k_1+0.66666667*k_2"},
        "A;A_004_1": {"p_1": "-1.66666667*k_1-1.33333333*k_2", "p_2": "1.33333333*k_1-0.33333334*k_2", "p_3": "0.33333334*k_1+1.66666667*k_2", "p_4": "-1.33333333*k_1-1.66666667*k_2", "p_5": "-0.33333334*k_1+1.33333333*k_2", "p_6": "1.66666667*k_1+0.33333334*k_2"},
        "A;A_005_1": {"p_1": "1.0*k_1+2.0*k_2", "p_2": "-2.0*k_1-1.0*k_2", "p_3": "1.0*k_1-1.0*k_2", "p_4": "-1.0*k_1-2.0*k_2", "p_5": "2.0*k_1+1.0*k_2", "p_6": "-1.0*k_1+1.0*k_2"},
        "A;A_006_1": {"p_1": "2.0*k_1", "p_2": "2.0*k_2", "p_3": "-2.0*k_1-2.0*k_2", "p_4": "-2.0*k_1", "p_5": "-2.0*k_2", "p_6": "2.0*k_1+2.0*k_2"},
        "A;A_007_1": {"p_1": "2.33333333*k_1+0.66666667*k_2", "p_2": "-0.66666667*k_1+1.66666666*k_2", "p_3": "-1.66666666*k_1-2.33333333*k_2", "p_4": "0.66666667*k_1+2.33333333*k_2", "p_5": "1.66666666*k_1-0.66666667*k_2", "p_6": "-2.33333333*k_1-1.66666666*k_2"},
        "A;A_008_1": {"p_1": "1.33333333*k_1+2.66666666*k_2", "p_2": "-2.66666666*k_1-1.33333333*k_2", "p_3": "1.33333333*k_1-1.33333333*k_2"},
        "A;A_009_1": {"p_1": "2.33333333*k_1+2.66666667*k_2", "p_2": "-2.66666667*k_1-0.33333334*k_2", "p_3": "0.33333334*k_1-2.33333333*k_2", "p_4": "2.66666667*k_1+2.33333333*k_2", "p_5": "-0.33333334*k_1-2.66666667*k_2", "p_6": "-2.33333333*k_1+0.33333334*k_2"},
        "A;A_010_1": {
            "p_1": "3.0*k_1+2.0*k_2",
            "p_2": "-2.0*k_1+1.0*k_2",
            "p_3": "-1.0*k_1-3.0*k_2",
            "p_4": "-3.0*k_1-2.0*k_2",
            "p_5": "2.0*k_1-1.0*k_2",
            "p_6": "1.0*k_1+3.0*k_2",
            "p_7": "2.0*k_1+3.0*k_2",
            "p_8": "1.0*k_1-2.0*k_2",
            "p_9": "-3.0*k_1-1.0*k_2",
            "p_10": "-2.0*k_1-3.0*k_2",
            "p_11": "-1.0*k_1+2.0*k_2",
            "p_12": "3.0*k_1+1.0*k_2",
        },
    },
    "k_multipole": {
        "2c": {("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[[sqrt(2)/2,sqrt(2)/2]]", "[1]"), ("Q", 3, "B1u", -1, -1, 0, 0, "q"): ("[[sqrt(2)/2,-sqrt(2)/2]]", "[sqrt(10)*y*(3*x**2-y**2)/4]")},
        "3a@3f": {
            ("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[sqrt(6)*cos(p_1)/3+sqrt(6)*cos(p_2)/3+sqrt(6)*cos(p_3)/3]", "[1]"),
            ("T", 3, "B1u", -1, -1, 0, 0, "q"): ("[sqrt(6)*sin(p_1)/3+sqrt(6)*sin(p_2)/3+sqrt(6)*sin(p_3)/3]", "[sqrt(10)*y*(3*x**2-y**2)/4]"),
            ("T", 1, "E1u", -1, -1, 0, 0, "q"): ("[-sin(p_2)+sin(p_3),2*sqrt(3)*sin(p_1)/3-sqrt(3)*sin(p_2)/3-sqrt(3)*sin(p_3)/3]", "[x,y]"),
            ("Q", 2, "E2g", -1, -1, 0, 0, "q"): ("[2*sqrt(3)*cos(p_1)/3-sqrt(3)*cos(p_2)/3-sqrt(3)*cos(p_3)/3,cos(p_2)-cos(p_3)]", "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]"),
        },
        "3b@1a": {
            ("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[sqrt(6)*cos(p_1)/3+sqrt(6)*cos(p_2)/3+sqrt(6)*cos(p_3)/3]", "[1]"),
            ("T", 3, "B1u", -1, -1, 0, 0, "q"): ("[sqrt(6)*sin(p_1)/3+sqrt(6)*sin(p_2)/3+sqrt(6)*sin(p_3)/3]", "[sqrt(10)*y*(3*x**2-y**2)/4]"),
            ("T", 1, "E1u", -1, -1, 0, 0, "q"): ("[-sin(p_2)+sin(p_3),2*sqrt(3)*sin(p_1)/3-sqrt(3)*sin(p_2)/3-sqrt(3)*sin(p_3)/3]", "[x,y]"),
            ("Q", 2, "E2g", -1, -1, 0, 0, "q"): ("[2*sqrt(3)*cos(p_1)/3-sqrt(3)*cos(p_2)/3-sqrt(3)*cos(p_3)/3,cos(p_2)-cos(p_3)]", "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]"),
        },
        "6b@6l": {
            ("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[sqrt(3)*cos(p_1)/3+sqrt(3)*cos(p_2)/3+sqrt(3)*cos(p_3)/3+sqrt(3)*cos(p_4)/3+sqrt(3)*cos(p_5)/3+sqrt(3)*cos(p_6)/3]", "[1]"),
            ("M", 1, "A2g", -1, -1, 0, 0, "q"): ("[sqrt(3)*sin(p_1)/3+sqrt(3)*sin(p_2)/3+sqrt(3)*sin(p_3)/3+sqrt(3)*sin(p_4)/3+sqrt(3)*sin(p_5)/3+sqrt(3)*sin(p_6)/3]", "[z]"),
            ("Q", 3, "B1u", -1, -1, 0, 0, "q"): ("[sqrt(3)*cos(p_1)/3+sqrt(3)*cos(p_2)/3+sqrt(3)*cos(p_3)/3-sqrt(3)*cos(p_4)/3-sqrt(3)*cos(p_5)/3-sqrt(3)*cos(p_6)/3]", "[sqrt(10)*y*(3*x**2-y**2)/4]"),
            ("T", 3, "B2u", -1, -1, 0, 0, "q"): ("[sqrt(3)*sin(p_1)/3+sqrt(3)*sin(p_2)/3+sqrt(3)*sin(p_3)/3-sqrt(3)*sin(p_4)/3-sqrt(3)*sin(p_5)/3-sqrt(3)*sin(p_6)/3]", "[sqrt(10)*x*(x**2-3*y**2)/4]"),
            ("Q", 1, "E1u", -1, -1, 0, 0, "q"): ("[-sqrt(2)*cos(p_2)/2+sqrt(2)*cos(p_3)/2+sqrt(2)*cos(p_5)/2-sqrt(2)*cos(p_6)/2,sqrt(6)*cos(p_1)/3-sqrt(6)*cos(p_2)/6-sqrt(6)*cos(p_3)/6-sqrt(6)*cos(p_4)/3+sqrt(6)*cos(p_5)/6+sqrt(6)*cos(p_6)/6]", "[x,y]"),
            ("T", 1, "E1u", -1, -1, 0, 0, "q"): ("[sqrt(6)*sin(p_1)/3-sqrt(6)*sin(p_2)/6-sqrt(6)*sin(p_3)/6-sqrt(6)*sin(p_4)/3+sqrt(6)*sin(p_5)/6+sqrt(6)*sin(p_6)/6,sqrt(2)*sin(p_2)/2-sqrt(2)*sin(p_3)/2-sqrt(2)*sin(p_5)/2+sqrt(2)*sin(p_6)/2]", "[x,y]"),
            ("Q", 2, "E2g", -1, -1, 0, 0, "q"): ("[sqrt(6)*cos(p_1)/3-sqrt(6)*cos(p_2)/6-sqrt(6)*cos(p_3)/6+sqrt(6)*cos(p_4)/3-sqrt(6)*cos(p_5)/6-sqrt(6)*cos(p_6)/6,sqrt(2)*cos(p_2)/2-sqrt(2)*cos(p_3)/2+sqrt(2)*cos(p_5)/2-sqrt(2)*cos(p_6)/2]", "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]"),
            ("T", 2, "E2g", -1, -1, 0, 0, "q"): ("[sqrt(2)*sin(p_2)/2-sqrt(2)*sin(p_3)/2+sqrt(2)*sin(p_5)/2-sqrt(2)*sin(p_6)/2,-sqrt(6)*sin(p_1)/3+sqrt(6)*sin(p_2)/6+sqrt(6)*sin(p_3)/6-sqrt(6)*sin(p_4)/3+sqrt(6)*sin(p_5)/6+sqrt(6)*sin(p_6)/6]", "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]"),
        },
        "6c@2c": {
            ("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[sqrt(3)*cos(p_1)/3+sqrt(3)*cos(p_2)/3+sqrt(3)*cos(p_3)/3+sqrt(3)*cos(p_4)/3+sqrt(3)*cos(p_5)/3+sqrt(3)*cos(p_6)/3]", "[1]"),
            ("M", 1, "A2g", -1, -1, 0, 0, "q"): ("[sqrt(3)*sin(p_1)/3+sqrt(3)*sin(p_2)/3+sqrt(3)*sin(p_3)/3+sqrt(3)*sin(p_4)/3+sqrt(3)*sin(p_5)/3+sqrt(3)*sin(p_6)/3]", "[z]"),
            ("Q", 3, "B1u", -1, -1, 0, 0, "q"): ("[sqrt(3)*cos(p_1)/3+sqrt(3)*cos(p_2)/3+sqrt(3)*cos(p_3)/3-sqrt(3)*cos(p_4)/3-sqrt(3)*cos(p_5)/3-sqrt(3)*cos(p_6)/3]", "[sqrt(10)*y*(3*x**2-y**2)/4]"),
            ("T", 3, "B2u", -1, -1, 0, 0, "q"): ("[sqrt(3)*sin(p_1)/3+sqrt(3)*sin(p_2)/3+sqrt(3)*sin(p_3)/3-sqrt(3)*sin(p_4)/3-sqrt(3)*sin(p_5)/3-sqrt(3)*sin(p_6)/3]", "[sqrt(10)*x*(x**2-3*y**2)/4]"),
            ("Q", 1, "E1u", -1, -1, 0, 0, "q"): ("[-sqrt(2)*cos(p_2)/2+sqrt(2)*cos(p_3)/2+sqrt(2)*cos(p_5)/2-sqrt(2)*cos(p_6)/2,sqrt(6)*cos(p_1)/3-sqrt(6)*cos(p_2)/6-sqrt(6)*cos(p_3)/6-sqrt(6)*cos(p_4)/3+sqrt(6)*cos(p_5)/6+sqrt(6)*cos(p_6)/6]", "[x,y]"),
            ("T", 1, "E1u", -1, -1, 0, 0, "q"): ("[sqrt(6)*sin(p_1)/3-sqrt(6)*sin(p_2)/6-sqrt(6)*sin(p_3)/6-sqrt(6)*sin(p_4)/3+sqrt(6)*sin(p_5)/6+sqrt(6)*sin(p_6)/6,sqrt(2)*sin(p_2)/2-sqrt(2)*sin(p_3)/2-sqrt(2)*sin(p_5)/2+sqrt(2)*sin(p_6)/2]", "[x,y]"),
            ("Q", 2, "E2g", -1, -1, 0, 0, "q"): ("[sqrt(6)*cos(p_1)/3-sqrt(6)*cos(p_2)/6-sqrt(6)*cos(p_3)/6+sqrt(6)*cos(p_4)/3-sqrt(6)*cos(p_5)/6-sqrt(6)*cos(p_6)/6,sqrt(2)*cos(p_2)/2-sqrt(2)*cos(p_3)/2+sqrt(2)*cos(p_5)/2-sqrt(2)*cos(p_6)/2]", "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]"),
            ("T", 2, "E2g", -1, -1, 0, 0, "q"): ("[sqrt(2)*sin(p_2)/2-sqrt(2)*sin(p_3)/2+sqrt(2)*sin(p_5)/2-sqrt(2)*sin(p_6)/2,-sqrt(6)*sin(p_1)/3+sqrt(6)*sin(p_2)/6+sqrt(6)*sin(p_3)/6-sqrt(6)*sin(p_4)/3+sqrt(6)*sin(p_5)/6+sqrt(6)*sin(p_6)/6]", "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]"),
        },
        "6a@6l": {
            ("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[sqrt(3)*cos(p_1)/3+sqrt(3)*cos(p_2)/3+sqrt(3)*cos(p_3)/3+sqrt(3)*cos(p_4)/3+sqrt(3)*cos(p_5)/3+sqrt(3)*cos(p_6)/3]", "[1]"),
            ("T", 0, "A1g", -1, -1, 0, 0, "q"): ("[sqrt(3)*sin(p_1)/3+sqrt(3)*sin(p_2)/3+sqrt(3)*sin(p_3)/3+sqrt(3)*sin(p_4)/3+sqrt(3)*sin(p_5)/3+sqrt(3)*sin(p_6)/3]", "[1]"),
            ("Q", 3, "B1u", -1, -1, 0, 0, "q"): ("[sqrt(3)*cos(p_1)/3+sqrt(3)*cos(p_2)/3+sqrt(3)*cos(p_3)/3-sqrt(3)*cos(p_4)/3-sqrt(3)*cos(p_5)/3-sqrt(3)*cos(p_6)/3]", "[sqrt(10)*y*(3*x**2-y**2)/4]"),
            ("T", 3, "B1u", -1, -1, 0, 0, "q"): ("[sqrt(3)*sin(p_1)/3+sqrt(3)*sin(p_2)/3+sqrt(3)*sin(p_3)/3-sqrt(3)*sin(p_4)/3-sqrt(3)*sin(p_5)/3-sqrt(3)*sin(p_6)/3]", "[sqrt(10)*y*(3*x**2-y**2)/4]"),
            ("Q", 1, "E1u", -1, -1, 0, 0, "q"): ("[-sqrt(2)*cos(p_2)/2+sqrt(2)*cos(p_3)/2+sqrt(2)*cos(p_5)/2-sqrt(2)*cos(p_6)/2,sqrt(6)*cos(p_1)/3-sqrt(6)*cos(p_2)/6-sqrt(6)*cos(p_3)/6-sqrt(6)*cos(p_4)/3+sqrt(6)*cos(p_5)/6+sqrt(6)*cos(p_6)/6]", "[x,y]"),
            ("T", 1, "E1u", -1, -1, 0, 0, "q"): ("[-sqrt(2)*sin(p_2)/2+sqrt(2)*sin(p_3)/2+sqrt(2)*sin(p_5)/2-sqrt(2)*sin(p_6)/2,sqrt(6)*sin(p_1)/3-sqrt(6)*sin(p_2)/6-sqrt(6)*sin(p_3)/6-sqrt(6)*sin(p_4)/3+sqrt(6)*sin(p_5)/6+sqrt(6)*sin(p_6)/6]", "[x,y]"),
            ("Q", 2, "E2g", -1, -1, 0, 0, "q"): ("[sqrt(6)*cos(p_1)/3-sqrt(6)*cos(p_2)/6-sqrt(6)*cos(p_3)/6+sqrt(6)*cos(p_4)/3-sqrt(6)*cos(p_5)/6-sqrt(6)*cos(p_6)/6,sqrt(2)*cos(p_2)/2-sqrt(2)*cos(p_3)/2+sqrt(2)*cos(p_5)/2-sqrt(2)*cos(p_6)/2]", "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]"),
            ("T", 2, "E2g", -1, -1, 0, 0, "q"): ("[sqrt(6)*sin(p_1)/3-sqrt(6)*sin(p_2)/6-sqrt(6)*sin(p_3)/6+sqrt(6)*sin(p_4)/3-sqrt(6)*sin(p_5)/6-sqrt(6)*sin(p_6)/6,sqrt(2)*sin(p_2)/2-sqrt(2)*sin(p_3)/2+sqrt(2)*sin(p_5)/2-sqrt(2)*sin(p_6)/2]", "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]"),
        },
        "6d@3f": {
            ("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[sqrt(3)*cos(p_1)/3+sqrt(3)*cos(p_2)/3+sqrt(3)*cos(p_3)/3+sqrt(3)*cos(p_4)/3+sqrt(3)*cos(p_5)/3+sqrt(3)*cos(p_6)/3]", "[1]"),
            ("Q", 6, "A2g", -1, -1, 0, 0, "q"): ("[sqrt(3)*cos(p_1)/3+sqrt(3)*cos(p_2)/3+sqrt(3)*cos(p_3)/3-sqrt(3)*cos(p_4)/3-sqrt(3)*cos(p_5)/3-sqrt(3)*cos(p_6)/3]", "[sqrt(462)*x*y*(x**2-3*y**2)*(3*x**2-y**2)/16]"),
            ("T", 3, "B1u", -1, -1, 0, 0, "q"): ("[sqrt(3)*sin(p_1)/3+sqrt(3)*sin(p_2)/3+sqrt(3)*sin(p_3)/3-sqrt(3)*sin(p_4)/3-sqrt(3)*sin(p_5)/3-sqrt(3)*sin(p_6)/3]", "[sqrt(10)*y*(3*x**2-y**2)/4]"),
            ("T", 3, "B2u", -1, -1, 0, 0, "q"): ("[sqrt(3)*sin(p_1)/3+sqrt(3)*sin(p_2)/3+sqrt(3)*sin(p_3)/3+sqrt(3)*sin(p_4)/3+sqrt(3)*sin(p_5)/3+sqrt(3)*sin(p_6)/3]", "[sqrt(10)*x*(x**2-3*y**2)/4]"),
            ("T", 1, "E1u", -1, 1, 0, 0, "q"): (
                "[5*sqrt(42)*sin(p_1)/42-2*sqrt(42)*sin(p_2)/21-sqrt(42)*sin(p_3)/42-sqrt(42)*sin(p_4)/42+5*sqrt(42)*sin(p_5)/42-2*sqrt(42)*sin(p_6)/21,sqrt(14)*sin(p_1)/14+sqrt(14)*sin(p_2)/7-3*sqrt(14)*sin(p_3)/14+3*sqrt(14)*sin(p_4)/14-sqrt(14)*sin(p_5)/14-sqrt(14)*sin(p_6)/7]",
                "[x,y]",
            ),
            ("T", 1, "E1u", -1, 2, 0, 0, "q"): (
                "[sqrt(14)*sin(p_1)/14+sqrt(14)*sin(p_2)/7-3*sqrt(14)*sin(p_3)/14-3*sqrt(14)*sin(p_4)/14+sqrt(14)*sin(p_5)/14+sqrt(14)*sin(p_6)/7,-5*sqrt(42)*sin(p_1)/42+2*sqrt(42)*sin(p_2)/21+sqrt(42)*sin(p_3)/42-sqrt(42)*sin(p_4)/42+5*sqrt(42)*sin(p_5)/42-2*sqrt(42)*sin(p_6)/21]",
                "[x,y]",
            ),
            ("Q", 2, "E2g", -1, -1, 0, 0, "q"): (
                "[11*sqrt(6)*cos(p_1)/42+sqrt(6)*cos(p_2)/21-13*sqrt(6)*cos(p_3)/42-13*sqrt(6)*cos(p_4)/42+11*sqrt(6)*cos(p_5)/42+sqrt(6)*cos(p_6)/21,-5*sqrt(2)*cos(p_1)/14+4*sqrt(2)*cos(p_2)/7-3*sqrt(2)*cos(p_3)/14+3*sqrt(2)*cos(p_4)/14+5*sqrt(2)*cos(p_5)/14-4*sqrt(2)*cos(p_6)/7]",
                "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]",
            ),
            ("Q", 4, "E2g", 1, -1, 0, 0, "q"): (
                "[5*sqrt(2)*cos(p_1)/14-4*sqrt(2)*cos(p_2)/7+3*sqrt(2)*cos(p_3)/14+3*sqrt(2)*cos(p_4)/14+5*sqrt(2)*cos(p_5)/14-4*sqrt(2)*cos(p_6)/7,11*sqrt(6)*cos(p_1)/42+sqrt(6)*cos(p_2)/21-13*sqrt(6)*cos(p_3)/42+13*sqrt(6)*cos(p_4)/42-11*sqrt(6)*cos(p_5)/42-sqrt(6)*cos(p_6)/21]",
                "[sqrt(35)*(x**2-2*x*y-y**2)*(x**2+2*x*y-y**2)/8,sqrt(35)*x*y*(x-y)*(x+y)/2]",
            ),
        },
        "12d@6l": {
            ("Q", 0, "A1g", -1, -1, 0, 0, "q"): ("[sqrt(6)*cos(p_1)/6+sqrt(6)*cos(p_10)/6+sqrt(6)*cos(p_11)/6+sqrt(6)*cos(p_12)/6+sqrt(6)*cos(p_2)/6+sqrt(6)*cos(p_3)/6+sqrt(6)*cos(p_4)/6+sqrt(6)*cos(p_5)/6+sqrt(6)*cos(p_6)/6+sqrt(6)*cos(p_7)/6+sqrt(6)*cos(p_8)/6+sqrt(6)*cos(p_9)/6]", "[1]"),
            ("T", 0, "A1g", -1, -1, 0, 0, "q"): ("[sqrt(6)*sin(p_1)/6+sqrt(6)*sin(p_10)/6+sqrt(6)*sin(p_11)/6+sqrt(6)*sin(p_12)/6+sqrt(6)*sin(p_2)/6+sqrt(6)*sin(p_3)/6+sqrt(6)*sin(p_4)/6+sqrt(6)*sin(p_5)/6+sqrt(6)*sin(p_6)/6+sqrt(6)*sin(p_7)/6+sqrt(6)*sin(p_8)/6+sqrt(6)*sin(p_9)/6]", "[1]"),
            ("M", 1, "A2g", -1, -1, 0, 0, "q"): ("[sqrt(6)*sin(p_1)/6-sqrt(6)*sin(p_10)/6-sqrt(6)*sin(p_11)/6-sqrt(6)*sin(p_12)/6+sqrt(6)*sin(p_2)/6+sqrt(6)*sin(p_3)/6+sqrt(6)*sin(p_4)/6+sqrt(6)*sin(p_5)/6+sqrt(6)*sin(p_6)/6-sqrt(6)*sin(p_7)/6-sqrt(6)*sin(p_8)/6-sqrt(6)*sin(p_9)/6]", "[z]"),
            ("Q", 6, "A2g", -1, -1, 0, 0, "q"): (
                "[sqrt(6)*cos(p_1)/6-sqrt(6)*cos(p_10)/6-sqrt(6)*cos(p_11)/6-sqrt(6)*cos(p_12)/6+sqrt(6)*cos(p_2)/6+sqrt(6)*cos(p_3)/6+sqrt(6)*cos(p_4)/6+sqrt(6)*cos(p_5)/6+sqrt(6)*cos(p_6)/6-sqrt(6)*cos(p_7)/6-sqrt(6)*cos(p_8)/6-sqrt(6)*cos(p_9)/6]",
                "[sqrt(462)*x*y*(x**2-3*y**2)*(3*x**2-y**2)/16]",
            ),
            ("Q", 3, "B1u", -1, -1, 0, 0, "q"): (
                "[sqrt(6)*cos(p_1)/6+sqrt(6)*cos(p_10)/6+sqrt(6)*cos(p_11)/6+sqrt(6)*cos(p_12)/6+sqrt(6)*cos(p_2)/6+sqrt(6)*cos(p_3)/6-sqrt(6)*cos(p_4)/6-sqrt(6)*cos(p_5)/6-sqrt(6)*cos(p_6)/6-sqrt(6)*cos(p_7)/6-sqrt(6)*cos(p_8)/6-sqrt(6)*cos(p_9)/6]",
                "[sqrt(10)*y*(3*x**2-y**2)/4]",
            ),
            ("T", 3, "B1u", -1, -1, 0, 0, "q"): (
                "[sqrt(6)*sin(p_1)/6+sqrt(6)*sin(p_10)/6+sqrt(6)*sin(p_11)/6+sqrt(6)*sin(p_12)/6+sqrt(6)*sin(p_2)/6+sqrt(6)*sin(p_3)/6-sqrt(6)*sin(p_4)/6-sqrt(6)*sin(p_5)/6-sqrt(6)*sin(p_6)/6-sqrt(6)*sin(p_7)/6-sqrt(6)*sin(p_8)/6-sqrt(6)*sin(p_9)/6]",
                "[sqrt(10)*y*(3*x**2-y**2)/4]",
            ),
            ("Q", 3, "B2u", -1, -1, 0, 0, "q"): (
                "[sqrt(6)*cos(p_1)/6-sqrt(6)*cos(p_10)/6-sqrt(6)*cos(p_11)/6-sqrt(6)*cos(p_12)/6+sqrt(6)*cos(p_2)/6+sqrt(6)*cos(p_3)/6-sqrt(6)*cos(p_4)/6-sqrt(6)*cos(p_5)/6-sqrt(6)*cos(p_6)/6+sqrt(6)*cos(p_7)/6+sqrt(6)*cos(p_8)/6+sqrt(6)*cos(p_9)/6]",
                "[sqrt(10)*x*(x**2-3*y**2)/4]",
            ),
            ("T", 3, "B2u", -1, -1, 0, 0, "q"): (
                "[sqrt(6)*sin(p_1)/6-sqrt(6)*sin(p_10)/6-sqrt(6)*sin(p_11)/6-sqrt(6)*sin(p_12)/6+sqrt(6)*sin(p_2)/6+sqrt(6)*sin(p_3)/6-sqrt(6)*sin(p_4)/6-sqrt(6)*sin(p_5)/6-sqrt(6)*sin(p_6)/6+sqrt(6)*sin(p_7)/6+sqrt(6)*sin(p_8)/6+sqrt(6)*sin(p_9)/6]",
                "[sqrt(10)*x*(x**2-3*y**2)/4]",
            ),
            ("Q", 1, "E1u", -1, -1, 0, 0, "q"): (
                "[5*sqrt(21)*cos(p_1)/42+sqrt(21)*cos(p_10)/42-5*sqrt(21)*cos(p_11)/42+2*sqrt(21)*cos(p_12)/21-2*sqrt(21)*cos(p_2)/21-sqrt(21)*cos(p_3)/42-5*sqrt(21)*cos(p_4)/42+2*sqrt(21)*cos(p_5)/21+sqrt(21)*cos(p_6)/42-sqrt(21)*cos(p_7)/42+5*sqrt(21)*cos(p_8)/42-2*sqrt(21)*cos(p_9)/21,sqrt(7)*cos(p_1)/14-3*sqrt(7)*cos(p_10)/14+sqrt(7)*cos(p_11)/14+sqrt(7)*cos(p_12)/7+sqrt(7)*cos(p_2)/7-3*sqrt(7)*cos(p_3)/14-sqrt(7)*cos(p_4)/14-sqrt(7)*cos(p_5)/7+3*sqrt(7)*cos(p_6)/14+3*sqrt(7)*cos(p_7)/14-sqrt(7)*cos(p_8)/14-sqrt(7)*cos(p_9)/7]",
                "[x,y]",
            ),
            ("T", 1, "E1u", -1, 1, 0, 0, "q"): (
                "[5*sqrt(21)*sin(p_1)/42+sqrt(21)*sin(p_10)/42-5*sqrt(21)*sin(p_11)/42+2*sqrt(21)*sin(p_12)/21-2*sqrt(21)*sin(p_2)/21-sqrt(21)*sin(p_3)/42-5*sqrt(21)*sin(p_4)/42+2*sqrt(21)*sin(p_5)/21+sqrt(21)*sin(p_6)/42-sqrt(21)*sin(p_7)/42+5*sqrt(21)*sin(p_8)/42-2*sqrt(21)*sin(p_9)/21,sqrt(7)*sin(p_1)/14-3*sqrt(7)*sin(p_10)/14+sqrt(7)*sin(p_11)/14+sqrt(7)*sin(p_12)/7+sqrt(7)*sin(p_2)/7-3*sqrt(7)*sin(p_3)/14-sqrt(7)*sin(p_4)/14-sqrt(7)*sin(p_5)/7+3*sqrt(7)*sin(p_6)/14+3*sqrt(7)*sin(p_7)/14-sqrt(7)*sin(p_8)/14-sqrt(7)*sin(p_9)/7]",
                "[x,y]",
            ),
            ("T", 1, "E1u", -1, 2, 0, 0, "q"): (
                "[sqrt(7)*sin(p_1)/14+3*sqrt(7)*sin(p_10)/14-sqrt(7)*sin(p_11)/14-sqrt(7)*sin(p_12)/7+sqrt(7)*sin(p_2)/7-3*sqrt(7)*sin(p_3)/14-sqrt(7)*sin(p_4)/14-sqrt(7)*sin(p_5)/7+3*sqrt(7)*sin(p_6)/14-3*sqrt(7)*sin(p_7)/14+sqrt(7)*sin(p_8)/14+sqrt(7)*sin(p_9)/7,-5*sqrt(21)*sin(p_1)/42+sqrt(21)*sin(p_10)/42-5*sqrt(21)*sin(p_11)/42+2*sqrt(21)*sin(p_12)/21+2*sqrt(21)*sin(p_2)/21+sqrt(21)*sin(p_3)/42+5*sqrt(21)*sin(p_4)/42-2*sqrt(21)*sin(p_5)/21-sqrt(21)*sin(p_6)/42-sqrt(21)*sin(p_7)/42+5*sqrt(21)*sin(p_8)/42-2*sqrt(21)*sin(p_9)/21]",
                "[x,y]",
            ),
            ("Q", 5, "E1u", 1, -1, 0, 0, "q"): (
                "[sqrt(7)*cos(p_1)/14+3*sqrt(7)*cos(p_10)/14-sqrt(7)*cos(p_11)/14-sqrt(7)*cos(p_12)/7+sqrt(7)*cos(p_2)/7-3*sqrt(7)*cos(p_3)/14-sqrt(7)*cos(p_4)/14-sqrt(7)*cos(p_5)/7+3*sqrt(7)*cos(p_6)/14-3*sqrt(7)*cos(p_7)/14+sqrt(7)*cos(p_8)/14+sqrt(7)*cos(p_9)/7,-5*sqrt(21)*cos(p_1)/42+sqrt(21)*cos(p_10)/42-5*sqrt(21)*cos(p_11)/42+2*sqrt(21)*cos(p_12)/21+2*sqrt(21)*cos(p_2)/21+sqrt(21)*cos(p_3)/42+5*sqrt(21)*cos(p_4)/42-2*sqrt(21)*cos(p_5)/21-sqrt(21)*cos(p_6)/42-sqrt(21)*cos(p_7)/42+5*sqrt(21)*cos(p_8)/42-2*sqrt(21)*cos(p_9)/21]",
                "[3*sqrt(14)*x*(x**4-10*x**2*y**2+5*y**4)/16,-3*sqrt(14)*y*(5*x**4-10*x**2*y**2+y**4)/16]",
            ),
            ("Q", 2, "E2g", -1, -1, 0, 0, "q"): (
                "[11*sqrt(3)*cos(p_1)/42-13*sqrt(3)*cos(p_10)/42+11*sqrt(3)*cos(p_11)/42+sqrt(3)*cos(p_12)/21+sqrt(3)*cos(p_2)/21-13*sqrt(3)*cos(p_3)/42+11*sqrt(3)*cos(p_4)/42+sqrt(3)*cos(p_5)/21-13*sqrt(3)*cos(p_6)/42-13*sqrt(3)*cos(p_7)/42+11*sqrt(3)*cos(p_8)/42+sqrt(3)*cos(p_9)/21,-5*cos(p_1)/14+3*cos(p_10)/14+5*cos(p_11)/14-4*cos(p_12)/7+4*cos(p_2)/7-3*cos(p_3)/14-5*cos(p_4)/14+4*cos(p_5)/7-3*cos(p_6)/14+3*cos(p_7)/14+5*cos(p_8)/14-4*cos(p_9)/7]",
                "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]",
            ),
            ("T", 2, "E2g", -1, 1, 0, 0, "q"): (
                "[11*sqrt(3)*sin(p_1)/42-13*sqrt(3)*sin(p_10)/42+11*sqrt(3)*sin(p_11)/42+sqrt(3)*sin(p_12)/21+sqrt(3)*sin(p_2)/21-13*sqrt(3)*sin(p_3)/42+11*sqrt(3)*sin(p_4)/42+sqrt(3)*sin(p_5)/21-13*sqrt(3)*sin(p_6)/42-13*sqrt(3)*sin(p_7)/42+11*sqrt(3)*sin(p_8)/42+sqrt(3)*sin(p_9)/21,-5*sin(p_1)/14+3*sin(p_10)/14+5*sin(p_11)/14-4*sin(p_12)/7+4*sin(p_2)/7-3*sin(p_3)/14-5*sin(p_4)/14+4*sin(p_5)/7-3*sin(p_6)/14+3*sin(p_7)/14+5*sin(p_8)/14-4*sin(p_9)/7]",
                "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]",
            ),
            ("T", 2, "E2g", -1, 2, 0, 0, "q"): (
                "[5*sin(p_1)/14+3*sin(p_10)/14+5*sin(p_11)/14-4*sin(p_12)/7-4*sin(p_2)/7+3*sin(p_3)/14+5*sin(p_4)/14-4*sin(p_5)/7+3*sin(p_6)/14+3*sin(p_7)/14+5*sin(p_8)/14-4*sin(p_9)/7,11*sqrt(3)*sin(p_1)/42+13*sqrt(3)*sin(p_10)/42-11*sqrt(3)*sin(p_11)/42-sqrt(3)*sin(p_12)/21+sqrt(3)*sin(p_2)/21-13*sqrt(3)*sin(p_3)/42+11*sqrt(3)*sin(p_4)/42+sqrt(3)*sin(p_5)/21-13*sqrt(3)*sin(p_6)/42+13*sqrt(3)*sin(p_7)/42-11*sqrt(3)*sin(p_8)/42-sqrt(3)*sin(p_9)/21]",
                "[sqrt(3)*(x-y)*(x+y)/2,-sqrt(3)*x*y]",
            ),
            ("Q", 4, "E2g", 1, -1, 0, 0, "q"): (
                "[5*cos(p_1)/14+3*cos(p_10)/14+5*cos(p_11)/14-4*cos(p_12)/7-4*cos(p_2)/7+3*cos(p_3)/14+5*cos(p_4)/14-4*cos(p_5)/7+3*cos(p_6)/14+3*cos(p_7)/14+5*cos(p_8)/14-4*cos(p_9)/7,11*sqrt(3)*cos(p_1)/42+13*sqrt(3)*cos(p_10)/42-11*sqrt(3)*cos(p_11)/42-sqrt(3)*cos(p_12)/21+sqrt(3)*cos(p_2)/21-13*sqrt(3)*cos(p_3)/42+11*sqrt(3)*cos(p_4)/42+sqrt(3)*cos(p_5)/21-13*sqrt(3)*cos(p_6)/42+13*sqrt(3)*cos(p_7)/42-11*sqrt(3)*cos(p_8)/42-sqrt(3)*cos(p_9)/21]",
                "[sqrt(35)*(x**2-2*x*y-y**2)*(x**2+2*x*y-y**2)/8,sqrt(35)*x*y*(x-y)*(x+y)/2]",
            ),
        },
    },
    "k_matrix": {
        "z1": ("A", "2c", {(0, 0): "sqrt(2)/2", (1, 1): "sqrt(2)/2"}),
        "z2": ("A;A_001_1", "3a@3f", {(1, 0): "sqrt(6)*(I*sin(p_1)+I*sin(p_2)+I*sin(p_3)+cos(p_1)+cos(p_2)+cos(p_3))/6", (0, 1): "sqrt(6)*(-I*sin(p_1)-I*sin(p_2)-I*sin(p_3)+cos(p_1)+cos(p_2)+cos(p_3))/6"}),
        "z3": ("A;A_002_1", "6b@6l", {(0, 0): "sqrt(3)*(cos(p_1)+cos(p_2)+cos(p_3))/3", (1, 1): "sqrt(3)*(cos(p_4)+cos(p_5)+cos(p_6))/3"}),
        "z4": ("A;A_003_1", "3b@1a", {(1, 0): "sqrt(6)*(I*sin(p_1)+I*sin(p_2)+I*sin(p_3)+cos(p_1)+cos(p_2)+cos(p_3))/6", (0, 1): "sqrt(6)*(-I*sin(p_1)-I*sin(p_2)-I*sin(p_3)+cos(p_1)+cos(p_2)+cos(p_3))/6"}),
        "z5": (
            "A;A_004_1",
            "6d@3f",
            {(1, 0): "sqrt(3)*(I*sin(p_1)+I*sin(p_2)+I*sin(p_3)-I*sin(p_4)-I*sin(p_5)-I*sin(p_6)+cos(p_1)+cos(p_2)+cos(p_3)+cos(p_4)+cos(p_5)+cos(p_6))/6", (0, 1): "sqrt(3)*(-I*sin(p_1)-I*sin(p_2)-I*sin(p_3)+I*sin(p_4)+I*sin(p_5)+I*sin(p_6)+cos(p_1)+cos(p_2)+cos(p_3)+cos(p_4)+cos(p_5)+cos(p_6))/6"},
        ),
        "z6": ("A;A_005_1", "6a@6l", {(0, 0): "sqrt(3)*(cos(p_1)+cos(p_2)+cos(p_3))/3", (1, 1): "sqrt(3)*(cos(p_4)+cos(p_5)+cos(p_6))/3"}),
        "z7": ("A;A_006_1", "6c@2c", {(0, 0): "sqrt(3)*(cos(p_1)+cos(p_2)+cos(p_3))/3", (1, 1): "sqrt(3)*(cos(p_4)+cos(p_5)+cos(p_6))/3"}),
        "z8": (
            "A;A_007_1",
            "6d@3f",
            {(1, 0): "sqrt(3)*(I*sin(p_1)+I*sin(p_2)+I*sin(p_3)-I*sin(p_4)-I*sin(p_5)-I*sin(p_6)+cos(p_1)+cos(p_2)+cos(p_3)+cos(p_4)+cos(p_5)+cos(p_6))/6", (0, 1): "sqrt(3)*(-I*sin(p_1)-I*sin(p_2)-I*sin(p_3)+I*sin(p_4)+I*sin(p_5)+I*sin(p_6)+cos(p_1)+cos(p_2)+cos(p_3)+cos(p_4)+cos(p_5)+cos(p_6))/6"},
        ),
        "z9": ("A;A_008_1", "3b@1a", {(1, 0): "sqrt(6)*(I*sin(p_1)+I*sin(p_2)+I*sin(p_3)+cos(p_1)+cos(p_2)+cos(p_3))/6", (0, 1): "sqrt(6)*(-I*sin(p_1)-I*sin(p_2)-I*sin(p_3)+cos(p_1)+cos(p_2)+cos(p_3))/6"}),
        "z10": (
            "A;A_009_1",
            "6d@3f",
            {(1, 0): "sqrt(3)*(I*sin(p_1)+I*sin(p_2)+I*sin(p_3)-I*sin(p_4)-I*sin(p_5)-I*sin(p_6)+cos(p_1)+cos(p_2)+cos(p_3)+cos(p_4)+cos(p_5)+cos(p_6))/6", (0, 1): "sqrt(3)*(-I*sin(p_1)-I*sin(p_2)-I*sin(p_3)+I*sin(p_4)+I*sin(p_5)+I*sin(p_6)+cos(p_1)+cos(p_2)+cos(p_3)+cos(p_4)+cos(p_5)+cos(p_6))/6"},
        ),
        "z11": ("A;A_010_1", "12d@6l", {(0, 0): "sqrt(6)*(cos(p_1)+cos(p_10)+cos(p_11)+cos(p_12)+cos(p_2)+cos(p_3))/6", (1, 1): "sqrt(6)*(cos(p_4)+cos(p_5)+cos(p_6)+cos(p_7)+cos(p_8)+cos(p_9))/6"}),
        "z12": (
            "A;A_004_1",
            "6d@3f",
            {(1, 0): "sqrt(3)*(I*sin(p_1)+I*sin(p_2)+I*sin(p_3)+I*sin(p_4)+I*sin(p_5)+I*sin(p_6)+cos(p_1)+cos(p_2)+cos(p_3)-cos(p_4)-cos(p_5)-cos(p_6))/6", (0, 1): "sqrt(3)*(-I*sin(p_1)-I*sin(p_2)-I*sin(p_3)-I*sin(p_4)-I*sin(p_5)-I*sin(p_6)+cos(p_1)+cos(p_2)+cos(p_3)-cos(p_4)-cos(p_5)-cos(p_6))/6"},
        ),
        "z13": (
            "A;A_007_1",
            "6d@3f",
            {(1, 0): "sqrt(3)*(I*sin(p_1)+I*sin(p_2)+I*sin(p_3)+I*sin(p_4)+I*sin(p_5)+I*sin(p_6)+cos(p_1)+cos(p_2)+cos(p_3)-cos(p_4)-cos(p_5)-cos(p_6))/6", (0, 1): "sqrt(3)*(-I*sin(p_1)-I*sin(p_2)-I*sin(p_3)-I*sin(p_4)-I*sin(p_5)-I*sin(p_6)+cos(p_1)+cos(p_2)+cos(p_3)-cos(p_4)-cos(p_5)-cos(p_6))/6"},
        ),
        "z14": (
            "A;A_009_1",
            "6d@3f",
            {(1, 0): "sqrt(3)*(I*sin(p_1)+I*sin(p_2)+I*sin(p_3)+I*sin(p_4)+I*sin(p_5)+I*sin(p_6)+cos(p_1)+cos(p_2)+cos(p_3)-cos(p_4)-cos(p_5)-cos(p_6))/6", (0, 1): "sqrt(3)*(-I*sin(p_1)-I*sin(p_2)-I*sin(p_3)-I*sin(p_4)-I*sin(p_5)-I*sin(p_6)+cos(p_1)+cos(p_2)+cos(p_3)-cos(p_4)-cos(p_5)-cos(p_6))/6"},
        ),
        "z15": ("A;A_010_1", "12d@6l", {(0, 0): "sqrt(6)*(cos(p_1)-cos(p_10)-cos(p_11)-cos(p_12)+cos(p_2)+cos(p_3))/6", (1, 1): "sqrt(6)*(cos(p_4)+cos(p_5)+cos(p_6)-cos(p_7)-cos(p_8)-cos(p_9))/6"}),
        "z16": ("A", "2c", {(0, 0): "sqrt(2)/2", (1, 1): "-sqrt(2)/2"}),
        "z17": ("A;A_002_1", "6b@6l", {(0, 0): "sqrt(3)*(cos(p_1)+cos(p_2)+cos(p_3))/3", (1, 1): "-sqrt(3)*(cos(p_4)+cos(p_5)+cos(p_6))/3"}),
        "z18": ("A;A_005_1", "6a@6l", {(0, 0): "sqrt(3)*(cos(p_1)+cos(p_2)+cos(p_3))/3", (1, 1): "-sqrt(3)*(cos(p_4)+cos(p_5)+cos(p_6))/3"}),
        "z19": ("A;A_006_1", "6c@2c", {(0, 0): "sqrt(3)*(cos(p_1)+cos(p_2)+cos(p_3))/3", (1, 1): "-sqrt(3)*(cos(p_4)+cos(p_5)+cos(p_6))/3"}),
        "z20": ("A;A_010_1", "12d@6l", {(0, 0): "sqrt(6)*(cos(p_1)+cos(p_10)+cos(p_11)+cos(p_12)+cos(p_2)+cos(p_3))/6", (1, 1): "-sqrt(6)*(cos(p_4)+cos(p_5)+cos(p_6)+cos(p_7)+cos(p_8)+cos(p_9))/6"}),
        "z21": ("A;A_010_1", "12d@6l", {(0, 0): "sqrt(6)*(cos(p_1)-cos(p_10)-cos(p_11)-cos(p_12)+cos(p_2)+cos(p_3))/6", (1, 1): "sqrt(6)*(-cos(p_4)-cos(p_5)-cos(p_6)+cos(p_7)+cos(p_8)+cos(p_9))/6"}),
        "z22": ("A;A_002_1", "6b@6l", {(0, 0): "sqrt(2)*(-cos(p_2)+cos(p_3))/2", (1, 1): "sqrt(2)*(cos(p_5)-cos(p_6))/2"}),
        "z23": ("A;A_002_1", "6b@6l", {(0, 0): "sqrt(6)*(2*cos(p_1)-cos(p_2)-cos(p_3))/6", (1, 1): "sqrt(6)*(-2*cos(p_4)+cos(p_5)+cos(p_6))/6"}),
        "z24": ("A;A_005_1", "6a@6l", {(0, 0): "sqrt(2)*(-cos(p_2)+cos(p_3))/2", (1, 1): "sqrt(2)*(cos(p_5)-cos(p_6))/2"}),
        "z25": ("A;A_005_1", "6a@6l", {(0, 0): "sqrt(6)*(2*cos(p_1)-cos(p_2)-cos(p_3))/6", (1, 1): "sqrt(6)*(-2*cos(p_4)+cos(p_5)+cos(p_6))/6"}),
        "z26": ("A;A_006_1", "6c@2c", {(0, 0): "sqrt(2)*(-cos(p_2)+cos(p_3))/2", (1, 1): "sqrt(2)*(cos(p_5)-cos(p_6))/2"}),
        "z27": ("A;A_006_1", "6c@2c", {(0, 0): "sqrt(6)*(2*cos(p_1)-cos(p_2)-cos(p_3))/6", (1, 1): "sqrt(6)*(-2*cos(p_4)+cos(p_5)+cos(p_6))/6"}),
        "z28": ("A;A_010_1", "12d@6l", {(0, 0): "sqrt(21)*(5*cos(p_1)+cos(p_10)-5*cos(p_11)+4*cos(p_12)-4*cos(p_2)-cos(p_3))/42", (1, 1): "sqrt(21)*(-5*cos(p_4)+4*cos(p_5)+cos(p_6)-cos(p_7)+5*cos(p_8)-4*cos(p_9))/42"}),
        "z29": ("A;A_010_1", "12d@6l", {(0, 0): "sqrt(7)*(cos(p_1)-3*cos(p_10)+cos(p_11)+2*cos(p_12)+2*cos(p_2)-3*cos(p_3))/14", (1, 1): "sqrt(7)*(-cos(p_4)-2*cos(p_5)+3*cos(p_6)+3*cos(p_7)-cos(p_8)-2*cos(p_9))/14"}),
        "z30": ("A;A_010_1", "12d@6l", {(0, 0): "sqrt(7)*(cos(p_1)+3*cos(p_10)-cos(p_11)-2*cos(p_12)+2*cos(p_2)-3*cos(p_3))/14", (1, 1): "sqrt(7)*(-cos(p_4)-2*cos(p_5)+3*cos(p_6)-3*cos(p_7)+cos(p_8)+2*cos(p_9))/14"}),
        "z31": ("A;A_010_1", "12d@6l", {(0, 0): "sqrt(21)*(-5*cos(p_1)+cos(p_10)-5*cos(p_11)+4*cos(p_12)+4*cos(p_2)+cos(p_3))/42", (1, 1): "sqrt(21)*(5*cos(p_4)-4*cos(p_5)-cos(p_6)-cos(p_7)+5*cos(p_8)-4*cos(p_9))/42"}),
        "z32": ("A;A_001_1", "3a@3f", {(1, 0): "sqrt(3)*(2*I*sin(p_1)-I*sin(p_2)-I*sin(p_3)+2*cos(p_1)-cos(p_2)-cos(p_3))/6", (0, 1): "sqrt(3)*(-2*I*sin(p_1)+I*sin(p_2)+I*sin(p_3)+2*cos(p_1)-cos(p_2)-cos(p_3))/6"}),
        "z33": ("A;A_001_1", "3a@3f", {(1, 0): "(I*sin(p_2)-I*sin(p_3)+cos(p_2)-cos(p_3))/2", (0, 1): "(-I*sin(p_2)+I*sin(p_3)+cos(p_2)-cos(p_3))/2"}),
        "z34": ("A;A_002_1", "6b@6l", {(0, 0): "sqrt(6)*(2*cos(p_1)-cos(p_2)-cos(p_3))/6", (1, 1): "sqrt(6)*(2*cos(p_4)-cos(p_5)-cos(p_6))/6"}),
        "z35": ("A;A_002_1", "6b@6l", {(0, 0): "sqrt(2)*(cos(p_2)-cos(p_3))/2", (1, 1): "sqrt(2)*(cos(p_5)-cos(p_6))/2"}),
        "z36": ("A;A_003_1", "3b@1a", {(1, 0): "sqrt(3)*(2*I*sin(p_1)-I*sin(p_2)-I*sin(p_3)+2*cos(p_1)-cos(p_2)-cos(p_3))/6", (0, 1): "sqrt(3)*(-2*I*sin(p_1)+I*sin(p_2)+I*sin(p_3)+2*cos(p_1)-cos(p_2)-cos(p_3))/6"}),
        "z37": ("A;A_003_1", "3b@1a", {(1, 0): "(I*sin(p_2)-I*sin(p_3)+cos(p_2)-cos(p_3))/2", (0, 1): "(-I*sin(p_2)+I*sin(p_3)+cos(p_2)-cos(p_3))/2"}),
        "z38": (
            "A;A_004_1",
            "6d@3f",
            {
                (1, 0): "sqrt(6)*(11*I*sin(p_1)+2*I*sin(p_2)-13*I*sin(p_3)+13*I*sin(p_4)-11*I*sin(p_5)-2*I*sin(p_6)+11*cos(p_1)+2*cos(p_2)-13*cos(p_3)-13*cos(p_4)+11*cos(p_5)+2*cos(p_6))/84",
                (0, 1): "sqrt(6)*(-11*I*sin(p_1)-2*I*sin(p_2)+13*I*sin(p_3)-13*I*sin(p_4)+11*I*sin(p_5)+2*I*sin(p_6)+11*cos(p_1)+2*cos(p_2)-13*cos(p_3)-13*cos(p_4)+11*cos(p_5)+2*cos(p_6))/84",
            },
        ),
        "z39": (
            "A;A_004_1",
            "6d@3f",
            {
                (1, 0): "sqrt(2)*(-5*I*sin(p_1)+8*I*sin(p_2)-3*I*sin(p_3)-3*I*sin(p_4)-5*I*sin(p_5)+8*I*sin(p_6)-5*cos(p_1)+8*cos(p_2)-3*cos(p_3)+3*cos(p_4)+5*cos(p_5)-8*cos(p_6))/28",
                (0, 1): "sqrt(2)*(5*I*sin(p_1)-8*I*sin(p_2)+3*I*sin(p_3)+3*I*sin(p_4)+5*I*sin(p_5)-8*I*sin(p_6)-5*cos(p_1)+8*cos(p_2)-3*cos(p_3)+3*cos(p_4)+5*cos(p_5)-8*cos(p_6))/28",
            },
        ),
        "z40": (
            "A;A_004_1",
            "6d@3f",
            {
                (1, 0): "sqrt(2)*(5*I*sin(p_1)-8*I*sin(p_2)+3*I*sin(p_3)-3*I*sin(p_4)-5*I*sin(p_5)+8*I*sin(p_6)+5*cos(p_1)-8*cos(p_2)+3*cos(p_3)+3*cos(p_4)+5*cos(p_5)-8*cos(p_6))/28",
                (0, 1): "sqrt(2)*(-5*I*sin(p_1)+8*I*sin(p_2)-3*I*sin(p_3)+3*I*sin(p_4)+5*I*sin(p_5)-8*I*sin(p_6)+5*cos(p_1)-8*cos(p_2)+3*cos(p_3)+3*cos(p_4)+5*cos(p_5)-8*cos(p_6))/28",
            },
        ),
        "z41": (
            "A;A_004_1",
            "6d@3f",
            {
                (1, 0): "sqrt(6)*(11*I*sin(p_1)+2*I*sin(p_2)-13*I*sin(p_3)-13*I*sin(p_4)+11*I*sin(p_5)+2*I*sin(p_6)+11*cos(p_1)+2*cos(p_2)-13*cos(p_3)+13*cos(p_4)-11*cos(p_5)-2*cos(p_6))/84",
                (0, 1): "sqrt(6)*(-11*I*sin(p_1)-2*I*sin(p_2)+13*I*sin(p_3)+13*I*sin(p_4)-11*I*sin(p_5)-2*I*sin(p_6)+11*cos(p_1)+2*cos(p_2)-13*cos(p_3)+13*cos(p_4)-11*cos(p_5)-2*cos(p_6))/84",
            },
        ),
        "z42": ("A;A_005_1", "6a@6l", {(0, 0): "sqrt(6)*(2*cos(p_1)-cos(p_2)-cos(p_3))/6", (1, 1): "sqrt(6)*(2*cos(p_4)-cos(p_5)-cos(p_6))/6"}),
        "z43": ("A;A_005_1", "6a@6l", {(0, 0): "sqrt(2)*(cos(p_2)-cos(p_3))/2", (1, 1): "sqrt(2)*(cos(p_5)-cos(p_6))/2"}),
        "z44": ("A;A_006_1", "6c@2c", {(0, 0): "sqrt(6)*(2*cos(p_1)-cos(p_2)-cos(p_3))/6", (1, 1): "sqrt(6)*(2*cos(p_4)-cos(p_5)-cos(p_6))/6"}),
        "z45": ("A;A_006_1", "6c@2c", {(0, 0): "sqrt(2)*(cos(p_2)-cos(p_3))/2", (1, 1): "sqrt(2)*(cos(p_5)-cos(p_6))/2"}),
        "z46": (
            "A;A_007_1",
            "6d@3f",
            {
                (1, 0): "sqrt(6)*(11*I*sin(p_1)+2*I*sin(p_2)-13*I*sin(p_3)+13*I*sin(p_4)-11*I*sin(p_5)-2*I*sin(p_6)+11*cos(p_1)+2*cos(p_2)-13*cos(p_3)-13*cos(p_4)+11*cos(p_5)+2*cos(p_6))/84",
                (0, 1): "sqrt(6)*(-11*I*sin(p_1)-2*I*sin(p_2)+13*I*sin(p_3)-13*I*sin(p_4)+11*I*sin(p_5)+2*I*sin(p_6)+11*cos(p_1)+2*cos(p_2)-13*cos(p_3)-13*cos(p_4)+11*cos(p_5)+2*cos(p_6))/84",
            },
        ),
        "z47": (
            "A;A_007_1",
            "6d@3f",
            {
                (1, 0): "sqrt(2)*(-5*I*sin(p_1)+8*I*sin(p_2)-3*I*sin(p_3)-3*I*sin(p_4)-5*I*sin(p_5)+8*I*sin(p_6)-5*cos(p_1)+8*cos(p_2)-3*cos(p_3)+3*cos(p_4)+5*cos(p_5)-8*cos(p_6))/28",
                (0, 1): "sqrt(2)*(5*I*sin(p_1)-8*I*sin(p_2)+3*I*sin(p_3)+3*I*sin(p_4)+5*I*sin(p_5)-8*I*sin(p_6)-5*cos(p_1)+8*cos(p_2)-3*cos(p_3)+3*cos(p_4)+5*cos(p_5)-8*cos(p_6))/28",
            },
        ),
        "z48": (
            "A;A_007_1",
            "6d@3f",
            {
                (1, 0): "sqrt(2)*(5*I*sin(p_1)-8*I*sin(p_2)+3*I*sin(p_3)-3*I*sin(p_4)-5*I*sin(p_5)+8*I*sin(p_6)+5*cos(p_1)-8*cos(p_2)+3*cos(p_3)+3*cos(p_4)+5*cos(p_5)-8*cos(p_6))/28",
                (0, 1): "sqrt(2)*(-5*I*sin(p_1)+8*I*sin(p_2)-3*I*sin(p_3)+3*I*sin(p_4)+5*I*sin(p_5)-8*I*sin(p_6)+5*cos(p_1)-8*cos(p_2)+3*cos(p_3)+3*cos(p_4)+5*cos(p_5)-8*cos(p_6))/28",
            },
        ),
        "z49": (
            "A;A_007_1",
            "6d@3f",
            {
                (1, 0): "sqrt(6)*(11*I*sin(p_1)+2*I*sin(p_2)-13*I*sin(p_3)-13*I*sin(p_4)+11*I*sin(p_5)+2*I*sin(p_6)+11*cos(p_1)+2*cos(p_2)-13*cos(p_3)+13*cos(p_4)-11*cos(p_5)-2*cos(p_6))/84",
                (0, 1): "sqrt(6)*(-11*I*sin(p_1)-2*I*sin(p_2)+13*I*sin(p_3)+13*I*sin(p_4)-11*I*sin(p_5)-2*I*sin(p_6)+11*cos(p_1)+2*cos(p_2)-13*cos(p_3)+13*cos(p_4)-11*cos(p_5)-2*cos(p_6))/84",
            },
        ),
        "z50": ("A;A_008_1", "3b@1a", {(1, 0): "sqrt(3)*(2*I*sin(p_1)-I*sin(p_2)-I*sin(p_3)+2*cos(p_1)-cos(p_2)-cos(p_3))/6", (0, 1): "sqrt(3)*(-2*I*sin(p_1)+I*sin(p_2)+I*sin(p_3)+2*cos(p_1)-cos(p_2)-cos(p_3))/6"}),
        "z51": ("A;A_008_1", "3b@1a", {(1, 0): "(I*sin(p_2)-I*sin(p_3)+cos(p_2)-cos(p_3))/2", (0, 1): "(-I*sin(p_2)+I*sin(p_3)+cos(p_2)-cos(p_3))/2"}),
        "z52": (
            "A;A_009_1",
            "6d@3f",
            {
                (1, 0): "sqrt(6)*(11*I*sin(p_1)+2*I*sin(p_2)-13*I*sin(p_3)+13*I*sin(p_4)-11*I*sin(p_5)-2*I*sin(p_6)+11*cos(p_1)+2*cos(p_2)-13*cos(p_3)-13*cos(p_4)+11*cos(p_5)+2*cos(p_6))/84",
                (0, 1): "sqrt(6)*(-11*I*sin(p_1)-2*I*sin(p_2)+13*I*sin(p_3)-13*I*sin(p_4)+11*I*sin(p_5)+2*I*sin(p_6)+11*cos(p_1)+2*cos(p_2)-13*cos(p_3)-13*cos(p_4)+11*cos(p_5)+2*cos(p_6))/84",
            },
        ),
        "z53": (
            "A;A_009_1",
            "6d@3f",
            {
                (1, 0): "sqrt(2)*(-5*I*sin(p_1)+8*I*sin(p_2)-3*I*sin(p_3)-3*I*sin(p_4)-5*I*sin(p_5)+8*I*sin(p_6)-5*cos(p_1)+8*cos(p_2)-3*cos(p_3)+3*cos(p_4)+5*cos(p_5)-8*cos(p_6))/28",
                (0, 1): "sqrt(2)*(5*I*sin(p_1)-8*I*sin(p_2)+3*I*sin(p_3)+3*I*sin(p_4)+5*I*sin(p_5)-8*I*sin(p_6)-5*cos(p_1)+8*cos(p_2)-3*cos(p_3)+3*cos(p_4)+5*cos(p_5)-8*cos(p_6))/28",
            },
        ),
        "z54": (
            "A;A_009_1",
            "6d@3f",
            {
                (1, 0): "sqrt(2)*(5*I*sin(p_1)-8*I*sin(p_2)+3*I*sin(p_3)-3*I*sin(p_4)-5*I*sin(p_5)+8*I*sin(p_6)+5*cos(p_1)-8*cos(p_2)+3*cos(p_3)+3*cos(p_4)+5*cos(p_5)-8*cos(p_6))/28",
                (0, 1): "sqrt(2)*(-5*I*sin(p_1)+8*I*sin(p_2)-3*I*sin(p_3)+3*I*sin(p_4)+5*I*sin(p_5)-8*I*sin(p_6)+5*cos(p_1)-8*cos(p_2)+3*cos(p_3)+3*cos(p_4)+5*cos(p_5)-8*cos(p_6))/28",
            },
        ),
        "z55": (
            "A;A_009_1",
            "6d@3f",
            {
                (1, 0): "sqrt(6)*(11*I*sin(p_1)+2*I*sin(p_2)-13*I*sin(p_3)-13*I*sin(p_4)+11*I*sin(p_5)+2*I*sin(p_6)+11*cos(p_1)+2*cos(p_2)-13*cos(p_3)+13*cos(p_4)-11*cos(p_5)-2*cos(p_6))/84",
                (0, 1): "sqrt(6)*(-11*I*sin(p_1)-2*I*sin(p_2)+13*I*sin(p_3)+13*I*sin(p_4)-11*I*sin(p_5)-2*I*sin(p_6)+11*cos(p_1)+2*cos(p_2)-13*cos(p_3)+13*cos(p_4)-11*cos(p_5)-2*cos(p_6))/84",
            },
        ),
        "z56": ("A;A_010_1", "12d@6l", {(0, 0): "sqrt(3)*(11*cos(p_1)-13*cos(p_10)+11*cos(p_11)+2*cos(p_12)+2*cos(p_2)-13*cos(p_3))/42", (1, 1): "sqrt(3)*(11*cos(p_4)+2*cos(p_5)-13*cos(p_6)-13*cos(p_7)+11*cos(p_8)+2*cos(p_9))/42"}),
        "z57": ("A;A_010_1", "12d@6l", {(0, 0): "(-5*cos(p_1)+3*cos(p_10)+5*cos(p_11)-8*cos(p_12)+8*cos(p_2)-3*cos(p_3))/14", (1, 1): "(-5*cos(p_4)+8*cos(p_5)-3*cos(p_6)+3*cos(p_7)+5*cos(p_8)-8*cos(p_9))/14"}),
        "z58": ("A;A_010_1", "12d@6l", {(0, 0): "(5*cos(p_1)+3*cos(p_10)+5*cos(p_11)-8*cos(p_12)-8*cos(p_2)+3*cos(p_3))/14", (1, 1): "(5*cos(p_4)-8*cos(p_5)+3*cos(p_6)+3*cos(p_7)+5*cos(p_8)-8*cos(p_9))/14"}),
        "z59": ("A;A_010_1", "12d@6l", {(0, 0): "sqrt(3)*(11*cos(p_1)+13*cos(p_10)-11*cos(p_11)-2*cos(p_12)+2*cos(p_2)-13*cos(p_3))/42", (1, 1): "sqrt(3)*(11*cos(p_4)+2*cos(p_5)-13*cos(p_6)+13*cos(p_7)-11*cos(p_8)-2*cos(p_9))/42"}),
    },
}
