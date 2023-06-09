"""
cartesian harmonics to cubic harmonics
    cartesian harmonics : [(coeff, cubic harmonics, idx)] (idx is the index of the _data_compatibility_relation["Oh"/"D6h"])
"""
_data_cartesian_to_cubic_harmonics = {
    # polar
    # rank 0
    "Q_{0}": [(1, "Q_{0}", 0)],
    # rank 1
    "Q_{x}": [(1, "Q_{x}", 14)],
    "Q_{y}": [(1, "Q_{y}", 15)],
    "Q_{z}": [(1, "Q_{z}", 16)],
    # rank 2
    "Q_{xx}": [(-1, "Q_{u}", 2), (1, "Q_{v}", 3)],
    "Q_{yy}": [(-1, "Q_{u}", 2), (-1, "Q_{v}", 3)],
    "Q_{zz}": [(2, "Q_{u}", 2)],
    "Q_{yz}": [(1, "Q_{yz}", 7)],
    "Q_{xz}": [(1, "Q_{zx}", 8)],
    "Q_{xy}": [(1, "Q_{xy}", 9)],
    # rank 3
    "Q_{xxx}": [(2, "Q_{x}^{\\alpha}", 14)],
    "Q_{yyy}": [(2, "Q_{y}^{\\alpha}", 15)],
    "Q_{zzz}": [(2, "Q_{z}^{\\alpha}", 16)],
    "Q_{xzz}": [(-1, "Q_{x}^{\\alpha}", 14), (-1, "Q_{x}^{\\beta}", 17)],
    "Q_{xxy}": [(-1, "Q_{y}^{\\alpha}", 15), (-1, "Q_{y}^{\\beta}", 18)],
    "Q_{yyz}": [(-1, "Q_{z}^{\\alpha}", 16), (-1, "Q_{z}^{\\beta}", 19)],
    "Q_{xyy}": [(-1, "Q_{x}^{\\alpha}", 14), (1, "Q_{x}^{\\beta}", 17)],
    "Q_{yzz}": [(-1, "Q_{y}^{\\alpha}", 15), (1, "Q_{y}^{\\beta}", 18)],
    "Q_{xxz}": [(-1, "Q_{z}^{\\alpha}", 16), (1, "Q_{z}^{\\beta}", 19)],
    "Q_{xyz}": [(1, "Q_{xyz}", 11)],
    # rank 4
    "Q_{xxxx}": [(2, "Q_{4}", 0), (-1, "Q_{4u}", 2), (1, "Q_{4v}", 3)],
    "Q_{yyyy}": [(2, "Q_{4}", 0), (-1, "Q_{4u}", 2), (-1, "Q_{4v}", 3)],
    "Q_{zzzz}": [(2, "Q_{4}", 0), (2, "Q_{4u}", 2)],
    "Q_{yyzz}": [(-1, "Q_{4}", 0), (-1, "Q_{4u}", 2), (1, "Q_{4v}", 3)],
    "Q_{xxzz}": [(-1, "Q_{4}", 0), (-1, "Q_{4u}", 2), (-1, "Q_{4v}", 3)],
    "Q_{xxyy}": [(-1, "Q_{4}", 0), (2, "Q_{4u}", 2)],
    "Q_{xxyz}": [(2, "Q_{4x}^{\\beta}", 7)],
    "Q_{xyyz}": [(2, "Q_{4y}^{\\beta}", 8)],
    "Q_{xyzz}": [(2, "Q_{4z}^{\\beta}", 9)],
    "Q_{yyyz}": [(1, "Q_{4x}^{\\alpha}", 4), (-1, "Q_{4x}^{\\beta}", 7)],
    "Q_{xzzz}": [(1, "Q_{4y}^{\\alpha}", 5), (-1, "Q_{4y}^{\\beta}", 8)],
    "Q_{xxxy}": [(1, "Q_{4z}^{\\alpha}", 6), (-1, "Q_{4z}^{\\beta}", 9)],
    "Q_{yzzz}": [(-1, "Q_{4x}^{\\alpha}", 4), (-1, "Q_{4x}^{\\beta}", 7)],
    "Q_{xxxz}": [(-1, "Q_{4y}^{\\alpha}", 5), (-1, "Q_{4y}^{\\beta}", 8)],
    "Q_{xyyy}": [(-1, "Q_{4z}^{\\alpha}", 6), (-1, "Q_{4z}^{\\beta}", 9)],
    # axial
    # rank 0
    "G_{0}": [(1, "G_{0}", 10)],
    # rank 1
    "G_{x}": [(1, "G_{x}", 4)],
    "G_{y}": [(1, "G_{y}", 5)],
    "G_{z}": [(1, "G_{z}", 6)],
    # rank 2
    "G_{xx}": [(-1, "G_{u}", 12), (1, "G_{v}", 13)],
    "G_{yy}": [(-1, "G_{u}", 12), (-1, "G_{v}", 13)],
    "G_{zz}": [(2, "G_{u}", 12)],
    "G_{yz}": [(1, "G_{yz}", 17)],
    "G_{xz}": [(1, "G_{zx}", 18)],
    "G_{xy}": [(1, "G_{xy}", 19)],
    # rank 3
    "G_{xxx}": [(2, "G_{x}^{\\alpha}", 4)],
    "G_{yyy}": [(2, "G_{y}^{\\alpha}", 5)],
    "G_{zzz}": [(2, "G_{z}^{\\alpha}", 6)],
    "G_{xzz}": [(-1, "G_{x}^{\\alpha}", 4), (-1, "G_{x}^{\\beta}", 7)],
    "G_{xxy}": [(-1, "G_{y}^{\\alpha}", 5), (-1, "G_{y}^{\\beta}", 8)],
    "G_{yyz}": [(-1, "G_{z}^{\\alpha}", 6), (-1, "G_{z}^{\\beta}", 9)],
    "G_{xyy}": [(-1, "G_{x}^{\\alpha}", 4), (1, "G_{x}^{\\beta}", 7)],
    "G_{yzz}": [(-1, "G_{y}^{\\alpha}", 5), (1, "G_{y}^{\\beta}", 8)],
    "G_{xxz}": [(-1, "G_{z}^{\\alpha}", 6), (1, "G_{z}^{\\beta}", 9)],
    "G_{xyz}": [(1, "G_{xyz}", 1)],
    # rank 4
    "G_{xxxx}": [(2, "G_{4}", 10), (-1, "G_{4u}", 12), (1, "G_{4v}", 13)],
    "G_{yyyy}": [(2, "G_{4}", 10), (-1, "G_{4u}", 12), (-1, "G_{4v}", 13)],
    "G_{zzzz}": [(2, "G_{4}", 10), (2, "G_{4u}", 12)],
    "G_{yyzz}": [(-1, "G_{4}", 10), (-1, "G_{4u}", 12), (1, "G_{4v}", 13)],
    "G_{xxzz}": [(-1, "G_{4}", 10), (-1, "G_{4u}", 12), (-1, "G_{4v}", 13)],
    "G_{xxyy}": [(-1, "G_{4}", 10), (2, "G_{4u}", 12)],
    "G_{xxyz}": [(2, "G_{4x}^{\\beta}", 17)],
    "G_{xyyz}": [(2, "G_{4y}^{\\beta}", 18)],
    "G_{xyzz}": [(2, "G_{4z}^{\\beta}", 19)],
    "G_{yyyz}": [(1, "G_{4x}^{\\alpha}", 14), (-1, "G_{4x}^{\\beta}", 17)],
    "G_{xzzz}": [(1, "G_{4y}^{\\alpha}", 15), (-1, "G_{4y}^{\\beta}", 18)],
    "G_{xxxy}": [(1, "G_{4z}^{\\alpha}", 16), (-1, "G_{4z}^{\\beta}", 19)],
    "G_{yzzz}": [(-1, "G_{4x}^{\\alpha}", 14), (-1, "G_{4x}^{\\beta}", 17)],
    "G_{xxxz}": [(-1, "G_{4y}^{\\alpha}", 15), (-1, "G_{4y}^{\\beta}", 18)],
    "G_{xyyy}": [(-1, "G_{4z}^{\\alpha}", 16), (-1, "G_{4z}^{\\beta}", 19)],
}

"""
cartesian harmonics to hexagonal harmonics
    cartesian harmonics : [(coeff, hexagonal harmonics, irrep)]
"""
_data_cartesian_to_hexagonal_harmonics = {
    # polar
    # rank 0
    "Q_{0}": [(1, "Q_{0}", 0)],
    # rank 1
    "Q_{x}": [(1, "Q_{x}", 12)],
    "Q_{y}": [(1, "Q_{y}", 13)],
    "Q_{z}": [(1, "Q_{z}", 9)],
    # rank 2
    "Q_{xx}": [(-1, "Q_{u}", 0), (1, "Q_{v}", 6)],
    "Q_{yy}": [(-1, "Q_{u}", 0), (-1, "Q_{v}", 6)],
    "Q_{zz}": [(2, "Q_{u}", 0)],
    "Q_{yz}": [(1, "Q_{yz}", 5)],
    "Q_{xz}": [(1, "Q_{zx}", 4)],
    "Q_{xy}": [(-1, "Q_{xy}", 7)],
    # rank 3
    "Q_{xxx}": [(-3, "Q_{3x}^{\\alpha}", 12), (1, "Q_{3}^{\\gamma}", 11)],
    "Q_{xxy}": [(-1, "Q_{3y}^{\\alpha}", 13), (1, "Q_{3}^{\\beta}", 10)],
    "Q_{xyy}": [(-1, "Q_{3x}^{\\alpha}", 12), (-1, "Q_{3}^{\\gamma}", 11)],
    "Q_{yyy}": [(-3, "Q_{3y}^{\\alpha}", 13), (-1, "Q_{3}^{\\beta}", 10)],
    "Q_{xxz}": [(-3, "Q_{3}^{\\alpha}", 9), (1, "Q_{3x}^{\\beta}", 14)],
    "Q_{yyz}": [(-1, "Q_{3}^{\\alpha}", 9), (-1, "Q_{3x}^{\\beta}", 14)],
    "Q_{xzz}": [(4, "Q_{3x}^{\\alpha}", 12)],
    "Q_{yzz}": [(4, "Q_{3y}^{\\alpha}", 13)],
    "Q_{zzz}": [(2, "Q_{3}^{\\alpha}", 9)],
    "Q_{xyz}": [(-1, "Q_{3y}^{\\beta}", 15)],
    # rank 4
    "Q_{xxxx}": [(3, "Q_{4}", 0), (1, "Q_{4x}^{\\beta}", 6), (-2, "Q_{4x}^{\\gamma}", 6)],
    "Q_{yyyy}": [(3, "Q_{4}", 0), (1, "Q_{4x}^{\\beta}", 6), (2, "Q_{4x}^{\\gamma}", 6)],
    "Q_{zzzz}": [(8, "Q_{4}", 0)],
    "Q_{xxxy}": [(1, "Q_{4y}^{\\beta}", 7), (1, "Q_{4y}^{\\gamma}", 7)],
    "Q_{xyyy}": [(-1, "Q_{4y}^{\\beta}", 7), (1, "Q_{4y}^{\\gamma}", 7)],
    "Q_{xxxz}": [(1, "Q_{4}^{\\alpha}", 2), (-3, "Q_{4x}^{\\alpha}", 4)],
    "Q_{xxyy}": [(1, "Q_{4}", 0), (-1, "Q_{4x}^{\\beta}", 6)],
    "Q_{xxyz}": [(-1, "Q_{4y}^{\\alpha}", 5), (1, "Q_{4}^{\\beta}", 3)],
    "Q_{xxzz}": [(-4, "Q_{4}", 0), (2, "Q_{4x}^{\\gamma}", 6)],
    "Q_{xyyz}": [(-1, "Q_{4}^{\\alpha}", 2), (-1, "Q_{4x}^{\\alpha}", 4)],
    "Q_{xyzz}": [(-2, "Q_{4y}^{\\gamma}", 7)],
    "Q_{xzzz}": [(4, "Q_{4x}^{\\alpha}", 4)],
    "Q_{yyyz}": [(-3, "Q_{4y}^{\\alpha}", 5), (-1, "Q_{4}^{\\beta}", 3)],
    "Q_{yyzz}": [(-4, "Q_{4}", 0), (-2, "Q_{4x}^{\\gamma}", 6)],
    "Q_{yzzz}": [(4, "Q_{4y}^{\\alpha}", 5)],
    # axial
    # rank 0
    "G_{0}": [(1, "G_{0}", 8)],
    # rank 1
    "G_{x}": [(1, "G_{x}", 4)],
    "G_{y}": [(1, "G_{y}", 5)],
    "G_{z}": [(1, "G_{z}", 1)],
    # rank 2
    "G_{xx}": [(-1, "G_{u}", 8), (1, "G_{v}", 14)],
    "G_{yy}": [(-1, "G_{u}", 8), (-1, "G_{v}", 14)],
    "G_{zz}": [(2, "G_{u}", 8)],
    "G_{yz}": [(1, "G_{yz}", 13)],
    "G_{xz}": [(1, "G_{zx}", 12)],
    "G_{xy}": [(-1, "G_{xy}", 15)],
    # rank 3
    "G_{xxx}": [(-3, "G_{3x}^{\\alpha}", 4), (1, "G_{3}^{\\gamma}", 3)],
    "G_{xxy}": [(-1, "G_{3y}^{\\alpha}", 5), (1, "G_{3}^{\\beta}", 2)],
    "G_{xyy}": [(-1, "G_{3x}^{\\alpha}", 4), (-1, "G_{3}^{\\gamma}", 3)],
    "G_{yyy}": [(-3, "G_{3y}^{\\alpha}", 5), (-1, "G_{3}^{\\beta}", 2)],
    "G_{xxz}": [(-3, "G_{3}^{\\alpha}", 1), (1, "G_{3x}^{\\beta}", 6)],
    "G_{yyz}": [(-1, "G_{3}^{\\alpha}", 1), (-1, "G_{3x}^{\\beta}", 6)],
    "G_{xzz}": [(4, "G_{3x}^{\\alpha}", 4)],
    "G_{yzz}": [(4, "G_{3y}^{\\alpha}", 5)],
    "G_{zzz}": [(2, "G_{3}^{\\alpha}", 1)],
    "G_{xyz}": [(-1, "G_{3y}^{\\beta}", 7)],
    # rank 4
    "G_{xxxx}": [(3, "G_{4}", 8), (1, "G_{4x}^{\\beta}", 14), (-2, "G_{4x}^{\\gamma}", 14)],
    "G_{yyyy}": [(3, "G_{4}", 8), (1, "G_{4x}^{\\beta}", 14), (2, "G_{4x}^{\\gamma}", 14)],
    "G_{zzzz}": [(8, "G_{4}", 8)],
    "G_{xxxy}": [(1, "G_{4y}^{\\beta}", 15), (1, "G_{4y}^{\\gamma}", 15)],
    "G_{xyyy}": [(-1, "G_{4y}^{\\beta}", 15), (1, "G_{4y}^{\\gamma}", 15)],
    "G_{xxxz}": [(1, "G_{4}^{\\alpha}", 10), (-3, "G_{4x}^{\\alpha}", 12)],
    "G_{xxyy}": [(1, "G_{4}", 8), (-1, "G_{4x}^{\\beta}", 14)],
    "G_{xxyz}": [(-1, "G_{4y}^{\\alpha}", 13), (1, "G_{4}^{\\beta}", 11)],
    "G_{xxzz}": [(-4, "G_{4}", 8), (2, "G_{4x}^{\\gamma}", 14)],
    "G_{xyyz}": [(-1, "G_{4}^{\\alpha}", 10), (-1, "G_{4x}^{\\alpha}", 12)],
    "G_{xyzz}": [(-2, "G_{4y}^{\\gamma}", 15)],
    "G_{xzzz}": [(4, "G_{4x}^{\\alpha}", 12)],
    "G_{yyyz}": [(-3, "G_{4y}^{\\alpha}", 13), (-1, "G_{4}^{\\beta}", 11)],
    "G_{yyzz}": [(-4, "G_{4}", 8), (-2, "G_{4x}^{\\gamma}", 14)],
    "G_{yzzz}": [(4, "G_{4y}^{\\alpha}", 13)],
}
