"""
sample input file for C3v molecule.
"""
C3v = {
    "model": "C3v",  # name of model.
    "group": "C3v-1",  # name of point group.
    "cell": {"c": 10},  # set large enough interlayer distance.
    #
    "site": {"A": ("[-1/6,-1/6,0]", "s"), "B": ("[-2/3,0,0]", "p")},  # positions of A and B sites and their orbitals.
    "bond": [("A", "A", 1), ("A", "B", 1)],  # nearest-neighbor A-A and B-B bonds.
    #
    "spinful": False,  # spinless.
}
