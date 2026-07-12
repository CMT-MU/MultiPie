# graphene_in.py
graphene_in = {
    "model": "graphene",  #  name of model.
    "group": 191,  # No. of space group.
    "cell": {"c": 4},  # set large enough interlayer distance.
    #
    "site": {"C": ("[1/3,2/3,0]", "pz")},  # positions of C site and its orbital.
    "bond": [("C", "C", [1, 2])],  # C-C bonds up to 2nd neighbors.
    #
    "spinful": False,  # spinless.
    #
    "SAMB_select": {  # select combined SAMB.
        "X": [],  # all Q/G/M/T.
    },
}
