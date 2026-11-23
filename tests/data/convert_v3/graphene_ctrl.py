# control of graphene (SAMB of MultiPie 3.1.3).
graphene_ctrl = {
    "samb": {"model": "graphene", "parameter": {"z1": 0.1, "z2": -1.0, "z3": -0.2, "z5": 0.3, "z15": 0.05}},
    "output": {
        "callback_info": ["Sz", "z15"],
        "dispersion": {"z": ["z2", "z16"]},
        "z_expectation": ["z1", "z5", "z15"],
    },
}
