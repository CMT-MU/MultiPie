# default model input.
default_model = {
    "model": "unknown",
    "group": "C1",
    "cell": {},
    "spinful": False,
    "site": {},
    "bond": [],
    "SAMB_select": {  # select combined SAMB.
        "X": ["Q", "G"],  # type, "Q/G/T/M", []=all.
        "l": [],  # rank, []=all.
        "Gamma": [],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
        "s": [],  # spin, 0/1, []=all.
    },
    "atomic_select": {  # select atomic SAMB.
        "X": [],  # type, "Q/G/T/M", []=all.
        "l": [],  # rank, []=all.
        "Gamma": [],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
        "s": [],  # spin, 0/1, []=all.
    },
    "site_select": {  # select site-cluster SAMB.
        "l": [],  # rank, []=all.
        "Gamma": [],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
    },
    "bond_select": {  # select bond-cluster SAMB.
        "X": [],  # type, "Q/T/M", []=all.
        "l": [],  # rank, []=all.
        "Gamma": [],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
    },
    "toroidal_priority": False,
    "max_neighbor": 10,  # for DFT fitting, e.g. =200.
    "search_cell_range": (-2, 3, -2, 3, -2, 3),  # (a1, a2, a3) range. for DFT fitting, e.g.=(-10, 10, -10, 10, -10, 10).
    "qtdraw": {
        "create": True,  # create QtDraw file ?
        "mode": "standard",  # mode, "standard/detail".
        "view": None,  # [a,b,c], if None, default=[6,5,1].
        #
        "rep_site": True,  # show representative site ?
        "rep_bond": True,  # show representative bond ?
        #
        "max_neighbor": 5,  # max. neighbor to draw.
        "cell_mode": None,  # "off/single/all", if None, single for SG, off for PG.
        "scale": None,  # default or magnification scale.
        "site_radius": 0.07,  # base site radius.
        "bond_width": 0.02,  # base bond width.
    },
    "pdf": {
        "create": True,  # create PDF file ?
        "common_samb": True,  # display common SAMB ?
        "harmonics": True,  # display harmonics ?
        "max_neighbor": 5,  # max. neighbor to write (None: all).
    },
}
