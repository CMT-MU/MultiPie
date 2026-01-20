"""
Default setting for MaterialModel.
"""

from multipie import __version__


# ==================================================
# default model input.
_default_model = {
    "model": "unknown",
    "group": "C1",
    "cell": {},
    "spinful": False,
    "site": {},
    "bond": [],
    "SAMB_select": {  # select combined SAMB.
        "X": ["Q", "G"],  # type, "Q/G/M/T", []=all.
        "l": [],  # rank, []=all.
        "Gamma": [],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
        "s": [],  # spin, 0/1, []=all.
    },
    "atomic_select": {  # select atomic SAMB.
        "X": [],  # type, "Q/G/M/T", []=all.
        "l": [],  # rank, []=all.
        "Gamma": [],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
        "s": [],  # spin, 0/1, []=all.
    },
    "site_select": {  # select site-cluster SAMB.
        "l": [],  # rank, []=all.
        "Gamma": [],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
    },
    "bond_select": {  # select bond-cluster SAMB.
        "X": [],  # type, "Q/G/M/T", []=all.
        "l": [],  # rank, []=all.
        "Gamma": [],  # "IR"=identity, []=all, "..."/["...","..."]=specified irreps.
    },
    "toroidal_priority": False,
    "max_neighbor": 10,  # for DFT fitting, e.g. =200.
    "search_cell_range": (-2, 3, -2, 3, -2, 3),  # (a1, a2, a3) range. for DFT fitting, e.g.=(-10, 10, -10, 10, -10, 10).
    "version": __version__,
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

# ==================================================
# site property for QtDraw.
_site_property = [  # (color, size, opacity)
    ("darkseagreen", 1.0, 1.0),  # 1st cluster
    ("lightblue", 0.9, 1.0),  # 2nd cluster
    ("sandybrown", 0.8, 1.0),  # 3rd cluster
    ("gold", 0.7, 1.0),  # 4th cluster
    ("darkkhaki", 0.6, 1.0),  # 5th cluster
    ("skyblue", 0.6, 1.0),  # 6th cluster
    ("thistle", 0.6, 1.0),  # 7th cluster
    ("darkgrey", 0.6, 1.0),  # 8th cluster
    ("burlywood", 0.6, 1.0),  # 9th cluster
    ("ghostwhite", 0.6, 1.0),  # other clusters
]

# ==================================================
# bond property for QtDraw.
_bond_property = [  # ((tail-color, head-color), width, opacity)
    (("snow", "silver"), 1.0, 1.0),  # 1st cluster
    (("lightcyan", "lightsteelblue"), 1.0, 1.0),  # 2nd cluster
    (("antiquewhite", "burlywood"), 1.0, 1.0),  # 3rd cluster
    (("palegoldenrod", "darkseagreen"), 1.0, 1.0),  # 4th cluster
    (("mistyrose", "lightpink"), 1.0, 1.0),  # 5th cluster
    (("aliceblue", "lightblue"), 1.0, 1.0),  # 6th cluster
    (("wheat", "sandybrown"), 1.0, 1.0),  # 7th cluster
    (("seashell", "thistle"), 1.0, 1.0),  # 8th cluster
    (("cornsilk", "peachpuff"), 1.0, 1.0),  # 9th cluster
    (("whitesmoke", "darkkhaki"), 1.0, 1.0),  # other clusters
]
