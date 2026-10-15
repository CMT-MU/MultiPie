# default control.
default_control = {
    "mode": "samb",  # samb (SAMB-based only), wannier (wannier-based only), symcw (SAMB and wannier based).
    "grid": (50, 50, 50),  # k-grid size (b1, b2, b3).
    "samb": {
        "model": None,  # model name for .pkl (optional in wannier mode).
        "select": {  # SAMB select condition, where S=site name, R=orbital rank, N=neighbor bond.
            # "site": [("A", [1, 2])],  # S or (S, [R]) in list.
            # "bond": [[0, 1, 2]],  # S1;S2, N, [N], (S1;S2, N), (S1;S2, [N]), (S1;S2, R1;R2), (S1;S2, R1;R2, N) or (S1;S2, R1;R2, [N]) in list.
            # "X": [],  # SAMB type Q/G/T/M in list. empty []=all.
            # "l": [],  # SAMB rank, 0, 1, ..., 11 in list. empty []=all.
            # "Gamma": ["A1g"],  # SAMB irreps. in list. "IR"=identity irrep., empty []=all.
            # "s": [],  # SAMB internal rank, 0,1 in list. empty []=all
        },
        "parameter": {  # SAMB with finite weight (float, or str of sympy const., e.g., "sqrt(2)"), or filename of z_j dict (.py file relative to model_name/info/).
            # "z1": 1.0,
        },
        "samb_figure": False,  # save SAMB QtDraw files ?
        "k_multipole": False,  # compute momentum multipole ?
        "NG_sum_rule": False,  # Nambu-Goldstone sum rule ?
    },
    "wannier": {  # Closest Wannier (CW) or SymWannier setting.
        "dir": "wannier",  # directory for wannier input, i.e., "./seedname/wannier".
        "seedname": None,  # read +(.win, .nnkp, _hr.dat, _wsvec.dat if exists, ...).
        "ket_wannier": [],  # correspondence between Wannier and MultiPie basis, [MultiPie ket name, "orbital@site(sublattice)"] in Wannier order. If empty, automatically determined.
    },
    "output": {  # physical quantity setting.
        "dir": "output",  #  output directory.
        "fourier": {"tb_gauge": True},  # TB gauge ?
        "dos": False,  # DOS: False (none), "auto" (partial DOS for each site cluster and orbital, spin summed), or [[ket name]] (sum of partial DOS of the kets in each list). total DOS is always included.
        "dos_e_num": 1000,  # number of energy points of DOS (range includes all eigenvalues).
        "dos_width": None,  # Gaussian broadening of DOS. if None, 1% of the band width.
        "dos_method": "gaussian",  # DOS method, "gaussian" or "tetrahedron" (optimized tetrahedron method by libtetrabz).
        "dos_spin": False,  # resolve spin, u/d, of spinful basis in "auto" DOS ?
        "n_electron": None,  # number of occupied states per primitive cell (one per band) to determine fermi energy. if None, 0 (samb mode) or that of seedname.win.
        "save_k_path": False,  # save Ek, Uk along k path in .npz ?
        "save_k_full": False,  # save Ek, Uk on full k grid in .npz ?
        "callback_k_path": None,  # file (.py, relative to topdir) defining callback_k_path(Ek, Uk, info), called with eigen system along k path.
        "callback_k_full": None,  # file (.py, relative to topdir) defining callback_k_full(Ek, Uk, info), called with eigen system on full k grid.
        "dispersion": {  # dispersion info.
            "power": None,  # plot as Ek^power. if None, power=1.
            "k_path": "",  # symmetry line, separated by "-". disconnected points by "|". if empty str, default values are used for k_path and k_point.
            "k_point": {"Γ": "[0,0,0]"},  # k-point (fractional, primitive cell of model, or of seedname.win without model).
            "local": [],  # local quantity, "Sx/Sy/Sz/Lx/Ly/Lz/Qu/Qv/Qyz/Qxz/Qxy".
            "weight": False,  # orbital weight (fat band): False (none), "auto" (each site cluster and orbital, spin summed), or [[ket name]].
            "z": [],  # expectation value of Z_j. [Under Construction]
        },
    },
}
