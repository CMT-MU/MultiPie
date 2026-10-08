"""
Regression tests for reading Wannier90 data in ModelAnalyzer (wannier/symcw mode).

The symcw tests write synthetic Wannier90 files (win, nnkp, hr.dat) from the H(R) of a MultiPie model,
expressed in another primitive cell, origin and orbital order, and check that the SAMB parameters are recovered.
"""

import os
import shutil

import numpy as np
import pytest

from multipie import Group, MaterialModel
from multipie.core.cmd import create_model
from multipie.core.model_analyzer import ModelAnalyzer
from multipie.util.util_model_analyzer import fourier_r_to_k
from multipie.util.util_wannier import (
    convert_hr_to_model,
    convert_to_conventional_frac,
    create_ket_wannier_multipie,
    find_lattice_transformation,
    get_or_add_vector,
    map_wannier_to_model,
    model_primitive_vector,
    read_hr,
    read_nnkp,
    read_win,
)

# Mn2Au-type, I4/mmm (Mn 4e, Au 2a).
A_LAT, C_LAT, Z_MN = 3.328, 8.539, 1 / 3
# primitive cell of Quantum ESPRESSO, ibrav=7.
A_QE = [[A_LAT / 2, -A_LAT / 2, C_LAT / 2], [A_LAT / 2, A_LAT / 2, C_LAT / 2], [-A_LAT / 2, -A_LAT / 2, C_LAT / 2]]
# Wannier90 (l, m) of orbitals.
W90_LM = {"s": (0, 1), "pz": (1, 1), "px": (1, 2), "py": (1, 3)}


# ==================================================
def wannier_info(A, cart):
    """
    Wannier info. for create_ket_wannier_multipie (one s projection per atom).
    """
    A = np.asarray(A, dtype=float)
    frac = {k: list(np.linalg.solve(A.T, np.asarray(v, dtype=float))) for k, v in cart.items()}
    n = list(range(len(cart)))
    return dict(
        A=A,
        atoms_frac=frac,
        atoms_cart=cart,
        fermi_energy=0.0,
        nw2n=n,
        nw2l=[0] * len(n),
        nw2m=[1] * len(n),
        nw2r=[1] * len(n),
        nw2s=[0] * len(n),
    )


# ==================================================
def test_get_or_add_vector_different_length():
    lst = [np.zeros(6)]
    assert get_or_add_vector(lst, np.zeros(12)) == 1
    assert get_or_add_vector(lst, np.zeros(6)) == 0


# ==================================================
def test_ket_different_multiplicity():
    # rutile TiO2 (P4_2/mnm), Ti 2a and O 4f.
    a, c, u = 4.594, 2.959, 0.305
    cart = {
        ("Ti", 1): [0, 0, 0],
        ("Ti", 2): [a / 2, a / 2, c / 2],
        ("O", 1): [u * a, u * a, 0],
        ("O", 2): [-u * a, -u * a, 0],
        ("O", 3): [(0.5 + u) * a, (0.5 - u) * a, c / 2],
        ("O", 4): [(0.5 - u) * a, (0.5 + u) * a, c / 2],
    }
    _, _, ket, _, _ = create_ket_wannier_multipie(wannier_info(np.diag([a, a, c]), cart))
    assert len(ket) == 6


# ==================================================
def test_wyckoff_centred_lattice_in_primitive_cell():
    cart = {("Mn", 1): [0, 0, Z_MN * C_LAT], ("Mn", 2): [0, 0, -Z_MN * C_LAT], ("Au", 1): [0, 0, 0]}
    info = wannier_info(A_QE, cart)
    no, conv = convert_to_conventional_frac(info["A"], info["atoms_frac"])
    assert no == 139
    assert Group(139).find_wyckoff_site(conv[("Mn", 1)])[0] == "4e"
    assert Group(139).find_wyckoff_site(conv[("Au", 1)])[0] == "2a"
    _, _, ket, _, _ = create_ket_wannier_multipie(info)
    assert len(ket) == 3


# ==================================================
def write_wannier(path, seed, A, centers, wann, HR):
    """
    Write seed.win, seed.nnkp and seed_hr.dat.

    Args:
        wann (list): (centre index, l, m) of each Wannier function.
        HR (dict): H(R), dict[((n1,n2,n3), a, b), value].
    """
    os.makedirs(path, exist_ok=True)
    nw = len(wann)
    A = np.asarray(A)
    B = 2 * np.pi * np.linalg.inv(A).T
    with open(os.path.join(path, f"{seed}.win"), "w") as f:
        f.write(f"num_wann = {nw}\nnum_bands = {nw}\nmp_grid = 1 1 1\n\nbegin unit_cell_cart\nAng\n")
        f.writelines("  %.10f %.10f %.10f\n" % tuple(v) for v in A)
        f.write("end unit_cell_cart\n\nbegin atoms_frac\n")
        f.writelines("X %.10f %.10f %.10f\n" % tuple(c) for c in centers)
        f.write("end atoms_frac\n")
    with open(os.path.join(path, f"{seed}.nnkp"), "w") as f:
        f.write("synthetic\n\nbegin real_lattice\n")
        f.writelines("  %.7f %.7f %.7f\n" % tuple(v) for v in A)
        f.write("end real_lattice\n\nbegin recip_lattice\n")
        f.writelines("  %.7f %.7f %.7f\n" % tuple(v) for v in B)
        f.write("end recip_lattice\n\nbegin kpoints\n  1\n  0.0 0.0 0.0\nend kpoints\n\n")
        f.write(f"begin projections\n {nw:5d}\n")
        for ic, l, m in wann:
            c = centers[ic]
            f.write(f"{c[0]:10.5f}{c[1]:11.5f}{c[2]:11.5f}{l:6d}{m:3d}{1:3d}\n")
            f.write("    0.0000000  0.0000000  1.0000000   1.0000000  0.0000000  0.0000000    1.00\n")
        f.write("end projections\n\nbegin nnkpts\n  6\n")
        f.writelines("  1  1  %d %d %d\n" % G for G in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)])
        f.write("end nnkpts\n")
    Rs = sorted({k[0] for k in HR})
    with open(os.path.join(path, f"{seed}_hr.dat"), "w") as f:
        f.write(f"synthetic\n{nw}\n{len(Rs)}\n" + " ".join(["1"] * len(Rs)) + "\n")
        for R in Rs:
            for b in range(nw):
                for a in range(nw):
                    v = HR.get((R, a, b), 0.0)
                    f.write(f"{R[0]:5d}{R[1]:5d}{R[2]:5d}{a+1:5d}{b+1:5d}{v.real:20.12f}{v.imag:20.12f}\n")


# ==================================================
def model_to_wannier(model, HR, A, t, site_order, translation=None):
    """
    Express H(R) of model in the cell A with origin shift t (x_model = x_wannier U + t).
    If translation is given, the projection centres are wrapped into the home cell of A, and then translated by translation[i].

    Returns:
        - (list) -- projection centres.
        - (list) -- (centre index, l, m) of each Wannier function.
        - (dict) -- H(R) in Wannier index.
        - (list) -- Wannier index of each MultiPie ket.
    """
    U = np.asarray(A) @ np.linalg.inv(model_primitive_vector(model))
    assert np.allclose(U, np.rint(U))
    Ui = np.linalg.inv(np.rint(U))
    ket = [tuple(k) for k in model["full_matrix"]["ket"]]
    pos = np.asarray(list(model.get_ket_site().values()), dtype=float)
    site_pos = {(k[0], k[1]): p for k, p in zip(ket, pos)}
    centers = [(site_pos[s] - t) @ Ui for s in site_order]
    if translation is not None:
        centers = [c - np.floor(c) + np.asarray(n) for c, n in zip(centers, translation)]

    wann, m2w = [], [None] * len(ket)
    for ic, s in enumerate(site_order):
        for lm, i in sorted((W90_LM[k[4]], i) for i, k in enumerate(ket) if (k[0], k[1]) == s):
            m2w[i] = len(wann)
            wann.append((ic, *lm))

    HR_w = {}
    for (n1, n2, n3, m, n), v in HR.items():
        a, b = m2w[m], m2w[n]
        R = (pos[n] + np.array([n1, n2, n3], dtype=float) - pos[m]) @ Ui - (centers[wann[b][0]] - centers[wann[a][0]])
        assert np.allclose(R, np.rint(R))
        key = (tuple(int(i) for i in np.rint(R)), a, b)
        HR_w[key] = HR_w.get(key, 0) + complex(v)

    return centers, wann, HR_w, m2w


# ==================================================
@pytest.fixture(scope="module")
def centred_model(tmp_path_factory):
    """
    Create I4/mmm model with H(R) from random SAMB parameters (once per module).
    """
    topdir = str(tmp_path_factory.mktemp("wannier"))
    name = "Mn2Au"
    model = {
        "model": name,
        "group": 139,
        "cell": {"a": A_LAT, "c": C_LAT},
        "site": {"Mn": (f"[0,0,{Z_MN}]", ["px", "py", "pz"]), "Au": ("[0,0,0]", "s")},
        "bond": [("Mn", "Mn", [1, 2]), ("Mn", "Au", [1])],
        "spinful": False,
        "pdf": {"create": False},
        "qtdraw": {"create": False},
    }
    create_model(model, topdir=topdir)
    mm = MaterialModel(topdir=topdir)
    mm.load(name)
    Zr = mm.get_samb_matrix({})["matrix"]
    rng = np.random.default_rng(1)
    parameter = {z: float(rng.normal()) for z in Zr}
    HR = {k: complex(v) for k, v in mm.get_hr(parameter, Zr).items()}

    return topdir, name, mm, parameter, HR


# ==================================================
def run_symcw(topdir, name, ket_wannier=None):
    control = {
        "mode": "symcw",
        "samb": {"model": name},
        "wannier": {"seedname": name, "ket_wannier": ket_wannier or []},
        "output": {"dispersion": {"k_path": None}},
    }
    cwd = os.getcwd()
    try:
        ma = ModelAnalyzer(topdir)
        ma.analyze(control)
    finally:
        os.chdir(cwd)
    return ma


# ==================================================
def test_find_lattice_transformation(centred_model):
    _, _, mm, _, _ = centred_model
    U, _ = find_lattice_transformation(A_QE, model_primitive_vector(mm))
    assert U is not None
    assert np.allclose(np.asarray(A_QE), U @ model_primitive_vector(mm))


# ==================================================
@pytest.mark.parametrize(
    "lattice, group, site",
    [("I", 139, "[0,0,0.3]"), ("F", 225, "[0.3,0,0]"), ("C", 65, "[0.3,0,0]"), ("A", 38, "[0,0.3,0.1]"), ("R", 166, "[0,0,0.3]")],
)
def test_model_primitive_vector(lattice, group, site):
    # Cartesian positions from conventional and primitive coordinates agree modulo primitive lattice vectors.
    mm = MaterialModel()
    mm.analyze({"model": "test", "group": group, "cell": {"a": 3.0, "b": 4.0, "c": 5.0}, "site": {"X": (site, "s")}})
    assert mm.group.info.lattice == lattice
    A = np.asarray(mm["unit_vector"], dtype=float)
    Ap = model_primitive_vector(mm)
    n_lattice = {"I": 2, "F": 4, "C": 2, "A": 2, "R": 3}[lattice]
    assert np.isclose(abs(np.linalg.det(Ap)) * n_lattice, abs(np.linalg.det(A)))
    for i in mm["site"]["cell"]["X"]:
        d = np.asarray(i.position, dtype=float) @ A - np.asarray(i.position_primitive, dtype=float) @ Ap
        n = np.linalg.solve(Ap.T, d)
        assert np.allclose(n, np.rint(n), atol=1e-8)


# ==================================================
SHEAR = np.array([[1, 0, 0], [3, 1, 0], [0, -2, 1]])
# large shear amplifies the rounding of projection centres in nnkp (5 decimals).
LARGE_SHEAR = np.array([[1, 0, 0], [40, 1, 0], [0, 40, 1]])


@pytest.mark.parametrize(
    "cell, shift, site_order, translation",
    [
        ("qe", [0.0, 0.0, 0.0], [("Mn", 1), ("Au", 1), ("Mn", 2)], None),
        ("qe", [0.13, -0.21, 0.07], [("Au", 1), ("Mn", 2), ("Mn", 1)], [(0, 0, 0), (1, -1, 0), (0, 0, 2)]),
        ("multipie", [0.0, 0.5, 0.25], [("Mn", 2), ("Au", 1), ("Mn", 1)], [(0, 0, 0), (0, 0, 0), (0, 0, 0)]),
        ("sheared", [0.31, 0.0, -0.17], [("Mn", 1), ("Mn", 2), ("Au", 1)], [(0, 1, 0), (-1, 0, 0), (0, 0, 0)]),
        ("large_shear", [0.0, 0.0, 0.0], [("Mn", 1), ("Mn", 2), ("Au", 1)], [(0, 0, 0), (0, 0, 0), (0, 0, 0)]),
    ],
)
def test_symcw_centred_lattice(centred_model, cell, shift, site_order, translation):
    topdir, name, mm, parameter, HR = centred_model
    Ap = model_primitive_vector(mm)
    A = {"qe": np.asarray(A_QE), "multipie": Ap, "sheared": SHEAR @ Ap, "large_shear": LARGE_SHEAR @ Ap}[cell]
    centers, wann, HR_w, m2w = model_to_wannier(mm, HR, A, np.asarray(shift), site_order, translation)
    path = os.path.join(topdir, name, "wannier")
    shutil.rmtree(path, ignore_errors=True)
    write_wannier(path, name, A, centers, wann, HR_w)

    # H(R) is recovered in the primitive cell of the model.
    win, nnkp = read_win(name, path), read_nnkp(name, path)
    hr_dict, _, _ = read_hr(f"{name}_hr.dat", path)
    HR_m = {
        k: v for k, (v, _) in convert_hr_to_model(hr_dict, map_wannier_to_model(nnkp, win["A"], mm)).items() if abs(v) > 1e-10
    }
    assert HR_m.keys() == {k for k, v in HR.items() if abs(v) > 1e-10}
    assert max(abs(HR_m[k] - HR[k]) for k in HR_m) < 1e-10

    # eigenvalues agree at the same (Cartesian) k points.
    k_w = np.random.default_rng(0).random((5, 3))
    k_m = k_w @ np.linalg.inv(np.rint(A @ np.linalg.inv(Ap))).T
    zero = np.zeros((len(m2w), 3))
    Hk_w = fourier_r_to_k({(R, a, b): v for (R, a, b), v in HR_w.items()}, zero, k_w, s=False)
    Hk_m = fourier_r_to_k({((int(n1), int(n2), int(n3)), m, n): v for (n1, n2, n3, m, n), v in HR.items()}, zero, k_m, s=False)
    assert np.allclose(np.linalg.eigvalsh(Hk_w), np.linalg.eigvalsh(Hk_m), atol=1e-8)

    ma = run_symcw(topdir, name)
    assert ma["wannier"]["multipie_to_wannier"] == m2w
    # primitive cell of analyzer is the same as that used for the mapping.
    assert np.allclose(ma["info"]["A"], Ap)
    assert np.allclose(np.asarray(ma["info"]["A"]) @ np.asarray(ma["info"]["B"]).T, 2 * np.pi * np.eye(3))
    assert max(abs(ma.parameter[z] - v) for z, v in parameter.items()) < 1e-8

    # user-given correspondence gives the same result.
    ket = [tuple(k) for k in mm["full_matrix"]["ket"]]
    w2m = [m2w.index(w) for w in range(len(m2w))]
    ket_wannier = [f"{ket[i][4]}@{ket[i][0]}({ket[i][1]})" for i in w2m]
    ma = run_symcw(topdir, name, ket_wannier)
    assert max(abs(ma.parameter[z] - v) for z, v in parameter.items()) < 1e-8

    # wrong correspondence is rejected.
    ket_wannier[0], ket_wannier[-1] = ket_wannier[-1], ket_wannier[0]
    with pytest.raises(ValueError):
        run_symcw(topdir, name, ket_wannier)


# ==================================================
def test_symcw_centre_not_on_site(centred_model):
    topdir, name, mm, _, HR = centred_model
    centers, wann, HR_w, _ = model_to_wannier(mm, HR, A_QE, np.zeros(3), [("Mn", 1), ("Au", 1), ("Mn", 2)])
    centers[0] = centers[0] + np.array([0.0, 0.0, 0.01])
    path = os.path.join(topdir, name, "wannier")
    shutil.rmtree(path, ignore_errors=True)
    write_wannier(path, name, A_QE, centers, wann, HR_w)

    with pytest.raises(ValueError, match="cannot be mapped"):
        run_symcw(topdir, name)


# ==================================================
def test_symcw_ambiguous_mapping(tmp_path):
    # P1 with two s sites related by a half translation: the mapping depends on the origin shift.
    topdir = str(tmp_path)
    name = "ambiguous"
    model = {
        "model": name,
        "group": 1,
        "cell": {"a": 4.0, "b": 3.0, "c": 5.0},
        "site": {"A": ("[0,0,0]", "s"), "B": ("[1/2,0,0]", "s")},
        "pdf": {"create": False},
        "qtdraw": {"create": False},
    }
    create_model(model, topdir=topdir)
    mm = MaterialModel(topdir=topdir)
    mm.load(name)
    Ap = model_primitive_vector(mm)
    path = os.path.join(topdir, name, "wannier")
    write_wannier(
        path,
        name,
        Ap,
        [np.zeros(3), np.array([0.5, 0, 0])],
        [(0, 0, 1), (1, 0, 1)],
        {((0, 0, 0), 0, 0): 1.0, ((0, 0, 0), 1, 1): 2.0},
    )
    win, nnkp = read_win(name, path), read_nnkp(name, path)

    with pytest.warns(UserWarning, match="several ways"):
        mapping = map_wannier_to_model(nnkp, win["A"], mm)
    assert np.allclose(mapping["t"], 0.0)

    # explicit correspondence: no warning.
    mapping = map_wannier_to_model(nnkp, win["A"], mm, ["s@B(1)", "s@A(1)"])
    assert [mm["full_matrix"]["ket"][i][0] for i in mapping["w2m"]] == ["B", "A"]


# ==================================================
def test_symcw_rounding_of_centres(tmp_path):
    # with a large shear, five-decimal rounding of two centres in nnkp gives opposite errors (about 2e-4 each).
    topdir = str(tmp_path)
    name = "rounding"
    U = np.array([[1, 0, 0], [40, 1, 0], [0, 0, 1]])
    c = [np.array([0.1, 0.1000049, 0.0]), np.array([0.2, 0.1999951, 0.0])]
    model = {
        "model": name,
        "group": 1,
        "cell": {"a": 4.0, "b": 3.0, "c": 5.0},
        "site": {"A": ("[0.100196,0.1000049,0]", "s"), "B": ("[0.199804,0.1999951,0]", "s")},  # c U modulo 1.
        "pdf": {"create": False},
        "qtdraw": {"create": False},
    }
    create_model(model, topdir=topdir)
    mm = MaterialModel(topdir=topdir)
    mm.load(name)
    A = U @ model_primitive_vector(mm)
    path = os.path.join(topdir, name, "wannier")
    HR_w = {((0, 0, 0), 0, 0): 1.0, ((0, 0, 0), 1, 1): 2.0, ((0, 0, 0), 0, 1): 0.5, ((0, 0, 0), 1, 0): 0.5}
    write_wannier(path, name, A, c, [(0, 0, 1), (1, 0, 1)], HR_w)
    win, nnkp = read_win(name, path), read_nnkp(name, path)
    hr_dict, _, _ = read_hr(f"{name}_hr.dat", path)

    # bond (c_B - c_A) U = r_B - r_A + (4, 0, 0).
    for ket_wannier in [None, ["s@A(1)", "s@B(1)"]]:
        mapping = map_wannier_to_model(nnkp, win["A"], mm, ket_wannier)
        assert mapping["U"].tolist() == U.tolist()
        m = mapping["w2m"]
        HR = convert_hr_to_model(hr_dict, mapping)
        assert set(HR) == {(0, 0, 0, m[0], m[0]), (0, 0, 0, m[1], m[1]), (4, 0, 0, m[0], m[1]), (-4, 0, 0, m[1], m[0])}


# ==================================================
def test_symcw_periodic_images_of_site(tmp_path):
    # s and px projections on the same site, but at centres differing by a lattice vector.
    topdir = str(tmp_path)
    name = "image"
    model = {
        "model": name,
        "group": 1,
        "cell": {"a": 4.0, "b": 3.0, "c": 5.0},
        "site": {"A": ("[0.1,0.2,0.3]", ["s", "px"])},
        "pdf": {"create": False},
        "qtdraw": {"create": False},
    }
    create_model(model, topdir=topdir)
    mm = MaterialModel(topdir=topdir)
    mm.load(name)
    path = os.path.join(topdir, name, "wannier")
    c = [np.array([0.1, 0.2, 0.3]), np.array([1.1, 0.2, 0.3])]
    HR_w = {((0, 0, 0), 0, 0): 1.0, ((0, 0, 0), 1, 1): 2.0, ((0, 0, 0), 0, 1): 0.5, ((0, 0, 0), 1, 0): 0.5}
    write_wannier(path, name, model_primitive_vector(mm), c, [(0, 0, 1), (1, 1, 2)], HR_w)
    win, nnkp = read_win(name, path), read_nnkp(name, path)
    hr_dict, _, _ = read_hr(f"{name}_hr.dat", path)

    for ket_wannier in [None, ["s@A(1)", "px@A(1)"], [["A", 1, "s"], ["A", 1, "px"]]]:
        mapping = map_wannier_to_model(nnkp, win["A"], mm, ket_wannier)
        m = mapping["w2m"]
        assert [mm["full_matrix"]["ket"][i][4] for i in m] == ["s", "px"]
        assert set(convert_hr_to_model(hr_dict, mapping)) == {
            (0, 0, 0, m[0], m[0]),
            (0, 0, 0, m[1], m[1]),
            (1, 0, 0, m[0], m[1]),
            (-1, 0, 0, m[1], m[0]),
        }


# ==================================================
def test_symcw_read_ks_not_implemented(centred_model):
    topdir, name, mm, _, HR = centred_model
    centers, wann, HR_w, _ = model_to_wannier(mm, HR, A_QE, np.zeros(3), [("Mn", 1), ("Au", 1), ("Mn", 2)])
    path = os.path.join(topdir, name, "wannier")
    shutil.rmtree(path, ignore_errors=True)
    write_wannier(path, name, A_QE, centers, wann, HR_w)
    control = {"mode": "symcw", "samb": {"model": name}, "wannier": {"seedname": name, "read_KS": True}}
    cwd = os.getcwd()
    try:
        with pytest.raises(NotImplementedError):
            ModelAnalyzer(topdir).analyze(control)
    finally:
        os.chdir(cwd)
