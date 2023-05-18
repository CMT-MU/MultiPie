"""
This file provides basis list of spin-orbital space.
"""
from sympy import S

from multipie.data.data_tag_harmonics_alias import _data_alias_oh, _data_alias_d6h


# ==================================================
""" max. of orbital angular momentum. """
_MAX_L = 3  # (s,p,d,f)


# ==================================================
""" spin basis """
_s_basis = ["U", "D"]
_s_dict = {"U": "1/2", "D": "-1/2"}


# ==================================================
""" LM basis """
_lm_basis = [f"({l},{m})" for l in range(_MAX_L + 1) for m in reversed(range(-l, l + 1))]


# ==================================================
""" cubic basis (s,p,d,f) """
_cubic_basis = list(_data_alias_oh.keys())


# ==================================================
""" hexagonal basis (s,p,d,f) """
_hexagonal_basis = list(_data_alias_d6h.keys())


# ==================================================
""" JM basis """
_jm_basis = [
    f"({L + S(s) / 2},{S(m) / 2},{L})"
    for L in range(_MAX_L + 1)
    for s in (-1, 1)
    for m in reversed(range(-(2 * L + s), 2 * L + s + 2, 2))
]


# ==================================================
""" LMS basis """
_lm_s_basis = [f"({l},{m},{s})" for l in range(_MAX_L + 1) for m in reversed(range(-l, l + 1)) for s in _s_basis]


# ==================================================
""" spinful cubic basis (s,p,d,f) """
_cubic_s_basis = [f"({o},{s})" for o in _data_alias_oh.keys() for s in _s_basis]


# ==================================================
""" spinful hexagonal basis (s,p,d,f) """
_hexagonal_s_basis = [f"({o},{s})" for o in _data_alias_d6h.keys() for s in _s_basis]


# ==================================================
""" standard spinless basis """
_standard_spinless_basis = {
    "lm": _lm_basis,
    "cubic": _cubic_basis,
    "hexagonal": _hexagonal_basis,
}


# ==================================================
""" standard spinful basis """
_standard_spinful_basis = {
    "jm": _jm_basis,
    "lm": _lm_s_basis,
    "cubic": _cubic_s_basis,
    "hexagonal": _hexagonal_s_basis,
}

# ==================================================
""" standard basis """
_standard_basis = [_standard_spinless_basis, _standard_spinful_basis]

# ==================================================
# ==================================================
# ==================================================
""" abbreviated spinless basis """
_abbrev_spinless_basis = {
    "lm": ["0", "1", "2", "3"],
    "cubic": ["p", "d", "f"],
    "hexagonal": ["p", "d", "f"],
}

# ==================================================
""" abbreviated spinful basis """
_abbrev_spinful_basis = {"jm": [f"({L + S(s) / 2},{L})" for L in range(_MAX_L + 1) for s in range(-1, 2, 2) if L + S(s) / 2 > 0]}


# ==================================================
"""
abbreviated basis list
"""
_abbrev_basis = [_abbrev_spinless_basis, _abbrev_spinful_basis]
