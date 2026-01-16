import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)

import sympy as sp
import numpy as np

from multipie.util.util import normalize_vector, read_dict_file


# ==================================================
def test_normalize_vector():
    print("=== test_normalize_vector (sympy) ===")
    v = np.array(
        [
            [sp.S(0), sp.S(1), sp.S(1), sp.S(0)],
            [sp.S(0), -sp.I, sp.I, sp.S(0)],
            [sp.S(1), sp.S(0), sp.S(0), sp.S(-1)],
            [sp.S(1), -sp.I, sp.I, sp.S(1)],
            [sp.S(2), sp.S(1), sp.S(1), sp.S(-2)],
            [sp.S(0), sp.S(0), sp.S(0), sp.S(0)],
        ]
    )
    # v = np.array([sp.S(0), sp.S(1), sp.S(1), sp.S(0)]) # single vector.
    vn = normalize_vector(v)
    print(vn.tolist())

    print("=== test_normalize_vector (float) ===")
    v = v.astype(complex)
    vn = normalize_vector(v)
    print(vn.tolist())


# ==================================================
def test_read_dict_file():
    d1 = {"a": 1}
    print(read_dict_file(d1))  # [{'a': 1}]

    d2 = [{"b": 2}]
    print(read_dict_file(d2))  # [{'b': 2}]

    d3 = read_dict_file(["model_example"], __top_dir__ + "others/tests", True)
    print([i["model"] for i in d3.values()])


# ==================================================
test_normalize_vector()
test_read_dict_file()
