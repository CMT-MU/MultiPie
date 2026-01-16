import sympy as sp
import numpy as np

from multipie.util.util_gram_schmidt import gram_schmidt


# ==================================================
def test_orthogonalize():
    print("=== test_orthogonalize (sympy) ===")
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
    vo, idx = gram_schmidt(v)
    print("nonzero idx =", idx)
    print(vo.tolist())

    print("=== test_orthogonalize (complex) ===")
    v = v.astype(complex)
    vo, idx = gram_schmidt(v)
    print("nonzero idx =", idx)
    print(vo.tolist())


# ==================================================
test_orthogonalize()
