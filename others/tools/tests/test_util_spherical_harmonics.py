import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)

import sympy as sp

from multipie.util.util import str_to_sympy

from others.tools.data.data_definition import max_s_rank
from others.tools.utils.util_spherical_harmonics import (
    Olm,
    transform_harmonics,
    _create_unitary_matrix_tesseral,
    is_transform_real,
)


# ==================================================
def test_Olm():
    rs = ["k_x", "k_y", "k_z"]
    rsp = [sp.S(1), sp.S(2), sp.S(3)]
    rip = [1, 2, 3]
    rfp = [1.0, 2.0, 3.0]
    print("=== test_Olm (monopole) ===")
    print("O(2,1)(x,y,z) =", Olm(2, 1))
    print("O(2,1)(kx,ky,kz) =", Olm(2, 1, rv=rs))
    print("O(2,1)(1,2,3) =", Olm(2, 1, rv=rsp))
    print("O(2,1)(1,2,3) =", Olm(2, 1, rv=rip))
    print("O(2,1)(1,2,3) =", Olm(2, 1, rv=rfp))

    print("=== test_Olm (dipole-hexadecapole) ===")
    for s in range(max_s_rank + 1):
        print(f"O(0,0;{s},{s})(x,y,z) =", Olm(0, 0, s, s))

    print("=== test_Olm (dipole at point) ===")
    print(f"O(1,1;{1},{1})(kx,ky,kz) =", Olm(1, 1, 1, 1, rv=rs))
    print(f"O(1,1;{1},{1})(1,2,3) =", Olm(1, 1, 1, 1, rv=rsp))
    print(f"O(1,1;{1},{1})(1,2,3) =", Olm(1, 1, 1, 1, rv=rip))
    print(f"O(1,1;{1},{1})(1,2,3) =", Olm(1, 1, 1, 1, rv=rfp))


# ==================================================
def test_transform():
    print("=== test_transform_harmonics ===")
    s = 0  # 1
    k = 0  # 1
    rfp = str_to_sympy("[[1/10,1/5,3/10], [1/5,3/10,2/5]]")
    for l in range(4):
        print(f"--- rank {l} ---")
        Olm_set = [[Olm(l, m, s, k, p) for p in rfp] for m in range(l, -l - 1, -1)]  # [m,(p,n)].
        U = _create_unitary_matrix_tesseral(l)  # [m,gamma].
        print("check real =", is_transform_real(U))
        Olg = transform_harmonics(Olm_set, U)
        print(f"[O({l};{s},{k})(ri)] =", Olg.tolist())


# ==================================================
test_Olm()
test_transform()
