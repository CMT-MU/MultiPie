import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)

import numpy as np
import sympy as sp

from others.tools.utils.util_symmetry_operation_spherical import (
    _rotation_matrix,
    _axis_to_euler_zyz,
    _wigner_D_matrix,
    representation_matrix,
)
from others.tools.utils.util_spherical_harmonics import _create_unitary_matrix_tesseral, create_harmonics_set, transform_harmonics
from others.tools.utils.util_symmetry_operation import symmetry_operation_matrix_multipole


# ==================================================
def test_symmetry_operation():
    print("=== Test for rotation_matrix ===")
    R = _rotation_matrix(2 * sp.pi / 3, [1, 1, 1])
    print("R[2pi/3, (1,1,1)] =", R)

    print("=== Test for axis_to_eular_zyz, wigner_D_matrix ===")
    angles = _axis_to_euler_zyz(2 * sp.pi / 3, [1, 1, 1])
    print("[2pi/3, (1,1,1)] => alpha, beta, gamma:", angles)

    print("=== Test for wigner_D_matrix ===")
    l = 1
    D = _wigner_D_matrix(l, 2 * sp.pi / 3, [1, 1, 1])
    print(f"D[{l}, 2pi/3, (1,1,1)] [m1,m2] =", D)
    U = _create_unitary_matrix_tesseral(l)
    Dc = (U.conjugate().T) @ D @ U
    print("U^dagger D U [cartesian] =", Dc)

    print("=== Test for representation_matrix ===")
    print("C3[111] =", representation_matrix(1, 3, [1, 1, 1], tesseral=True))
    print("C3[111] for rank 2 =", representation_matrix(2, 3, [1, 1, 1], tesseral=True))
    print("C3[001] =", representation_matrix(1, 3, [0, 0, 1], tesseral=True))
    print("-C3[001] =", representation_matrix(1, 3, [0, 0, 1], inversion=True, tesseral=True))

    print("=== Test for harmonics transformation ===")
    l, s, k = 2, 0, 0
    olm = create_harmonics_set(l, s, k, sv=True)
    print(f"O{l};{s},{k} =", olm)
    du = np.array([0, 0, 1, 0, 0])
    print("du =", du)
    print("C3[111] * du =", transform_harmonics(olm, du))
    g = representation_matrix(2, 3, [1, 1, 1], tesseral=False)
    dup = g @ du
    print("du' =", transform_harmonics(olm, dup))

    print("m[001] =", symmetry_operation_matrix_multipole(1, "m[001]", axial=False, tesseral=True))


# ==================================================
test_symmetry_operation()
