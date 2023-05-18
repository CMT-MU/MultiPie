import numpy as np
import sympy as sp
from gcoreutils.nsarray import NSArray
from multipie.data.data_transform_matrix import _data_trans_Ah, _data_trans_lattice_p, _data_trans_lattice_t


# 3x3 matrix to convert from reduced to cartesian coordinate.
Ah = (NSArray(_data_trans_Ah), NSArray(_data_trans_Ah).inverse())

# 4x4 matrix to convert from conventioanl to primitive coordinate.
latticeP = {lat: (NSArray(d), NSArray(d).inverse()) for lat, d in _data_trans_lattice_p.items()}

# set of translational vectors.
latticeT = {lat: NSArray(d) for lat, d in _data_trans_lattice_t.items()}


# ==================================================
def to_cartesian(crystal, vm):
    """
    convert from reduced to cartesian coordinates.

    Args:
        vm (NSArray): vector/matrix/bond to convert.

    Returns:
        NSArray: converted array.
    """
    if crystal in ["hexagonal", "trigonal"]:
        return vm.transform(Ah[0], Ah[1])
    else:
        return vm


# ==================================================
def to_reduced(crystal, vm):
    """
    convert from cartesian to reduced coordinates.

    Args:
        vm (NSArray): vector/matrix/bond to convert.

    Returns:
        NSArray: converted array.
    """
    if crystal in ["hexagonal", "trigonal"]:
        return vm.transform(Ah[1], Ah[0])
    else:
        return vm


# ==================================================
def to_conventional(lattice, vm, plus_set=False):
    """
    convert from primitive to conventional cells in reduced coordinate.

    Args:
        lattice (str): crystal lattice (A/B/C/P/I/F/R).
        vm (NSArray): vector/matrix/bond to convert.
        plus_set (bool, optional): when True, add partial translations.

    Returns:
        NSArray: converted array. when plus_set is True, return [set(t0),set(t1),...]
    """
    if lattice == "P":
        return vm
    else:
        P = latticeP[lattice]
        vmP = vm.transform(P[0], P[1])
        if plus_set:
            tl = latticeT[lattice]
            lst = []
            if vm.style == "matrix":
                for tli in tl:
                    A = vmP.copy()
                    A[0:3, 3] = tli
                    lst.append(A)
            else:
                A = NSArray("[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]")
                for tli in tl:
                    A[0:3, 3] = tli
                    lst.append(vmP.transform(A))
            return NSArray.concat(lst)
        else:
            return vmP


# ==================================================
def to_primitive(lattice, vm):
    """
    convert from conventional to primitive cells in reduced coordinate.

    Args:
        lattice (str): crystal lattice (A/B/C/P/I/F/R).
        vm (NSArray): vector/matrix/bond to convert.

    Returns:
        NSArray: converted array.
    """
    if lattice == "P":
        return vm
    else:
        P = latticeP[lattice]
        return vm.transform(P[1], P[0])


# ==================================================
def to_axial(mp):
    """
    convert from matrix for polar vector to one for axial vector (cartesian).

    Args:
        mp (NSArray): transform matrix for polar vector (3x3).

    Returns:
        NSArray: transform matrix for axial vector (3x3).
    """

    def axial(m):
        return np.array(
            [
                [
                    m[1][1] * m[2][2] - m[1][2] * m[2][1],
                    m[1][2] * m[2][0] - m[1][0] * m[2][2],
                    m[1][0] * m[2][1] - m[1][1] * m[2][0],
                ],
                [
                    m[0][2] * m[2][1] - m[0][1] * m[2][2],
                    m[0][0] * m[2][2] - m[0][2] * m[2][0],
                    m[0][1] * m[2][0] - m[0][0] * m[2][1],
                ],
                [
                    m[0][1] * m[1][2] - m[0][2] * m[1][1],
                    m[0][2] * m[1][0] - m[0][0] * m[1][2],
                    m[0][0] * m[1][1] - m[0][1] * m[1][0],
                ],
            ]
        )

    return NSArray(NSArray._apply(axial, mp, "matrix"), "matrix")


# ==================================================
def rotation_matrix(n, axis):
    """
    rotation matrix (3x3).

    Args:
        n (int): n-fold rotation (1, 2, 3, 4, 6, -1, -3, -4, -6).
        axis (NSArray): rotation axis (normalization free).

    Returns:
        NSArray: rotatoin matrix (3x3).
    """
    th = 2 * sp.pi / n
    a = axis.normalize()

    c = sp.cos(th)
    s = sp.sin(th)
    c1 = 1 - c

    m3 = sp.Matrix(
        [
            [c + a[0] * a[0] * c1, a[0] * a[1] * c1 - a[2] * s, a[2] * a[0] * c1 + a[1] * s],
            [a[0] * a[1] * c1 + a[2] * s, c + a[1] * a[1] * c1, a[1] * a[2] * c1 - a[0] * s],
            [a[2] * a[0] * c1 - a[1] * s, a[1] * a[2] * c1 + a[0] * s, c + a[2] * a[2] * c1],
        ]
    )
    m3 = NSArray(m3.tolist(), "matrix")

    return m3


# ==================================================
def reflection_matrix(axis):
    """
    reflection matrix (3x3).

    Args:
        axis (NSArray): reflection axis (normalization free).

    Returns:
        NSArray: reflection matrix (3x3).
    """
    return -rotation_matrix(2, axis)


# ==================================================
def rotoinversion_matrix(n, axis):
    """
    rotoinversion matrix (3x3).

    Args:
        n (int): n-fold rotation.
        axis (NSArray): rotaton axis (normalization free).

    Returns:
        NSArray: rotoinversion matrix (3x3).
    """
    return -rotation_matrix(n, axis)
