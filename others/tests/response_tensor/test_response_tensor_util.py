from multipie.response_tensor.util.response_tensor_util import (
    S0,
    S1,
    S12,
    A12,
    S123,
    A123,
    S1234,
    Sb1234,
    A1234,
    Ab1234,
    M1234,
    Mb1234,
)

_dv = {1: "x", 2: "y", 3: "z"}
_dq = {(1, 1): "1", (2, 2): "2", (3, 3): "3", (2, 3): "4", (3, 1): "5", (1, 2): "6", (3, 2): "-4", (1, 3): "-5", (2, 1): "-6"}
_dq1 = {(1, 1): "1", (2, 2): "2", (3, 3): "3", (2, 3): "4", (3, 1): "5", (1, 2): "6", (3, 2): "4", (1, 3): "5", (2, 1): "6"}


# ==================================================
def print_S0():
    print(S0())


# ==================================================
def print_S1():
    d = []
    for i1 in range(1, 4):
        s = "S_{" + _dv[i1] + "} = " + str(S1(i1))
        d.append(s)
    d = sorted(d)
    for i in d:
        print(i)


# ==================================================
def print_S12():
    d = []
    for i1 in range(1, 4):
        for i2 in range(1, 4):
            s = "S_{" + _dv[i1] + _dv[i2] + "} = " + str(S12(i1, i2))
            d.append(s)
    d = sorted(d)
    for i in d:
        print(i)


# ==================================================
def print_A12():
    d = []
    for i1 in range(1, 4):
        for i2 in range(1, 4):
            s = "A_{" + _dv[i1] + _dv[i2] + "} = " + str(A12(i1, i2))
            d.append(s)
    d = sorted(d)
    for i in d:
        print(i)


# ==================================================
def print_S123():
    d = []
    for i1 in range(1, 4):
        for i2 in range(1, 4):
            for i3 in range(1, 4):
                s = "S_{" + _dq1[(i1, i2)] + _dv[i3] + "} = " + str(S123(i1, i2, i3))
                d.append(s)
    d = sorted(d)
    for i in d:
        print(i)


# ==================================================
def print_A123():
    d = []
    for i1 in range(1, 4):
        for i2 in range(1, 4):
            for i3 in range(1, 4):
                s = "A_{" + _dq[(i1, i2)] + _dv[i3] + "} = " + str(A123(i1, i2, i3))
                d.append(s)
    d = sorted(d)
    for i in d:
        print(i)


# ==================================================
def print_S1234():
    d = []
    for i1 in range(1, 4):
        for i2 in range(1, 4):
            for i3 in range(1, 4):
                for i4 in range(1, 4):
                    s = "S_{" + _dq[(i1, i2)] + _dq[(i3, i4)] + "} = " + str(S1234(i1, i2, i3, i4))
                    d.append(s)
    d = sorted(d)
    for i in d:
        print(i)


# ==================================================
def print_Sb1234():
    d = []
    for i1 in range(1, 4):
        for i2 in range(1, 4):
            for i3 in range(1, 4):
                for i4 in range(1, 4):
                    s = "\\bar{S}_{" + _dq[(i1, i2)] + _dq[(i3, i4)] + "} = " + str(Sb1234(i1, i2, i3, i4))
                    d.append(s)
    d = sorted(d)
    for i in d:
        print(i)


# ==================================================
def print_A1234():
    d = []
    for i1 in range(1, 4):
        for i2 in range(1, 4):
            for i3 in range(1, 4):
                for i4 in range(1, 4):
                    s = "A_{" + _dq[(i1, i2)] + _dq[(i3, i4)] + "} = " + str(A1234(i1, i2, i3, i4))
                    d.append(s)
    d = sorted(d)
    for i in d:
        print(i)


# ==================================================
def print_Ab1234():
    d = []
    for i1 in range(1, 4):
        for i2 in range(1, 4):
            for i3 in range(1, 4):
                for i4 in range(1, 4):
                    s = "\\bar{A}_{" + _dq[(i1, i2)] + _dq[(i3, i4)] + "} = " + str(Ab1234(i1, i2, i3, i4))
                    d.append(s)
    d = sorted(d)
    for i in d:
        print(i)


# ==================================================
def print_M1234():
    d = []
    for i1 in range(1, 4):
        for i2 in range(1, 4):
            for i3 in range(1, 4):
                for i4 in range(1, 4):
                    s = "M_{" + _dq[(i1, i2)] + _dq[(i3, i4)] + "} = " + str(M1234(i1, i2, i3, i4))
                    d.append(s)
    d = sorted(d)
    for i in d:
        print(i)


# ==================================================
def print_Mb1234():
    d = []
    for i1 in range(1, 4):
        for i2 in range(1, 4):
            for i3 in range(1, 4):
                for i4 in range(1, 4):
                    s = "\\bar{M}_{" + _dq[(i1, i2)] + _dq[(i3, i4)] + "} = " + str(Mb1234(i1, i2, i3, i4))
                    d.append(s)
    d = sorted(d)
    for i in d:
        print(i)


# ==================================================
print_S0()

print_S1()

print_S12()
print_A12()

print_S123()
print_A123()

print_S1234()
print_Sb1234()
print_A1234()
print_Ab1234()
print_M1234()
print_Mb1234()
