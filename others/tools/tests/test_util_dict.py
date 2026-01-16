import pickle
from collections import namedtuple


from multipie.util.util_dict import Dict


# ==================================================
def test_Dict():
    print("=== test_Dict ===")
    tp = namedtuple("SphericalHarmonics", ["head", "s", "k", "l", "m"])
    data = Dict(tp)

    data[("Q", 0, 0, 1, 1)] = "a"
    data[("Q", 0, 0, 1, 0)] = "b"
    data[("T", 0, 0, 0, -1)] = "c"
    data[("G", 1, 1, 1, 1)] = "d"
    data[("T", 1, 1, 2, 0)] = "e"
    data[("M", 1, 1, 1, -1)] = "f"
    data[("Q", 1, 0, 3, 1)] = "g"
    data[("M", 1, 0, 1, 0)] = "h"
    data[("Q", 1, 0, 1, -1)] = "i"

    print(data.name)

    print("--- original ---")
    print(*data.field)
    for k, v in data.items():
        print(tuple(k), v)

    print("--- select s=1 ---")
    print(*data.field)
    for k, v in data.select(s=1).items():
        print(tuple(k), v)

    print("--- sort ---")
    print(*data.field)
    for k, v in data.sort().items():
        print(tuple(k), v)

    print("--- sort by s, m ---")
    print(*data.field)
    for k, v in data.sort("s", "m").items():
        print(tuple(k), v)

    print("--- sort by head(Q/G/T/M), s, m(descending) ---")
    print(*data.field)
    for k, v in data.sort(("head", ["Q", "G", "T", "M"]), "s", ("m", False)).items():
        print(tuple(k), v)

    print("--- pickle (wrote dict.pkl) ---")
    with open("dict.pkl", "wb") as f:
        pickle.dump(data, f)

    with open("dict.pkl", "rb") as f:
        loaded_dict = pickle.load(f)

    print(loaded_dict == data)
    print(loaded_dict.field)
    for i in loaded_dict.named_keys():
        print(i)


# ==================================================
test_Dict()
