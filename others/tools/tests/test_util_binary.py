from multipie.util.util_binary import BinaryManager, convert_binary_to_text

# ==================================================
header1 = "* Root cluster site: [[x,y,z]]"
data1 = {
    "integer": 42,
    "float": 3.14159,
    "boolean": True,
    "complex_number": complex(1, 2),
    "string": "hello",
    "list": [1, 2, 3, 4],
    "nested_dict": {"key1": "value1", "key2": 12345},
}

header2 = "* Scalar basis on root cluster: { (cs,l,m) : [float] } c00, c10, c11, s11, c20, c21, c22, s21, ..."
data2 = {
    "integer": 24,
    "float": 1.72,
    "boolean": False,
    "complex_number": complex(2, 1),
    "string": "goodbye",
    "list": [4, 3, 2, 1],
    "nested_dict": {"key2": "value2", "key1": 54321},
}


# ==================================================
def test_util_binary():
    bm = BinaryManager(verbose=True)

    bm["data1"] = data1
    bm.add_comment(header1)
    bm["data2"] = data2
    bm.add_comment(header2)

    bm.save_binary("sample")
    bm.clear()
    bm.load_binary("sample")

    print("=====")
    print(repr(bm))
    print("=====")

    convert_binary_to_text("sample")


# ==================================================
if __name__ == "__main__":
    test_util_binary()
