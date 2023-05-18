from multipie.tag.tag_symmetry_operation import TagSymmetryOperation
from multipie.tag.tag_list import TagList


# ==================================================
def test_tag_symmetry_operation():
    print("=== tag_symmetry_operation ===")

    lst = [
        "1",
        "2[001]",
        "2[010]",
        "2[100]",
        "2[-101]",
        "2[01-1]",
        "2[011]",
        "2[1-10]",
        "2[101]",
        "2[110]",
        "3+[-1-11]",
        "3+[-11-1]",
        "3+[1-1-1]",
        "3+[111]",
        "3-[-1-11]",
        "3-[-11-1]",
        "3-[1-1-1]",
        "3-[111]",
        "4+[001]",
        "4+[010]",
        "4+[100]",
        "4-[001]",
        "4-[010]",
        "4-[100]",
        "-1",
        "m[001]",
        "m[010]",
        "m[100]",
        "m[-101]",
        "m[01-1]",
        "m[011]",
        "m[1-10]",
        "m[101]",
        "m[110]",
        "-3+[-1-11]",
        "-3+[-11-1]",
        "-3+[1-1-1]",
        "-3+[111]",
        "-3-[-1-11]",
        "-3-[-11-1]",
        "-3-[1-1-1]",
        "-3-[111]",
        "-4+[001]",
        "-4+[010]",
        "-4+[100]",
        "-4-[001]",
        "-4-[010]",
        "-4-[100]",
        "2[010]:[0,0,0]",
        "2[010]:[0,1/2,1/2]",
        "m[010]:[0,1/2,1/2]",
    ]

    tags = TagList.from_str(TagSymmetryOperation, lst)

    print("--- repr ---")
    print(tags[:3])

    print("--- latex ---")
    print(tags[:3].latex())

    print("--- symbol ---")
    print(tags[:3].symbol())

    print("--- str ---")
    print(tags[:3].str_list())

    print("--- first and last 3 ---")
    print(*(tags[:3] + tags[-3:]).latex())

    print("--- mirror ---")
    print(*tags.select(mirror=True).latex())


# ================================================== main
test_tag_symmetry_operation()
