from multipie.tag.tag_response_tensor import TagResponseTensor
from multipie.tag.tag_list import TagList


# ==================================================
def test_tag_response_tensor():
    print("=== tag_response_tensor ===")

    head_lst = ["Q", "G", "T", "M"]
    comp_dict = {0: ("s"), 1: ("s"), 2: ("s", "a"), 3: ("s", "a"), 4: ("sss", "ssa", "aas", "aaa", "sa", "as")}

    tags = []
    for head in head_lst:
        for rank, comp_lst in comp_dict.items():
            for comp in comp_lst:
                tags.append(TagResponseTensor.create(head, rank, comp))
    tags = TagList(tags)

    print("--- repr ---")
    print(tags[:3])

    print("--- latex ---")
    print(tags[:3].latex())

    print("--- symbol ---")
    print(tags[:3].symbol())

    print("--- str ---")
    print(tags[:3].str_list())

    print("--- first and last 3 ---")
    print(*(tags[:3] + tags[-3:]))

    print("--- rank 1 ---")
    print(*tags.select(rank=1))

    print("--- rank 1 and G or M ---")
    print(*(tags.select(rank=1, head="G") + tags.select(rank=1, head="M")))

    print("--- asymmetric ---")
    print(*tags.select(comp="a"))

    print("--- tensor_head ('C') ---")
    print(tags[0].tensor_head)


# ================================================== main
test_tag_response_tensor()
