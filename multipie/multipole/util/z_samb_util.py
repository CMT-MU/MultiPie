"""
This file provides utility functions for combined multipole set.
"""
from itertools import product

from gcoreutils.nsarray import NSArray

from multipie.tag.tag_list import TagList
from multipie.const import is_electric


# ==================================================
def create_z_samb(cg, hs, tag1_list, tag2_list, toroidal_priority=False, **kwargs):
    """
    initialize combined multipoles.

    Args:
        cg (clebsch_gordan_pg): Clebsh-Gordan coefficients.
        hs (HarmonicsPG): the class to manage harmonics.
        tag1_list (TagList): multipole/harmonics tag list.
        tag2_list (TagList): multipole/harmonics tag list.
        toroidal_priority (bool): create toroidal multipoles (G,T) in priority? else prioritize conventional multipoles (Q,M).
        kwargs (dict, optional): select conditions for multipoles,
                                 keywords in ["head", "rank", "irrep", "mul", "comp", "s", "k"].
    Returns:
        dict: {TagMultipole: [(coefficient, TagMultipole(atomic), TagMultipole(site/bond)] }.
    """
    key_list = ["head", "rank", "irrep", "mul", "comp", "s", "k"]
    for key in kwargs.keys():
        if key not in key_list:
            raise Exception(f"{key} cannot be specified in kwargs.")

    kwargs = {k: [v] if type(v) != list else v for k, v in kwargs.items()}
    keys, values = list(kwargs.keys()), list(kwargs.values())
    idx_pairs = product(*[range(len(v_list)) for v_list in values])
    kwargs_list = [{keys[i]: values[i][j] for i, j in enumerate(idx_pair)} for idx_pair in idx_pairs]

    head_list = ["G", "Q", "T", "M"] if toroidal_priority else ["Q", "G", "M", "T"]
    Ztag_list = [tag.replace(m_type="") for tag in hs.key_list()]
    Ztag_list += [tag.replace(head="M") for tag in TagList(Ztag_list).select(head="G")]
    Ztag_list += [tag.replace(head="T") for tag in TagList(Ztag_list).select(head="Q")]

    Ztag_list = [tag for kwargs in kwargs_list for tag in TagList(Ztag_list).select(**kwargs)]
    Ztag_list = sorted(set(Ztag_list), key=Ztag_list.index)
    Ztag_list = [tag for l in range(12) for Z in head_list for tag in Ztag_list if tag.rank == l and tag.head == Z]
    Ztag_list = TagList(Ztag_list)

    is_inversion = True if "g" in Ztag_list[0].irrep or "u" in Ztag_list[0].irrep else False
    irrep_comp_list = [(tag.irrep, tag.comp) for tag in Ztag_list] if is_inversion else [tag.irrep for tag in Ztag_list]
    irrep_comp_list = sorted(set(irrep_comp_list), key=irrep_comp_list.index)

    head_s_k_l_list1 = [(tag1.head, tag1.s, tag1.k, tag1.l) for tag1 in tag1_list]
    head_s_k_l_list1 = sorted(set(head_s_k_l_list1), key=head_s_k_l_list1.index)
    head_l_list2 = [(tag2.head, tag2.l) for tag2 in tag2_list]
    head_l_list2 = sorted(set(head_l_list2), key=head_l_list2.index)

    z_samb = {}
    for X, s, k, l1 in head_s_k_l_list1:
        tag1_list_ = tag1_list.select(head=X, s=s, k=k, rank=l1)
        for Y, l2 in head_l_list2:
            tag2_list_ = tag2_list.select(head=Y, rank=l2)
            tag1_tag2_list = [(tag1, tag2) for tag1 in tag1_list_ for tag2 in tag2_list_]
            dim = len(tag1_list_) * len(tag2_list_)
            min_l, max_l = abs(l1 - l2), l1 + l2
            t1 = +1 if is_electric(X) else -1
            t2 = +1 if is_electric(Y) else -1
            t_type = "electric" if t1 * t2 == 1 else "magnetic"

            for irrep_comp in irrep_comp_list:
                if is_inversion:
                    irrep, comp = irrep_comp
                    Ztag_list_ = Ztag_list.select(irrep=irrep, comp=comp)
                else:
                    irrep = irrep_comp
                    Ztag_list_ = Ztag_list.select(irrep=irrep)

                Ztag_list_ = [tag for tag in Ztag_list_ if tag.t_type == t_type and min_l <= tag.l and tag.l <= max_l]

                tags, vecs = [], []
                for tag in Ztag_list_:
                    lst = []
                    for tag1, tag2 in tag1_tag2_list:
                        c = cg.cg(tag1, tag2, tag)
                        lst.append(c)

                    if lst != [0] * dim:
                        vecs.append(NSArray(lst, style="vector"))
                        tags.append(tag.replace(s=s, k=k))

                # orthogonalize
                if vecs != []:
                    vecs = NSArray(vecs, style="vector")
                    vecs, idx = NSArray.orthogonalize(vecs, dim)
                    vecs = vecs[idx]
                    tags = [tags[i] for i in idx]

                    for tag, vec in zip(tags, vecs):
                        lst = []
                        for i, c in enumerate(vec):
                            if c != 0:
                                tag1, tag2 = tag1_tag2_list[i]
                                lst.append((c, tag1, tag2))
                        _, tag1, tag2 = lst[0]
                        z_samb[(tag, tag1, tag2)] = lst

    return z_samb
