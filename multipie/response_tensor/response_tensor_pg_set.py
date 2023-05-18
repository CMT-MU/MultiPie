"""
ResponseTensorPGSet manages response tensors upto rank 4 for all point groups.
"""
from multipie.tag.tag_group import TagGroup
from multipie.tag.tag_list import TagList
from multipie.response_tensor.response_tensor_pg import ResponseTensorPG


# ==================================================
class ResponseTensorPGSet(dict):  # dict of (point group, response tensors), {TagGroup: ResponseTensorPG}.
    """
    a set of ResponseTensorPG.

    Attributes:
        tag (str): class name tag.
    """

    # ==================================================
    def __init__(self):
        """
        initialize the class.
        """
        self.tag = __class__.__name__
        for tag in TagGroup.create():
            self[tag] = ResponseTensorPG(tag)

    # ==================================================
    def __str__(self):
        return self.tag

    # ==================================================
    def __repr__(self):
        return self.tag

    # ==================================================
    def latex(self):
        return self.tag

    # ==================================================
    def __getitem__(self, tag):
        if type(tag) == str:
            return self.get(TagGroup(tag))
        else:
            return self.get(tag)

    # ==================================================
    def key_list(self):
        return TagList(self.keys())
