"""
CharacterPGSet manages a set of point-group character tables.
"""
from multipie.tag.tag_group import TagGroup
from multipie.tag.tag_irrep import TagIrrep
from multipie.tag.tag_list import TagList
from multipie.character.character_pg import CharacterPG
from multipie.data.data_compatibility_relation import _data_compatibility_relation


# ==================================================
class CharacterPGSet(dict):  # dict of (point-group, characters), {TagGroup: CharacterPG}.
    """
    a set of point-group characters.

    Attributes:
        tag (str): class name tag.
    """

    # ==================================================
    def __init__(self):
        """
        initialize the class.
        """
        self.tag = __class__.__name__
        """class name tag."""
        for tag in TagGroup.create():
            self[tag] = CharacterPG(tag)

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

    # ==================================================
    @classmethod
    def compatibility_relation(cls):
        """
        compatibility relation.

        Returns:
            {TagGroup:[TagIrrep]}: compatibility relation table, (group tag, list of irreps.).
        """
        comp_rel = {}
        cr = _data_compatibility_relation

        for tag, ir in cr.items():
            ir = TagList.from_str(TagIrrep, ir)
            comp_rel[TagGroup(tag)] = ir

        return comp_rel
