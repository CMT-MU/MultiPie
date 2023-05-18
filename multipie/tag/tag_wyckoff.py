"""
TagWyckoff manages tag of wyckoff position.
"""
from dataclasses import dataclass
from multipie.tag.tag import Tag
from multipie.tag.tag_list import TagList
from multipie.data.data_tag_point_group import _data_tag_point_group
from multipie.data.data_wyckoff_pg import _data_wyckoff_pg
from multipie.data.data_wyckoff_sg import _data_wyckoff_sg


# ==================================================
@dataclass(frozen=True, order=True)
class TagWyckoff(Tag):
    """
    tag of wyckoff position.
    """

    letter: str = "a"  # Wyckoff letter.
    n: int = 1  # number of sites.

    # ==================================================
    @classmethod
    def _from_str(cls, tag):
        letter = tag[-1]
        n = int(tag[:-1])
        d = {"letter": letter, "n": n}
        return d

    # ==================================================
    def __str__(self):
        s = str(self.n) + self.letter
        return s

    # ==================================================
    def latex(self):
        return str(self)

    # ==================================================
    @classmethod
    def create(cls, tag):
        """
        point-group tags.

        Args:
            tag (str): point-group/space-group tag.

        Returns:
            [TagWyckoff]: tags.
        """
        if tag in _data_tag_point_group.keys():
            tags = TagList.from_str(TagWyckoff, [i[0] for i in _data_wyckoff_pg[tag]])
        else:
            tags = TagList.from_str(TagWyckoff, [i[0] for i in _data_wyckoff_sg[tag]])

        return tags
