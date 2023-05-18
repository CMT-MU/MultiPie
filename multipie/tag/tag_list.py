"""
This file provides utility functions for tag.
"""
import copy
from dataclasses import asdict


# ==================================================
class TagList(list):
    """
    list of tags.
    """

    # ==================================================
    def __init__(self, lst=None):
        """
        create tags.

        Args:
            lst ([Tag]): list of tag class.
        """
        if lst is None:
            lst = []
        elif not hasattr(lst, "__iter__"):
            lst = [lst]
        super().__init__(lst)

    # ==================================================
    def select(self, **kwargs):
        """
        select tags with given keywords.

        Args:
            kwargs (dict): field-value pairs to select.

        Returns:
            TagList: selected tag list.
        """
        key = kwargs.items()
        x = [i for i in self if (asdict(i).items() & key) == set(key)]
        return TagList(x)

    # ==================================================
    def latex(self):
        """
        convert to list of LaTeX.

        Returns:
            [str] or str: latex string.
        """
        return [i.latex() for i in self]

    # ==================================================
    def symbol(tags, commutative=False):
        """
        convert to list of non-commutative symbol.

        Returns:
            [sympy]: non-commutative symbol.
        """
        return [i.symbol(commutative=commutative) for i in tags]

    # ==================================================
    def str_list(tags):
        """
        convert to list of string.

        Returns:
            [str] or str: tag string.
        """
        return [str(i) for i in tags]

    # ==================================================
    @classmethod
    def from_str(self, tag_type, tag_list):
        """
        create tags.

        Args:
            tag_type (Tag): Tag class.
            tag_list ([str]): tag list in str.
        """
        return TagList([tag_type(tag) for tag in tag_list])

    # ==================================================
    def __getslice__(self, start, stop):
        lst = super().__getslice__(start, stop)
        return self.__class__(lst)

    # ==================================================
    def __getitem__(self, key):
        if isinstance(key, slice):
            lst = super().__getitem__(key)
            return self.__class__(lst)
        else:
            return super().__getitem__(key)

    # ==================================================
    def __iadd__(self, other):
        if isinstance(other, TagList):
            for i in other:
                self.append(i)
        return self

    # ==================================================
    def __add__(self, other):
        s = copy.copy(self)
        s += other
        return s
