from dataclasses import dataclass, InitVar
import sympy as sp
from copy import copy


# ==================================================
@dataclass(frozen=True, order=True)
class Tag:
    """
    default Tag class.
    """

    tag_str: InitVar[str]

    # ==================================================
    def __post_init__(self, tag_str):
        if tag_str == "":
            return
        d = self._from_str(tag_str)
        for k, v in d.items():
            object.__setattr__(self, k, v)

    # ==================================================
    @classmethod
    def _from_str(cls, tag):
        """
        create tag from str. impliment explicitly in each subclass.

        Args:
            tag (str): str to be parsed.

        Returns:
            dict: field-value pairs.
        """
        raise NotImplementedError()

    # ==================================================
    def __str__(self):
        """
        default function to print. impliment explicitly in each subclass.
        """
        raise NotImplementedError()

    # ==================================================
    def latex(self):
        """
        default function to convert latex format. impliment explicitly in each subclass.
        """
        raise NotImplementedError()

    # ==================================================
    def symbol(self, commutative=False):
        return sp.Symbol(self.latex(), commutative=commutative)

    # ==================================================
    def replace(self, **kwargs):
        """
        replace field-value.

        Args:
            kwargs (dict): field-value pairs to replace.

        Retruns:
            Tag: replaced tag.
        """
        ret = copy(self)
        for k, v in kwargs.items():
            object.__setattr__(ret, k, v)

        return ret
