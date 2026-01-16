"""
For customized dict.
"""

from collections.abc import MutableMapping
from collections import namedtuple


# ==================================================
class Dict(MutableMapping):
    # ==================================================
    def __init__(self, name_field, *args, **kwargs):
        """
        Named tuple dict.

        Args:
            name_field (namedtuple): Dict name, field.
        """
        self._key_type = name_field
        if len(args) == 1 and type(args[0]) == dict:
            self._data = {self._key_type(*tuple(k)): v for k, v in args[0].items()}
        else:
            self._data = dict(*args, **kwargs)

    # ==================================================
    def __getitem__(self, key):
        key = tuple(key) if isinstance(key, (list, tuple)) else (key,)  # for single value.
        return self._data[key]

    # ==================================================
    def __setitem__(self, key, value):
        key = tuple(key) if isinstance(key, (list, tuple)) else (key,)  # for single value.
        self._data[self._key_type(*key)] = value

    # ==================================================
    def __delitem__(self, key):
        key = tuple(key) if isinstance(key, (list, tuple)) else (key,)  # for single value.
        del self._data[key]

    # ==================================================
    def __iter__(self):
        return iter(self._data)

    # ==================================================
    def __len__(self):
        return len(self._data)

    # ==================================================
    def get(self, key, default=None):
        return self._data.get(key, default)

    # ==================================================
    def named_keys(self):
        """
        Keys as named tuple.

        Returns:
            - (dict_keys) -- named tuple keys.
        """
        return self._data.keys()

    # ==================================================
    def keys(self):
        if len(self.field) == 1:
            return dict.fromkeys(k[0] for k in self._data.keys()).keys()
        else:
            return dict.fromkeys(tuple(k) for k in self._data.keys()).keys()

    # ==================================================
    def values(self):
        return self._data.values()

    # ==================================================
    def named_items(self):
        """
        Items with named tuple keys.

        Returns:
            - (dict_items) -- items with named tuple keys.
        """
        return self._data.items()

    # ==================================================
    def items(self):
        return dict(zip(self.keys(), self.values())).items()

    # ==================================================
    @property
    def name(self):
        """
        Key field name.

        Returns:
            - (str) -- key field name.
        """
        return self._key_type.__name__

    # ==================================================
    @property
    def field(self):
        """
        Field names.

        Returns:
            - (list) -- field names.
        """
        return self._key_type._fields

    # ==================================================
    @property
    def key_type(self):
        """
        Key type.

        Returns:
            - (namedtuple) -- key type.
        """
        return self._key_type

    # ==================================================
    @staticmethod
    def select_key(keys, **conditions):
        """
        Select dict by conditions.

        Args:
            keys (dict_keys): named keys.
            conditions (dict): key field name and value(s) to match.
                - single value: full match.
                - list/tuple: any of given values.

        Returns:
            - (list) -- list of keys.
        """
        return [
            key
            for key in keys
            if all(
                getattr(key, attr) in value if isinstance(value, (list, tuple)) else getattr(key, attr) == value
                for attr, value in conditions.items()
            )
        ]

    # ==================================================
    def select(self, **conditions):
        """
        Select dict by conditions.

        Args:
            **conditions: key field name and value(s) to match.
                - single value: full match.
                - list/tuple: any of given values.

        Returns:
            - (Dict) -- selected Dict.
        """
        selected = Dict.select_key(self._data.keys(), **conditions)
        return Dict(self._key_type, {key: self._data[key] for key in selected})

    # ==================================================
    @staticmethod
    def sort_key(keys, key_type=None, *attributes):
        """
        Sort key.

        Args:
            keys (dict_keys): named_keys.
            key_type (namedtuple, optional): key type.

        Returns:
            - (list) -- sorted keys.

        Note:
            - attributes: tuple for sort property.
                - ("key_name", custum order list, ascending?)
                - ("key_name", custum order list)
                - ("key_name", ascending?)
                - "key_name"
        """
        if not attributes:
            if key_type is None:
                raise ValueError("key_type is required when attributes is empty")
            attributes = key_type._fields

        def sort_key(key):
            values = []
            for attr in attributes:
                if isinstance(attr, tuple):
                    if len(attr) == 3:
                        attr_name, order, asc = attr
                    elif len(attr) == 2:
                        if isinstance(attr[1], list):
                            attr_name, order = attr
                            asc = True
                        else:
                            attr_name, asc = attr
                            order = None
                    else:
                        raise ValueError("Invalid attribute tuple format.")
                else:
                    attr_name, order, asc = attr, None, True

                value = getattr(key, attr_name)

                if order is not None:
                    idx = order.index(value) if value in order else float("inf")
                    values.append(idx if asc else -idx)
                else:
                    values.append(value if asc else -value)

            return tuple(values)

        return sorted(keys, key=sort_key)

    # ==================================================
    def sort(self, *attributes):
        """
        Sort dict.

        Returns:
            - (Dict) -- sorted Dict.

        Note:
            - attributes: tuple for sort property.
                - ("key_name", custum order list, ascending?)
                - ("key_name", custum order list)
                - ("key_name", ascending?)
                - "key_name"
        """
        sorted_keys = self.sort_key(self._data.keys(), self._key_type, *attributes)
        return Dict(self._key_type, {key: self[key] for key in sorted_keys})

    # ==================================================
    def __repr__(self):
        items = (f"{tuple(key)}: {value}" for key, value in self._data.items())
        return "{" + ", ".join(items) + "}"

    # ==================================================
    def __getstate__(self):
        return {
            "name": self._key_type.__name__,
            "field": self._key_type._fields,
            "data": {tuple(k): v for k, v in self._data.items()},
        }

    # ==================================================
    def __setstate__(self, state):
        self._key_type = namedtuple(state["name"], state["field"])
        self._data = {self._key_type(*k): v for k, v in state["data"].items()}
