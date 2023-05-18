"""
compatibility relation
    pg_tag: [ irrep ] correspondences in each column
"""
_data_compatibility_relation_cubic = {
    "Oh": ["A1g", "A2g", "Eg", "Eg", "T1g", "T1g", "T1g", "T2g", "T2g", "T2g", "A1u", "A2u", "Eu", "Eu", "T1u", "T1u", "T1u", "T2u", "T2u", "T2u"],
    "Td": ["A1", "A2", "E", "E", "T1", "T1", "T1", "T2", "T2", "T2", "A2", "A1", "E", "E", "T2", "T2", "T2", "T1", "T1", "T1"],
    "O": ["A1", "A2", "E", "E", "T1", "T1", "T1", "T2", "T2", "T2", "A1", "A2", "E", "E", "T1", "T1", "T1", "T2", "T2", "T2"],
    "Th": ["Ag", "Ag", "Ega", "Egb", "Tg", "Tg", "Tg", "Tg", "Tg", "Tg", "Au", "Au", "Eua", "Eub", "Tu", "Tu", "Tu", "Tu", "Tu", "Tu"],
    "T": ["A", "A", "Ea", "Eb", "T", "T", "T", "T", "T", "T", "A", "A", "Ea", "Eb", "T", "T", "T", "T", "T", "T"],
    "D4h": ["A1g", "B1g", "A1g", "B1g", "Eg", "Eg", "A2g", "Eg", "Eg", "B2g", "A1u", "B1u", "A1u", "B1u", "Eu", "Eu", "A2u", "Eu", "Eu", "B2u"],
    "D2d": ["A1", "B1", "A1", "B1", "E", "E", "A2", "E", "E", "B2", "B1", "A1", "B1", "A1", "E", "E", "B2", "E", "E", "A2"],
    "D2d-1": ["A1", "B2", "A1", "B2", "E", "E", "A2", "E", "E", "B1", "B1", "A2", "B1", "A2", "E", "E", "B2", "E", "E", "A1"],
    "C4v": ["A1", "B1", "A1", "B1", "E", "E", "A2", "E", "E", "B2", "A2", "B2", "A2", "B2", "E", "E", "A1", "E", "E", "B1"],
    "D4": ["A1", "B1", "A1", "B1", "E", "E", "A2", "E", "E", "B2", "A1", "B1", "A1", "B1", "E", "E", "A2", "E", "E", "B2"],
    "C4h": ["Ag", "Bg", "Ag", "Bg", "Ega", "Egb", "Ag", "Egb", "Ega", "Bg", "Au", "Bu", "Au", "Bu", "Eua", "Eub", "Au", "Eub", "Eua", "Bu"],
    "S4": ["A", "B", "A", "B", "Eb", "Ea", "A", "Ea", "Eb", "B", "B", "A", "B", "A", "Ea", "Eb", "B", "Eb", "Ea", "A"],
    "C4": ["A", "B", "A", "B", "Ea", "Eb", "A", "Eb", "Ea", "B", "A", "B", "A", "B", "Ea", "Eb", "A", "Eb", "Ea", "B"],
    "D2h": ["Ag", "Ag", "Ag", "Ag", "B3g", "B2g", "B1g", "B3g", "B2g", "B1g", "Au", "Au", "Au", "Au", "B3u", "B2u", "B1u", "B3u", "B2u", "B1u"],
    "C2v": ["A1", "A1", "A1", "A1", "B2", "B1", "A2", "B2", "B1", "A2", "A2", "A2", "A2", "A2", "B1", "B2", "A1", "B1", "B2", "A1"],
    "D2": ["A", "A", "A", "A", "B3", "B2", "B1", "B3", "B2", "B1", "A", "A", "A", "A", "B3", "B2", "B1", "B3", "B2", "B1"],
    "C2h": ["Ag", "Ag", "Ag", "Ag", "Bg", "Ag", "Bg", "Bg", "Ag", "Bg", "Au", "Au", "Au", "Au", "Bu", "Au", "Bu", "Bu", "Au", "Bu"],
    "Cs": ["A'", "A'", "A'", "A'", "A''", "A'", "A''", "A''", "A'", "A''", "A''", "A''", "A''", "A''", "A'", "A''", "A'", "A'", "A''", "A'"],
    "C2": ["A", "A", "A", "A", "B", "A", "B", "B", "A", "B", "A", "A", "A", "A", "B", "A", "B", "B", "A", "B"],
    "Ci": ["Ag", "Ag", "Ag", "Ag", "Ag", "Ag", "Ag", "Ag", "Ag", "Ag", "Au", "Au", "Au", "Au", "Au", "Au", "Au", "Au", "Au", "Au"],
    "C1": ["A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"],
}

_data_compatibility_relation_hexagonal = {
    "D6h": ["A1g", "A2g", "B1g", "B2g", "E1g", "E1g", "E2g", "E2g", "A1u", "A2u", "B1u", "B2u", "E1u", "E1u", "E2u", "E2u"],
    "D3h": ["A1'", "A2'", "A1''", "A2''", "E''", "E''", "E'", "E'", "A1''", "A2''", "A1'", "A2'", "E'", "E'", "E''", "E''"],
    "D3h-1": ["A1'", "A2'", "A2''", "A1''", "E''", "E''", "E'", "E'", "A1''", "A2''", "A2'", "A1'", "E'", "E'", "E''", "E''"],
    "C6v": ["A1", "A2", "B2", "B1", "E1", "E1", "E2", "E2", "A2", "A1", "B1", "B2", "E1", "E1", "E2", "E2"],
    "D6": ["A1", "A2", "B1", "B2", "E1", "E1", "E2", "E2", "A1", "A2", "B1", "B2", "E1", "E1", "E2", "E2"],
    "C6h": ["Ag", "Ag", "Bg", "Bg", "E1ga", "E1gb", "E2ga", "E2gb", "Au", "Au", "Bu", "Bu", "E1ua", "E1ub", "E2ua", "E2ub"],
    "C3h": ["A'", "A'", "A''", "A''", "Ea''", "Eb''", "Ea'", "Eb'", "A''", "A''", "A'", "A'", "Ea'", "Eb'", "Ea''", "Eb''"],
    "C6": ["A", "A", "B", "B", "E1a", "E1b", "E2a", "E2b", "A", "A", "B", "B", "E1a", "E1b", "E2a", "E2b"],
    "D3d": ["A1g", "A2g", "A1g", "A2g", "Eg", "Eg", "Eg", "Eg", "A1u", "A2u", "A1u", "A2u", "Eu", "Eu", "Eu", "Eu"],
    "D3d-1": ["A1g", "A2g", "A2g", "A1g", "Eg", "Eg", "Eg", "Eg", "A1u", "A2u", "A2u", "A1u", "Eu", "Eu", "Eu", "Eu"],
    "C3v": ["A1", "A2", "A2", "A1", "E", "E", "E", "E", "A2", "A1", "A1", "A2", "E", "E", "E", "E"],
    "C3v-1": ["A1", "A2", "A1", "A2", "E", "E", "E", "E", "A2", "A1", "A2", "A1", "E", "E", "E", "E"],
    "D3": ["A1", "A2", "A1", "A2", "E", "E", "E", "E", "A1", "A2", "A1", "A2", "E", "E", "E", "E"],
    "D3-1": ["A1", "A2", "A2", "A1", "E", "E", "E", "E", "A1", "A2", "A2", "A1", "E", "E", "E", "E"],
    "C3i": ["Ag", "Ag", "Ag", "Ag", "Ega", "Egb", "Ega", "Egb", "Au", "Au", "Au", "Au", "Eua", "Eub", "Eua", "Eub"],
    "C3": ["A", "A", "A", "A", "Ea", "Eb", "Ea", "Eb", "A", "A", "A", "A", "Ea", "Eb", "Ea", "Eb"],
}

# python 3.9 = dict1 | dict2
_data_compatibility_relation = {**_data_compatibility_relation_cubic, **_data_compatibility_relation_hexagonal}
