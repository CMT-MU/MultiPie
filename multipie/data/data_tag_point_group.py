"""
point-group information.
    point-group tag : (No, Schoenflies symbol, international symbol(short), setting, crystal)
"""
_data_tag_point_group = {
    "C1": (1, "C_{1}", "1", "", "triclinic"),
    "Ci": (2, "C_{i}", "-1", "", "triclinic"),
    "C2": (3, "C_{2}", "2", "b-axis", "monoclinic"),
    "Cs": (4, "C_{s}", "m", "b-axis", "monoclinic"),
    "C2h": (5, "C_{2h}", "2/m", "b-axis", "monoclinic"),
    "D2": (6, "D_{2}", "222", "", "orthorhombic"),
    "C2v": (7, "C_{2v}", "mm2", "", "orthorhombic"),
    "D2h": (8, "D_{2h}", "mmm", "", "orthorhombic"),
    "C4": (9, "C_{4}", "4", "", "tetragonal"),
    "S4": (10, "S_{4}", "-4", "", "tetragonal"),
    "C4h": (11, "C_{4h}", "4/m", "", "tetragonal"),
    "D4": (12, "D_{4}", "422", "", "tetragonal"),
    "C4v": (13, "C_{4v}", "4mm", "", "tetragonal"),
    "D2d": (14, "D_{2d}", "-42m", "-42m", "tetragonal"),
    "D2d-1": (14, "D_{2d}-1", "-4m2", "-4m2", "tetragonal"),
    "D4h": (15, "D_{4h}", "4/mmm", "", "tetragonal"),
    "C3": (16, "C_{3}", "3", "", "trigonal"),
    "C3i": (17, "C_{3i}", "-3", "", "trigonal"),
    "D3": (18, "D_{3}", "312", "312", "trigonal"),
    "D3-1": (18, "D_{3}-1", "321", "321", "trigonal"),
    "C3v": (19, "C_{3v}", "3m1", "3m1", "trigonal"),
    "C3v-1": (19, "C_{3v}-1", "31m", "31m", "trigonal"),
    "D3d": (20, "D_{3d}", "-31m", "-31m", "trigonal"),
    "D3d-1": (20, "D_{3d}-1", "-3m1", "-3m1", "trigonal"),
    "C6": (21, "C_{6}", "6", "", "hexagonal"),
    "C3h": (22, "C_{3h}", "-6", "", "hexagonal"),
    "C6h": (23, "C_{6h}", "6/m", "", "hexagonal"),
    "D6": (24, "D_{6}", "622", "", "hexagonal"),
    "C6v": (25, "C_{6v}", "6mm", "", "hexagonal"),
    "D3h": (26, "D_{3h}", "-6m2", "-6m2", "hexagonal"),
    "D3h-1": (26, "D_{3h}-1", "-62m", "-62m", "hexagonal"),
    "D6h": (27, "D_{6h}", "6/mmm", "", "hexagonal"),
    "T": (28, "T", "23", "", "cubic"),
    "Th": (29, "T_{h}", "m-3", "", "cubic"),
    "O": (30, "O", "432", "", "cubic"),
    "Td": (31, "T_{d}", "-43m", "", "cubic"),
    "Oh": (32, "O_{h}", "m-3m", "", "cubic"),
}


"""
point groups with complex character.
"""
_data_tag_point_group_complex = ("C4", "S4", "C4h", "C3", "C3i", "C6", "C3h", "C6h", "T", "Th")


"""
point groups with complex character.
"""
_data_tag_point_group_real = tuple([i for i in _data_tag_point_group if i not in _data_tag_point_group_complex])


"""
point groups categorized by crystal.
"""
_data_tag_point_group_crystal = {val[4]: [] for val in _data_tag_point_group.values()}
for tag, val in _data_tag_point_group.items():
    _data_tag_point_group_crystal[val[4]].append(tag)
_data_tag_point_group_crystal = {c: tuple(v) for c, v in _data_tag_point_group_crystal.items()}


"""
point groups (cubic series).
"""
_data_tag_point_group_cubic_series = (
    _data_tag_point_group_crystal["cubic"]
    + _data_tag_point_group_crystal["tetragonal"]
    + _data_tag_point_group_crystal["orthorhombic"]
    + _data_tag_point_group_crystal["monoclinic"]
    + _data_tag_point_group_crystal["triclinic"]
)


"""
point groups (hexagonal series).
"""
_data_tag_point_group_hexagonal_series = _data_tag_point_group_crystal["hexagonal"] + _data_tag_point_group_crystal["trigonal"]
