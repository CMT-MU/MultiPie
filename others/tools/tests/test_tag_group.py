# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from multipie.core.group import Group


# ==================================================
def test_tag_group():
    print("=== tag_group ===")

    tag = ["C3v", "D3-1", "D6h^3", "Oh^4"]

    for t in tag:
        g = Group(t)
        t1 = str(g)
        print(f"'{t}' = '{t1}'")
        print(g.latex())
        print(g.latex(detail=True))


# ================================================== main
test_tag_group()
