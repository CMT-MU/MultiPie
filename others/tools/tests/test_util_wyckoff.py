# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np

from multipie.util.util_wyckoff import find_vector


# ==================================================
def test_find_vector():
    print("=== test_find_vector (numpy) ===")
    vs = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    print(find_vector([4, 5, 6], vs))
    print(find_vector([1, 2, 4], vs))
    print(find_vector([0, 0, 0], [[0, 0, 0], [1, 1, 1]]))
    vs = [np.random.randint(0, 10, 3) for _ in range(1000)]
    print(find_vector([3, 5, 7], vs))


# ==================================================
test_find_vector()
