# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

__version__ = "2.4.0"

import os
import importlib.resources as res
import multipie

__top_dir__ = str(res.files(multipie).parent) + os.sep
__bin_dir__ = os.path.join(__top_dir__, "multipie", "binary_data") + os.sep
