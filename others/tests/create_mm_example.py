# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import os
from multipie.util.util import setup_logging
from multipie.core.cmd import create_model
from model_example import models

MODEL_DIR = os.getcwd() + "/others/tests/model/"

# ==================================================
setup_logging()
create_model(models, topdir=MODEL_DIR, verbose=True)
