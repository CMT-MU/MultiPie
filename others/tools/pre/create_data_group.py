# SPDX-License-Identifier: MPL-2.0
# Copyright (C) 2021- The MultiPie developers
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
Create group data.
Created "output/data_group_pg/sg/mpg/msg.py" should be moved to "multipie/data".
"""

from create_data_group_msg import create_msg

create_msg()

from create_data_group_mpg import create_mpg

create_mpg()

from create_data_group_sg import create_sg

create_sg()

from create_data_group_pg import create_pg

create_pg()
