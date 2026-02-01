__version__ = "2.0.9"

import os
import importlib.resources as res
import multipie

__top_dir__ = str(res.files(multipie).parent) + os.sep
__bin_dir__ = os.path.join(__top_dir__, "multipie", "binary_data") + os.sep
