import os
from multipie.util.util import setup_logging
from multipie.core.cmd import create_samb
from model_example import models

MODEL_DIR = os.getcwd() + "/others/tests/model/"

# ==================================================
setup_logging()
create_samb(models, topdir=MODEL_DIR, verbose=True)
