"""
Model analyzer class.

This module provides model analyzer.
"""

import os
from multipie.core.material_model import MaterialModel


# ==================================================
class ModelAnalyzer:
    # ==================================================
    def __init__(self, model=None, topdir=None, verbose=False):
        """
        Model analyzer.

        Args:
            model (str, optional): model name.
            topdir (str, optional): top directory. [default: cwd]
            verbose (bool, optional): verbose comment ?
        """
        if topdir is None:
            topdir = os.getcwd()

        self.mm = MaterialModel(topdir, verbose=verbose)
        if model:
            self.mm.load(model)
            self.loaded = True
        else:
            self.loaded = False

    # ==================================================
    def load(self, model):
        """
        Load model file.

        Args:
            model (str): model name.
        """
        self.mm.load(model)
