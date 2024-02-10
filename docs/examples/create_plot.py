"""
example of SAMB construction and plot for C3v molecule and graphene.
"""
import numpy as np
from multipie.model.create_model import create_model
from utility import plot_model


# ==================================================
# C3v molecule parameters.
N = 100
C3v_params = {
    "x_name": "$z_8$",
    "x": np.linspace(0.0, 1.0, N),
    "z_001": [0.0] * N,
    "z_002": [1.0] * N,
    "z_003": [0.0] * N,
    "z_004": [0.0] * N,
    "z_005": [0.0] * N,
    "z_006": [0.2] * N,
    "z_007": [0.0] * N,
    "z_008": np.linspace(0.0, 1.0, N),
    "z_009": [0.0] * N,
}

# ==================================================
# graphene parameters.
graphene_params = {
    "z_001": -0.163,
    "z_002": -7.274,
    "z_003": 0.880,
    "z_004": -0.693,
    "z_005": 0.0761,
    "z_006": 0.202,
    "z_007": -0.080,
}


# ==================================================
# model parameter.
models = ["C3v.py", "graphene.py"]
model_param = [("C3v", C3v_params), ("graphene", graphene_params)]

# ==================================================
# create and models.
data_dir = __file__[: __file__.rfind("/")] + "/"
create_model(models, data_dir, verbose=True)
plot_model(model_param, data_dir, verbose=True)
