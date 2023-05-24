from multipie import __top_dir__
from model_example import models
from multipie.model.create_model import create_model

create_dir = __top_dir__ + "others/tests/model/material_model"


# ==================================================
create_model(
    models,
    topdir=create_dir,
    symbolic=True,
    parallel=True,
    verbose=True,
    pdf=True,
    formatter=True,
    view_mode="debug",
    qtdraw=True,
)
