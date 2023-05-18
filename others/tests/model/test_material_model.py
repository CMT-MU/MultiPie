from multipie import __top_dir__
from multipie.model.tests.model_example import models
from multipie.model.create_model import create_model

create_dir = __top_dir__ + "multipie/model/tests/material_model"


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
