from model_example import models
from multipie.model.create_model import create_model

create_dir = __file__[: __file__.rfind("/")] + "/material_model"


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
