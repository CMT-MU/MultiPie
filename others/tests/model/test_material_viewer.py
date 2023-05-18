from qtdraw.qt_draw import QtDraw
from multipie import __top_dir__
from multipie.model.tests.model_example import models

create_dir = __top_dir__ + "multipie/model/tests/material_model/"


# ==================================================
def test_material_viewer():
    for dic in models:
        model = dic["model"]
        print("=== " + model + " ===")
        fname = create_dir + model + "/" + model + "_view.qtdw"
        QtDraw(filename=fname).show()


# ==================================================
test_material_viewer()
