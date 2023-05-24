from qtdraw.qt_draw import QtDraw
from model_example import models

data_dir = __file__[: __file__.rfind("/")] + "/material_model/"


# ==================================================
def test_material_viewer():
    for dic in models:
        model = dic["model"]
        print("=== " + model + " ===")
        fname = data_dir + model + "/" + model + "_view.qtdw"
        QtDraw(filename=fname).show()


# ==================================================
test_material_viewer()
