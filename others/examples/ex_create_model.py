from multipie import __top_dir__
from multipie.model.create_model import create_model

from models.cm import CM
from models.graphene import graphene
from models.te import Te

topdir = __top_dir__ + "multipie/examples/material_model"


models = [CM, graphene, Te]
create_model(
    models, topdir=topdir, symbolic=True, parallel=True, verbose=True, pdf=True, formatter=True, view_mode=None, qtdraw=True
)
