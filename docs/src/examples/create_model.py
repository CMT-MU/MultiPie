# create_model.py
from multipie import create_samb, create_samb_qtdraw, create_samb_matrix

from C3v import C3v  # C3v model input.
from C3v_param import C3v_sel_par  # C3v SAMB selection and parameters.
from graphene import graphene  # graphene model input.
from graphene_param import graphene_sel_par  # graphene SAMB selection and parameters.

verbose = True

# ==================================================
models = [C3v, graphene]  # model input files.
select_inputs = [C3v_sel_par, graphene_sel_par]  # selection and paramter files.

create_samb(models, verbose=verbose)  # create model.
create_samb_matrix(select_inputs, verbose=verbose)  # create full matrix.
create_samb_qtdraw(models, verbose=verbose)  # create full matrix with parameters.
