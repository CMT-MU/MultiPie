from multipie.response_tensor.response_tensor_pg import ResponseTensorPG


# ==================================================
def test_response_tensor_pg():
    rt = ResponseTensorPG("D3")
    rt._info()
    rt._dump()


# ================================================== main
test_response_tensor_pg()
