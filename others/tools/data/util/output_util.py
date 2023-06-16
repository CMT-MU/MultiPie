# ==================================================
def write_data(filename, dic, header=None, var=None, mode="w"):
    with open(filename, mode=mode, encoding="utf-8") as f:
        print('"""' + header + '"""', file=f)
        print(f"{var} =", dic, file=f, end="\n\n")
