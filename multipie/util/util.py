"""
For versatile utility.
"""

import os
import re
import sys
import subprocess
import shutil
import ast
import time
import logging
import copy
import importlib.util
import numpy as np
import sympy as sp
from datetime import datetime
from sympy import SympifyError
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, rationalize
from functools import wraps

TOL = 1e-11
_FORMATTER_cmd = "black --line-length=300"


# ==================================================
def _check_shape(a, shape):
    """
    Check array shape.

    Args:
        a (ndarray): array.
        shape (tuple): shape, (), (n,), (n,m), ...

    Returns:
        - (bool) -- if a is given shape, return True otherwise False.

    Note:
        - "0" in shape means any size.
    """
    if shape is None:
        return True
    return a.ndim == len(shape) and all(s == 0 or x == s for x, s in zip(a.shape, shape))


# ==================================================
def str_to_sympy(s, check_var=None, check_shape=None, rational=True, subs=None, **assumptions):
    """
    Convert a string to a sympy.

    Args:
        s (str): a string.
        check_var (list, optional): variables to accept, None (all).
        check_shape (tuple, optional): shape, (), (n,), (n,m), ...
        rational (bool, optional): use rational number ?
        subs (dict, optional): replace dict for local variables.
        **assumptions (dict, optional): common assumptions for all variables.

    Returns:
        - (ndarray) -- (list of) sympy.

    Notes:
        - if format error occurs, raise ValueError.
        - if s cannot be converted to a sympy, raise ValueError.
    """
    # reserved words in sympy (functions, constants, etc.).
    reserved = set(sp.__all__) | {"pi", "E", "I", "oo", "zoo"}

    # extract candidate variable names.
    var = sorted(set(re.findall(r"[A-Za-z_]\w*", s)))
    # remove reserved ones.
    var = [v for v in var if v not in reserved]

    # check var validation.
    if (check_var is not None) and not (set(var) <= set(check_var)):
        raise ValueError(f"not found variable '{var}' in '{check_var}'.")

    # set up local symbol environment.
    local_dict = {v: sp.Symbol(v, **assumptions) for v in var}
    if subs:
        local_dict.update(subs)

    # setup parser transformations.
    transformations = standard_transformations + (implicit_multiplication,)
    if rational:
        transformations += (rationalize,)

    # parse string.
    try:
        s = re.sub(r",\s*]", "]", s)
        expression = parse_expr(s, transformations=transformations, local_dict=local_dict)
    except (SympifyError, SyntaxError, TypeError):
        raise ValueError(f"invalid string '{s}'.")

    expression = np.asarray(expression, dtype=object)

    if not _check_shape(expression, check_shape):
        raise ValueError(f"invalid shape, {expression.shape}!={check_shape}.")

    if expression.ndim == 0:
        return expression.item()
    return expression


# ==================================================
def to_latex(a, style="scalar"):
    """
    convert list to latex list.

    Args:
        a (array-like): list of sympy.
        style (str, optional): style, "scalar/vector/matrix".

    Returns:
        - (ndarray or str) -- (list of) LaTeX string without "$".
    """
    a = np.array(a, dtype=object)

    if style == "scalar":
        if a.ndim == 0:
            return sp.latex(a.item())
        else:
            return np.vectorize(lambda x: sp.latex(x))(a).astype(object)

    elif style == "vector":

        def vec_latex(v):
            return r"\left[ " + r",\, ".join(sp.latex(x) for x in v) + r" \right]"

        if a.ndim == 1:
            return vec_latex(a)
        elif a.ndim > 1:
            s = a.shape
            sz, v = s[:-1], s[-1]
            return np.asarray([vec_latex(i) for i in a.reshape(-1, v)], dtype=object).reshape(sz)
        else:
            raise ValueError(f"invalid array shape, {a.shape}.")

    elif style == "matrix":

        def mat_latex(m):
            rows = [" & ".join(sp.latex(x) for x in row) for row in m]
            return r"\begin{bmatrix} " + r" \\ ".join(rows) + r" \end{bmatrix}"

        if a.ndim == 2:
            return mat_latex(a)
        elif a.ndim > 2:
            s = a.shape
            sz, v = s[:-2], s[-2:]
            return np.asarray([mat_latex(i) for i in a.reshape(-1, *v)], dtype=object).reshape(sz)
        else:
            raise ValueError(f"invalid array shape, {a.shape}.")

    raise ValueError(f"unknown style, {style}.")


# ==================================================
def replace(a, s):
    """
    Replace expression (exchange among variables is ok).

    Args:
        a (ndarray): array.
        s (dict): dict for substitution.

    Returns:
        - (ndarray) -- replaced array.
    """
    return np.vectorize(lambda i: i.subs(s, simultaneous=True))(a)


# ==================================================
def timer(name=None, verbose=True):
    def decorator(func):
        label = name if isinstance(name, str) else func.__name__

        @wraps(func)
        def wrapper(*args, **kwargs):
            start = time.time()
            if verbose:
                logging.info(f"=== ({label}) begin ===")
            result = func(*args, **kwargs)
            end = time.time()
            if verbose:
                logging.info(f"=== ({label}) end ({end - start:.7f} [s] elapsed) ===")
            return result

        return wrapper

    if callable(name):  # in case without arg.
        return decorator(name)
    else:
        return decorator


# ==================================================
def normalize_vector(vec, tol=TOL):
    """
    Normalize vector (sympy or complex/float).

    Args:
        vec (array-like): list of vectors.
        tol (float, optional): absolute norm tolerance for float.

    Returns:
        - (ndaray) -- normalized vector.
    """
    vec = np.asarray(vec)
    norm = np.linalg.norm if vec.dtype in [float, complex] else lambda x: sp.sqrt(np.dot(x.conjugate(), x))

    if vec.ndim == 1:
        n_vec = norm(vec)
        if n_vec > tol:
            vec /= n_vec
        return vec

    n_vec = np.apply_along_axis(norm, 1, vec)
    if vec.dtype in [float, complex]:
        n_vec[np.isclose(n_vec, tol)] = 1.0
    else:
        n_vec[n_vec == 0] = 1

    vec = vec / n_vec[:, np.newaxis]

    return vec


# ==================================================
def deep_update(d, u):
    """
    Update dict with deepcopy.

    Args:
        d (dict): dict to update (inplace).
        u (dict): additional dict.
    """
    for k, v in u.items():
        if isinstance(v, dict) and isinstance(d.get(k), dict):
            deep_update(d[k], v)
        else:
            d[k] = copy.deepcopy(v)


# ==================================================
def read_dict(filename, r_dir=None):
    """
    Read dict.

    Args:
        filename (str): file name.
        r_dir (str, optional): directory to read, if None, cwd is used.

    Returns:
        - (dict) -- read dict.
    """
    if r_dir is None:
        r_dir = os.getcwd()
    filename = os.path.join(r_dir, filename)

    with open(filename, mode="r", encoding="utf-8") as f:
        s = f.read()

    if s[: s.find("{")].count("=") > 0:
        s = s.split("=")[-1].strip(" ")

    c = ast.get_docstring(ast.parse(s))
    if c is not None:
        s = s.replace(c, "").replace('"""', "")
    dic = ast.literal_eval(s)

    return dic


# ==================================================
def read_dict_file(data, topdir=None, verbose=False):
    """
    Read dict file or dict.

    Args:
        data (str or dict): filename or dict.
        topdir (str, optional): top directory.
        verbose (bool, optional): verbose ?

    Returns:
        - (dict) -- dict[var, data_dict].

    Note:
        - if topdir is None, current directory is used.
    """
    if topdir is None:
        topdir = os.getcwd()

    cwd = os.getcwd()
    os.chdir(topdir)

    # dict or [dict.
    if isinstance(data, dict):
        return {data["model"]: data}
    if isinstance(data, (tuple, list)) and data and isinstance(data[0], dict):
        return {i["model"]: i for i in data}

    if not isinstance(data, str) and not (isinstance(data, (tuple, list)) and data and isinstance(data[0], str)):
        raise TypeError(f"invalid type, {type(data)}.")

    # read files with dict data.
    if isinstance(data, str):
        data = [data]

    data = [i + ".py" if not i.endswith(".py") else i for i in data]

    dict_vars = {}
    for i in data:
        if not os.path.exists(i):
            print(f"cannot open '{topdir}/{i}'.")
            exit(1)
        if verbose:
            print(f"read '{topdir}/{i}'.")
        # create module name.
        module_name = os.path.splitext(os.path.basename(i))[0]

        # import.
        spec = importlib.util.spec_from_file_location(module_name, i)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        dict_vars.update({k: v for k, v in vars(module).items() if isinstance(v, dict) and not k.startswith("__")})

    os.chdir(cwd)

    if not dict_vars:
        raise ValueError(f"no dict is found in {data}.")

    return dict_vars


# ==================================================
def write_dict(dic, filename, var=None, comment="", w_dir=None):
    """
    Write dict.

    Args:
        dic (dict): dict to write.
        filename (str): file name.
        var (str, optional): dict variable, if None, filename is used.
        comment (str, optional): comment.
        w_dir (str, optional): directory to write, if None, cwd is used.
    """
    filename = os.path.basename(filename)
    base, ext = os.path.splitext(filename)

    if var is None:
        var = base
    if w_dir is None:
        w_dir = os.getcwd()
    filename = os.path.join(w_dir, filename)
    if comment != "":
        comment = '"""\n' + comment + '"""\n'

    with open(filename, mode="w", encoding="utf-8") as f:
        s = comment + f"{var} = " + str(dic)
        print(s, file=f)

    if ext in [".py", ".qtdw"]:
        do_black(w_dir, ext)


# ==================================================
def setup_logging(level=logging.INFO):
    """
    Setup logging.

    Args:
        level (int, optional): log level.
    """
    logging.basicConfig(format="%(message)s", level=level, force=True, stream=sys.stdout)


# ==================================================
def time_stamp():
    """
    Get current time stamp.

    Returns:
        - (str) -- time stamp.
    """
    now = datetime.now()
    formatted = now.strftime("%Y-%m-%d %H:%M:%S")
    return formatted


# ==================================================
def progress_bar_step(length=50, label=""):
    """
    Show progress bar.

    Args:
        length (int, optional): width.
        label (str, optional): prefix.
    """
    pos = 0
    while True:
        bar = "█" * pos + "-" * (length - pos)
        sys.stdout.write(f"\r{label} |{bar}|")
        sys.stdout.flush()
        pos = (pos + 1) % (length + 1)
        yield


# ==================================================
def progress_bar_done(length=50, label=""):
    """
    Finalize progress bar.

    Args:
        length (int, optional): width.
        label (str, optional): prefix.
    """
    bar = "█" * length
    sys.stdout.write(f"\r{label} |{bar}| Done!\n")
    sys.stdout.flush()


# ==================================================
def check_qtdraw():
    """
    Check if qtdraw is installed or not.

    Returns:
        - (bool) -- installed ?
    """
    try:
        import qtdraw

        return True
    except ImportError:
        return False


# ==================================================
def check_black():
    """
    Check if black is installed or not.

    Returns:
        - (bool) -- installed ?
    """
    try:
        import black

        return True
    except ImportError:
        return False


# ==================================================
def check_latex():
    """
    Check if LaTeX is installed or not.

    Returns:
        - (bool) -- installed ?
    """
    if shutil.which("latex") is None:
        return False

    try:
        result = subprocess.run(["latex", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.returncode == 0
    except Exception:
        return False


# ==================================================
def do_black(w_dir, ext=".py"):
    """
    Execute black for Python file.

    Args:
        w_dir (str): directory.
        ext (str, optional): extention.
    """
    if check_black():
        cmd = _FORMATTER_cmd + " *" + ext
        subprocess.run(cmd, shell=True, capture_output=True, cwd=w_dir, text=True)
