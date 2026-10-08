"""
For versatile utility.
"""

import os
import re
import glob
import sys
import subprocess
import shutil
import ast
import io
import tokenize
import time
import logging
import copy
import difflib
import numpy as np
import sympy as sp
from datetime import datetime
from sympy import SympifyError
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, rationalize
from functools import wraps

TOL = 1e-11
_FORMATTER_opt = ["--line-length=300"]
_FORMATTER_max_length = 8000  # max. total length of file names in one black command.


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
def check_dict_keys(d, ref, name="input", allowed=None, free=None):
    """
    Check if all keys in dict are known, and raise error for unknown keys.

    Args:
        d (dict): dict to check.
        ref (dict): reference dict (default values).
        name (str, optional): name of dict used in error message.
        allowed (dict, optional): dict[path, [key]], acceptable keys at path instead of those in ref.
        free (list, optional): paths whose keys are not checked.

    Raises:
        ValueError: if unknown keys are found.

    Note:
        - path is keys joined by "/", e.g., "output/dispersion".
        - sub dicts are checked recursively when both values in d and ref are dict.
    """
    if allowed is None:
        allowed = {}
    if free is None:
        free = []

    unknown = []

    def check(d, ref, path):
        if path in free:
            return
        known = allowed[path] if path in allowed else ref.keys()
        for k, v in d.items():
            p = f"{path}/{k}" if path else str(k)
            if k not in known:
                close = difflib.get_close_matches(str(k), [str(i) for i in known], n=1)
                hint = f" (did you mean '{close[0]}'?)" if close else ""
                unknown.append(f"'{p}'{hint}")
            elif isinstance(v, dict) and isinstance(ref.get(k), dict):
                check(v, ref[k], p)

    check(d, ref, "")
    if unknown:
        raise ValueError(f"unknown key(s) in {name}: {', '.join(unknown)}.")


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
def _is_assignable(target):
    """
    Is target a valid python assignment target ?
    """
    try:
        ast.parse(f"{target} = 0")
        return True
    except SyntaxError:
        return False


# ==================================================
def _rename_targets(src):
    """
    Replace assignment targets which are not identifiers, e.g., "my-model_z = {...}" written by write_dict.

    Args:
        src (str): source.

    Returns:
        - (str) -- source with temporary identifiers.
        - (dict) -- original names, dict[temporary identifier, name].

    Note:
        - only a target at the beginning of a logical line, consisting of names, numbers and operators without spaces, followed by "= {", and not a valid python assignment target is replaced. Strings and comments are not changed.
    """
    try:
        tokens = list(tokenize.generate_tokens(io.StringIO(src).readline))
    except (tokenize.TokenError, SyntaxError):
        return src, {}

    skip = (tokenize.NL, tokenize.COMMENT, tokenize.INDENT, tokenize.DEDENT, tokenize.ENCODING)
    targets = []
    line = []
    for tok in tokens:
        if tok.type in skip:
            continue
        if tok.type not in (tokenize.NEWLINE, tokenize.ENDMARKER):
            line.append(tok)
            continue
        k = next((i for i, t in enumerate(line) if t.type == tokenize.OP and t.string == "="), None)
        if k and k + 1 < len(line) and line[k + 1].string == "{":
            target = line[:k]
            contiguous = all(a.end == b.start for a, b in zip(target, target[1:]))
            if (
                target[0].start[1] == 0
                and contiguous
                and all(t.type in (tokenize.NAME, tokenize.NUMBER, tokenize.OP) for t in target)
            ):
                name = "".join(t.string for t in target)
                if not _is_assignable(name):
                    targets.append((target[0].start, target[-1].end, name))
        line = []

    names = {}
    # split lines in the same way as tokenize.
    lines = io.StringIO(src).readlines()
    for (row, col0), (_, col1), name in reversed(targets):
        n = len(names)
        tmp = f"_multipie_var{n}"
        while tmp in src or tmp in names:
            n += 1
            tmp = f"_multipie_var{n}"
        names[tmp] = name
        lines[row - 1] = lines[row - 1][:col0] + tmp + lines[row - 1][col1:]

    return "".join(lines), names


# ==================================================
def _parse_dict_source(src, filename="<string>", strict=False):
    """
    Parse dicts in python source without executing it.

    Args:
        src (str): source.
        filename (str, optional): file name for error message.
        strict (bool, optional): raise ValueError for statements other than dicts and docstrings ?

    Returns:
        - (list) -- list of (variable name, dict), "dict" is used for a dict without assignment.

    Note:
        - only literals are allowed (ast.literal_eval).
        - "var = {...}" and a bare "{...}" are read, other statements are ignored (strict=False).
        - var need not be a python identifier (write_dict uses file name, e.g., "my-model_z = {...}").
    """
    names = {}
    try:
        tree = ast.parse(src, filename)
    except SyntaxError:
        src, names = _rename_targets(src)
        if not names:
            raise
        tree = ast.parse(src, filename)

    lst = []
    for node in tree.body:
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Dict):
            lst.append(("dict", ast.literal_eval(node.value)))
        elif (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
            and isinstance(node.value, ast.Dict)
        ):
            name = node.targets[0].id
            lst.append((names.get(name, name), ast.literal_eval(node.value)))
        elif strict and not (
            isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str)
        ):
            raise ValueError(f"'{filename}' (line {node.lineno}): only a dict (and docstrings) is allowed.")

    return lst


# ==================================================
def read_dict(filename, r_dir=None):
    """
    Read dict.

    Args:
        filename (str): file name.
        r_dir (str, optional): directory to read, if None, cwd is used.

    Returns:
        - (dict) -- read dict.

    Note:
        - the file must contain exactly one dict, "var = {...}" or "{...}", and no other statements except for docstrings and comments.
    """
    if r_dir is None:
        r_dir = os.getcwd()
    filename = os.path.join(r_dir, filename)

    with open(filename, mode="r", encoding="utf-8") as f:
        src = f.read()

    lst = _parse_dict_source(src, filename, strict=True)
    if len(lst) != 1:
        raise ValueError(f"'{filename}' must contain exactly one dict, but {len(lst)} found.")

    return lst[0][1]


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
        - relative file names are read from topdir. The current directory is not changed.
    """
    if topdir is None:
        topdir = os.getcwd()

    # dict or [dict].
    if isinstance(data, dict):
        return {"dict1": data}
    if isinstance(data, (tuple, list)) and data and isinstance(data[0], dict):
        return {f"dict{no+1}": i for no, i in enumerate(data)}

    if not isinstance(data, str) and not (isinstance(data, (tuple, list)) and data and isinstance(data[0], str)):
        raise TypeError(f"invalid type, {type(data)}.")

    # read files with dict data.
    if isinstance(data, str):
        data = [data]

    data = [i + ".py" if not i.endswith(".py") else i for i in data]

    result = {}
    counts = {}

    def add(name, value):
        n = counts.get(name, 0) + 1
        counts[name] = n
        if n == 1:
            result[name] = value
        else:
            result[f"{name}{n}"] = value

    for filename in data:
        with open(os.path.join(topdir, filename), "r", encoding="utf-8") as f:
            src = f.read()

        for name, value in _parse_dict_source(src, filename):
            add(name, value)

    if not result:
        raise ValueError(f"no dict is found in {data}.")

    return result


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
        do_black(os.path.dirname(filename), base + ext)


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
def do_black(w_dir, pattern="*.py"):
    """
    Execute black for Python file.

    Args:
        w_dir (str): directory.
        pattern (str, optional): pattern for black.
    """
    if check_black():
        w_dir = os.path.abspath(w_dir if w_dir else ".")
        if os.path.exists(os.path.join(w_dir, pattern)):  # file name.
            files = [os.path.join(w_dir, pattern)]
        else:  # glob pattern, expanded here without shell.
            files = sorted(glob.glob(os.path.join(glob.escape(w_dir), pattern)))
        # run with absolute paths and without changing cwd, and with -P (cwd is not added to sys.path),
        # so that a file such as w_dir/black.py is never imported instead of black.
        cmd = [sys.executable, "-P", "-m", "black"] + _FORMATTER_opt + ["--"]
        # split into batches to keep the command line short (e.g., 32767 characters on Windows).
        batch, length = [], 0
        for f in files:
            if batch and length + len(f) + 1 > _FORMATTER_max_length:
                subprocess.run(cmd + batch, capture_output=True, text=True)
                batch, length = [], 0
            batch.append(f)
            length += len(f) + 1
        if batch:
            subprocess.run(cmd + batch, capture_output=True, text=True)


# ==================================================
def simplify_ex(ex, full_factor=False):
    """
    Simplify expression.

    Args:
        ex (sympy): sympy expression.
        full_factor (bool, optional): use factor at last?

    Returns
        - (sympy) -- simplified expression.
    """
    ex = sp.together(ex)
    ex = sp.cancel(ex)
    ex = sp.factor_terms(ex, radical=True)
    ex = sp.radsimp(ex)

    if full_factor:
        ex = sp.factor(ex)

    return ex


# ==================================================
def simplify(obj, full_factor=False):
    """
    Simplify expressions.

    Args:
        obj (sympy or ndarray): expressions.
        full_factor (bool, optional): use factor at last?

    Returns:
        - (sympy or ndarray) -- simplied expressions.
    """
    if isinstance(obj, sp.Expr):
        return simplify_ex(obj, full_factor)

    if isinstance(obj, np.ndarray):
        f = np.vectorize(lambda x: simplify_ex(x, full_factor) if isinstance(x, sp.Expr) else x, otypes=[object])
        return f(obj)

    raise TypeError(f"Unsupported type: {type(obj)}")


# ==================================================
def get_n_jobs():
    """
    Get number of parallel jobs.

    Returns:
        - (int) -- number of jobs for joblib, given by environment variable MULTIPIE_N_JOBS. [default: -1 (all cores)]
    """
    n = os.environ.get("MULTIPIE_N_JOBS", "").strip()
    if n == "":
        return -1
    try:
        n = int(n)
    except ValueError:
        raise ValueError(f"MULTIPIE_N_JOBS must be integer, '{n}' is given.")
    if n == 0:
        raise ValueError("MULTIPIE_N_JOBS must not be 0.")

    return n
