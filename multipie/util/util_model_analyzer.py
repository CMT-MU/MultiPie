"""
Utility for ModelAnalyzer calss.
"""

import numpy as np


# ==================================================
def grid_path(pts, gpath, N1=100, A=None):
    """
    Create grid points along path.

    Args:
        pts (dict): definitions of points, e.g., {"A":[1,2,3],"B":[3,4,5]}.
        gpath (str): path, e.g., "A-B|C-D-E".
        N1 (int, optional): number of divisions.
        A (ndarray, optional): conversion matrix.

    Returns:
        - (ndarray) -- grid points.
        - (ndarray) -- linear positions.
        - (dict) -- dict[disconnected linear position, label].
    """
    if A is None:
        d = np.array(list(pts.values())[0]).shape[0]
        A = np.eye(d)

    gpath = gpath.replace(" ", "").replace("\t", "").replace("\n", "")
    glabel = gpath.split("-")
    gpath = gpath.split("|")
    gpath = [i.split("-") for i in gpath]
    gpath = [[(i1, i2) for i1, i2 in zip(i[:-1], i[1:])] for i in gpath]
    gpath = sum(gpath, [])

    grid = []
    glin = []
    gdis = [0]
    x = 0
    for s, e in gpath:
        s = np.array(pts[s])
        e = np.array(pts[e])
        dv = (e - s) / N1
        d = np.linalg.norm(dv @ A.T)
        for j in range(N1):
            glin.append(x)
            grid.append(s + j * dv)
            x += d
        grid.append(s + N1 * dv)
        glin.append(x)
        gdis.append(x)

    grid = np.asarray(grid)
    glin = np.asarray(glin)
    gdis = dict(zip(gdis, glabel))

    return grid, glin, gdis
