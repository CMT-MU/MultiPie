"""
This class mangaes PDF creation via LaTeX.
"""

import os
import subprocess
import math
import numpy as np


# ==================================================
def wrap_string(a, left="", right=""):
    """
    Wrap array elementwise with left and right strings.

    Args:
        a (array-like): array.
        left (str, optional): left string to wrap.
        right (str, optional): right string to wrap.

    Returns:
        - (array-like) -- wrapped array.
    """
    return np.frompyfunc(lambda i: left + i + right, 1, 1)(np.array(a, dtype=str))


# ==================================================
def regularize_table(lst2d, padding=None):
    """
    Regularize an irregular table.

    Args:
        lst2d (list): 2d list with different number of columns.
        padding (any, optional): padding value.

    Returns:
        - (list) -- regularized 2d list.
    """
    if type(lst2d) != list:
        raise KeyError(f"non list type ({type(lst2d)}) is given.")

    col = max(map(len, lst2d))
    tbl = []
    for r in lst2d:
        p = [padding] * col
        p[: len(r)] = r
        tbl.append(p)
    return tbl


# ==================================================
def list_to_table(lst1d, col, p=None):
    """
    Convert from list to table.

    Args:
        lst1d (list): 1d list.
        col (int): number of columns.
        p (any, optional): padding value (no padding for None).

    Returns:
        - (list) -- 2d list.
    """
    if type(lst1d) != list:
        raise KeyError(f"non list type ({type(lst1d)}) is given.")

    n = len(lst1d)
    row = int(math.ceil(n / col))
    tbl = []
    for i in range(row):
        tbl.append(lst1d[col * i : col * (i + 1)])
    d = row * col - n
    if p is not None and d != 0:
        tbl[-1].extend([p] * d)

    return tbl


# ==================================================
class PDFviaLaTeX:
    """
    PDF creation via LaTeX.
    """

    # Attributes:
    #    __content (list): text in each line [str]
    #    __package (list): package [(name, option)]
    #    __twice (bool): compile twice ?
    #    __fname (str): file name without ".tex"
    #    __dir (str): output directory
    #    __opt (dict): option
    #    __replace (list): replace string at creation
    #    __pos (dict): position string
    #
    # ==================================================
    def __init__(self, fname, package=[], style="normal", pt=10, no_page=False, landscape=False, english=False, dir=""):
        """
        PDF creator class.

        Args:
            fname (str): file name without ".tex".
            package (list, optional): package [(name,option)].
            style (str, optional): margin style (narrowest/narrow/normal/wide/widest).
            pt (int, optional): default point (8/9/10/11/12/14).
            no_page (bool, optional): no page number ?
            landscape (bool, optional): landscape ?
            english (bool, optional): english mode ?
            dir (str, optional): output directory.
        """
        self.__content = []
        self.__twice = False
        self.__fname = fname
        self.__dir = dir if dir else "."
        self.__opt = {"style": style, "pt": pt, "no_page": no_page, "landscape": landscape, "english": english}

        self.__package = [
            ("color", "usenames"),
            ("amsmath", ""),
            ("breqn", ""),
            ("graphicx", ""),
            ("amssymb", ""),
            ("bm", ""),
            ("enumitem", ""),
            ("braket", ""),
            ("multirow", ""),
            ("booktabs", ""),
            ("longtable", ""),
            ("caption", "justification=raggedright, singlelinecheck=off"),
        ]  #: LaTeX package (name, option)
        self.__package += package

        self.__replace = [
            (r"\left[\begin{matrix}", r"\begin{pmatrix}"),
            (r"\end{matrix}\right]", r"\end{pmatrix}"),
        ]  #: replace string (target, replace)

        self.__pos = {"c": "center", "l": "flushleft", "r": "flushright"}

    # ==================================================
    def build(self):
        """
        Create PDF.
        """
        txt = self._create_source(content=self.__content, package=self.__package, **self.__opt)
        for x in self.__replace:
            txt = txt.replace(*x)

        pdfdir = self.__dir
        cwd = os.getcwd()
        os.chdir(pdfdir)
        f = open(self.__fname + ".tex", mode="wt")
        f.write(txt)
        f.close()

        cmd = f"ptex2pdf -l -ot -synctex=0 -halt-on-error {self.__fname}.tex"  # for texlive.
        rm_file = [self.__fname + ext for ext in [".aux", ".log"]]

        try:
            subprocess.run(cmd.split(), capture_output=True, check=True, cwd=pdfdir)
            if self.__twice:
                subprocess.run(cmd.split(), capture_output=True, check=True, cwd=pdfdir)
        except subprocess.CalledProcessError:
            raise Exception(f"LaTeX compile error. See, {self.__fname}.log")

        for rm in rm_file:
            if os.path.exists(rm):
                os.remove(rm)

        if pdfdir:
            os.chdir(cwd)

    # ==================================================
    def _create_source(self, content, package=[], style="normal", pt=10, no_page=False, landscape=False, english=False):
        """
        Create tex source.

        Args:
            content (str or list): content of document.
            package (list, optional): package [(name,option)].
            style (str, optional): margin style (narrowest/narrow/normal/wide/widest).
            pt (int, optional): default point size (8/9/10/11/12/14).
            no_page (bool, optional): no page number ?
            landscape (bool, optional): landscape ?
            english (bool, optional): english mode ? (use US letter size).

        Returns:
            - (str) -- tex source.
        """
        px, py = (210, 297)  # A4
        st = "jsarticle"
        foot = 7
        layout_dict = {
            "narrowest": (-15, -15),
            "narrow": (-5, -5),
            "normal": (0, 0),
            "wide": (5, 5),
            "widest": (10, 10),
        }  # margin (width,height)

        scale_dict = {8: 9.6, 9: 9.82, 10: 10.0, 11: 10.0, 12: 10.0, 14: 9.7}
        mx, my = layout_dict[style]
        if landscape:
            px, py = py, px
            px *= scale_dict[pt] / pt
            py *= scale_dict[pt] / pt
        tx, ty = (px - 25.4 * 2 - 2 * mx, py - 25.4 * 2 - 2 * my - foot)

        package_txt = []
        if package:
            package_txt += ["\n%%% package\n"]
        for p, o in package:
            txt = r"\usepackage"
            if o:
                txt += "[" + o + "]"
            txt += "{" + p + "}"
            package_txt.append(txt)

        ls = ",landscape" if landscape else ""
        head = [
            r"\documentclass[fleqn," + str(pt) + "pt" + ls + "]{" + st + "}",
        ]
        if landscape:
            head.append(r"\special{papersize=\the\paperwidth,\the\paperheight}")

        preamble = [
            "\n%%% layout setting\n",
            f"\\paperwidth={px}truemm",
            f"\\paperheight={py}truemm",
            f"\\textwidth={tx}truemm",
            f"\\textheight={ty}truemm",
            f"\\oddsidemargin={mx}truemm",
            f"\\topmargin={my}truemm",
            f"\\footskip={foot}truemm",
            "\n%%% fixed\n",
            r"\hoffset=0truemm",
            r"\voffset=0truemm",
            r"\headheight=0truemm",
            r"\headsep=0truemm",
            r"\marginparsep=0truemm",
            r"\marginparwidth=0truemm",
            "\n%%% main text\n",
            r"\begin{document}" + "\n",
            r"\setcounter{MaxMatrixCols}{32}" + "\n",
        ]
        if english:
            preamble += [
                r"\renewcommand{\tablename}{Table~}",
                r"\renewcommand{\figurename}{Figure~}",
            ]
        if no_page:
            preamble.append(r"\pagestyle{empty}" + "\n")
        postamble = ["\n" + r"\end{document}" + "\n"]

        txt = "\n".join(head + package_txt + preamble + content + postamble)

        return txt

    # ==================================================
    def _create_table(self, tbl, cols=None, math=False):
        """
        Create table content.

        Args:
            tbl (list): content.
            cols (int, optional): number of columns.
            math (bool, optional): math mode ?

        Returns:
            - (list) -- 2d table with given cols [[str]].
            - (list) -- number of divisions in each rows [int].
        """
        if type(tbl) != list:
            raise KeyError("non-list type is given for table.")
        if len(tbl) == 0:
            raise KeyError("empty list is given for table.")
        if type(tbl[0]) != list:  # in case of 1d list
            tbl = [tbl]
        nmax = max(map(len, tbl))
        if cols is None or nmax <= cols:
            n = [1] * len(tbl)
            t = tbl
        else:
            n = []
            t = []
            for i in tbl:
                t1 = list_to_table(i, cols, "")
                t.extend(t1)
                n.append(len(t1))
        t = regularize_table(t, "")

        if math:
            t = wrap_string(t, "$ ", " $")

        return t, n

    # ==================================================
    def title(self, title, size="LARGE", pos="c"):
        """
        Add title.

        Args:
            title (str): title.
            size (str, optional): font size (large/Large/LARGE/Huge).
            pos (str, optional): position (c/l/r) (center/left/right).
        """
        pos = self.__pos[pos]
        txt = [r"\begin{" + pos + "}", f"\\{size}", title, r"\end{" + pos + "}"]
        self.text(txt)

    # ==================================================
    def text(self, content=""):
        """
        Add text.

        Args:
            content (str or list, optional): text.
        """
        if type(content) == list:
            self.__content += content
        else:
            self.__content += [content]

    # ==================================================
    def equation(self, eqs, num=False, long=False):
        """
        Add equation.

        Args:
            eqs (str or list): equation or list of eqs.
            num (bool, optional): with eq. number ?
            long (bool, optional): long eq. ? (single eq. only).
        """
        n = "" if num else "*"
        if type(eqs) == list:
            if long:
                raise KeyError("long equation is unsupported.")
            if not eqs:
                return
            txt = (
                [r"\begin{align" + n + "}"] + ["& " + i + r" \\" for i in eqs[:-1]] + ["& " + eqs[-1]] + [r"\end{align" + n + "}"]
            )
        else:
            env = "dmath" if long else "align"
            txt = [r"\begin{" + env + n + "}", eqs, r"\end{" + env + n + "}"]
        self.text(txt)

    # ==================================================
    def figure(self, fname, width=0.8, pos="c", caption=""):
        """
        Add figure (not in the figure environment, if caption is empty).

        Args:
            fname (str): file name of figure.
            width (float, optional): figure width relative to page width.
            pos (str, optional): position (c/l/r) (center/left/right).
            caption (str, optional): caption of figure.
        """
        pos = self.__pos[pos]
        txt = []
        if caption:
            txt.append(r"\begin{figure}[ht!]")
        if pos:
            txt.append(r"\begin{" + pos + "}")
        if caption:
            txt.append(r"\caption{" + caption + "}")
        txt.append(r"\includegraphics[width=" + str(width) + r"\textwidth]{" + fname + "}")
        if pos:
            txt.append(r"\end{" + pos + "}")
        if caption:
            txt.append(r"\end{figure}")
        self.text(txt)

    # ==================================================
    def simple_table(self, tbl, cols=None, cpos=None, tmath=False):
        """
        Add simple table.

        Args:
            tbl (list): content of table (1d or 2d).
            cols (int, optional): number of columns.
            cpos (str, optional): position of each cols. default is center.
            tmath (bool, optional): math mode ?
        """
        t, _ = self._create_table(tbl, cols, tmath)

        if cpos is None:
            cpos = "c" * len(t[0])

        txt = [r"\begin{tabular}{" + cpos + "}"]
        for i in t[:-1]:
            txt.append(" & ".join(i) + r" \\")
        txt.append(" & ".join(t[-1]))
        txt.append(r"\end{tabular}")

        self.text(txt)

    # ==================================================
    def _table(
        self,
        tbl,
        row,
        col,
        rc="",
        hl=[],
        caption="",
        stretch=1.0,
        cpos=None,
        rcmath=False,
        rmath=False,
        cmath=False,
        tmath=False,
        center=True,
    ):
        """
        Add table (reshaped by using number of col).

        Args:
            same as self.table() except "long".
        """
        cr = r" \\"
        hr = r" \hline"
        if rcmath and rc != "":
            rc = "$ " + rc + " $"
        if rmath:
            row = wrap_string(row, "$ ", " $")
        if cmath:
            col = wrap_string(col, "$ ", " $")

        txt = []
        if caption:
            txt.append(r"\begin{table}[ht!]")
            if center:
                txt.append(r"\begin{center}")
            txt.append(r"\caption{" + caption + "}")
            txt.append(r"\renewcommand{\arraystretch}{" + str(stretch) + "}")

        # tabular
        t, n = self._create_table(tbl, len(col), tmath)
        if cpos is None:
            cpos = "c" * (len(col) + 1)

        head = rc + " & " + " & ".join(col) + cr + hr
        txt.append(r"\begin{tabular}{" + cpos + "}" + hr + hr)
        txt.append(head)

        ti = 0
        for i in range(len(row)):
            last = hr if i in hl else ""
            if n[i] == 1:
                txt.append(str(row[i]) + " & " + " & ".join(t[ti]) + cr + last)
            else:
                txt.append(str(row[i]) + " & " + " & ".join(t[ti]) + cr)
                for j in range(1, n[i] - 1):
                    txt.append("& " + " & ".join(t[ti + j]) + cr)
                txt.append("& " + " & ".join(t[ti + n[i] - 1]) + cr + last)
            ti += n[i]
        txt.append(hr + hr)

        txt.append(r"\end{tabular}")
        # end tabular

        if caption:
            if center:
                txt.append(r"\end{center}")
            txt.append(r"\end{table}")

        self.text(txt)

    # ==================================================
    def _long_table(
        self,
        tbl,
        row,
        col,
        rc="",
        hl=[],
        caption="",
        stretch=1.0,
        cpos=None,
        rcmath=False,
        rmath=False,
        cmath=False,
        tmath=False,
        center=True,
    ):
        """
        Add table (reshaped by using number of col).

        Args:
            same as self.table() except "long".
        """
        cr = r" \\"
        hr = r" \hline"
        if rcmath and rc != "":
            rc = "$ " + rc + " $"
        if rmath:
            row = wrap_string(row, "$ ", " $")
        if cmath:
            col = wrap_string(col, "$ ", " $")
        t, n = self._create_table(tbl, len(col), tmath)
        if cpos is None:
            cpos = "c" * (len(col) + 1)

        mc = r"\multicolumn{" + str(len(col)) + "}"
        head = rc + " & " + " & ".join(col) + cr + hr

        txt = []
        cnt = "" if center else "[l]"

        txt.append(r"\renewcommand{\arraystretch}{" + str(stretch) + "}")
        txt.append(r"\begin{longtable}" + cnt + "{" + cpos + "}")
        if caption:
            txt.append(r"\caption{" + caption + "}")
            txt.append(cr)

        txt.append(hr + hr)
        txt.append(head + r" \endfirsthead" + "\n")

        txt.append(mc + r"{l}{\tablename\ \thetable{}}" + cr)
        txt.append(hr + hr)
        txt.append(head + r" \endhead" + "\n")

        txt.append(hr + hr)
        txt.append(mc + r"{r}{\footnotesize\it continued ...}" + cr + r" \endfoot" + "\n")

        txt.append(hr + hr)
        txt.append(mc + "{r}{}" + cr + r" \endlastfoot" + "\n")

        ti = 0
        for i in range(len(row)):
            last = hr if i in hl else ""
            if n[i] == 1:
                txt.append(str(row[i]) + " & " + " & ".join(t[ti]) + cr + last)
            else:
                txt.append(str(row[i]) + " & " + " & ".join(t[ti]) + cr)
                for j in range(1, n[i] - 1):
                    txt.append("& " + " & ".join(t[ti + j]) + cr)
                txt.append("& " + " & ".join(t[ti + n[i] - 1]) + cr + last)
            ti += n[i]

        txt.append(r"\end{longtable}")

        self.text(txt)
        self.__twice = True

    # ==================================================
    def table(
        self,
        tbl,
        row,
        col,
        rc="",
        hl=[],
        caption="",
        stretch=1.0,
        cpos=None,
        rcmath=False,
        rmath=False,
        cmath=False,
        tmath=False,
        long=False,
        center=True,
    ):
        """
        Add table (reshaped by using number of col).

        Args:
            tbl (list): content of table [[str]].
            row (list): row name [str].
            col (list): column name [str].
            rc (str, optional): row-col name.
            hl (list, optional): position with horizontal line [int]. if True, use horizontal lines for all lines.
            caption (str, optional): table caption.
            stretch (float, optional): baseline stretch (ratio).
            cpos (str, optional): column position. default is center.
            rcmath (bool, optional): math mode for row-col name ?
            rmath (bool, optional): math mode for row name ?
            cmath (bool, optional): math mode for column name ?
            tmath (bool, optional): math mode for table ?
            long (bool, optional): long table ?
            center (bool, optional): centering ?
        """
        if hl is True:
            hl = list(range(len(row) - 1))
        if long:
            self._long_table(tbl, row, col, rc, hl, caption, stretch, cpos, rcmath, rmath, cmath, tmath, center)
        else:
            self._table(tbl, row, col, rc, hl, caption, stretch, cpos, rcmath, rmath, cmath, tmath, center)
