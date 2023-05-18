"""
FileManager manages current directory and formatter.
"""
import os
import time
import subprocess
from gcoreutils.io_util import read_dict, write_dict
from gcoreutils.latex_util import check_latex_installed
from multipie import get_binary
from multipie.group.point_group import PointGroup
from multipie.group.space_group import SpaceGroup


# ==================================================
class MultiPieManager:
    """
    manage directory, read, write etc.
    """

    # ==================================================
    def __init__(self, topdir=None, verbose=False, symbolic=True, formatter=True, pdf=True, parallel=True, qtdraw=True):
        """
        initialize the class.

        Args:
            topdir (str): top directory for output, ends without slash.
            verbose (bool, optional): verbose parallel info. ?
            symbolic (bool, optional): output in sympy format ?
            formatter (bool, optional): format by using black ?
            pdf (bool, optional): create pdf/tex file ?
            parallel (bool, optional): use parallel code ?
            qtdraw (bool, optional): use QtDraw ?
        """
        self._start = time.time()
        self._lap = self._start

        # initialize.
        if topdir is not None:
            os.makedirs(os.path.abspath(topdir), exist_ok=True)
            os.chdir(topdir)

        self._topdir = os.getcwd().replace(os.sep, "/")
        self._dirname = self._topdir

        self._core = get_binary()

        self._verbose = verbose
        self._symbolic = symbolic
        self._formatter = formatter
        self._pdf = pdf
        if not check_latex_installed():
            self._pdf = False
        self._parallel = parallel
        self._qtdraw = qtdraw
        try:
            import qtdraw
        except ImportError:
            self._qtdraw = False

        self._group = PointGroup("C1", core=self._core)
        self._point_group = self._group
        self._molecule = True

    # ==================================================
    @property
    def dirname(self):
        return self._dirname

    # ==================================================
    @property
    def core(self):
        return self._core

    # ==================================================
    @property
    def verbose(self):
        return self._verbose

    # ==================================================
    @property
    def symbolic(self):
        return self._symbolic

    # ==================================================
    @property
    def pdf(self):
        return self._pdf

    # ==================================================
    @property
    def parallel(self):
        return self._parallel

    # ==================================================
    @property
    def qtdraw(self):
        return self._qtdraw

    # ==================================================
    @property
    def group(self):
        return self._group

    # ==================================================
    @property
    def point_group(self):
        return self._point_group

    # ==================================================
    @property
    def molecule(self):
        return self._molecule

    # ==================================================
    def elapsed(self, from_stamp=True):
        """
        elapsed time.

        Args:
            from_stamp (bool, optional): if True, return difference from lap.

        Returns:
            float: elapsed time.
        """
        if from_stamp:
            return time.time() - self._lap
        else:
            return time.time() - self._start

    # ==================================================
    def set_group(self, group):
        """
        set group (and point group).

        Args:
            group (str or int or TagGroup): group.

        Notes:
            - space/point group can be accessed by self.group property.
            - associated point group can be accessed by self.point_group property.
        """
        if type(group) not in [str, int]:
            group = str(group)
        self._molecule = not (type(group) == int or "^" in group)
        if self._molecule:
            self._group = PointGroup(group, core=self._core)
            self._point_group = self._group
        else:
            self._group = SpaceGroup(group, core=self._core)
            self._point_group = self._group.pg

    # ==================================================
    def set_stamp(self):
        """
        set current time.
        """
        self._lap = time.time()

    # ==================================================
    def log(self, text, stamp="", end=None, file=None):
        """
        write log if verbose is True.

        Args:
            text (str): text to write.
            stamp (str or None, optional): attach elapsed time if "start" or "", otherwise no stamp is attached.
            end (str, optional): same as print end option.
            file (str, optional): same as print file option.
        """
        if self._verbose:
            if stamp is not None:
                text += " ( " + str(round(self.elapsed(stamp != "start"), 3)) + " [sec] )."
            print(text, end=end, file=file)

    # ==================================================
    def filename(self, filename, full=True):
        """
        get (full) file name.

        Args:
            filename (str): file name.
            full (bool, optional): with full path ?

        Returns:
            str: full file name (with full path).
        """
        if full:
            return self.dirname + "/" + filename
        else:
            return os.path.split(filename)[1]

    # ==================================================
    def read(self, file_dict):
        """
        read dict file or dict itself.

        Args:
            file_dict (str or dict): filename of dict. or dict.

        Returns:
            dict: read dict.
        """
        if type(file_dict) == str:
            full = self._topdir + "/" + file_dict
            if os.path.isfile(full):
                dic = read_dict(full)
                self.log(f"  * read '{full}'.", None)
            else:
                raise Exception(f"cannot open {full}.")
        else:
            dic = file_dict

        return dic

    # ==================================================
    def write(self, filename, dic, header=None, var=None):
        """
        write dict to file.

        Args:
            filename (str): file name.
            dic (dict): dict to write.
            header (str, optional): header of dict.
            var (str, optional): variable name for dict.
        """
        full = self.filename(filename)
        write_dict(full, dic, header, var)
        self.log(f"  * wrote '{filename}'.", None)

    # ==================================================
    def create_subdir(self, subdir):
        """
        create sub directory under topdir.
        """
        self._dirname = os.path.abspath(self._topdir + "/" + subdir)
        os.makedirs(self._dirname, exist_ok=True)

    # ==================================================
    def formatter(self):
        if self._formatter:
            cmd = "black --line-length=130 ."
            if self._qtdraw:
                cmd += " *.qtdw"
            try:
                subprocess.run(cmd, capture_output=True, check=True, cwd=self._dirname, shell=True)
            except subprocess.CalledProcessError:
                raise Exception("Formatting by black is failed.")
