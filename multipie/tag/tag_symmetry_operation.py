"""
TagSymmetryOperation manages tag of symmetry operation.
"""
import sympy as sp
from dataclasses import dataclass
from multipie.tag.tag import Tag


# ==================================================
@dataclass(frozen=True, order=True)
class TagSymmetryOperation(Tag):
    """
    tag of symmetry operation.
    """

    mirror: bool = False  # mirror ?
    inversion: bool = False  # inversion ?
    n: int = 1  # rotation.
    axis: str = "[0,0,1]"  # rotation axis.
    t: str = "[0,0,0]"  # partial translation.

    t0: bool = True  # zero translation ?
    pg: bool = True  # point group ?

    # ==================================================
    @classmethod
    def _from_str(cls, tag):
        """
        format of symmetry-operation string (translation is omittable).
            inversion(-) n-fold rotation(2,3+/-,4+/-,6+/-)mirror(m) axis[n1n2n3] :translation[t1,t2,t3].
        """
        s = tag
        s = s.replace("[", "").replace("]", "")  # remove []

        st = s.split(":")
        s = st[0]

        inversion = False
        mirror = False

        pos = 0
        c = s[pos]
        if c == "-":
            inversion = True
            pos += 1
        elif c == "m":
            mirror = True
            pos += 1

        n = 1
        if not mirror:
            n = int(s[pos])
            pos += 1

        if n == 3 or n == 4 or n == 6:
            d = s[pos]
            pos += 1
            if d == "-":
                n = -n

        axis = "[0,0,1]"
        if n != 1 or mirror:
            axis = TagSymmetryOperation._parse_axis(s[pos:])

        # partial translation.
        if len(st) > 1:
            pg = False
            t = st[1].split(",")
        else:
            pg = True
            t = ["0"] * 3
        t0 = t == ["0", "0", "0"]
        t = f"[{t[0]},{t[1]},{t[2]}]"

        d = {"mirror": mirror, "inversion": inversion, "n": n, "axis": axis, "t": t, "t0": t0, "pg": pg}
        return d

    # ==================================================
    def __str__(self):
        mirror = self.mirror
        inversion = self.inversion
        n = self.n
        axis = self.axis
        t = self.t
        pg = self.pg

        s = ""
        if mirror:
            s += "m"
        else:
            if inversion:
                s += "-"
            s += str(abs(n))
        if abs(n) > 2:
            s += "+" if n > 0 else "-"
        if n != 1 or mirror:
            s += axis.replace(",", "")
        if not pg:
            s += ":" + t

        return s

    # ==================================================
    def latex(self):
        astr = self.axis.replace("[", "").replace("]", "").replace(",", "")
        s = ""
        if not self.pg:
            s += "\\{"
        if self.mirror:
            s += r"{\rm m}_{" + astr + "}"
        else:
            if self.inversion:
                s += "-"
            an = abs(self.n)
            s += str(an)
            us = False
            if an > 2:
                s += "^{+}" if self.n > 0 else "^{-}"
                us = True
            if an != 1:
                if us:
                    s += "_{\\,\\," + astr + "}"
                else:
                    s += "{}_{" + astr + "}"

        if not self.pg:
            s += "|"
            if self.t0:
                s += "0"
            else:
                t = self.t.strip("[]").split(",")
                t = [sp.latex(sp.S(i)) for i in t]
                s += " ".join(t)
            s += "\\}"

        return s

    # ==================================================
    def is_point_group(self):
        return self.pg

    # ==================================================
    def is_no_translation(self):
        return self.t0

    # ==================================================
    @classmethod
    def _parse_axis(cls, tag):
        """
        parse axis string.

        Args:
            tag (str): axis string (without [],comma,space).

        Returns:
            str: axis string as [x,y,z].
        """
        s = tag.replace("[", "").replace("]", "").replace(",", "")

        pos = 0
        axis = [0, 0, 1]
        for i in range(3):
            if s[pos] == "-":
                axis[i] = int(s[pos : pos + 2])
                pos += 2
            else:
                axis[i] = int(s[pos])
                pos += 1
        return f"[{axis[0]},{axis[1]},{axis[2]}]"
