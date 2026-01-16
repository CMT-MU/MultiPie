"""
For tag.
"""

import sympy as sp

from multipie.util.util import str_to_sympy


# ==================================================
class TagSymmetryOperation:
    # ==================================================
    @staticmethod
    def parse_axis(tag):
        """
        Parse axis string.

        Args:
            tag (str): axis string (without [],comma,space).

        Returns:
            - (str) axis, [a,b,c].
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
        axis = str(axis).replace(" ", "")

        return axis

    # ==================================================
    @staticmethod
    def parse(tag, hexagonal=False):
        """
        Parse symmetry operation.

        Args:
            tag (str): symmetry operation.
            hexagonal (bool, optional): convert to cartesian coordinate for hexagonal systems ?

        Returns:
            - (dict): info, (mirror:bool, inversion:bool, n:int, axis:str, det:int, partial_trans:str, zero_trans:bool, point_group:bool, time_reversal:bool).

        Note:
            - format of symmetry-operation string (translation is omittable).
            - inversion(-) n-fold rotation(2,3+/-,4+/-,6+/-)mirror(m) axis[n1n2n3]['] :translation[t1,t2,t3].
        """
        s = tag
        tr = tag.count("'") > 0
        s = s.replace("[", "").replace("]", "")  # remove []
        s = s.replace("'", "")  # remove '

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
            axis = TagSymmetryOperation.parse_axis(s[pos:])

        # partial translation.
        if len(st) > 1:
            point_group = False
            t = st[1].split(",")
            t = ("[" + ",".join(t) + "]").replace(" ", "")
        else:
            point_group = True
            t = "[0,0,0]"
        zero_trans = t == "[0,0,0]"

        if inversion or mirror:
            det = -1
        else:
            det = 1

        if hexagonal:
            axis = str_to_sympy(axis)
            axis = str([sp.simplify(axis[0] - axis[1] / 2), sp.simplify(sp.sqrt(3) * axis[1] / 2), axis[2]]).replace(" ", "")
            t = str_to_sympy(t)
            t = str([sp.simplify(t[0] - t[1] / 2), sp.simplify(sp.sqrt(3) * t[1] / 2), t[2]]).replace(" ", "")

        d = {
            "mirror": mirror,
            "inversion": inversion,
            "n": n,
            "axis": axis,
            "det": det,
            "partial_trans": t,
            "zero_trans": zero_trans,
            "point_group": point_group,
            "time_reversal": tr,
        }

        return d

    # ==================================================
    @staticmethod
    def str(info, hexagonal=False):
        """
        Convert to string.

        Args:
            info (dict): symmetry operation info.
            hexagonal (bool, optional): convert from cartesian coordinate for hexagonal systems ?

        Returns:
            - (str) -- in string format.
        """
        d = {  # default.
            "mirror": False,
            "inversion": False,
            "n": 1,
            "axis": "[0,0,1]",
            "det": 1,
            "partial_trans": "[0,0,0]",
            "zero_trans": True,
            "point_group": True,
            "time_reversal": False,
        }
        d.update(info)

        mirror = d["mirror"]
        inversion = d["inversion"]
        n = d["n"]
        axis = d["axis"]
        t = d["partial_trans"]
        pg = d["point_group"]
        tr = d["time_reversal"]

        if hexagonal:
            axis = str_to_sympy(axis)
            axis = str([sp.simplify(axis[0] + axis[1] / sp.sqrt(3)), sp.simplify(2 * axis[1] / sp.sqrt(3)), axis[2]]).replace(
                " ", ""
            )
            t = str_to_sympy(t)
            t = str([sp.simplify(t[0] + t[1] / sp.sqrt(3)), sp.simplify(2 * t[1] / sp.sqrt(3)), t[2]]).replace(" ", "")

        s = ""
        if mirror:
            s += "m"
        else:
            if inversion:
                s += "-"
            s += str(abs(n))
        if tr:
            s += "'"
        if abs(n) > 2:
            s += "+" if n > 0 else "-"
        if n != 1 or mirror:
            s += axis.replace(",", "")
        if not pg:
            s += ":" + t

        return s

    # ==================================================
    @staticmethod
    def latex(tag):
        """
        Convert to LaTeX.

        Args:
            tag (str): symmetry operation.

        Returns:
            - (str) -- in LaTeX format, without $.
        """
        d = TagSymmetryOperation.parse(tag)

        mirror = d["mirror"]
        inversion = d["inversion"]
        n = d["n"]
        axis = d["axis"]
        t = d["partial_trans"]
        pg = d["point_group"]
        t0 = d["zero_trans"]
        tr = d["time_reversal"]

        astr = axis.replace("[", "").replace("]", "").replace(",", "")
        s = ""
        if not pg:
            s += "\\{"
        if mirror:
            s += r"{\rm m}_{" + astr + "}"
        else:
            if inversion:
                s += "-"
            an = abs(n)
            s += str(an)
            us = False
            if an > 2:
                s += "^{+}" if n > 0 else "^{-}"
                us = True
            if an != 1:
                if us:
                    s += "_{\\,\\," + astr + "}"
                else:
                    s += "{}_{" + astr + "}"

        if tr:
            s += "{}'"

        if not pg:
            s += "|"
            if t0:
                s += "0"
            else:
                t = t.strip("[]").split(",")
                t = [sp.latex(sp.S(i)) for i in t]
                s += " ".join(t)
            s += "\\}"

        return s


# ==================================================
class TagIrrep:
    # ==================================================
    @staticmethod
    def parse(tag):
        """
        Parse irrep.

        Args:
            tag (str): irrep.

        Returns:
            - (dict) -- info, (tag:str, head:str, dimension:int, parity:str, complex:str, subscript:str, superscript:str).
        """
        dim_dict = {"A": 1, "B": 1, "E": 2, "T": 3}

        head = tag[0]
        sub = tag[1:]
        sub = sub.replace("'", "")
        np = tag.count("'")
        sup = "'" * np
        d = dim_dict[head]

        if tag.count("a"):
            c = "a"
            d = 1
        elif tag.count("b"):
            c = "b"
            d = 1
        else:
            c = ""
        p = ""
        if tag.count("g"):
            p = "even"
        elif tag.count("u"):
            p = "odd"

        d = {
            "tag": tag,
            "dimension": d,
            "parity": p,
            "complex": c,
            "head": head,
            "subscript": sub,
            "superscript": sup,
        }

        return d

    # ==================================================
    @staticmethod
    def str(info):
        """
        Convert to string.

        Args:
            info (dict): irrep. info.

        Returns:
            - (str) -- in string format.
        """
        d = {  # default.
            "tag": "A",
            "dimension": 1,  # dimension.
            "parity": "",  # null/even/odd.
            "complex": "",  # null/a/b.
            "head": "A",
            "subscript": "",
            "superscript": "",
        }
        d.update(info)

        s = d["head"] + d["subscript"] + d["superscript"]

        return s

    # ==================================================
    @staticmethod
    def latex(tag):
        """
        Convert to LaTeX.

        Args:
            tag (str): irrep.

        Returns:
            - (str) -- in LaTeX format, without $.
        """
        info = TagIrrep.parse(tag)

        sup = info["superscript"].replace("'", r"\prime")
        sub = info["subscript"]
        if sub.count("a"):
            sup += "(a)"
            sub = sub.replace("a", "")
        elif sub.count("b"):
            sup += "(b)"
            sub = sub.replace("b", "")

        s = info["head"]
        if sub:
            s += "_{" + sub + "}"
        if sup:
            s += "^{" + sup + "}"

        return s


# ==================================================
class TagMultipole:
    s_dict = {0: "s", 1: "p", 2: "d", 3: "f"}
    t_dict = {"Q": "electric", "G": "electric", "T": "magnetic", "M": "magnetic"}
    p_dict = {"Q": "polar", "G": "axial", "T": "polar", "M": "axial"}
    d_default = {
        "X": "Q",
        "l": 0,
        "s": 0,
        "k": 0,
        "x": "q",
        "Gamma": "",
        "n": -1,
        "p": -1,
        "component": None,
        "t_type": t_dict["Q"],
        "p_type": p_dict["Q"],
        "mp_type": "spherical",
    }

    # ==================================================
    @classmethod
    def parse(cls, idx, mp_type, component=None):
        """
        Parse multipole.

        Args:
            idx (tuple): multipole index.
            mp_type (str): multipole type, "spherical/point_group".
            component (int, optional): component.

        Returns:
            - (dict) -- info, (X:str, l:int, s:int, k:int, x:str, Gamma:str, n:int, p:int, component:int, t_type:str, p_type:str, mp_type: str).

        Notes:
            - spherical: (X, l, s, k, x), [type, rank, internal rank, internal component, internal type].
            - point_group: (X, l, Gamma, n, p, s, k, x), [type, rank, irrep., multiplicity, additional multiplicity, internal rank, internal component, internal type].
        """
        d = TagMultipole.d_default

        n_idx = len(idx)
        if mp_type == "spherical":
            if component is None:
                component = 0
            if n_idx == 2:  # (X,l).
                X, l = idx
                s, k, x = 0, 0, "q"
            elif n_idx == 4:  # (X,l,s,k).
                X, l, s, k = idx
                x = "q"
            elif n_idx == 5:
                X, l, s, k, x = idx
            else:
                raise Exception(f"invalid index (spherical), {idx}.")
            d.update({"X": X, "l": l, "s": s, "k": k, "x": x})
        elif mp_type == "point_group":
            if component is None:
                component = -1
            if n_idx == 5:  # (X,l,Gamma,n,p).
                X, l, Gamma, n, p = idx
                s, k, x = 0, 0, "q"
            elif n_idx == 7:  # (X,l,Gamma,n,p,s,k).
                X, l, Gamma, n, p, s, k = idx
                x = "q"
            elif n_idx == 8:
                X, l, Gamma, n, p, s, k, x = idx
            else:
                raise Exception(f"invalid index (point group), {idx}.")
            if Gamma[0] in ["A", "B"]:
                component = -1
            d.update({"X": X, "l": l, "Gamma": Gamma, "n": n, "p": p, "s": s, "k": k, "x": x})
        else:
            raise Exception(f"invalid multipole type, {mp_type}.")

        X = d["X"]
        d.update({"t_type": cls.t_dict[X], "p_type": cls.p_dict[X], "component": component, "mp_type": mp_type})

        return d

    # ==================================================
    @staticmethod
    def str(info):
        """
        Convert to string.

        Args:
            info (dict): multipole info.

        Returns:
            - (str) -- multipole in str.
        """
        d = TagMultipole.d_default
        d.update(info)

        X, l, s, k, x, Gamma, n, p, component, t_type, p_type, mp_type = list(info.values())

        txt = ""
        txt += X
        txt += "(" + f"{l}"

        if mp_type == "spherical":
            txt += f",{component}" + ")"

            if s > 0:
                txt += "^("
                if s > 0:
                    txt += f"{s},{k}"
                txt += ")"
        else:
            if component != -1:
                txt += f",{component}"
            txt += ")"

            if s > 0:
                txt += "^("
                if s > 0:
                    txt += f"{s},{k}"
                txt += ")"

            txt += "(" + Gamma

            if not (n == -1 and p == -1):
                p = "" if p == -1 else chr(ord("a") + p - 1)
                if n != -1:
                    txt += f",{n}{p}"
                else:
                    txt += f",{p}"

            txt += ")"

        return txt

    # ==================================================
    @staticmethod
    def latex(idx, mp_type, component=None, tag="", vector=True, internal=False):
        """
        Convert to LaTeX.

        Args:
            idx (tuple): multipole index.
            mp_type (str): multipole type, "spherical/point_group".
            component (int, optional): component.
            tag (str, optional): superscript.
            vector (bool, optional): with vector arrow ?
            internal (bool, optional): add internal info. ?

        Returns:
            - (str) -- in LaTeX format, without $.
        """
        info = TagMultipole.parse(idx, mp_type, component)

        X, l, s, k, x, Gamma, n, p, component, t_type, p_type, mp_type = list(info.values())

        latex = ""
        if vector and s > 0:
            latex += r"\vec{"
        latex += r"\mathbb{" + X + "}"
        if vector and s > 0:
            latex += "}"
        latex += "_{" + f"{l}"

        if mp_type == "spherical":
            latex += f",{component}" + "}"

            if tag != "" or s > 0:
                latex += "^{("
                if s > 0:
                    latex += f"{s},{k}"
                if tag != "":
                    if s > 0:
                        latex += ";"
                    latex += tag
                latex += ")}"

            if s > 0 and internal:
                latex += f"[{x}]"
        else:
            if component != -1:
                latex += f",{component}"
            latex += "}"

            if tag != "" or s > 0:
                latex += "^{("
                if s > 0:
                    latex += f"{s},{k}"
                if tag != "":
                    if s > 0:
                        latex += ";"
                    latex += tag
                latex += ")}"

            if s > 0 and internal:
                latex += f"[{x}]"

            latex += "(" + TagIrrep.latex(Gamma)

            if mp_type == "point_group":
                if not (n == -1 and p == -1):
                    p = "" if p == -1 else chr(ord("a") + p - 1)
                    if n != -1:
                        latex += f",{n}{p}"
                    else:
                        latex += f",{p}"

            latex += ")"

        return latex


# ==================================================
class TagBasis:
    # ==================================================
    @staticmethod
    def parse(tag):
        """
        Parse basis.

        Args:
            idx (str): tag.

        Returns:
            - (dict) -- info, (spinless:bool, spherical:bool, orbital(M/gamma/J):str, spin(spin/J_M):str).
        """
        spinless = tag.count("(") == 0
        if spinless:
            spherical = not (tag[0] in ["s", "p", "d", "f"])
            if spherical:
                orbital = tag
            else:
                orbital = tag[1:]
            spin = None
        else:
            v = tag.strip("()").split(",")
            spherical = not (v[0][0] in ["s", "p", "d", "f"])
            spin_block = v[1][0] in ["u", "d"]
            if spin_block:
                spin = v[1]
                if spherical:
                    orbital = v[0]
                else:
                    orbital = v[0][1:]
            else:
                orbital = v[0]
                spin = v[1]

        d = {"spinless": spinless, "spherical": spherical, "orbital": orbital, "spin": spin}

        return d

    # ==================================================
    @staticmethod
    def str(info, rank):
        """
        Convert to string.

        Args:
            info (dict): basis info.
            rank (int): rank.

        Returns:
            - (str) -- basis tag.
        """
        s_rank = {0: "s", 1: "p", 2: "d", 3: "f"}
        if info["spinless"]:
            if info["spherical"]:
                s = info["orbital"]
            else:
                s = s_rank[rank] + info["orbital"]
        else:
            if info["spherical"]:
                s = f"({info["orbital"]},{info["spin"]})"
            else:
                orb = s_rank[rank] + info["orbital"]
                s = f"({orb},{info["spin"]})"

        return s

    # ==================================================
    @staticmethod
    def latex(tag, rank, ket=True):
        """
        Convert to LaTeX.

        Args:
            tag (str): basis tag.
            rank (int): rank.
            ket (bool, optional): ket ? otherwise bra.

        Returns:
            - (str) -- in LaTeX format, without $.
        """
        info = TagBasis.parse(tag)

        s_rank = {0: "s", 1: "p", 2: "d", 3: "f"}
        if info["spinless"]:
            if info["spherical"]:
                s = f"{s_rank[rank]},{info["orbital"]}"
            else:
                if rank == 0:
                    s = f"{s_rank[rank]}"
                else:
                    s = f"{s_rank[rank]}" + "_{" + f"{info["orbital"]}" + "}"
        else:
            spin_block = info["spin"] in ["u", "d"]
            if spin_block:
                spin = r"\uparrow" if info["spin"] == "u" else r"\downarrow"
            else:
                spin = info["spin"]
            if info["spherical"]:
                if spin_block:
                    s = f"{s_rank[rank]},{info["orbital"]},{spin}"
                else:
                    orbital = sp.latex(sp.sympify(info["orbital"]))
                    spin = sp.latex(sp.sympify(info["spin"]))
                    s = f"{orbital},{spin};{s_rank[rank]}"
            else:
                orb = f"{s_rank[rank]}" + "_{" + f"{info["orbital"]}" + "}"
                s = f"{orb},{spin}"

        if ket:
            s = r"\ket{" + s + "}"
        else:
            s = r"\bra{" + s + "}"

        return s
