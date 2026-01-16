"""
Create binary data for each group and global.
"""

import sys
from multipie.core.multipie_info import __top_dir__

sys.path.append(__top_dir__)
BIN_DIR = __top_dir__ + "others/tools/binary_data"

import os
import logging
import re

from multipie.core.multipie_info import __bin_dir__
from multipie.util.util_binary import BinaryManager
from multipie.util.util import timer

# ==================================================
h_global_info = """
* Information of all groups.
- "info" (str): info.
  - "id_set" (str): (dict) info. of group's IDs.
    - "PG/SG/MPG/MSG" (str):  (dict) ID set.
      - PG: all/crystal/complex/irrep.
      - SG: all/crystal/PG.
      - MPG: all/crystal/PG/type.
      - MSG: all/crystal/PG/SG/MPG/type.
  - "tag" (str): (dict) from id to tag.
  - "id" (str):  (dict) from tag to id.
- "character" (str): (dict) info. of character.
  - "alias" (str): (dict) alias.
    - "wp/wm" (str): (sympy) definition of wp and wm, exp(2 pi i/3) or exp(-2 pi i/3).
  - "compatibility_relation" (str): (dict) compatibility relation.
    - (PG_tag1,PG_tag2) (str,str): ([(str,str)]) list of compatibility relation, (irrep1,irrep2).
- "harmonics" (str): (dict) info. of harmonics.
  - "max_rank" (str): (int) max rank of harmonics.
  - "max_s_rank" (str): (int) max internal rank.
  - "variable" (str): (ndarray(3,sympy)) coordinate variable, [x,y,z].
  - "internal_basis" (str): (Dict) internal basis.
    - (s,x) (int,str): (ndarray(2s+1,sympy)) list of basis.
  - "atomic_basis" (str): (dict) atomic basis.
    - "spinless" (str): (dict) spinless basis.
      - "lm/cubic/hexagonal" (str): (dict) basis.
        - L (int): ([str]) list of basis, M or gamma.
    - "spinful" (str): (dict) spinful basis.
      - "jml/lms/cubic/hexagonal" (str): (dict) basis.
        - L (int): ([str]) list of basis, (J,M) or (M,s) or (gamma,s).
  - "basis_function" (str): (dict) orbital basis function.
    - name (str): ((sympy,ndarray(2L+1,sympy))) cartesian expression and u-matrix, <m|gamma>.T.
  - "wannier90" (str): (dict) alias for Wannier90.
    - tag (str): (str) name corresponding to Wannier90 tag.
  - "tesseral" (str): (dict) tesseral tag.
    - name (str): ((int,int,str)) tesseral tag (L,M,"c/s").
  - "rank_name" (str): (dict) conversion between rank and name.
    - name/rank (str/int): (int/str) rank/name.
  - "orbital_decomposition" (str): (dict) spherical harmonics irrep. decomposition.
    - PG_tag (str): (dict) spherical harmonics irrep. decomposition for each PG.
      - L (int): ([(int,str)]) list of (n,irrep).
  NOTE:
    - see "multipie/data/data_definition.py" for detail.
- "response_tensor" (str): (dict) info. of response tensor.
  - "multipole" (str): (dict) multipole info.
    - "cubic/hexagonal" (str): ([str]) list of multipoles.
  - "cartesian_multipole" (str): (dict) linear combination info.
    - "cubic/hexagonal" (str): { comp (str): ([(int,str)]) }, relation between cartesian component and multipole.
- "root_cluster" (str): (dict) info. of root cluster.
  - "position" (str): (ndarray) list of root-cluster sites (cartesian).
  - "op_mapping" (str): (dict) operation-site mapping.
    - "Oh/D6h" (str): (dict) data.
      - SO_tag (str): (int) rc_site_index.
- "comment_group" (str): (dict) comment for each group.
  - key (str): (str) comment.
- "comment_group_opt" (str): (dict) optional comment for each group.
  - key (str): (str) optional comment.
"""
# ==================================================
h_each_group_pg = """
* Data for each PG.
- "info" (str): (namedtupole) info. of group.
  - tag: Schoenflies in text.
  - international: international short symbol in LaTeX.
  - schoenflies: Schoenflies symbol in LaTeX.
  - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
  - setting: setting.
  - PG: = id.
  - SG: associated SG, 1st SG in the same PGs.
  - MPG: associated MPG, 1st MPG with type II in the same PGs.
  - MSG: associated MSG, PG -> MPG -> MSG.
  - lattice: "0".
  - hexagonal_g: trigonal or hexagonal ?
  - SO: symmetry operations.
  - SO_gen: generator of SOs.
- "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
  - "tag" (str): ([str]) list of SO tag.
  - "generator" (str): ([str]) list of generator SO tag.
  - "cartesian" (str): (ndarray(n,3,3,sympy)) list of SO matrix (cartesian).
  - "fractional" (str): (ndarray(n,3,3,sympy)) list of SO matrix (fractional).
  - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
  - "so_matrix" (str): (Dict) dict of SO matrix (X="Q/G", s=0,1, cartesian).
    - (X,s) (str,int): ([ndarray(2s+1,2s+1,sympy)]) list of SO matrix.
  - "product" (str): (dict) product table.
    - (SO_tag,SO_tag) (str,str): (str) SO_tag of product.
  - "mapping" (str): (dict) mapping of SO to parent PG.
    - SO_tag (str): ([str]) list of SO_tag.
  NOTE:
    - SOs are in the same order of ITA (BCS).
    - SOs of SG and PG (without translational part) are in the same order.
    - SO matrices are in the same order of "tag".
    - monopole (s=0) and dipole (s=1) basis are [1] and [x,y,z].
- "wyckoff" (str): (dict) Wyckoff site/bond of group.
  - "site" (str): (dict) info. of Wyckoff site.
    - sw_tag (str): (dict) site Wyckoff tag.
      - "symmetry" (str): (str) site symmetry.
      - "conventional" (str):(ndarray(n,3,sympy)) list of Wyckoff site (fractional).
      - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional).
      - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
      - "bond" (str): ([str]) list of bw_tag.
  - "bond" (str): (dict) info. of Wyckoff bond.
    - bw_tag (str): (dict) bond Wyckoff tag.
      - "conventional" (str): (ndarray(n,6,sympy)) list of Wyckoff bond (vector+center) (fractional).
      - "primitive" (str): (ndarray(n,6,sympy)) list of Wyckoff bond (vector+center) (fractional).
      - "reference" (str): (ndarray(n,6,sympy)) list of representative bond (vector+center) (fractional).
      - "mapping" (str): ([[int]]) list of SO number (from 1) for each bond (opposite direction with negative sign).
  NOTE:
    - only for standard setting, (unique b axes, origin choice 2, hexagonal axes, abc).
    - Wyckoff sites are in the same order of ITA (BCS).
- "character" (str): (dict) character table of group.
  - "conjugacy" (str): ([[str]]) list of SO_tag in each conjugacy class.
  - "dimension" (str): (dict) dimension of irrep.
    - irrep_tag (str): (int) dimension.
  - "table" (str): (dict) character table (first SO in each conjugacy class).
    - irrep_tag (str): (ndarray(n,sympy) list of character.
  - "table_full" (str): (dict) character table for all SO.
    - irrep_tag (str): (ndarray(n,sympy) list of character.
  - "polar_axial_conversion" (str): (dict) conversion between polar and axial.
    - irrep_tag (str): (str) converted irrep_tag.
  - "symmetric_product" (str): (dict) symmetric product decomposition.
    - (irrep_tag,irrep_tag) (str,str): [(int,str)] list of (n, irrep_tag).
  - "anti_symmetric_product" (str): (dict) anti-symmetric product decomposition.
    - irrep_tag (str): ([(int,str)]) list of (n,irrep_tag).
  NOTE:
    - wp = exp(+2pi i/3), wm = exp(-2pi i/3).
    - for complex character, use PG_tag with "-c", otherwise without "-c", real version is given.
- "harmonics" (str): (Dict) polar (Q) and axial (G) harmonics (p=-1,s=k=0,x="q").
  - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,sympy), ndarray(2l+1,dim,sympy), ndarray(dim,sympy)) [cartesian ex.], u-matrix, <m|gamma>, and [tesseral ex.] for each component.
  NOTE:
    - u-matrix is in order of m=[l, l-1, ..., -l][gamma].
    - "-1" in multiplicity and componet represents no multiplicity and single component, respectively.
    - for complex expresson, use point group tag with "-c".
    - multipoles are sorted in order of [Gamma, l, Q/G, n, p].
- "atomic_samb" (str): (dict) atomic SAMB.
  - "jml/lgs/lg" (str): (dict) matrix of atomic multipole in (J,M;L)/(L,gama,s)/(L,gamma) basis.
    - (L1,L2) (int,int): (Dict) atomic multipole data for each bra-ket block (x="q).
      - (X,l,Gamma,n,p,s,k,x) (str,int,str,int,int,int,int,str): (ndarray(dim,2(2L1+1),2(2L2+1),sympy), ndarray(dim,sympy)) [matrix] and [cartesian ex.] for each component.
  NOTE:
    - half size of matrix for spinless "lg".
    - multipoles are sorted in order of [Q/G/T/M, s, k, l, Gamma, n, p].
- "cluster_samb" (str): (dict) site/bond cluster SAMB data.
  - "site" (str): (dict) site-cluster SAMB data.
    - sw_tag (str): (Dict) SAMB data at site Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,ns,sympy), ndarray(dim,sympy)) [SAMB] and [cartesian ex.] for each component.
  - "bond_s" (str): (dict) symmetric-bond-cluster SAMB data.
    - bw_tag (str): (Dict) SAMB data at bond Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,ns,sympy), ndarray(dim,sympy)) [SAMB] and [cartesian ex.] for each component.
  - "bond_a" (str): (dict) anti-symmetric-bond-cluster SAMB data.
    - bw_tag (str): (Dict) SAMB data at bond Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,ns,sympy), ndarray(dim,sympy)) [SAMB] and [cartesian ex.] for each component.
  - "vector" (str): (dict) bond-cluster vector SAMB data.
    - bw_tag (str): (Dict) vector SAMB data at bond Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,ns,sympy), ndarray(dim,sympy), ndarray(dim,2s+1,sympy)) [SAMB], [cartesian ex.], and [cartesian vector ex.] for each component.
  NOTE:
    - "-1" in multiplicity and componet represents no multiplicity and single component, respectively.
    - (X,l,Gamma,n) in PGMultipoleType indicates index of PG harmonics (p=-1,s=k=0,x="q").
    - multipoles are sorted in order of [Gamma, l, k, Q/G, n, p].
"""
# ==================================================
h_each_group_sg = """
* Data for each SG.
- "info" (str): (namedtuple) info.
  - tag: Schoenflies in text.
  - international: international short symbol in LaTeX.
  - schoenflies: Schoenflies symbol in LaTeX.
  - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
  - setting: setting comment.
  - PG: associated PG, unique.
  - SG: = id.
  - MPG: associated MPG, SG -> PG -> MPG.
  - MSG: associated MSG, 1st MSG with type II in the same SGs.
  - lattice: A, B, C, P, I, F, R.
  - hexagonal_g: trigonal or hexagonal ?
  - SO: symmetry operations.
  - SO_gen: generator of SOs.
- "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
  - "tag" (str): ([str]) list of SO tag.
  - "generator" (str): ([str]) list of generator SO tag.
  - "cartesian" (str): (ndarray(n,4,4,sympy)) list of SO matrix (cartesian, conventional).
  - "fractional" (str): (ndarray(n,4,4,sympy)) list of SO matrix (fractional, conventional).
  - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
  - "cartesian_primitive" (str): (ndarray(n,4,4,sympy)) list of SO matrix (cartesian, primitive).
  - "fractional_primitive" (str): (ndarray(n,4,4,sympy)) list of SO matrix (fractional, primitive).
  - "plus_set" (str): (ndarray(n,3,sympy)) list of plus set vector (fractional, conventional).
  - "so_matrix" (str): (Dict) dict of SO matrix (X="Q/G", s=0,1, cartesian).
    - (X,s) (str,int): ([ndarray(2s+1,2s+1,sympy)]) list of SO matrix.
  NOTE:
    - SOs are in the same order of ITA (BCS).
    - SOs of SG and PG (without translational part) are in the same order.
    - SO matrices are in the same order of "tag".
- "wyckoff" (str): (dict) Wyckoff site/bond of group.
  - "site" (str): (dict) info. of Wyckoff site.
    - sw_tag (str): (dict) site Wyckoff tag.
      - "symmetry" (str): (str) site symmetry.
      - "conventional" (str):(ndarray(n,3,sympy)) list of Wyckoff site (fractional, conventional, no plus set).
      - "primitive" (str): (ndarray(n,3,sympy)) list of Wyckoff site (fractional, primitive).
      - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional, conventional, plus set).
      - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
      - "bond" (str): ([str]) list of bw_tag.
  - "bond" (str): (dict) info. of Wyckoff bond.
    - bw_tag (str): (dict) bond Wyckoff tag.
      - "conventional" (str): (ndarray(n,6,sympy)) list of Wyckoff bond (vector+center) (fractional, conventional, no plus set).
      - "primitive" (str): (ndarray(n,6,sympy)) list of Wyckoff bond (vector+center) (fractional, primitive).
      - "reference" (str): (ndarray(n,6,sympy)) list of representative bond (vector+center) (fractional, conventional, plus set).
      - "mapping" (str): ([[int]]) list of SO number (from 1) for each bond (opposite direction with negative sign).
  NOTE:
    - only for standard setting, (unique b axes, origin choice 2, hexagonal axes, abc).
    - Wyckoff sites are in the same order of ITA (BCS).
- "cluster_samb" (str): (dict) site/bond cluster SAMB data.
  - "site" (str): (dict) site-cluster SAMB data.
    - sw_tag (str): (Dict) SAMB data at site Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,ns,sympy), ndarray(dim,sympy)) [SAMB] and [cartesian ex.] for each component.
  - "bond_s" (str): (dict) symmetric-bond-cluster SAMB data.
    - bw_tag (str): (Dict) SAMB data at bond Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,ns,sympy), ndarray(dim,sympy)) [SAMB] and [cartesian ex.] for each component.
  - "bond_a" (str): (dict) anti-symmetric-bond-cluster SAMB data.
    - bw_tag (str): (Dict) SAMB data at bond Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,ns,sympy), ndarray(dim,sympy)) [SAMB] and [cartesian ex.] for each component.
  - "vector" (str): (dict) bond-cluster vector SAMB data.
    - bw_tag (str): (Dict) vector SAMB data at bond Wyckoff.
      - (X,l,Gamma,n,p,s,k,x) (PGMultipoleType): (ndarray(dim,ns,sympy), ndarray(dim,sympy), ndarray(dim,2s+1,sympy)) [SAMB], [cartesian ex.], and [cartesian vector ex.] for each component.
  NOTE:
    - "-1" in multiplicity and componet represents no multiplicity and single component, respectively.
    - (X,l,Gamma,n) in PGMultipoleType indicates index of PG harmonics (p=-1,s=k=0,x="q").
    - multipoles are sorted in order of [Gamma, l, k, Q/G, n, p].
"""
# ==================================================
h_each_group_mpg = """
* Data for each MPG.
- "info" (str): (namedtuple) info.
  - tag: MPG_id (PG.no.ID).
  - international: international short symbol in LaTeX.
  - schoenflies: Schoenflies symbol in LaTeX.
  - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
  - setting: setting.
  - PG: associated PG, unique.
  - SG: associated SG, MPG -> MSG -> SG.
  - MPG: = id.
  - MSG: associated MSG, 1st MSG in the same MPGs.
  - lattice: "0".
  - hexagonal_g: trigonal or hexagonal ?
  - type: type I, II, III.
  - SO: symmetry operations.
- "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
  - "tag" (str): ([str]) list of SO tag.
  - "cartesian" (str): (ndarray(n,3,3,sympy)) list of SO matrix (cartesian).
  - "fractional" (str): (ndarray(n,3,3,sympy)) list of SO matrix (fractional).
  - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
  - "tr_sign" (str): (ndarray(n,sympy)) list of time-reversal sign of SO matrix.
  NOTE:
    - SOs are in the same order of ITA (BCS).
    - SO matrices are in the same order of "tag".
- "wyckoff" (str): (dict) Wyckoff site of group.
  - "site" (str): (dict) info. of Wyckoff site.
    - sw_tag (str): (dict) site Wyckoff tag.
      - "symmetry" (str): (str) site symmetry.
      - "conventional" (str):(ndarray(n,3,sympy)) list of Wyckoff site (fractional).
      - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional).
      - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
- "active_multipole" (list): active MPs list, [str].
"""
# ==================================================
h_each_group_msg = """
* Data for each MSG.
- "info" (str): (namedtuple) info.
  - tag: BNS_id (SG.no).
  - BNS: Belov-Neronova-Smirnova (international) notation in LaTeX.
  - OG: Opechowski-Guccione notation in LaTeX.
  - schoenflies: Schoenflies symbol in LaTeX.
  - crystal: triclinic, monoclinic, orthorhombic, tetragonal, trigonal, hexagonal, cubic.
  - setting: setting.
  - PG: associated PG, unique.
  - SG: associated SG, unique.
  - MPG: associated MPG, unique.
  - MSG: = id.
  - lattice: A, B, C, P, I, F, R.
  - hexagonal_g: trigonal or hexagonal ?
  - type: type I, II, III, IV.
  - SO: symmetry operations.
- "symmetry_operation" (str): (dict) symmetry operation (SO) of group.
  - "tag" (str): ([str]) list of SO tag.
  - "cartesian" (str): (ndarray(n,4,4,sympy)) list of SO matrix (cartesian, conventional).
  - "fractional" (str): (ndarray(n,4,4,sympy)) list of SO matrix (fractional, conventional).
  - "det" (str): (ndarray(n,sympy)) list of determinant of SO matrix.
  - "tr_sign" (str): (ndarray(n,sympy)) list of time-reversal sign of SO matrix.
  NOTE:
    - SOs are in the same order of ITA (BCS).
    - SO matrices are in the same order of "tag".
- "wyckoff" (str): (dict) Wyckoff site of group.
  - "site" (str): (dict) info. of Wyckoff site.
    - sw_tag (str): (dict) site Wyckoff tag.
      - "symmetry" (str): (str) site symmetry.
      - "conventional" (str):(ndarray(n,3,sympy)) list of Wyckoff site (fractional).
      - "reference" (str): (ndarray(n,3,sympy)) list of representative site (fractional).
      - "mapping" (str): ([[int]]) list of SO number (from 1) for each site.
"""
# ==================================================
h_each_group_opt = """
- "harmonics_multipole" (str): (Dict) multipolar harmonics (s=1,2,3) with polar (q) and axial (g) internal basis.
  - (X,l,Gamma,n,p,s,k,x) (str,int,str,int,int,int,int,str): (ndarray(dim,sympy), ndarray(dim,2s+1,sympy), ndarray(dim,sympy)) [cartesian ex.], [multipole basis ex.], and [multipole ex.] for each component.
  NOTE:
    - "-1" in multiplicity and componet represents no multiplicity and single component, respectively.
    - multipoles are sorted in order of [q/g, s, Gamma, l, Q/G, k, n, p].
- "harmonics_irrep" (str): (Dict) Q and G harmonics grouped by irrep.
  - (X,Gamma) (str,str): ([(int,int)]) list of (l,n).
- "representation_matrix" (str): (dict) representation matrix for each group.
  - "tag" (str): ([str]) list of SO_tag.
  - "matrix" (str): (Dict) representation matrix.
    - Gamma (str): (ndarray(n,dim,dim,sympy)) list of rep. matrix.
  NOTE:
    - rep_matrix is in the same order of SO.
"""


# ==================================================
def extract_first_quoted_name(line):
    match = re.search(r'"([^"]+)"', line)
    return match.group(1) if match else None


# ==================================================
def group_by_top_level_dashes_from_text(text):
    lines = text.split("\n")
    groups = []
    current_group = []

    for line in lines:
        if line.startswith("- "):
            if current_group:
                groups.append("\n".join(current_group))
            current_group = [line]
        elif current_group:
            current_group.append(line)
    if current_group:
        groups.append("\n".join(current_group))

    dic = {extract_first_quoted_name(g): g for g in groups}

    return dic


# ==================================================
@timer
def create_global_info():
    info = BinaryManager("info", topdir=BIN_DIR)
    root_cluster = BinaryManager("root_cluster", topdir=BIN_DIR)

    print(f"creating global info.")
    global_info = BinaryManager()
    global_info.add_comment(h_global_info)

    for k in info.keys():
        global_info[k] = info[k]
    global_info["root_cluster"] = dict(root_cluster.items())

    global_info["comment_group"] = {
        "PG": group_by_top_level_dashes_from_text(h_each_group_pg),
        "SG": group_by_top_level_dashes_from_text(h_each_group_sg),
        "MPG": group_by_top_level_dashes_from_text(h_each_group_mpg),
        "MSG": group_by_top_level_dashes_from_text(h_each_group_msg),
    }
    global_info["comment_group_opt"] = group_by_top_level_dashes_from_text(h_each_group_opt)

    global_info.save_binary("info")


# ==================================================
@timer
def create_each_group():
    info = BinaryManager("info", topdir=BIN_DIR)
    group = BinaryManager("group", topdir=BIN_DIR)
    atomic_multipole_group = BinaryManager("atomic_multipole_group", topdir=BIN_DIR)
    samb = BinaryManager("cluster_samb", topdir=BIN_DIR)
    harmonics = BinaryManager("harmonics", topdir=BIN_DIR)

    for id_s in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        tag = info["tag"][id_s]
        print(f"creating group data for '{tag}'.")
        each_group = BinaryManager(subdir="PG")
        each_group.add_comment(f"Data for {id_s}={tag}.")
        each_group["info"] = group[id_s]["info"]
        each_group["symmetry_operation"] = group[id_s]["symmetry_operation"]
        each_group["character"] = group[id_s]["character"]
        each_group["wyckoff"] = group[id_s]["wyckoff"]
        each_group["harmonics"] = harmonics[id_s]["harmonics"]
        if tag.count("-c") == 0:
            each_group["cluster_samb"] = samb[id_s]
            each_group["atomic_samb"] = atomic_multipole_group[id_s]
        each_group.save_binary(id_s)

    for id_s in info["id_set"]["SG"]["all"]:
        tag = info["tag"][id_s]
        print(f"creating group data for '{tag}'.")
        each_group = BinaryManager(subdir="SG")
        each_group.add_comment(f"Data for {id_s}={tag}.")
        each_group["info"] = group[id_s]["info"]
        each_group["symmetry_operation"] = group[id_s]["symmetry_operation"]
        each_group["wyckoff"] = group[id_s]["wyckoff"]
        each_group["cluster_samb"] = samb[id_s]
        each_group.save_binary(id_s)

    for id_s in info["id_set"]["MPG"]["all"]:
        tag = info["tag"][id_s]
        print(f"creating group data for '{tag}'.")
        each_group = BinaryManager(subdir="MPG")
        each_group.add_comment(f"Data for {id_s}={tag}.")
        each_group["info"] = group[id_s]["info"]
        each_group["symmetry_operation"] = group[id_s]["symmetry_operation"]
        each_group["wyckoff"] = group[id_s]["wyckoff"]
        each_group["active_multipole"] = group[id_s]["active_multipole"]
        each_group.save_binary(id_s)

    for id_s in info["id_set"]["MSG"]["all"]:
        tag = info["tag"][id_s]
        print(f"creating group data for '{tag}'.")
        each_group = BinaryManager(subdir="MSG")
        each_group.add_comment(f"Data for {id_s}={tag}.")
        each_group["info"] = group[id_s]["info"]
        each_group["symmetry_operation"] = group[id_s]["symmetry_operation"]
        each_group["wyckoff"] = group[id_s]["wyckoff"]
        each_group.save_binary(id_s)


# ==================================================
@timer
def create_each_group_opt():
    info = BinaryManager("info", topdir=BIN_DIR)
    harmonics_multipole = BinaryManager("harmonics_multipole", topdir=BIN_DIR)
    harmonics = BinaryManager("harmonics", topdir=BIN_DIR)
    representation_matrix = BinaryManager("representation_matrix", topdir=BIN_DIR)

    for id_s in info["id_set"]["PG"]["all"] + info["id_set"]["PG"]["complex"]:
        tag = info["tag"][id_s]
        print(f"creating group data (optional) for '{tag}'.")
        each_group = BinaryManager(subdir="PG")
        each_group.add_comment(f"Data for {id_s}={tag}.")

        each_group["harmonics_irrep"] = harmonics[id_s]["irrep"]
        each_group["harmonics_multipole"] = harmonics_multipole[id_s]
        each_group["representation_matrix"] = representation_matrix[id_s]

        each_group.save_binary(f"{id_s}_opt")


# ================================================== main
if __name__ == "__main__":
    logging.basicConfig(format="%(message)s", level=logging.INFO)
    os.makedirs(__bin_dir__, exist_ok=True)

    os.makedirs(__bin_dir__ + "PG", exist_ok=True)
    os.makedirs(__bin_dir__ + "SG", exist_ok=True)
    os.makedirs(__bin_dir__ + "MPG", exist_ok=True)
    os.makedirs(__bin_dir__ + "MSG", exist_ok=True)

    create_global_info()
    create_each_group()
    create_each_group_opt()
