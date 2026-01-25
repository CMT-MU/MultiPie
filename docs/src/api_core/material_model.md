# Material Model

## Dict data

This class has the following dict data.

- General info.
  - **model**: (str) model name.
  - **version**: (str) MultiPie version.
  - **created**: (str) created time stamp.

- Structural info.
  - **group**: (str) group tag.
  - **crystal**: (str) crystal class.
  - **cell**: (dict) unit-cell parameters.
    - a: (float) a.
    - b: (float) b.
    - c: (float) c.
    - alpha: (float) alpha (deg.).
    - beta: (float) beta (deg.).
    - gamma: (float) gamma (deg.).
  - **cell_info**: (dict) unit-cell information
    - cell: (dict) = cell.
    - volume: (float) volume.
    - A: (list) [a1,a2,a3,t] (4x4, each column).
    - G: (list) metric matrix (4x4).
  - **unit_vector**: (list) unit vectors, a1, a2, a3 (conventional cell).
  - **unit_vector_primitive**: (list) unit vectors, a1p, a2p, a3p (primitive cell).

- Generation condition.
  - **SAMB_select**: (dict) SAMB selection.
  - **atomic_select**: (dict) atomic SAMB selection.
  - **site_select**: (dict) site-cluster SAMB selection.
  - **bond_select**: (dict) bond-cluster SAMB selection.
  - **max_neighbor**: (int) max. neighbor in search of bonds.
  - **search_cell_range**: (tuple) cell range in search of bonds, (a1,a2,b1,b2,c1,c2).
  - **toroidal_priority**: (bool) toroidal priority, false=Q-G or true=G-Q.

- Site and bond info.
  - **site_grid**: (dict) site grid.
  - **site**: (dict) site info.
    - representative: (dict) representative site, dict[name, RepSiteType].
    - cell: (dict) conventional cell site, dict[name, [CellSiteType]].
  - **bond**: bond info.
    - representative: (dict) representative bond, dict[name, RepBondType].
    - cell: (dict) conventional cell bond, dict[name, [CellBondType]].
    - info: (list) bond properties, [BondInfoType].

- SAMB.
  - **SAMB_number**: (int) total number of SAMBs.
  - **SAMB_number_min**: (int) total number of SAMBs excluding duplicated ones.
  - **atomic_samb**: (dict) atomic SAMB, dict[BraketInfoType, SAMB Dict].
  - **atomic_id**: (dict) atomic SAMB ID, dict[tag, (index, component)].
  - **cluster_samb**: (dict) site/bond cluster SAMB, dict[wyckoff, SAMB Dict].
  - **cluster_id**: (dict) site/bond cluster SAMB ID, dict[tag, (wyckoff, index, component)].
  - **combined_samb**: (dict) combined SAMB, dict[SAMBType, SAMB Dict].
  - **combined_id**: (dict) combined SAMB ID, dict[tag, (symbol, UniqueSAMBType, SAMB index, component)].
  - **common_id**: (dict) common ID of combined SAMB, dict[SAMBType, ([Ztag_list], [site/bond-tag])].
  - **cluster_info**: (dict) cluster info, dict[site/bond_name, dict[(bra_rank,ket_rank), (wyckoff,z_list)] ].

- Misc.
  - **basis_type**: (str) atomic basis type.
  - **full_matrix**: (dict) full matrix info.
    - ket: (list) ket info., [(atom,sublattice,rank,orbital)].
    - index: (dict) index info., dict[(atom,SL,rank), (top,size)].
  - **braket**: (dict) braket info., dict[site/bond-tag, [BraketInfoType]].
  - **wyckoff**: (dict) wyckoff info., dict[site/bond-tag, wyckoff].
  - **qtdraw_prop**: (dict) QtDraw properties.
  - **pdf_ctrl**: (dict) PDF creation control.

## SAMB Dict data

Combined, atomic, and site/bond-cluster SAMBs have the following dict data.

- Combined SAMB, $Z_i$, Dict[index, SAMB].
  - index: (PGMultipoleType) (X,l,$\Gamma$,s,k,n,p,x)
  - SAMB: (tuple) SAMB data for set of components.
    - (list) [[(coeff, X_index, X_component, Y_index, Y_component)]],
    - (list) [symmetry].
- Atomic SAMB, $X_j$, Dict[index, SAMB].
  - index: (PGMultipoleType) (X,l,$\Gamma$,s,k,n,p,x)
  - SAMB: (tuple) SAMB data for set of components.
    - (list) [matrix],
    - (list) [symmetry]
- Site/bond-cluster SAMB, $Y_k$, Dict[index, SAMB].
  - index: (PGMultipoleType) (X,l,$\Gamma$,s,k,n,p,x)
  - SAMB: (tuple) SAMB data for set of components.
    - (list) [vector],
    - (list) [symmetry]

## MaterialModel Class

```{eval-rst}
.. automodule:: multipie.core.material_model
```
