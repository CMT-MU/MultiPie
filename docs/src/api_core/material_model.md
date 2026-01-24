# Material Model

This class has the following dict data.

- General info.
  - **model**: model name.
  - **version**: MultiPie version.
  - **created**: created time stamp.

- Structural info.
  - **group**: group tag.
  - **crystal**: crystal class.
  - **cell**: unit-cell parameters.
  - **cell_info**: unit-cell information (parameters, volume, A, G matrices).
  - **unit_vector**: unit vectors (conventional cell).
  - **unit_vector_primitive**: unit vectors (primitive cell).

- Generation condition.
  - **SAMB_select**: SAMB selection.
  - **atomic_select**: atomic SAMB selection.
  - **site_select**: site-cluster SAMB selection.
  - **bond_select**: bond-cluster SAMB selection.
  - **max_neighbor**: max. neighbor in search of bonds.
  - **search_cell_range**: cell range in search of bonds.
  - **toroidal_priority**: toroidal priority, false=Q-G or true=G-Q.

- Site and bond info.
  - **site_grid**: site grid.
  - **site**: represetative and conventional unit-cell sites.
  - **bond**: represetative and conventional unit-cell bonds.

- SAMB.
  - **SAMB_number**: total number of SAMBs.
  - **SAMB_number_min**: total number of SAMBs excluding duplicated ones.
  - **atomic_samb**: atomic SAMB.
  - **atomic_id**: atomic SAMB ID.
  - **cluster_samb**: site/bond cluster SAMB.
  - **cluster_id**: site/bond cluster SAMB ID.
  - **combined_samb**: combined SAMB.
  - **combined_id**: combined SAMB ID.
  - **common_id**: common ID of combined SAMB.
  - **cluster_representative**: cluster-name, ID, and representative for each combined SAMB.

- Misc.
  - **basis_type**: atomic basis type.
  - **full_matrix**: "ket": ket=[(atom,SL,rank,orbital)], "index"={(atom,SL,rank): (top,size)}.
  - **braket**: dict from site_bond tag to list of BraketInfoType.
  - **wyckoff**: dict from site_bond tag to wyckoff.
  - **qtdraw_prop**: QtDraw properties.
  - **pdf_ctrl**: PDF creation control.

## MaterialModel Class

```{eval-rst}
.. automodule:: multipie.core.material_model
```
