<img width="128" src="multipie_logo.png">

# MultiPie

## Molecule or Crystal Model (* only for crystal)
- **info**
    - **model** : model name.
    - **molecule** : molecule or crystal ?
    - **group** : (tag, detailed str)
    - **crystal** : crystal class
    - **cell*** : {a, b, c, alpha, beta, gamma}
    - **volume*** : unit cell volume
    - **a1*** : unit cell a1 vector
    - **a2*** : unit cell a2 vector
    - **a3*** : unit cell a3 vector
    - **option** :
        - **view** : view index
        - **view_mode** : QtDraw mode, standard/arrow/debug
        - **output** : output directory.
        - **minimal_samb** : output minimal SAMB ?
    - **generate** :
        - **model_type** : tight_binding/phonon
        - **time_reversal_type** : electric/magnetic/both
        - **irrep** : irrep list
    - **k_point*** : representative k points
    - **k_path*** : high-symmetry line in k space
    - **dimension** : dimension of full matrix
    - **spinful** : spinful or not
    - **orbital** : list of orbitals (U,D: up/down spin)
    - **ket** : ket basis list, orbital@site
    - **ket_site** : list of sites
    - **site** : input for "site" { name: (position, orbitals) }
    - **rep_site** : representative site { name: (position, wp, orbitals, site-symmetry) }
    - **cell_site** : { name_idx(pset): (position, SOs) }
    - **bond** : input for "bond" [ (tail, head, neighbors) ]
    - **rep_bond** : representative bond { name: (vector@center, wp, directional, neighbor, site-symmetry) }
    - **cell_bond** : { name_idx(pset): (vector@center, SOs) }

- **name**
    - **alias** : { cluster_tag: name or name: cluster_tag }
    - **site** : { site_tag: (name, no) }
    - **site_name** : { position : (site_tag, pset) }
    - **bond** : { bond_tag: (tail:head:n:multiplicity, no) }
    - **bond_name** : { vector@center : (bond_tag, pset) }

- **data**
    - **plus_set*** : [ plus_set list ]
    - **cluster_site** : { cluster_tag: site_list }
    - **cluster_bond** : { cluster_tag: bond_list }
    - **site** : { site_tag: (position, SO, (bra_site_no, ket_site_no)) }
    - **bond** : { bond_tag: (vector@center, SO, (bra_site_no, ket_site_no), vector, tail;head) }
    - **cluster_atomic** : { (bra_site_no, ket_site_no): [(bra_no, ket_no, matrix_tag)] }
    - **atomic_braket** : { matrix_tag : (bra_list, ket_list) }

- **detail**
    - **rep_bond_all** : { tail_head: [rep_bond in 0-9th neighbors] }
    - **cell_range*** : search range for bonds
    - **max_neighbor** : max. of neighbors to search
    - **A*** : transform matrix, [a1,a2,a3]
    - **version** : MultiPie version
