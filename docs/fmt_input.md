<img width="128" src="multipie_logo.png">

# MultiPie

## Input for SAMB construction (* optional [default])
- **model** : model name (str).
- **group** : group name (Schoenflies notation) or group number (space group) (str/int).
- **site** : site position, orbital info. (dict) { name: ("position", orbital or orbital list) }.
- **bond*** : bond (list) [ ("tail", "head", (list of) neighbors) ], [[]].
- **spinful*** : spinful basis ? (bool), [False].
- **cell*** : unit-cell constants (dict) { "a", "b", "c", "alpha", "beta", "gamma" }, [a=b=c=1,alpha=beta=gamma=90].
- **option***
  - **view*** : view point (int list), [None].
  - **view_mode*** : mode for QtDraw file (str) ("standard"/"arrow"/"debug"), ["standard"].
  - **output*** : output folder (str), [model name].
  - **minimal_samb*** (bool) : minimal output in _samb.pdf ? [True].
- **generate***
  - **model_type*** : model type (str), ("tight_binding"/"phonon"), ["tight_binding"].
  - **time_reversal_type*** : time-reversal type (str), ("electric"/"magnetic"/"both"), ["electric"].
  - **irrep*** : irrep. (str list), [identity irrep.].
- **k_point*** : k-point (dict) {name: "position"}, [{ "Γ": "[0,0,0]", "X": "[1/2,0,0]" }].
- **k_path*** : k-path (str) (concatenate by "-" or "\|"), ["Γ-X"].
