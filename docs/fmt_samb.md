<img width="128" src="multipie_logo.png">

# MultiPie

## SAMB (* only for crystal with fourier transform)
- **info**
    - **atomic** : { "M_#" : ["amp_#"] }
    - **site_cluster** : { "S_#" : ["smp_#"] }
    - **bond_cluster** : { "B_#" : ["bmp_#"] }
    - **uniform*** : { "S_#"/"B_#" : ["ump_#"] }
    - **structure*** : { "B_#" : ["kmp_#"] }
    - **Z** : { ("M_#", "S_#"/"B_#") : ["z_#"] }
    - **version** : MultiPie version
    - **harmonics** : { head : [TagMultipole] }

- **data**
    - **atomic** : { "amp_#" : ( TagMultipole, shape, [(i, j, matrix element)] ) }
    - **site_cluster** : { "smp_#" : ( TagMultipole, [vector component] ) }
    - **bond_cluster** : { "bmp_#" : ( TagMultipole, [vector component] ) }
    - **uniform** : { "ump_#" : ( TagMultipole, shape, [(i, j, matrix element)] ) }
    - **structure*** : { "kmp_#" : (TagMultipole, "structure factor") }
    - **Z** : {"z_#" : ( TagMultipole, [(coeff, "amp_#", "smp_#"/"bmp_#/ump_#")] ) }
    - **Zk*** : {"z_#" : ( TagMultipole, [(coeff, "amp_#", "ump_#", "kmp_#")] ) }
