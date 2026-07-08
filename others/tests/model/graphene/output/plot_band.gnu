unset key 
unset grid 
lwidth = 3 
set xrange [:3.966982270564949] 
set yrange [-1.4696938456699067:1.4696938456699067] 
set tics font 'Times New Roman, 30' 

set size ratio 0.7 

set arrow from  0,  -1.4696938456699067 to 0, 1.4696938456699067 nohead 
set arrow from  1.2901817879218855,  -1.4696938456699067 to 1.2901817879218855, 1.4696938456699067 nohead 
set arrow from  2.2855900490007226,  -1.4696938456699067 to 2.2855900490007226, 1.4696938456699067 nohead 
set arrow from  3.966982270564949,  -1.4696938456699067 to 3.966982270564949, 1.4696938456699067 nohead 
set xtics ('Γ' 0,'M' 1.2901817879218855,'K' 2.2855900490007226,'Γ' 3.966982270564949,) 

ef = 0.0 
set terminal postscript eps color enhanced 

set output 'graphene.eps' 

plot 'graphene.txt' u 1:2 w l lw lwidth dt (3,1) lc 'salmon', 0.0 lw 0.5 dt (2,1) lc 'black' 

set terminal pdf 

set output 'graphene.pdf' 

plot 'graphene.txt' u 1:2 w l lw lwidth dt (3,1) lc 'salmon', 0.0 lw 0.5 dt (2,1) lc 'black'