unset key 
unset grid 
lwidth = 2 
set xrange [:9.659601854375268] 
set yrange [-1.4696938456699067:1.4696938456699067] 
set tics font 'Times New Roman, 14' 

set size ratio 0.7 

set palette defined ( -1.0 "royalblue", 0 "gray90", 1.0 "salmon")
set cbrange [-1.0:1.0]

set arrow from  0,  -1.4696938456699067 to 0, 1.4696938456699067 nohead 
set arrow from  3.141592653589796,  -1.4696938456699067 to 3.141592653589796, 1.4696938456699067 nohead 
set arrow from  5.565411782325234,  -1.4696938456699067 to 5.565411782325234, 1.4696938456699067 nohead 
set arrow from  9.659601854375268,  -1.4696938456699067 to 9.659601854375268, 1.4696938456699067 nohead 
set xtics ("Γ" 0,"M" 3.141592653589796,"K" 5.565411782325234,"Γ" 9.659601854375268,) 

ef = 0.0 
set terminal postscript eps color enhanced 

set output 'graphene_dispersion.eps' 

plot 'graphene_dispersion.txt' u 1:2 w l lw lwidth lc 'salmon', 0.0 lw 0.5 lc 'black' 

set terminal pdf 

set output 'graphene_dispersion.pdf' 

plot 'graphene_dispersion.txt' u 1:2 w l lw lwidth lc 'salmon', 0.0 lw 0.5 lc 'black'