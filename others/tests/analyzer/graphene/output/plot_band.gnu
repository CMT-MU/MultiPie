unset key 
unset grid 
lwidth = 2 
set xrange [:3.9669822705649658] 
set yrange [-9.947911407510002:14.139327908730003] 
set tics font 'Times New Roman, 14' 

set size ratio 0.7 

set palette defined ( -1.0 "royalblue", 0 "gray90", 1.0 "salmon")
set cbrange [-1.0:1.0]

set arrow from  0,  -9.947911407510002 to 0, 14.139327908730003 nohead 
set arrow from  1.2901817879218849,  -9.947911407510002 to 1.2901817879218849, 14.139327908730003 nohead 
set arrow from  2.2855900490007173,  -9.947911407510002 to 2.2855900490007173, 14.139327908730003 nohead 
set arrow from  3.9669822705649658,  -9.947911407510002 to 3.9669822705649658, 14.139327908730003 nohead 
set xtics ("Γ" 0,"M" 1.2901817879218849,"K" 2.2855900490007173,"Γ" 3.9669822705649658,) 

ef = -0.4603 
set terminal postscript eps color enhanced 

set output 'graphene_dispersion.eps' 

plot 'graphene_dispersion.txt' u 1:2 w l lw lwidth lc 'salmon', 0.0 lw 0.5 lc 'black' 

set terminal pdf 

set output 'graphene_dispersion.pdf' 

plot 'graphene_dispersion.txt' u 1:2 w l lw lwidth lc 'salmon', 0.0 lw 0.5 lc 'black'