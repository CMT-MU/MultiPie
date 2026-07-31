unset key
unset grid
lwidth = 2
set xrange [:1.0604532164016165]
set yrange [-0.5685369160876159:0.5685369160876159]
set tics font 'Times New Roman, 14'

set size ratio 0.7

set palette defined ( -1.0 "royalblue", 0 "gray90", 1.0 "salmon")
set cbrange [-1.0:1.0]

set arrow from  0,  -0.5685369160876159 to 0, 0.5685369160876159 nohead
set arrow from  0.5302266082008091,  -0.5685369160876159 to 0.5302266082008091, 0.5685369160876159 nohead
set arrow from  1.0604532164016165,  -0.5685369160876159 to 1.0604532164016165, 0.5685369160876159 nohead
set xtics ("A'" 0,"Γ" 0.5302266082008091,"A" 1.0604532164016165,)

Ef = 0.0
set terminal postscript eps color enhanced

set output 'Te_dispersion.eps'

plot 'Te_dispersion.txt' u 1:2:3 w l lw lwidth lc palette,0.0 lw 0.5 lc 'black'

set terminal pdf

set output 'Te_dispersion.pdf'

plot 'Te_dispersion.txt' u 1:2:3 w l lw lwidth lc palette,0.0 lw 0.5 lc 'black'