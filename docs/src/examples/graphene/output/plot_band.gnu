Ef = 0.0  # fermi energy (data are already shifted by Ef).
Ef_line = 0.0  # position of fermi level in plot.
xmin = 0.0
xmax = 9.910784060203705
ymin = -1.4696938456699067
ymax = 1.4696938456699067
lwidth = 2

unset key
unset grid
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set tics font 'Times New Roman, 14'

set size ratio 0.7

set palette defined ( -1.0 "royalblue", 0 "gray90", 1.0 "salmon")
set cbrange [-1.0:1.0]

set arrow from 0, ymin to 0, ymax nohead
set arrow from 3.6275987324314865, ymin to 3.6275987324314865, ymax nohead
set arrow from 5.721993841688887, ymin to 5.721993841688887, ymax nohead
set arrow from 9.910784060203705, ymin to 9.910784060203705, ymax nohead
set xtics ("{/Symbol G}" 0, "M" 3.6275987324314865, "K" 5.721993841688887, "{/Symbol G}" 9.910784060203705)

set terminal postscript eps color enhanced
set output 'graphene_dispersion.eps'
plot 'graphene_dispersion.txt' u 1:2 w l lw lwidth lc 'salmon', Ef_line lw 0.5 lc 'black'

set terminal pdf enhanced
set output 'graphene_dispersion.pdf'
plot 'graphene_dispersion.txt' u 1:2 w l lw lwidth lc 'salmon', Ef_line lw 0.5 lc 'black'

