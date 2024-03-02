set size ratio -1
set xrange [0:4]
set yrange [0:4]
set xlabel 'Variable 1'
set ylabel 'Variable 2'
plot 'datos.dat' using 1:2:(0.1):3 with circles lc var notitle
