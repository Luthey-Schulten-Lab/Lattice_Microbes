set y2tics 
set xlabel "Time (s)"
set ylabel "E, ES, P, Enzyme, Product Count"
set y2label "RNA Count"
set title  "Coupled RDME/CME Solution"
set key below

set terminal png 
set output 'Coupled.png'

plot 'T3.2-CMETraces.dat' u 1:2 w l lw 4 t 'E', \
     'T3.2-CMETraces.dat' u 1:4 w l lw 4 t 'ES', \
     'T3.2-CMETraces.dat' u 1:5 w l lw 4 t 'P', \
     'T3.2-RDMETraces.dat' u 1:2 axis x1y2 w l lw 4 t 'RNA', \
     'T3.2-RDMETraces.dat' u 1:3 w l lw 4 t 'Enzyme_{tot}', \
     'T3.2-RDMETraces.dat' u 1:4 w l lw 4 t 'Product'
