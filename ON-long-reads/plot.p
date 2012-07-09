set term png
set output 'homoref.png'
set xlabel "Accuracy (1 - Specificity)"
set ylabel "Sensitivity"
set style line 1 lt 1 pt 3 lw 3
plot 'homoref.txt' w l ls 1

set term png
set output 'het.png'
set xlabel "Accuracy (1 - Specificity)"
set ylabel "Sensitivity"
set style line 1 lt 1 pt 3 lw 3
plot 'het.txt' w l ls 1

set term png
set output 'homonull.png'
set xlabel "Accuracy (1 - Specificity)"
set ylabel "Sensitivity"
set style line 1 lt 1 pt 3 lw 3
plot 'homonull.txt' w l ls 1

