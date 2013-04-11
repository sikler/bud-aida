set term post eps enh color dashed "Helvetica" 25
set bar small
set tics scale 2,1
set style data p
set key top right Left reverse samplen 2 width +1 box noauto

set ylabel offset 1,0

set size 1,1.2

set pointsize 2

set mxtics 5
set mytics 5

set lmargin 8
set rmargin 1
set tmargin 1
set bmargin 3.5

set style line 1 lt 1 lw 3 lc  1 pt 6
set style line 2 lt 2 lw 3 lc  rgb "dark-green" pt 8
set style line 3 lt 3 lw 3 lc  3 pt 4
set style line 4 lt 4 lw 3 lc  4 pt 12

set fit errorvariables
