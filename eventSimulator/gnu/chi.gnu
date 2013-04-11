load "myStyleEps.gnu"

set style data his
set log y

set output "../eps/chi.eps"

set xlabel "{/Symbol c}"
set log y ; set format y "10^{%T}"

set xrange [0:15]

set label 1 "" at graph 0.7,0.7 center

set output "../eps/chi.eps"
set label 1 "p < 0.4 GeV/c"
plot \
  "<awk '{if($7<0.4 && ($8== 211 || $8== -211)) print}' ../out/chi.dat | \
    histogram 2 0. 15. 75 -" t "{/Symbol p}" ls 1, \
  "<awk '{if($7<0.4 && ($8== 321 || $8== -321)) print}' ../out/chi.dat | \
    histogram 2 0. 15. 75 -" t "K" ls 2, \
  "<awk '{if($7<0.4 && ($8==2212 || $8==-2212)) print}' ../out/chi.dat | \
    histogram 2 0. 15. 75 -" t "p" ls 3

set output "../eps/chi_1.eps"
set label 1 "0.4 < p < 0.8 GeV/c"
plot \
  "<awk '{if($7>0.4 && $7<0.8 && ($8== 211 || $8== -211)) print}' ../out/chi.dat | \
    histogram 2 0. 15. 75 -" t "{/Symbol p}" ls 1, \
  "<awk '{if($7>0.4 && $7<0.8 && ($8== 321 || $8== -321)) print}' ../out/chi.dat | \
    histogram 2 0. 15. 75 -" t "K" ls 2, \
  "<awk '{if($7>0.4 && $7<0.8 && ($8==2212 || $8==-2212)) print}' ../out/chi.dat | \
    histogram 2 0. 15. 75 -" t "p" ls 3

set output "../eps/chi_2.eps"
set label 1 "p > 0.8 GeV/c"
plot \
  "<awk '{if($7>0.8 && ($8== 211 || $8== -211)) print}' ../out/chi.dat | \
    histogram 2 0. 15. 75 -" t "{/Symbol p}" ls 1, \
  "<awk '{if($7>0.8 && ($8== 321 || $8== -321)) print}' ../out/chi.dat | \
    histogram 2 0. 15. 75 -" t "K" ls 2, \
  "<awk '{if($7>0.8 && ($8==2212 || $8==-2212)) print}' ../out/chi.dat | \
    histogram 2 0. 15. 75 -" t "p" ls 3
