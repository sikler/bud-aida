load "myStyleEps.gnu"

set key width -2

set style data his

set output "../eps/residual_charge_pixels.eps"
set label 99 "Pixels" at graph 0.05,0.9 left
set xlabel "{/Symbol d}y_{cluster} [keV]"

# removed overflows
plot \
  "<awk '{if($1>=1 && $16!=0. && $14!=0.) print $14*1e+3}' \
    ../out/result_pixels.dat | \
    histogram 1 -20. 20. 100 -" t "Sum"    ls 1, \
  "<awk '{if($1>=1 && $16!=0. && $15!=0.) print $15*1e+3}' \
    ../out/result_pixels.dat | \
    histogram 1 -20. 20. 100 -" t "Fitter" ls 3

set output "../eps/residual_charge_strips.eps"
set label 99 "Strips"

# removed overflows
plot \
  "<awk '{if($1>=1 && $16!=0. && $14!=0.) print $14*1e+3}' \
    ../out/result_strips.dat | \
    histogram 1 -20. 20. 100 -" t "Sum"    ls 1, \
  "<awk '{if($1>=1 && $16!=0. && $15!=0.) print $15*1e+3}' \
    ../out/result_strips.dat | \
    histogram 1 -20. 20. 100 -" t "Fitter" ls 3
