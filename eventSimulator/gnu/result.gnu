load "myStyleEps.gnu"

set style data his

set key width -3

set label 99 "Pixels" at graph 0.05,0.9 left

set output "../eps/residual_x.eps"
set xlabel "{/Symbol d}x [{/Symbol m}m]"
if(0) plot [-100:100] \
  "<awk '{if($1==2) print}' ../out/result_pixels.dat | \
    histogram  2 -100. 100. 200 -" t "Weighted"   ls 1, \
  "<awk '{if($1==2) print}' ../out/result_pixels.dat | \
    histogram  6 -100. 100. 200 -" u ($1+1):2 t "First-last" ls 2, \
  "<awk '{if($1==2) print}' ../out/result_pixels.dat | \
    histogram 10 -100. 100. 200 -" t "Fitter"     ls 3

plot \
  "<histogram  2 -100. 100. 200 ../out/result_pixels.dat" t "Weighted"   ls 1, \
  "<histogram  6 -100. 100. 200 ../out/result_pixels.dat" t "First-last" ls 2, \
  "<histogram 10 -100. 100. 200 ../out/result_pixels.dat" t "Fitter"     ls 3

set label 99 "Strips"
set output "../eps/residual_x_strips.eps"
set xlabel "{/Symbol d}x [{/Symbol m}m]"
if(0) plot [-100:100] \
  "<awk '{if($1==3) print}' ../out/result_strips.dat | \
    histogram  2 -100. 100. 200 -" t "Weighted"   ls 1, \
  "<awk '{if($1==3) print}' ../out/result_strips.dat | \
    histogram  6 -100. 100. 200 -" u ($1+1):2 t "First-last" ls 2, \
  "<awk '{if($1==3) print}' ../out/result_strips.dat | \
    histogram 10 -100. 100. 200 -" t "Fitter"     ls 3

plot \
  "<histogram  2 -100. 100. 200 ../out/result_strips.dat" t "Weighted"   ls 1, \
  "<histogram  6 -100. 100. 200 ../out/result_strips.dat" t "First-last" ls 2, \
  "<histogram 10 -100. 100. 200 ../out/result_strips.dat" t "Fitter"     ls 3

set label 99 "Pixels"
set output "../eps/residual_y.eps"
set xlabel "{/Symbol d}y [{/Symbol m}m]"
plot \
  "<histogram  3 -100. 100. 200 ../out/result_pixels.dat" t "Weighted"   ls 1, \
  "<histogram  7 -100. 100. 200 ../out/result_pixels.dat" t "First-last" ls 2, \
  "<histogram 11 -100. 100. 200 ../out/result_pixels.dat" t "Fitter"     ls 3

set output "../eps/residual_para.eps"
set xlabel "{/Symbol d}_{para} [{/Symbol m}m]"
plot \
  "<histogram  4 -100. 100. 200 ../out/result_pixels.dat" t "Weighted"   ls 1, \
  "<histogram  8 -100. 100. 200 ../out/result_pixels.dat" t "First-last" ls 2, \
  "<histogram 12 -100. 100. 200 ../out/result_pixels.dat" t "Fitter"     ls 3

set output "../eps/residual_perp.eps"
set xlabel "{/Symbol d}_{perp} [{/Symbol m}m]"
plot \
  "<histogram  5 -100. 100. 200 ../out/result_pixels.dat" t "Weighted"   ls 1, \
  "<histogram  9 -100. 100. 200 ../out/result_pixels.dat" t "First-last" ls 2, \
  "<histogram 13 -100. 100. 200 ../out/result_pixels.dat" t "Fitter"     ls 3
