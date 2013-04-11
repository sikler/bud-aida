load "myStyleEps.gnu"

#set size 0.66,1
#set lmargin 5
#set rmargin 0

set style data his

set key width -1

set xtics -100,50

set output "../eps/resolution_x_pixels.eps"
set label 99 "Pixels" at graph 0.05,0.9 left
set xlabel "{/Symbol d}x [{/Symbol m}m]"
plot [][0:] \
   "<awk '{if($1>=2 && $17!=\"nan\") print}' \
     ../out/result_pixels.dat | grep -v nan | \
     histogram 10 -100. 100. 200 -" ls 3 t "Fitter", \
   "<awk '{if($1>=2 && $17!=\"nan\") \
     print $17 * sqrt(-2*log(rand())) * cos(3.1415926*rand()) * (2*sqrt(2))}' \
     ../out/result_pixels.dat | grep -v nan | \
     histogram 1 -100. 100. 200 -" ls 4 t "Predicted"

set output "../eps/resolution_y_pixels.eps"
set xlabel "{/Symbol d}z [{/Symbol m}m]"
plot [][0:] \
   "<awk '{if($1>=2 && $18!=\"nan\") print}' \
     ../out/result_pixels.dat | grep -v nan | \
     histogram 11 -100. 100. 200 -" ls 3 t "Fitter", \
   "<awk '{if($1>=2 && $18!=\"nan\") \
     print $18 * sqrt(-2*log(rand())) * cos(3.1415926*rand()) * (2*sqrt(2))}' \
     ../out/result_pixels.dat | grep -v nan | \
     histogram 1 -100. 100. 200 -" ls 4 t "Predicted"

set output "../eps/resolution_x_strips.eps"
set label 99 "Strips" at graph 0.05,0.9 left
set xlabel "{/Symbol d}x [{/Symbol m}m]"
plot [][0:] \
   "<awk '{if($1>=2 && $17!=\"nan\") print}' \
     ../out/result_strips.dat | \
     histogram 10 -100. 100. 200 -" ls 3 t "Fitter", \
   "<awk '{if($1>=2 && $17!=\"nan\") \
     print  $17 * sqrt(-2*log(rand())) * cos(3.1415926*rand()) * (2*sqrt(2))}' \
     ../out/result_strips.dat | grep -v nan | \
     histogram 1 -100. 100. 200 -" ls 4 t "Predicted"
