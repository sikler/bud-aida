load "myStyleEps.gnu"

set key width -2

set style data e

set xrange [1.5:10.5]
set yrange [0:75]

w(sig) = 5e-1/sig**2

set xlabel "n_{pixel}"
set xtics 1,1

fwhm = 2 * sqrt(2*log(2))
t(x) = fwhm * x

set output "../eps/residual_fit_para.eps"
set label 99 "Pixels" at graph 0.05,0.9 left
set key bottom right
#set ylabel "$\\sigma_\\parallel$ [$\\mu$m]"
set ylabel "FWHM of {/Symbol d}_{para} [{/Symbol m}m]"

plot \
  "../out/fitRes_pixel.dat" \
     u 1:(t($2)):(t( $8)) ls 1 t "Weighted", \
  "" u 1:(t($4)):(t($10)) ls 2 t "First-last", \
  "" u 1:(t($6)):(t($12)) ls 3 t "Fitter", \
  "<awk '{if($1>=3) print}' ../out/fitRes_pixel.dat" \
     u 1:(t($2)):(w(t( $8))) sm acs ls 1, \
  "" u 1:(t($4)):(w(t($10))) sm acs ls 2, \
  "" u 1:(t($6)):(w(t($12))) sm acs ls 3

set output "../eps/residual_fit_perp.eps"
set key top right
set label 99 at graph 0.05,0.9 left
#set ylabel "$\\sigma_\\perp$ [$\\mu$m]"
set ylabel "FWHM of {/Symbol d}_{perp} [{/Symbol m}m]"
plot \
  "../out/fitRes_pixel.dat" \
     u 1:(t($3)):(t( $9)) ls 1 t "Weighted", \
  "" u 1:(t($5)):(t($11)) ls 2 t "First-last", \
  "" u 1:(t($7)):(t($13)) ls 3 t "Fitter" , \
  "<awk '{if($1>=3) print}' ../out/fitRes_pixel.dat" \
     u 1:(t($3)):(w(t( $9))) sm acs ls 1, \
  "" u 1:(t($5)):(w(t($11))) sm acs ls 2, \
  "" u 1:(t($7)):(w(t($13))) sm acs ls 3

set output "../eps/residual_fit_strips.eps"
set key bottom right
set xlabel "n_{strip}"
set ylabel "FWHM of {/Symbol d}x [{/Symbol m}m]"

set label 99 "Strips" # at graph 0.95,0.9 right

plot \
  "../out/fitRes_strip.dat" \
     u 1:(t($2)):(t( $8)) ls 1 t "Weighted", \
  "" u 1:(t($4)):(t($10)) ls 2 t "First-last", \
  "" u 1:(t($6)):(t($12)) ls 3 t "Fitter", \
  "<awk '{if($1>=3) print}' ../out/fitRes_strip.dat" \
     u 1:(t($2)):(w(t( $8))) sm acs ls 1, \
  "" u 1:(t($4)):(w(t($10))) sm acs ls 2, \
  "" u 1:(t($6)):(w(t($12))) sm acs ls 3

##
set key width -2
set auto y
dx = 0.1

t(x) = fwhm/2 * x

set output "../eps/residual_fit_charge_pixels.eps"
#set key top left
set xlabel "$n_\\text{pixel}$"
set ylabel "$\\delta y_\\text{cluster}$ [keV]"

set label 99 "Pixels" # at graph 0.95,0.9 right

`./ellipse.awk i=0 ../out/fitRes_pixel.dat`
`./ellipse.awk i=1 ../out/fitRes_pixel.dat`

plot [][-12:6] \
  "../out/fitRes_pixel.dat" \
     u ($1-dx):(($14)):(t($18)) ls 1 t "Sum", \
  "" u ($1+dx):(($16)):(t($20)) ls 3 t "Fitter", \
  "<awk '{if($1>=3) print}' ../out/fitRes_pixel.dat" \
     u ($1-dx):(($14)):(w(t($15))) sm acs ls 1, \
  "" u ($1+dx):(($16)):(w(t($17))) sm acs ls 3

set output "../eps/residual_fit_charge_strips.eps"
#set key top left
set xlabel "$n_\\text{strip}$"
set ylabel "$\\delta y_\\text{cluster}$ [keV]"

set label 99 "Strips" # at graph 0.95,0.9 right

`./ellipse.awk i=0 ../out/fitRes_strip.dat`
`./ellipse.awk i=1 ../out/fitRes_strip.dat`

plot [][-12:6]\
  "../out/fitRes_strip.dat" \
     u ($1-dx):(($14)):(t($18)) ls 1 t "Sum", \
  "" u ($1+dx):(($16)):(t($20)) ls 3 t "Fitter", \
  "<awk '{if($1>=3) print}' ../out/fitRes_strip.dat" \
     u ($1-dx):(($14)):(w(t($15))) sm acs ls 1, \
  "" u ($1+dx):(($16)):(w(t($17))) sm acs ls 3

