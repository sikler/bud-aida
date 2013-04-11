load "myStyleEps.gnu"

hmin = 5

gmin = 0.0
gmax = 2.5

r(x) = rand(0)-0.5

set key off

###
set output "../eps/hits.eps"
set xlabel "Number of hits"
set ylabel "Counts"

set log y ; set format y "10^{%T}"
plot "<grep -v int ../out/gains.dat" u 0:4 w p ps 0.25
unset log y ; set format y

###
set output "../eps/gain.eps"

set xlabel "Gain correction"
set ylabel "Counts"

set arrow from 1   ,graph 0.95 to 1   ,graph 0.85 lt 3 lw 5
set arrow from 0.99,graph 0.95 to 0.99,graph 0.85 lt 5 lw 5

plot [gmin:gmax] \
  "<grep -v inf ../out/gains.dat | ./histogram 5 0. 2. 2000 -" w his

unset arrow

###
set output "../eps/gain_vs_nhits.eps"

set xlabel "Number of hits"
set ylabel "Gain correction"

c = 0.8

set log x
plot [hmin:][gmin:gmax] \
  "<grep -v inf ../out/gains.dat" u ($4+c*r(0)):5 pt 6 ps 0.25, \
  0.99 lt 5 lw 5

###
set output "../eps/gainUncertainty_vs_nhits.eps"

set ylabel "Uncertainty ({/Symbol s}) of the gain correction"

f(x) = a/sqrt(x)
a = 0.1

set fit quiet

fit f(x) "< grep -v inf ../out/gains.dat" u 4:6 via a

set label \
  "{/Symbol s} {/Symbol \273} %.3f/{/Symbol @{\140\140\140}\326}n_{hits}",a \
  at graph 0.7,0.7 center

plot [hmin:][0:2*f(hmin)] "" u ($4+c*r(0)):6 pt 6 ps 0.25, \
     f(x) lt 3 lw 5

unset log x
unset label

#####
set output "../eps/gainCompare.eps"

set xlabel "True gain factor"
set ylabel "Reconstructed gain correction"

set datafile missing "0"
set size square

f(x) = a*(x-1)+b

fit [:1] f(x) "<./compare.awk ../out/gains.orig ../out/gains.dat" u 1:2:3 via a,b
set label "g_{rec} = (%.2f",a,"{/Symbol \327} g_{true}{/Symbol -}1) + %.2f",b \
    at graph 0.5,0.1 center tc lt 3

plot [gmin:gmax][gmin:gmax] \
 "<./compare.awk ../out/gains.orig ../out/gains.dat" w e pt 6 ps 1, \
     f(x) lt 3 lw 3

! rm fit.log
