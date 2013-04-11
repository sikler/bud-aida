set output "../eps/chi_$0_$3.eps"

set label 1 ""            at graph 0.5 ,0.9  center
set label 2 "$\\eta = 0$" at graph 0.50,0.50 left
set label 3 ""            at graph 0.50,0.40 left

if("$3" eq "low") set label 1 "$4"
set label 2 ""
set label 3 ""

#if("$0" eq "CMS") set label 1 "$4"
if("$0" eq "CMS") set label 2 "$$\\eta = 0$$"
if("$0" eq "CMS") set label 3 "$$p_T = $2~\\mathrm{GeV}/c$$"

#set style data his

# central chi distribution
f_chi(x,r) = 2**(1-r/2.) * x**(r-1.) * exp(-x**2/2.)/gamma(r/2.)
# use it with rescaler
g(x) = a/b * f_chi(x/b,r)

if("$0" eq "Atlas"     ) r =  6.
if("$0" eq "Atlas_full") r = 45.
if("$0" eq "Alice"     ) r =  8.

if(("$0" eq "CMS") && ("$3" eq "low" )) r = 10.
if(("$0" eq "CMS") && ("$3" eq "high")) r = 16.

set multiplot
plot \
  "<awk '{print}' ../out/$1_$2_pi.spec | \
    histogram 2 0. 40. 200 -" ls 1 t pi_, \
  "<awk '{print}' ../out/$1_$2_ka.spec | \
    histogram 2 0. 40. 200 -" ls 2 t ka_, \
  "<awk '{print}' ../out/$1_$2_pr.spec | \
    histogram 2 0. 40. 200 -" ls 3 t pr_, \
  "<awk '{print}' ../out/$1_$2_el.spec | \
    histogram 2 0. 40. 200 -" ls 4 t el_#, \
  "<awk '{print}' ../out/$1_$2_??.spec | \
    histogram 2 0. 40. 200 -" ls 0, \

unset xtics
unset ytics
unset label
set xlabel

w(x) = sqrt(x)

a = 0.7*1e+4
b = 1.
fit g(x) "<awk '{print}' ../out/$1_$2_pi.spec | \
    histogram 2 0. 40. 200 -" u 1:2:($2 > 1 ? sqrt($2) : 1.) via a,b,r
plot g(x) ls 9

a = 0.02 * 1e+4
b = 1.
fit [][1:] g(x) "<awk '{print}' ../out/$1_$2_el.spec | \
    histogram 2 0. 40. 200 -" u 1:2:(w($2)) via a,b,r
plot g(x) ls 9

a = 0.10 * 1e+4
if("$3" eq "low" ) b = 1.5
if("$3" eq "high") b = 1.1
fit [][1:] g(x) "<awk '{print}' ../out/$1_$2_ka.spec | \
    histogram 2 0. 40. 200 -" u 1:2:(w($2)) via a,b,r
plot g(x) ls 9

a = 0.18 * 1e+4
if("$3" eq "low" ) b = 3.0
if("$3" eq "high") b = 1.5

fit [][1:] g(x) "<awk '{print}' ../out/$1_$2_pr.spec | \
    histogram 2 0. 40. 200 -" u 1:2:(w($2)) via a,b,r
plot g(x) ls 9

set nomultiplot

set xtics

! rm fit.log
