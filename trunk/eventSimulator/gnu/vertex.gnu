load "myStyleEps.gnu"

set fit quiet

set macros
#div = "\"baseline\""
kme = "\"kMeans\""
gau = "\"GaussM\""

opt = "\"fPNN\""

fpn = "\"fPNN\""

#############################

set mxtics 2

f(x) = a * exp(-0.5 * (x/s)**2)
a = 1e+0
s = 0.1

w(x) = (x > 1 ? sqrt(x) : 1)

c = 10.

# Macro
#call "fitMacro.gnu" "div"
call "fitMacro.gnu" "kme"
call "fitMacro.gnu" "gau"
call "fitMacro.gnu" "opt"

set xlabel "Vertex multiplicity"
set ylabel "{/Symbol s}_{{/Symbol D}z} [mm]"

f(x) = a * x ** n
a = 0.2
n = -1.

print "../eps/sigmaZ_macro.eps"
set output "../eps/sigmaZ_macro.eps"
#plot [0:11][0:0.15]
plot [*:*][0:1] \
  "../out/sigmaZ_kme.dat" u 1:($2*c):($3*c) t @kme w e ls 2, \
  "../out/sigmaZ_gau.dat" u 1:($2*c):($3*c) t @gau w e ls 3, \
  "../out/sigmaZ_opt.dat" u 1:($2*c):($3*c) t @opt w e ls 4

#  "../out/sigmaZ_div.dat" u 1:($2*c):($3*c) t @div w e ls 1

unset label

w = 0.1

print "../eps/lost_macro.eps"
set output "../eps/lost_macro.eps"
set ylabel "Average fraction of lost vertex tracks"
plot [0:11][0:0.1] \
 "<awk '{a[$1,$2]+=$5; n[$1,$2]+=$2} END { \
   for(k=1; k<=10; k+=1) print k,a[\"kme\",k]/n[\"kme\",$2]}' \
   ../out/macro.dat" ls 2, "" u 1:2:(w) sm acs ls 2, \
 "<awk '{a[$1,$2]+=$5; n[$1,$2]+=$2} END { \
   for(k=1; k<=10; k+=1) print k,a[\"gau\",k]/n[\"gau\",$2]}' \
   ../out/macro.dat" ls 3, "" u 1:2:(w) sm acs ls 3, \
 "<awk '{a[$1,$2]+=$5; n[$1,$2]+=$2} END { \
   for(k=1; k<=10; k+=1) print k,a[\"opt\",k]/n[\"opt\",$2]}' \
   ../out/macro.dat" ls 4, "" u 1:2:(w) sm acs ls 4, \
  "-" t @kme w lp ls 2, \
  "-" t @gau w lp ls 3, \
  "-" t @opt w lp ls 4
1 -1
e
1 -1
e
1 -1
e

unset label
set xtics auto

unset label
unset log ; set format
set xtics auto

set key bottom right

#############################
print "../eps/multiDepEfficiency.eps"
set output "../eps/multiDepEfficiency.eps"
set style data p
set xlabel "Track multiplicity"
set ylabel "Vertexing efficiency" offset 0,0

set mxtics 2

w = 1e+2

q = 0.8

# n-bol k jo
g(n,k) = exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)) * q**k * (1-q)**(n-k)
f(n) = (n >= 2 ? 1 - g(n,0) - g(n,1) - g(n,2) : 1/0)

set samples 1000
h(n) = (n < 1.49 ? 0 : (n < 1.51 ? 1/0 : 1))

#fit [2:10] f(x) "../out/multi.dat" u 1:2:(1/$1**2) via q

#set label "$q$ = %.2f", q at graph 0.9,0.7 right
#plot [0:10.99][0:1.1]
plot [2.9:100.1][0:1.1] \
 "<awk '{if($1>1) print}' ../out/multi.dat" \
                    u 1:3 ls 2, \
                 "" u 1:4 ls 3, \
                 "" u 1:8 ls 4, \
  (x > 1.5 ? 1 : 1/0) ls 2, \
  (x > 1.5 ? 1 : 1/0) ls 3, \
  (x > 1.5 ? 1 : 1/0) ls 4, \
  "-" t @kme w lp ls 2, \
  "-" t @gau w lp ls 3, \
  "-" t @opt w lp ls 4
1 -1
e
1 -1
e
1 -1
e

#u 1:2 ls 1, f(x) ls 1

unset label
set samples 100

#############################
set key top right

print "../eps/multiDepUnefficiency.eps"
set output "../eps/multiDepUnefficiency.eps"
set style data lp
set xlabel "Track multiplicity"
set ylabel "Vertexing unefficiency"
set log y ; set format y "10^{%T}"
plot [0:110.99][1e-4:1.1] \
 "../out/multi.dat" u 1:(1-$3) t @kme w p ls 2, \
                 "" u 1:(1-$4) t @fpn w p ls 3

unset log y ; set format y

set style data p

set xrange [0:11]
set yrange [0:1]

#set xlabel "Vertex multiplicity (pile-up)"

w = 0.1

set xlabel ""

set ylabel "" offset 2,0

#############################
print "../eps/sim_0.eps"
set output "../eps/sim_0.eps"
set ylabel "Fraction of lost vertices"
set key top left
plot [][0:1] \
  "<awk '{if($1==1 && $3==0) print}' ../out/result.dat" u 2:4 \
    ls 2, \
  "<awk '{if($1==2 && $3==0) print}' ../out/result.dat" u 2:4 \
    ls 3, \
  "<awk '{if($1==3 && $3==0) print}' ../out/result.dat" u 2:4 \
    ls 4, \
  "-" t @kme w lp ls 2, \
  "-" t @gau w lp ls 3, \
  "-" t @opt w lp ls 4
1 -1
e
1 -1
e
1 -1
e


#############################
print "../eps/sim_1.eps"
set output "../eps/sim_1.eps"
set ylabel "Fraction of singly reconstructed vertices"
set key bottom left
plot [][0:1] \
  "<awk '{if($1==1 && $3==1) print}' ../out/result.dat" u 2:4 \
    ls 2, \
  "<awk '{if($1==2 && $3==1) print}' ../out/result.dat" u 2:4 \
    ls 3, \
  "<awk '{if($1==3 && $3==1) print}' ../out/result.dat" u 2:4 \
    ls 4, \
  "-" t @kme w lp ls 2, \
  "-" t @gau w lp ls 3, \
  "-" t @opt w lp ls 4
1 -1
e
1 -1
e
1 -1
e

set key top left

#############################
print "../eps/sim_2.eps"
set output "../eps/sim_2.eps"
set ylabel "Fraction of split vertices"
set xlabel "Vertex multiplicity (pile-up)"
plot [][0:1] \
  "<awk '{if($1==1 && $3==2) print}' ../out/result.dat" u 2:4 \
    ls 2, \
  "<awk '{if($1==2 && $3==2) print}' ../out/result.dat" u 2:4 \
    ls 3, \
  "<awk '{if($1==3 && $3==2) print}' ../out/result.dat" u 2:4 \
    ls 4, \
  "-" t @kme w lp ls 2, \
  "-" t @gau w lp ls 3, \
  "-" t @opt w lp ls 4
1 -1
e
1 -1
e
1 -1
e

set xlabel ""

#############################
print "../eps/rec_0.eps"
set output "../eps/rec_0.eps"
set ylabel "Fraction of fake vertices"
set key top left
plot [][0:1] \
  "<awk '{if($1==1 && $3==0) print}' ../out/result.dat" u 2:5 \
    ls 2, \
  "<awk '{if($1==2 && $3==0) print}' ../out/result.dat" u 2:5 \
    ls 3, \
  "<awk '{if($1==3 && $3==0) print}' ../out/result.dat" u 2:5 \
    ls 4, \
  "-" t @kme w lp ls 2, \
  "-" t @gau w lp ls 3, \
  "-" t @opt w lp ls 4
1 -1
e
1 -1
e
1 -1
e


#############################
print "../eps/rec_1.eps"
set output "../eps/rec_1.eps"
set ylabel "Fraction of correctly found vertices"
set key bottom right
plot [][0:1] \
  "<awk '{if($1==1 && $3==1) print}' ../out/result.dat" u 2:5 \
    ls 2, \
  "<awk '{if($1==2 && $3==1) print}' ../out/result.dat" u 2:5 \
    ls 3, \
  "<awk '{if($1==3 && $3==1) print}' ../out/result.dat" u 2:5 \
    ls 4, \
  "-" t @kme w lp ls 2, \
  "-" t @gau w lp ls 3, \
  "-" t @opt w lp ls 4
1 -1
e
1 -1
e
1 -1
e

set key top left

#############################
print "../eps/rec_2.eps"
set output "../eps/rec_2.eps"
set ylabel "Fraction of merged vertices"
set xlabel "Vertex multiplicity (pile-up)"
plot [][0:1] \
  "<awk '{if($1==1 && $3==2) print}' ../out/result.dat" u 2:5 \
    ls 2, \
  "<awk '{if($1==2 && $3==2) print}' ../out/result.dat" u 2:5 \
    ls 3, \
  "<awk '{if($1==2 && $3==2) print}' ../out/result.dat" u 2:5 \
    ls 4, \
  "-" t @kme w lp ls 2, \
  "-" t @gau w lp ls 3, \
  "-" t @opt w lp ls 4
1 -1
e
1 -1
e
1 -1
e


#############################
print "../eps/sigma.eps"
set output "../eps/sigma.eps"
f(eta,pt) = sqrt( ( 50e-4)**2 + \
                  (100e-4/pt)**2 * cosh(eta)**3) 

set xlabel "$\\eta$"
set ylabel "$p_T$ [GeV/$c$]"
set cblabel "$\\sigma_z$ [cm]"

set view map
set samples 1000
set iso 101,101

set lmargin screen 0.15
set rmargin screen 0.90
set tmargin screen 0.95
set bmargin screen 0.15

set log z  ; set format z  "$10^{%T}$"

set xrange [-2.5:2.5]
set yrange [0:2]
set contour base
set cntrparam levels incremental 1e-3,10

unset surf
splot f(x,y)

! rm fit.log
