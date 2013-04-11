load "myStyleEps.gnu"

#ww = 2.4
#set term epslatex color dashed dl 2 size 5.0/ww, 1.2*3.5/ww "" 10

#set output "../eps/chi_mass.eps"

set xlabel "$\\chi$"

set style data l

set key width -2

#set label 1 ""            at graph 0.5 ,0.9  center
#set label 2 "$\\eta = 0$" at graph 0.50,0.50 left
#set label 3 ""            at graph 0.50,0.40 left

#set lmargin 2
#set rmargin 0
#set tmargin 0
#set bmargin 2

set tmargin 0
set bmargin 1
set lmargin 3
set rmargin 0

#set mytics 2

# pt = 0.4
set xlabel
set xrange [0:21]
set yrange [0:9000]
#set ytics ("" 0, "" 500, "" 1000, "" 1500, "" 2000, "" 2500, "" 3000)
set ytics ("" 0, "" 2000, "" 4000, "" 6000, "" 8000, "" 10000)
#set mytics 2
set key off
#call "mass_.gnu" "Atlas"      "atlas_0.0"      "0.4" "low" "Exp A'"
set key on
call "mass_.gnu" "CMS"        "cms_0.0"        "0.4" "low" "Exp C"
set key off
call "mass_.gnu" "Alice"      "alice_0.0"      "0.4" "low" "Exp B"
set ytics 0,2000
call "mass_.gnu" "Atlas_full" "atlas_full_0.0" "0.4" "low" "Exp A"

# pt = 0.8
set xlabel "$\\chi$"

unset label 1

set xrange [0:11]
set yrange [0:9000]
set ytics ("" 0, "" 2000, "" 4000, "" 6000, "" 8000, "" 10000)
set key off
#set xlabel "$\\chi$"
#call "mass_.gnu" "Atlas"      "atlas_0.0"      "0.8" "high" "Exp A'"
set key on
set xlabel "$\\chi$"
call "mass_.gnu" "CMS"        "cms_0.0"        "0.8" "high" "Exp C"
set key off
set xlabel "$\\chi$"
call "mass_.gnu" "Alice"      "alice_0.0"      "0.8" "high" "Exp B"
set ytics 0,2000
set xlabel "$\\chi$"
call "mass_.gnu" "Atlas_full" "atlas_full_0.0" "0.8" "high" "Exp A"
