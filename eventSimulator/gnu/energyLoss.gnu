load "myStyleEps.gnu"

set log x
set log y

set xlabel "p [GeV/c]"
set ylabel "log({/Symbol e}/[MeV/cm]) [for 450 {/Symbol m}m Si]" offset 1,0

set xrange [exp(-4):exp(4)]
set yrange [exp(0):exp(3.5)]

load "epsilon.gnu"

set key noopaque

set colorbox front
color1(gray) = sqrt(gray)     + exp(-gray/0.01)
color2(gray) = gray**3        + exp(-gray/0.01)
color3(gray) = sin(2*pi*gray) + exp(-gray/0.01)
set palette model RGB functions color1(gray), color2(gray), color3(gray)

set lmargin 1
set rmargin 0
set bmargin 1
set tmargin 0

set pm3d corners2color c1 explicit map

set parametric
set urange [-4:4]

set samples 1000

lw = 1

set output "../eps/energyLoss_noOverflow.eps"
splot "<awk '{if($9==0) print log($2),log($3)}' \
        ../out/estimates.dat | \
        histogram 1 -4. 4. 100 2 0. 3.5 70 -" \
  u (exp($1)):(exp($2)):($3) w pm3d, \
     exp(u),epsilon(exp(u)/mpi),1 w l t "{/Symbol p}" lt 1 lw lw, \
     exp(u),epsilon(exp(u)/mka),1 w l t "K"           lt 2 lc rgb "dark-green" lw lw, \
     exp(u),epsilon(exp(u)/mpr),1 w l t "p"           lt 3 lw lw, \
     exp(u),epsilon(exp(u)/mel),1 w l t "e"           lt 4 lw lw

set output "../eps/energyLoss_Overflow.eps"
splot "<awk '{if($9==1) print log($2),log($3)}' \
        ../out/estimates.dat | \
        histogram 1 -4. 4. 100 2 0. 3.5 70 -" \
  u (exp($1)):(exp($2)):($3) w pm3d, \
     exp(u),epsilon(exp(u)/mpi),1 w l t "{/Symbol p}" lt 1 lw lw, \
     exp(u),epsilon(exp(u)/mka),1 w l t "K"           lt 2 lc rgb "dark-green" lw lw, \
     exp(u),epsilon(exp(u)/mpr),1 w l t "p"           lt 3 lw lw, \
     exp(u),epsilon(exp(u)/mel),1 w l t "e"           lt 4 lw lw

set output "../eps/energyLoss.eps"
splot "<awk '{print log($2),log($3)}' \
        ../out/estimates.dat | \
        histogram 1 -4. 4. 100 2 0. 3.5 70 -" \
  u (exp($1)):(exp($2)):($3) w pm3d, \
     exp(u),epsilon(exp(u)/mpi),1 w l t "{/Symbol p}" lt 1 lw lw, \
     exp(u),epsilon(exp(u)/mka),1 w l t "K"           lt 2 lc rgb "dark-green" lw lw, \
     exp(u),epsilon(exp(u)/mpr),1 w l t "p"           lt 3 lw lw, \
     exp(u),epsilon(exp(u)/mel),1 w l t "e"           lt 4 lw lw 
