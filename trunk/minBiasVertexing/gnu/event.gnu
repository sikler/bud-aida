set term post eps enh color solid "Helvetica" 25
set output "event.eps"

set bar small

unset arrow

load "../out/z.gnu"
load "../out/z0.gnu"

set xlabel "z [cm]"
unset key

plot [-20:20] "../out/z.dat" u 1:(rand(0)):2 w xe lw 2 pt 4, \
 "<./histogram 1 -20. 20. 200 ../out/z.dat" u 1:($2/$3) w his lw 2 lt -1
