print " fitMacro ",$0

set print "../out/sigmaZ_$0.dat"
fit f(x) "<awk '{if($$1==\"$0\" && $$2== 1) print}' ../out/macro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 1,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2== 2) print}' ../out/macro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 2,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2== 3) print}' ../out/macro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 3,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2== 4) print}' ../out/macro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 4,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2== 5) print}' ../out/macro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 5,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2== 6) print}' ../out/macro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 6,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2== 7) print}' ../out/macro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 7,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2== 8) print}' ../out/macro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 8,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2== 9) print}' ../out/macro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 9,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2==10) print}' ../out/macro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 10,s,s_err

set print
