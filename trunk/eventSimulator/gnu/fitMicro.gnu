set print "../out/sigmaZ_$0.dat"
fit f(x) "<awk '{if($$1==\"$0\" && $$2== 5) print}' ../out/micro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print  5,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2==15) print}' ../out/micro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 15,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2==25) print}' ../out/micro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 25,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2==35) print}' ../out/micro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 35,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2==45) print}' ../out/micro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 45,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2==55) print}' ../out/micro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 55,s,s_err
fit f(x) "<awk '{if($$1==\"$0\" && $$2==65) print}' ../out/micro.dat | \
           histogram 3 -0.2 0.2 200 -" u 1:2:(w($$2)) via a,s
print 65,s,s_err
