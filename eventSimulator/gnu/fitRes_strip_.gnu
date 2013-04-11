##############
# x direction
s = 20.
fit f(x) \
  "<awk '{if($$1==$0) print}' ../out/result_strips.dat | \
   histogram  2 -100. 100. 200 -" via a,s
sx1 = abs(s); sx1_e = s_err

s = 20.
fit f(x) \
  "<awk '{if($$1==$0) print}' ../out/result_strips.dat | \
   histogram  6 -100. 100. 200 -" via a,s
sx2 = abs(s); sx2_e = s_err

s = 20.
fit f(x) \
  "<awk '{if($$1==$0) print}' ../out/result_strips.dat | \
   histogram 10 -100. 100. 200 -" via a,s
sx3 = abs(s); sx3_e = s_err

##############
# x direction
s = 20.
fit f(x) \
  "<awk '{if($$1==$0) print}' ../out/result_strips.dat | \
   histogram  3 -100. 100. 200 -" via a,s
sy1 = abs(s); sy1_e = s_err

s = 20.
fit f(x) \
  "<awk '{if($$1==$0) print}' ../out/result_strips.dat | \
   histogram  7 -100. 100. 200 -" via a,s
sy2 = abs(s); sy2_e = s_err

s = 20.
fit f(x) \
  "<awk '{if($$1==$0) print}' ../out/result_strips.dat | \
   histogram 11 -100. 100. 200 -" via a,s
sy3 = abs(s); sy3_e = s_err

##############
# charge
m = -5.
s = 5.
fit g(x) \
  "<awk '{if($$1==$0 && $$16!=0.) print}' ../out/result_strips.dat | \
   histogram 14 -20. 20. 100 -" via a,m,s
mSu =     m ; mSu_e = m_err
sSu = abs(s); sSu_e = s_err

m = -1.
s = 5.
fit [-10:10] g(x) \
  "<awk '{if($$1==$0 && $$16!=0.) print}' ../out/result_strips.dat | \
   histogram 15 -20. 20. 100 -" via a,m,s
mFi =     m ; mFi_e = m_err
sFi = abs(s); sFi_e = s_err


print $0, sx1,sy1,     sx2,sy2,     sx3,sy3,    \
          sx1_e,sy1_e, sx2_e,sy2_e, sx3_e,sy3_e, \
          mSu, mSu_e, mFi, mFi_e, \
          sSu, sSu_e, sFi, sFi_e

