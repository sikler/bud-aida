#ifndef _MostProbable_h_
#define _MostProbable_h_

class MostProbable
{
 public:
  MostProbable();

  double value(const double & bg);
  double dEdx (const double & bg);
  double dpdx (const double & bg);

 private:
  double Z_A, I, rho, depth, K, me, C, x0,x1,a,k;
  double d0;
};

#endif
