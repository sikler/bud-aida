#ifndef _MostProbable_h_
#define _MostProbable_h_

class MostProbable
{
 public:
  MostProbable();

  double value(const double & bg);
  double getValue(double p, int pid);

  double dEdx (const double & bg);
  double dpdx (const double & bg);

  int guessPid(double p, double y);
  int  surePid(double p, double y);

 private:
  double depth, I, me, C, x0,x1,a,k;
  double d0;
};

#endif
