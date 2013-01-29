#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstdio>

#include <vector>

using namespace std;

#define Overflow 300000.

#define Shift 2.

/*****************************************************************************/
double getFlatRandom()
{
  return drand48();
}

/*****************************************************************************/
double getGaussRandom()
{
  return sqrt(-2*log(getFlatRandom())) * cos(M_PI*getFlatRandom());
}

/*****************************************************************************/
double getChi2(double y, double Delta)
{
/*
  double nu = 0.75;
  double s  = 1.9 + 0.102 * y;

  double t = Overflow;

  if(y < t)
  {
    if(Delta < y - nu * s) return -2*nu*(Delta - y)/s - nu*nu;
                      else return (Delta - y)*(Delta - y)/(s*s);
  }
  else
  {
    if(Delta < t + s) return -(Delta - t)/s + 1;
                 else return 0.;
  }
*/
//  return 2*fabs(y-Delta);

double s = 10.;

  double x = y - Delta;

  if(x > 0) x += Shift;
       else x -= Shift;

  return x*x / (s*s);
}

/*****************************************************************************/
double getProb(double y, double Delta)
{
/*
  double nu = 0.75;
  double s  = 1.9 + 0.102 * y;

  double chi2;
  if(Delta < y - nu * s) chi2 =  -2*nu*(Delta - y)/s - nu*nu;
                    else chi2 = (Delta - y)*(Delta - y)/(s*s);

  return exp(-chi2/2 - 0.8*log(y/Delta));
*/
double s = 10.;

  double x = y - Delta;

  if(x > 0) x += Shift;
       else x -= Shift;
  
  return exp(-x*x/2/(s*s));
}

/*****************************************************************************/
double generate(double Delta)
{
 double r;
 do
 {
   r = getGaussRandom();
 }
 while(fabs(r) < Shift);

 if(r > 0) r -= Shift;
      else r += Shift;

 double s = 10.;

 r *= s;

 r += Delta;
 

/*
  double y;

  do
  {
    y = getFlatRandom() * Delta * 10;
  }
  while(getProb(y,Delta) < getFlatRandom());

	  if(y > Overflow) y = Overflow + 1e-3;
*/
/*
  double r = -log(getFlatRandom());

  if(getFlatRandom() < 0.5) r = -r;
*/
  return r;
}

/*****************************************************************************/
double getChi2(const vector<double> & v, double x)
{
  double chi2 = 0.;

  for(vector<double>::const_iterator y = v.begin(); y!= v.end(); y++)
    chi2 += getChi2(*y, x);

  return chi2;
}

/*****************************************************************************/
int main()
{
  ofstream file("../out/distribution.dat");

  int nHits[10] = { 3, 4, 6, 8,10,
                   12,16,20,24,28 };

  for(int ii = 0; ii < 10; ii++)
  {
    int n = nHits[ii];

    for(double Delta = 100.; Delta < 101; Delta += 50)
    {
    cerr << " n = " << n << " Delta = " << Delta << " ";

    for(int i = 0; i < 1e+4; i++)
    {
      if((i % (10000/40)) == 0) cerr << ".";

      vector<double> v;

      for(int j = 0; j < n; j++)
      {
        double x = generate(Delta);
        v.push_back(x);
      }

      double minChi2 = 1e+30;
      double minx = 0.;

      for(double x = Delta/5.; x <= 5*Delta; x += 0.5)
      {
        double chi2 = getChi2(v,x);

        if(chi2 < minChi2)
         { minChi2 = chi2; minx = x; }
      }

      double sec = (getChi2(v,minx + 0.001) + getChi2(v,minx - 0.001)
              - 2 * getChi2(v,minx))/1e-3/1e-3;

      file << " " << n << " " << Delta 
           << " " << minx << " " << sqrt(2/sec)*getGaussRandom() << " " << sec
<< endl;
    }
    cerr << endl;
    }

    file << endl << endl;
  }

  file.close();

  return 0;
}
