#include "../interface/MostProbable.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

#define Sqr(x) ((x) * (x))

/*****************************************************************************/
MostProbable::MostProbable()
{
  // MeV, cm
  Z_A = 0.49848;
  I   = 173e-6;     // MeV
  rho = 2.33;       // g/cm3

  depth = 450e-4;   // cm
  K     = 0.307075; // MeV cm^2 / g
  me    = 0.511;    // MeV

  // Read density correction
  ifstream file("../../siEnergyLoss/data/density_Sternheimer84.par");

  int i = 0;

  while(!file.eof())
  {
    string s;
    file >> s; file >> s;

    float value;
    file >> value;

    switch(i)
    {
      case 0 : C  = value; break;
      case 1 : x0 = value; break;
      case 2 : x1 = value; break;
      case 3 : a  = value; break;
      case 4 : k  = value; break;
    }

    i++;
  }

  d0 = 2*M_LN10*x0 - C + a*pow(x1 - x0, k);
}

/*****************************************************************************/
double MostProbable::value(const double & bg)
{
  // Density correction
  double delta;
  double x = log10(bg);

  if(x >= x0)
    delta = 2*M_LN10*x - C + (x1 > x ? a*pow(x1 - x , k) : 0);
  else
    delta = d0 * pow(10., 2*(x - x0));

  // Beta
  double beta = bg/sqrt(Sqr(bg) + 1);

  // Xi
  double xi = K/2 * Z_A * (depth * rho)/Sqr(beta);

  // Most probable
  double mp = xi/depth * (log(2*me*xi * Sqr(bg/I)) + 0.200 - Sqr(beta) - delta);

  return mp; // MeV/cm
}

/*****************************************************************************/
double MostProbable::dEdx(const double & bg) // dE/dx
{
  // Density correction
  double delta;
  double x = log10(bg);

  if(x >= x0)
    delta = 2*M_LN10*x - C + (x1 > x ? a*pow(x1 - x , k) : 0);
  else
    delta = d0 * pow(10., 2*(x - x0));

  // Beta
  double beta = bg/sqrt(Sqr(bg) + 1);

  // Tmax
  double Tmax = 2*me*Sqr(bg); // approximate me/M terms neglected

  // Continuous loss
  double dEdx = K * Z_A * rho / Sqr(beta) *
               (1./2 * log(2*me*Sqr(bg/I)*Tmax) - Sqr(beta) - delta/2);

  return dEdx;
}

/*****************************************************************************/
double MostProbable::dpdx(const double & bg) // dp/dx
{
  // Beta
  double beta = bg/sqrt(Sqr(bg) + 1);

  return dEdx(bg)/beta; // CHECK
}
