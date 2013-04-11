#include "../interface/MostProbable.h"
#include "../interface/Levels.h"
#include "../../DataFormats/interface/Particles.h"

#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

#define sqr(x) ((x) * (x))

/*****************************************************************************/
MostProbable::MostProbable()
{
  // Silicon, MeV, cm
  I  = 16 * pow(Z,0.9) * 1e-6; // MeV

  depth = l0; // cm;

  me = mass[elec] * 1e+3;      // GeV -> MeV

  cerr << "  [MostProbable] reading density correction (SternHeimer84)..";

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

  cerr << " [done]" << endl;
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
  double beta = bg/sqrt(sqr(bg) + 1);

  // Xi
  double xi = (K*1e+3)/2 * Z/A * (depth * rho)/sqr(beta);

  // Most probable
  double mp = xi/depth * (log(2*me*xi * sqr(bg/I)) + 0.200 - sqr(beta) - delta);

  return mp; // MeV/cm
}

/*****************************************************************************/
double MostProbable::getValue(double p, int pid)
{
  return value(p / mass[pid]);
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
  double beta = bg/sqrt(sqr(bg) + 1);

  // Tmax
  double Tmax = 2*me*sqr(bg); // approximate me/M terms neglected

  // Continuous loss
  double dEdx = (K*1e+3) * Z/A * rho / sqr(beta) *
               (1./2 * log(2*me*sqr(bg/I)*Tmax) - sqr(beta) - delta/2);

  return dEdx;
}

/*****************************************************************************/
double MostProbable::dpdx(const double & bg) // dp/dx
{
  // Beta
  double beta = bg/sqrt(sqr(bg) + 1);

  return dEdx(bg)/beta;
}

/*****************************************************************************/
// for depositMap and calibGain
int MostProbable::guessPid(double p, double y)
{
  int pid = unknown;

  double eps[nParts];

  double ylim = 3.20;

  for(int i = elec; i <= prot; i++)
    eps[i] = value(p/mass[i]);

  if(    y < (eps[elec] + eps[pion])/2 && p < 0.16)       pid = elec;
  else
  {
    if(  y < (eps[pion] + eps[kaon])/2 || y < ylim || p > 2)
    {
      if(p > mass[pion]) pid = pion; // only if p/m > 1
    }
    else
    {
      if(y < (eps[kaon] + eps[prot])/2) { if(p < 0.70) pid = kaon; }
                                   else { if(p < 1.40) pid = prot; }
    }
  }

  return pid;
}

/*****************************************************************************/
// for V0 selection
int MostProbable::surePid(double p, double y)
{
  int pid = unknown;

  double eps[nParts];

  for(int i = elec; i <= prot; i++)
    eps[i] = value(p/mass[i]);

  if(y < (eps[pion] + eps[kaon])/2)
  { 
    if(p > 0.16 && p < 0.70) pid = pion;
  }
  else
  {
    if(y < (eps[kaon] + eps[prot])/2) { if(p < 0.70) pid = kaon; }
                                 else { if(p < 1.40) pid = prot; }
  }

  return pid;
}

