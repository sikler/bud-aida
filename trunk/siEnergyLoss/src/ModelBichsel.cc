#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

#include "../interface/ModelBichsel.h"
#include "../interface/CubicSpline.h"

#include "../../DataFormats/interface/RandomGenerator.h"

/*****************************************************************************/
ModelBichsel::ModelBichsel(const char * elem)
{
  readMeanNumberOfCollisions(elem);
  readCumulativeCrossSection(elem);

  theRandom = new RandomGenerator();
}

/*****************************************************************************/
ModelBichsel::~ModelBichsel()
{
  delete theRandom;
}

/*****************************************************************************/
void ModelBichsel::readMeanNumberOfCollisions(const char * elem)
{
  cerr << "  [ModelBichsel] reading mean number of collisions ("<<elem<<")..";

  ifstream file;
  char fileName[256];
  sprintf(fileName,
     "../../siEnergyLoss/data/meanNumberOfCollisions_%s.dat", elem);
  file.open(fileName);

  // - betaGamma M_0 M_1

  int d; double f, betaGamma, sigmaT;
  for(int i = 0; i < mncRows; i++)
  {
    file >> d;
    file >> betaGamma; mnc[0][i] = log(betaGamma); // xa
    file >> sigmaT   ; mnc[1][i] = log(sigmaT);    // ya
    file >> f;
  }
  file.close();

  CubicSpline theSpline;
  theSpline.prepare(mnc[0],mnc[1], mncRows, 1e+30, 1e+30, mnc[2]);

  cerr << " [done]" << endl;
}

/*****************************************************************************/
void ModelBichsel::readCumulativeCrossSection(const char * elem)
{
  cerr << "  [ModelBichsel] reading cumulative collision cross-section ("<<elem<<")..";

  ifstream file;
  char fileName[256];
  sprintf(fileName, 
     "../../siEnergyLoss/data/cumulativeCollisionCrossSection_%s.dat", elem);
  file.open(fileName);

  // - - energies [eV]
  // ccs [betagamma] [0=cum prob / 1=log(delta) / 2=der] [row]

  int k0 = 0;
  for(int j = 0 ; j < 2 ; j++)
  {
    for(int k = 0; k < ccsCols/2; k++)
      file >> cbg[k0 + k];

    int d; double delta, cumulativeProb;

    for(int i = 0; i < ccsRows; i++)
    {
      file >> d; file >> cumulativeProb;

      for(int k = 0; k < ccsCols/2; k++)
      {
        file >> delta;
        ccs[k0+k][0][i] = cumulativeProb; // xa
        ccs[k0+k][1][i] = log(delta);     // ya
      }
    }

    k0 += ccsCols/2;
  }

  file.close();

  CubicSpline theSpline;
  double der[ccsCols] = {1e+30, 6e+4, 7e+4, 7e+4, 7e+4,
                          7e+4, 7e+4, 7e+4, 7e+4, 7e+4};
  for(int k = 0; k < ccsCols; k++)
    theSpline.prepare(ccs[k][0],ccs[k][1], ccsRows, 1e+30, der[k], ccs[k][2]);

  cerr << " [done]" << endl;
}

/*****************************************************************************/
// thickness [um], return in MeV
double ModelBichsel::generate(double thickness)
{
  double x     = 0.; // cumulative length of segment
  double delta = 0.; // the total energy loss per segment

  CubicSpline theSpline;

  while(1)
  {
    // distance to next collision
    double Dx = -log(theRandom->getFlatRandom()) / sigmaT;

    x = x + Dx; // total length of segment so far

    if(x > thickness) break; 

    double r = theRandom->getFlatRandom();

    double f, logEloss[2];

    for(int j=0; j<2; j++)
      theSpline.interpolate(ccs[ic+j][0],ccs[ic+j][1],ccs[ic+j][2], ccsRows,
                            r, &logEloss[j],&f,&f);

    delta += exp(logEloss[0]*(1-xc) +
                 logEloss[1]*   xc);  // the energy loss in segment
  }

  return delta * 1e-6; // MeV
}

/*****************************************************************************/
// thickness [um], return in MeV
double ModelBichsel::generateDemo(double thickness, ofstream & file)
{
  double x     = 0.; // cumulative length of segment
  double delta = 0.; // the total energy loss per segment

  CubicSpline theSpline;

  while(1)
  {
    // distance to next collision
    double Dx = -log(theRandom->getFlatRandom()) / sigmaT;

    x = x + Dx; // total length of segment so far

    if(x > thickness) break;

    double r = theRandom->getFlatRandom();

    double f, logEloss[2];
    for(int j=0; j<2; j++)
      theSpline.interpolate(ccs[ic+j][0],ccs[ic+j][1],ccs[ic+j][2], ccsRows,
                            r, &logEloss[j],&f,&f);

    delta += exp(logEloss[0]*(1-xc) +
                 logEloss[1]*   xc);  // the energy loss in segment

    file << " " << x << " " << exp(logEloss[0]*(1-xc) +
                                   logEloss[1]*   xc) << endl;
  }

  return delta * 1e-6; // MeV
}

/*****************************************************************************/
void ModelBichsel::prepare(double betaGamma)
{
  // Get sigmaT
  CubicSpline theSpline;
  double f;
  theSpline.interpolate(mnc[0],mnc[1],mnc[2],mncRows,
                        log(betaGamma), &sigmaT,&f,&f);
  sigmaT = exp(sigmaT);

  // Look for row
  for(ic = 0 ; ic < ccsCols-1; ic++)
    if(betaGamma >= cbg[ic] && betaGamma < cbg[ic+1])
      break;

  // Over or under
  if(betaGamma < cbg[0])         ic = 0; 
  if(betaGamma > cbg[ccsCols-1]) ic = ccsCols-2;
  
  xc = log(betaGamma / cbg[ic]) / log(cbg[ic+1] / cbg[ic]);
}

