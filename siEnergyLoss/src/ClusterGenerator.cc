#include "../interface/ClusterGenerator.h"
#include "../interface/TouchedChannels.h"
#include "../interface/Levels.h"
#include "../interface/ModelBichsel.h"

#include "../../DataFormats/interface/RandomGenerator.h"
#include "../../DataFormats/interface/TLayer.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <algorithm>

#define sqr(x) ((x) * (x))

using namespace std;

/*****************************************************************************/
ClusterGenerator::ClusterGenerator()
{
  theModel  = new ModelBichsel("Si");
  theRandom = new RandomGenerator();
}

/*****************************************************************************/
ClusterGenerator::~ClusterGenerator()
{
  delete theModel;
  delete theRandom;
}

/*****************************************************************************/
struct sortByPos
{
  bool operator() (const Pixel & a, const Pixel & b) const
  {
    return (a.meas.x[0] < b.meas.x[0]);
  }
};

/*****************************************************************************/
void ClusterGenerator::addCoupling(Hit & hit, TLayer * unit)
{
  Pixel empty;
 
  // Sort hit.allPixels according to meas.x[0] value
  sort(hit.allPixels.begin(), hit.allPixels.end(), sortByPos());

  // First 
  empty.meas.x[0] = hit.allPixels.front().meas.x[0] - 1;
  empty.meas.x[1] = hit.allPixels.front().meas.x[1];
  empty.meas.y    = 0.;
  hit.allPixels.insert(hit.allPixels.begin(), empty);

  // Last
  empty.meas.x[0] = hit.allPixels.back().meas.x[0] + 1;
  empty.meas.x[1] = hit.allPixels.back().meas.x[1];
  empty.meas.y    = 0.;
  hit.allPixels.push_back(empty);

  double c = unit->coupling;
  vector<double> coup(hit.allPixels.size(), 0.);

  int i = 1;
  for(vector<Pixel>::iterator pixel = hit.allPixels.begin() + 1;
                              pixel!= hit.allPixels.end()   - 1; pixel++)
  {
    coup[i-1] +=      c  * pixel->meas.y;
    coup[i  ] += (1-2*c) * pixel->meas.y;
    coup[i+1] +=      c  * pixel->meas.y;

    i++;
  }

  i = 0;
  for(vector<Pixel>::iterator pixel = hit.allPixels.begin();
                              pixel!= hit.allPixels.end(); pixel++ )
  {
    pixel->meas.y = coup[i];
    i++;
  }
}

/*****************************************************************************/
double ClusterGenerator::getIntegral(double x, int x0, int x1, double sigma)
{
  return fabs( erf((x-x0)/(sigma*M_SQRT2)) -
               erf((x-x1)/(sigma*M_SQRT2)) )/2;
}

/*****************************************************************************/
void ClusterGenerator::matchOrAdd(vector<Pixel> & pixels,
                                  int x, int y, double Delta)
{
  bool isFound = false;

  for(vector<Pixel>::iterator pixel = pixels.begin();
                              pixel!= pixels.end(); pixel++)
    if(pixel->meas.x[0] == x &&
       pixel->meas.x[1] == y)
    {
      pixel->meas.y += Delta;
      isFound = true;

      break;
    }

  if(!isFound)
  { 
    Pixel pixel;
    pixel.meas.x[0] = x;
    pixel.meas.x[1] = y;
    pixel.meas.y = Delta;

    Pixel::Calc calc;
    calc.isTouched = false;
    pixel.calc = calc;

    pixels.push_back(pixel);
  }
}

/*****************************************************************************/
Hit ClusterGenerator::create
  (double betaGamma, double theta, double phi, 
   TChipId & chipId, TLayer * mat, int ilayer)
{
  Hit hit;

  hit.ilayer = ilayer;
  hit.chipId = chipId;

  theModel->prepare(betaGamma);

  double dpitch[2];
  dpitch[0] = mat->pitch[2] / mat->pitch[0] * fabs(tan(phi));
  dpitch[1] = mat->pitch[2] / mat->pitch[1] * fabs(tan(theta));

  // ExB FIXME EXB
  if(theRandom->getFlatRandom() < 0.5) dpitch[0] += 0.5; 
                                  else dpitch[0] -= 0.5;

  dpitch[0] = fabs(dpitch[0]); // FIXME????

  ///////////////
  double endpoint[2][2]; // ordered

  for(int k = 0; k < 2; k++)
  {
    if(k == 0) 
      endpoint[0][k] = mat->nrows/2    + theRandom->getFlatRandom();
    else
      endpoint[0][k] = mat->ncolumns/2 + theRandom->getFlatRandom();

    endpoint[1][k] = endpoint[0][k] + dpitch[k];
  }

  ///////////////
  {
    double dpos[3];
    for(int k = 0; k < 3; k++)
      if(k < 2)
        dpos[k] = (endpoint[1][k] - endpoint[0][k]) * mat->pitch[k];
      else
        dpos[k] =                                     mat->pitch[k];

    hit.length = sqrt(sqr(dpos[0]) + sqr(dpos[1]) + sqr(dpos[2])); 
    hit.lambda = sqrt(sqr(dpos[0]) + sqr(dpos[1]));

    for(int k = 0; k < 2; k++)
      hit.ulambda[k] = dpos[k] / hit.lambda;
  }

  ///////////////
  TouchedChannels theTouchedChannels(hit, mat);
  hit.allPixels.reserve(100);
  hit.allPixels = theTouchedChannels.findChannels(endpoint);

  double simDelta = 0, recDelta = 0.;

  // drift to which side? FIXME
  bool shift = (theRandom->getFlatRandom() < 0.5);

  int nj = hit.allPixels.size();
  for(int jj = 0; jj < nj; jj++)
  {
    Pixel * pixel = &hit.allPixels[jj];
    double l = pixel->calc.l;

    // With diffusion; try to have at most 1um sections
    int n = int(l*1e+4);
    if(n == 0) n = 1;

    double dl = l / n;

    for(int i = 0; i < n; i++)
    {
      double dDelta = theModel->generate(dl*1e+4);

      double pos[3];

      pos[0] = (pixel->meas.c[0].pos[0] * (n-i - 0.5) +
                pixel->meas.c[1].pos[0] * (  i + 0.5)) / n;

      pos[1] = (pixel->meas.c[0].pos[1] * (n-i - 0.5) +
                pixel->meas.c[1].pos[1] * (  i + 0.5)) / n;

      if(shift)
      pos[2] = (mat->pitch[2] * (n-i - 0.5) +
                            0 * (  i + 0.5)) / n;
      else
      pos[2] = (            0 * (n-i - 0.5) +
                mat->pitch[2] * (  i + 0.5)) / n;

      // In pitch units
      double sigmaX = mat->diffSigma / mat->pitch[0]
                       * sqrt(pos[2] / mat->pitch[2]);
      // * 1.01379;// with ExB FIXME

      double sigmaY = mat->diffSigma / mat->pitch[1]
                       * sqrt(pos[2] / mat->pitch[2]);

      for(int x0 = pixel->meas.x[0] - 1; x0 <= pixel->meas.x[0] + 1; x0++)
      for(int y0 = pixel->meas.x[1] - 1; y0 <= pixel->meas.x[1] + 1; y0++)
      {
        double integral = getIntegral(pos[0], x0,x0+1, sigmaX) *
          (mat->isPixel ? getIntegral(pos[1], y0,y0+1, sigmaY)
                         : (y0 == pixel->meas.x[1] ? 1 : 0 ) );

        double Delta = integral * dDelta;

        if(Delta > 0)
          matchOrAdd(hit.allPixels, x0,y0, Delta);
      }
    }
  }

  // add coupling
  if(!mat->isPixel)
    addCoupling(hit, mat);

  // add noise
  for(vector<Pixel>::iterator pixel = hit.allPixels.begin();
                              pixel!= hit.allPixels.end(); pixel++)
  {
    pixel->meas.y *= mat->gain[chipId.code];
    pixel->meas.y += theRandom->getGaussRandom() * mat->noise;

    simDelta += pixel->meas.y;
  }

  // overflow
  for(vector<Pixel>::iterator pixel = hit.allPixels.begin();
                              pixel!= hit.allPixels.end(); pixel++)
    if(pixel->meas.y > mat->overflow)
       pixel->meas.y = mat->overflow + 1e-4; // add 0.1 keV

  // Threshold, remove if below
  for(vector<Pixel>::iterator pixel = hit.allPixels.begin();
                              pixel!= hit.allPixels.end(); /* empty */)
  {
    if(pixel->meas.y < mat->threshold)
      hit.allPixels.erase(pixel);
    else
    {
      recDelta += pixel->meas.y;

      pixel++;
    }
  }

  hit.dpos[0]     = (endpoint[1][0] - endpoint[0][0]);
  hit.dpos[1]     = (endpoint[1][1] - endpoint[0][1]);

  hit.pos_orig[0] = (endpoint[0][0] + endpoint[1][0])/2;
  hit.pos_orig[1] = (endpoint[0][1] + endpoint[1][1])/2;

  hit.charge_orig = simDelta;

  return hit;
}

