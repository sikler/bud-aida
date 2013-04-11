#include "../interface/ClusterReco.h"
#include "../interface/Levels.h"
#include "../interface/ElossModel.h"
#include "../interface/NewtonMinimizer.h"
#include "../interface/FitStripCluster.h"

#include "../../DataFormats/interface/TTrack.h"
#include "../../DataFormats/interface/TLayer.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cstring>

using namespace std;

double csim, cdef, psim[2], betaGamma;

/*****************************************************************************/
ClusterReco::ClusterReco(const vector<TLayer> & materials_) :
  materials(materials_) 
{
  filePixel.open("../out/result_pixels.dat");
  fileStrip.open("../out/result_strips.dat");

  fileFit.push_back(&filePixel);
  fileFit.push_back(&fileStrip);
}

/*****************************************************************************/
ClusterReco::~ClusterReco()
{
  fileFit.clear();

  filePixel.close();
  fileStrip.close();
}

/*****************************************************************************/
void ClusterReco::fixHit(Hit & hit)
{
  // FIXME !!!!!! helycsere
  hit.filledPixels = hit.allPixels;
  hit.allPixels.clear();
}

/*****************************************************************************/
void ClusterReco::estimatePositionLF(Hit & hit, double & epsilon)
{
  int xfirs[2] = {1e+3,1e+3};
  int xlast[2] = {-1,-1};

  double qfirs[2] = {0,0};
  double qlast[2] = {0,0};

  // Find extremes
  for(vector<Pixel>::const_iterator pixel = hit.filledPixels.begin();
                                    pixel!= hit.filledPixels.end(); pixel++)
    for(int k = 0; k < 2; k++)
    {
      if(pixel->meas.x[k] < xfirs[k]) xfirs[k] = pixel->meas.x[k];
      if(pixel->meas.x[k] > xlast[k]) xlast[k] = pixel->meas.x[k];
    }

  // Sum charge
  for(vector<Pixel>::const_iterator pixel = hit.filledPixels.begin();
                                    pixel!= hit.filledPixels.end(); pixel++)
    for(int k = 0; k < 2; k++)
    { 
      if(pixel->meas.x[k] == xfirs[k]) qfirs[k] += pixel->meas.y;
      if(pixel->meas.x[k] == xlast[k]) qlast[k] += pixel->meas.y;
    }

  for(int k = 0; k < 2; k++)
  {  
    // get W, from hit.dpos[] // FIXME
    double W_eff = fabs(hit.dpos[k]) - (xlast[k] - xfirs[k] - 1);

    double delta = 1./2 * (qlast[k] - qfirs[k]) / (qlast[k] + qfirs[k]) * W_eff;

    hit.pos[k] = (xfirs[k] + xlast[k] + 1.)/2 + delta;
  }

  double q = 0.;
  for(vector<Pixel>::const_iterator pixel = hit.filledPixels.begin();
                                    pixel!= hit.filledPixels.end(); pixel++)
    q += pixel->meas.y;

  epsilon = q / hit.length;
}

/*****************************************************************************/
void ClusterReco::estimatePositionWeighted(Hit & hit, double & epsilon)
{
  double s = 0., sx = 0., sy = 0.;

  for(vector<Pixel>::const_iterator pixel = hit.filledPixels.begin();
                                    pixel!= hit.filledPixels.end(); pixel++)
  {
    Pixel::Meas meas = pixel->meas;

    double w = meas.y;

    s  += w;
    sx += w * (meas.x[0] + 0.5); 
    sy += w * (meas.x[1] + 0.5); 
  }

  hit.pos[0] = sx/s;
  hit.pos[1] = sy/s;

  double q = 0.;
  for(vector<Pixel>::const_iterator pixel = hit.filledPixels.begin();
                                    pixel!= hit.filledPixels.end(); pixel++)
    q += pixel->meas.y;

  epsilon = q / hit.length;
}

/*****************************************************************************/
void ClusterReco::estimateDeposit(Hit & hit, TLayer * unit,
                                  double epsilon, double * charge)
{
  for(int k = 0; k < 2; k++)
    charge[k] = 0.;

  // Plain sum
  for(vector<Pixel>::const_iterator pixel = hit.filledPixels.begin();
                                    pixel!= hit.filledPixels.end(); pixel++)
    charge[0] += pixel->meas.y;

  double t = unit->threshold;
  double o = unit->overflow;

  // With model
  for(vector<Pixel>::const_iterator pixel = hit.allPixels.begin();
                                    pixel!= hit.allPixels.end(); pixel++)
    if(pixel->calc.isTouched)
    {
      if(pixel->meas.y > 0)
      {
        if(pixel->meas.y < o)
          charge[1] += pixel->meas.y;
        else
        {
          charge[1] += max(o + cForOverflow, epsilon * pixel->calc.l);
          hit.hasOverflow = true;
        }
      }
      else
        charge[1] += t * tanh(epsilon * pixel->calc.l / t);
    }
    else
    {
      if(pixel->meas.y > o)
        hit.hasOverflow = true;
    }
}

/*****************************************************************************/
pair<double,double> ClusterReco::processHit(Hit & hit)
{
  TLayer * unit = &materials[hit.ilayer];

  // FIXME
  for(int k = 0; k < 2; k++)
    psim[k] = hit.pos_orig[k];

  csim = hit.charge_orig;

  // Undo coupling, strips only
  vector<Pixel> unfoldedPixels, solvedPixels;

  int n = hit.filledPixels.size();

  vector<Pixel> origStrips = hit.filledPixels;

  if(!unit->isPixel)
  {
    double a = unit->coupling;

    // Simple unfolding
    vector<double> ups(n,0);
 
    vector<double> c(n);
    c[0] = 1 / (1 - 2*a);

    for(int i = 1; i < n; i++)
      c[i] = c[i-1] / (2 - 1/a);

    for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      ups[i] += c[abs(i-j)] * hit.filledPixels[j].meas.y;

    unfoldedPixels = hit.filledPixels;
    for(int i = 0; i < n; i++)
      unfoldedPixels[i].meas.y = ups[i];

    {
      enum { isNormal, isBelow, isOver };

      vector<pair<double,int> > b;

      b.push_back(pair<double,int>(unit->threshold, isBelow));

      for(vector<Pixel>::const_iterator m = hit.filledPixels.begin(); 
                                        m!= hit.filledPixels.end(); m++)
      {
        if(m->meas.y > unit->overflow)
          b.push_back(pair<double,int>(m->meas.y, isOver  ));
        else
          b.push_back(pair<double,int>(m->meas.y, isNormal));
      }

      b.push_back(pair<double,int>(unit->threshold, isBelow));

      FitStripCluster theFitter(b, unit->coupling, unit->noise);
      double chi2; pair<double,double> result; vector<double> x;
      /*int width =*/ theFitter.run(chi2, result, x);

      solvedPixels = hit.filledPixels;
      for(int i = 0; i < n; i++)
        solvedPixels[i].meas.y = x[i+1];
    }
  }

  // Estimate position
  double epsilon;
  double prec[3][2];

  if(!unit->isPixel)
    hit.filledPixels = unfoldedPixels;

  estimatePositionWeighted(hit, epsilon);
  for(int k = 0; k < 2; k++)
    prec[0][k] = hit.pos[k];

  estimatePositionLF      (hit, epsilon);
  for(int k = 0; k < 2; k++)
    prec[1][k] = hit.pos[k];

  // Initial values
  vector<double> pars(3);
  pars[0] = hit.pos[0];
  pars[1] = hit.pos[1];
  pars[2] = log(epsilon);        // log(3 MeV/cm)

  // Get model
  if(!unit->isPixel)
    hit.filledPixels = solvedPixels;

  ElossModel theModel(hit, unit);

  cForOverflow = theModel.getSigma(unit->overflow)/Nu;

  TMatrixD errs(3,3);
  int nstep;
  
  vector<bool> isFix(3, false);
  if(!unit->isPixel)
    isFix[1] = true;

  NewtonMinimizer theMinimizer;
  theMinimizer.minimize(theModel, isFix, pars, errs, nstep);

  for(int k = 0; k < 2; k++)
    prec[2][k] = pars[k];

  // Write out
  int type;
  if(unit->isPixel) type = 0;
               else type = 1;

  *fileFit[type] << " " << n;

  // Go though the three methods
  for(int j = 0; j < 3; j++)
  {
    // Cartesian residuals
    for(int k = 0; k < 2; k++)
      *fileFit[type] << " " << (prec[j][k] - psim[k]) * unit->pitch[k] * 1e+4;

    // Projected residuals
    double p[2];
    {
      double b[2], v[2];
      for(int k = 0; k < 2; k++)
        b[k] = (prec[j][k] - psim[k]) * unit->pitch[k];

      for(int k = 0; k < 2; k++)
        v[k] = hit.ulambda[k];

      p[0] = (b[0]*v[0] + b[1]*v[1]);
      p[1] = (b[0]*v[1] - b[1]*v[0]);
    }

    for(int k = 0; k < 2; k++)
      *fileFit[type] << " " << p[k] * 1e+4; 
  }

  // Deposit
  hit = theModel.getHit(pars);
  double crec[2] = {0,0};

  
  if(unit->isPixel)
  {
    estimateDeposit(hit, unit, exp(pars[2]), crec);
  }
  else
  {
    for(vector<Pixel>::const_iterator pixel = origStrips.begin();
                                      pixel!= origStrips.end(); pixel++)
      crec[0] += pixel->meas.y;
 
    for(vector<Pixel>::const_iterator pixel = solvedPixels.begin();
                                      pixel!= solvedPixels.end(); pixel++)
      crec[1] += pixel->meas.y;
  }

/*
  // FIXME
  if(hit.hasOverflow || crec[1] < crec[0]
    crec[1] = crec[0];
*/

  for(int k = 0; k < 2; k++)
    *fileFit[type] << " " << crec[k] - csim;
  *fileFit[type] << " " << cdef - csim;

  // add diagonal elements
  for(int k = 0; k < 3; k++)
    errs[k][k] += 1e-6;

  errs.Invert();

  // in um
  for(int k = 0; k < 3; k++)
    *fileFit[type] << " " << sqrt(errs[k][k]) * unit->pitch[k] * 1e+4;

  *fileFit[type] << " " << sqrt(errs[2][2]);

  *fileFit[type] << " " << betaGamma
          << " " << csim    / hit.lambda
          << " " << crec[1] / hit.lambda
          << " " << (unit->isPixel ? "pixel" : "strip")
          << endl; 

  return pair<double,double>(crec[0], crec[1]);
}

/*****************************************************************************/
void ClusterReco::run(TTrack & recTrack)
{
  if(recTrack.hits.size() > 0)
  for(vector<Hit>::iterator hit = recTrack.hits.begin();
                            hit!= recTrack.hits.end(); hit++)
  {
    fixHit(*hit);

    pair<double,double> Delta = processHit(*hit);

    hit->charge = Delta.second;
  }
}

