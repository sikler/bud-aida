#include "../interface/RecoClusters.h"

#include "../interface/ClusterGenerator.h"
#include "../interface/Levels.h"
#include "../interface/TouchedChannels.h"
#include "../interface/ElossModel.h"
#include "../interface/NewtonMinimizer.h"

#include "../interface/LinearProgramming.h"
#include "../interface/FitStripCluster.h"

#include "../../DataFormats/interface/TTrack.h"

#include "CLHEP/Matrix/Matrix.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cstring>

#define sqr(x) ((x) * (x))

#undef Debug
//#define Debug

using namespace std;
using namespace CLHEP;

double csim, cdef, psim[2], betaGamma;

/*****************************************************************************/
RecoClusters::RecoClusters()
{
//  if(strcmp(arc[3],"-pixel") == 0)
// FIXME
//    fileFit.open("../out/result_pixels.dat");
/*
  else
    fileFit.open("../out/result_strips.dat");
*/
}

/*****************************************************************************/
RecoClusters::~RecoClusters()
{
//  fileFit.close();
}

/*****************************************************************************/
void RecoClusters::fixHit(Hit & hit)
{
  hit.threshold = Threshold;
  hit.overflow  = Overflow;

  // FIXME !!!!!! helycsere
  hit.filledPixels = hit.allPixels;
  hit.allPixels.clear();
}

/*****************************************************************************/
void RecoClusters::estimatePositionLF(Hit & hit, double & epsilon)
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
    // get W, from hit.dpos[], FIXME
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
void RecoClusters::estimatePositionWeighted(Hit & hit, double & epsilon)
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
void RecoClusters::estimateDeposit(Hit & hit, double epsilon, double * charge)
{
  for(int k = 0; k < 2; k++)
    charge[k] = 0.;

  for(vector<Pixel>::const_iterator pixel = hit.filledPixels.begin();
                                    pixel!= hit.filledPixels.end(); pixel++)
    charge[0] += pixel->meas.y;

  double t = hit.threshold;
  double o = hit.overflow;

  for(vector<Pixel>::const_iterator pixel = hit.allPixels.begin();
                                    pixel!= hit.allPixels.end(); pixel++)
    if(pixel->calc.isTouched)
    {
      if(pixel->meas.y > 0)
      {
        if(pixel->meas.y < o)
          charge[1] += pixel->meas.y;
        else
          charge[1] += max(o + 25, epsilon * pixel->calc.l);
      }
      else
        charge[1] += t * tanh(epsilon * pixel->calc.l / t);
    }
}

/*****************************************************************************/
pair<double,double> RecoClusters::processHit(Hit & hit)
{
  // FIXME
  bool isPixel = true;
  
  // Undo coupling, strips only
  if(isPixel) hit.coupling = 0.;
         else hit.coupling = Coupling;

  vector<Pixel> unfoldedPixels, solvedPixels;

  int n = hit.filledPixels.size();

  vector<Pixel> origStrips = hit.filledPixels;

  if(hit.coupling > 0)
  {
    double a = hit.coupling;

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

      b.push_back(pair<double,int>(hit.threshold, isBelow));

      for(vector<Pixel>::const_iterator m = hit.filledPixels.begin(); 
                                        m!= hit.filledPixels.end(); m++)
      {
        if(m->meas.y > hit.overflow)
          b.push_back(pair<double,int>(m->meas.y, isOver  ));
        else
          b.push_back(pair<double,int>(m->meas.y, isNormal));
      }

      b.push_back(pair<double,int>(hit.threshold, isBelow));

      FitStripCluster theFitter(b, hit.coupling, Noise);
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

  if(hit.coupling > 0)
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
  pars[2] = log(epsilon);        // log(3000 keV/cm)

  // Get model
  if(hit.coupling > 0)
  {
    hit.filledPixels = solvedPixels;
  }

  ElossModel theModel(hit);

  HepMatrix errs(3,3);
  int nstep;
  
  vector<bool> isFix(3, false);
  if(!isPixel) // FIXME
    isFix[1] = true;

  NewtonMinimizer theMinimizer;
  theMinimizer.minimize(theModel, isFix, pars, errs, nstep);
  cerr << " nstep=" << nstep << endl;

  for(int k = 0; k < 2; k++)
    prec[2][k] = pars[k];

  // Write out
//  fileFit << " " << n;

  // Go though the three methods
  for(int j = 0; j < 3; j++)
  {
    // Cartesian residuals
//    for(int k = 0; k < 2; k++)
//      fileFit << " " << (prec[j][k] - psim[k]) * hit.unit.pitch[k] * 1e+4;

    // Projected residuals
    double p[2];
    {
      double b[2], v[2];
      for(int k = 0; k < 2; k++)
        b[k] = (prec[j][k] - psim[k]) * hit.unit.pitch[k];

      for(int k = 0; k < 2; k++)
        v[k] = hit.ulambda[k];

      p[0] = (b[0]*v[0] + b[1]*v[1]);
      p[1] = (b[0]*v[1] - b[1]*v[0]);
    }

//    for(int k = 0; k < 2; k++)
//      fileFit << " " << p[k] * 1e+4; 
  }

  // Deposit
  hit = theModel.getHit(pars);
  double crec[2] = {0,0};

  
  if(hit.coupling == 0.)
  {
//    Hit hit = theModel.getHit(pars);
    estimateDeposit(hit, exp(pars[2]), crec);
  }
  else
  {
    for(vector<Pixel>::const_iterator pixel = origStrips.begin();
                                      pixel!= origStrips.end(); pixel++)
      crec[0] += pixel->meas.y;
 
    for(vector<Pixel>::const_iterator pixel = hit.filledPixels.begin();
                                      pixel!= hit.filledPixels.end(); pixel++)
      crec[1] += pixel->meas.y;
  }

//  for(int k = 0; k < 2; k++)
//    fileFit << " " << crec[k] - csim;
//  fileFit << " " << cdef - csim;

  int ierr;
  errs.invert(ierr);

  // in um
//  for(int k = 0; k < 3; k++)
//    fileFit << " " << sqrt(errs[k][k]) * hit.unit.pitch[k] * 1e+4;
//
//  fileFit << " " << sqrt(errs[2][2]);
//
//  fileFit << " " << betaGamma
//          << " " << csim    / hit.lambda
//          << " " << crec[1] / hit.lambda << endl; 

  cerr << "  pos Weighted   : " << prec[0][0] << " " << prec[0][1] << endl;
  cerr << "  pos Last-First : " << prec[1][0] << " " << prec[1][1] << endl;
  cerr << "  pos Fitter     : " << prec[2][0] << " " << prec[2][1] << endl;

  return pair<double,double>(crec[0], crec[1]);
}

/*****************************************************************************/
void RecoClusters::run(TTrack & recTrack)
{
  cerr << " doing clusters " << recTrack.lhits.size() << endl;

  if(recTrack.lhits.size() > 0)
  for(vector<Hit>::iterator hit = recTrack.lhits.begin();
                            hit!= recTrack.lhits.end(); hit++)
  {
    cerr << "  doing cluster " << int(hit - recTrack.lhits.begin()) + 1;
//    hit->print();
//    cerr << " -----------------------------------------" << endl;
    fixHit(*hit);
    cerr << " (" << hit->filledPixels.size() << ")";
//    hit->print();
//    hit->unit.print();
//    while(getchar()==0); 

    pair<double,double> Delta = processHit(*hit);
    cerr << " Delta = " << Delta.first << " " << Delta.second
         << " l=" << hit->length << endl;

    hit->charge = Delta.second;
  }
}

