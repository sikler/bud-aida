#include "../interface/ClusterGenerator.h"
#include "../interface/Levels.h"
#include "../interface/TouchedChannels.h"
#include "../interface/ElossModel.h"
#include "../interface/NewtonMinimizer.h"

#include "../interface/LinearProgramming.h"
#include "../interface/FitStripCluster.h"

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
Hit readHit(ifstream & file)
{
  Hit hit;

  hit.threshold = Threshold;
  hit.overflow  = Overflow;

  /////////////
  {
    file >> hit.unit.nrows;
    file >> hit.unit.ncolumns;
  
    file >> hit.unit.pitch[0];
    file >> hit.unit.pitch[1];
    file >> hit.unit.pitch[2];
  }

  /////////////
  int nPixels;
  file >> nPixels; 

  file >> hit.dpos[0]; // projected path, in pitch units
  file >> hit.dpos[1];

  double dpos[3];
  for(int k = 0; k < 2; k++)
    dpos[k] = hit.dpos[k] * hit.unit.pitch[k];
  dpos[2] = hit.unit.pitch[2];

  /////////////
  hit.lambda = sqrt( sqr(dpos[0]) + sqr(dpos[1]) );
  hit.length = sqrt( sqr(dpos[0]) + sqr(dpos[1]) + sqr(dpos[2]) );

  for(int k = 0; k < 2; k++)
    hit.ulambda[k] = dpos[k] / hit.lambda;

  for(int i = 0; i < nPixels; i++)
  {
    Pixel pixel;
    Pixel::Meas meas;

    file >> meas.x[0]; // position [pitch units]
    file >> meas.x[1];

    file >> meas.y;    // deposited energy, read in e-

    pixel.meas = meas;

    hit.filledPixels.push_back(pixel);
  }

  // # orig
  char s[256];
  file >> s;
  file >> s;

  // Simulated values
  for(int i = 0; i < 2; i++)
    file >> psim[i];

  // Simulated and default reconstructed deposits
  file >> csim;
  file >> cdef;

  file >> betaGamma;

  return hit;
} 

/*****************************************************************************/
void estimatePositionLF(Hit & hit, double & epsilon)
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
void estimatePositionWeighted(Hit & hit, double & epsilon)
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
void estimateDeposit(Hit & hit, double epsilon, double * charge)
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
//          charge[1] += 2 * o - o * tanh(2 - epsilon * pixel->calc.l / o);
      }
      else
        charge[1] += t * tanh(epsilon * pixel->calc.l / t);
    }
}

/*****************************************************************************/
void mapChi2(ElossModel & theModel, const vector<double> & pars_)
{
  vector<double> pars(3);
  pars[2] = pars_[2];

  ofstream file("../out/map.dat"); 

  double Dp = 0.5;
  double dp = 0.5e-2;

  for(pars[0] = pars_[0] - Dp; pars[0] < pars_[0] + Dp; pars[0] += dp)
  {
  for(pars[1] = pars_[1] - Dp; pars[1] < pars_[1] + Dp; pars[1] += dp)
    {
      file << " " << pars[0]
           << " " << pars[1]
           << " " << theModel.getValue(pars)
           << endl;
  }
  file << endl;
  }
}

/*****************************************************************************/
// Not used FIXME
void checkDerivatives(ElossModel & theModel, const vector<double> & pars_)
{
  vector<double> pars(pars_);

  double dp = 1e-5;

    HepVector beta(3);
    HepMatrix alpha(3,3);
    theModel.getDerivatives(pars, beta, alpha);

  // FIXME
  for(int k = 0; k < 3; k++)
  {
    double val[3];

    pars[k] -= dp; val[0] = theModel.getValue(pars);
    pars[k] += dp; val[1] = theModel.getValue(pars);
    pars[k] += dp; val[2] = theModel.getValue(pars);
    pars[k] -= dp;

    cerr << " der[" << k << "] " << (val[2] - val[0])/(2 * dp)
                << " -beta=" << -beta[k] << endl;

    cerr << " alpha[ " << k << "] " << (val[2] + val[0] - 2*val[1])/dp/dp
                << "   " << alpha[k][k] << endl;
  }

   for(int j = 0; j < 3; j++)
   for(int k = 0; k < 3; k++)
   if(j != k)
   {
     double val[2][2];

     val[0][0] = theModel.getValue(pars); pars[j] += dp;
     val[1][0] = theModel.getValue(pars); pars[k] += dp;
     val[1][1] = theModel.getValue(pars); pars[j] -= dp;
     val[0][1] = theModel.getValue(pars); pars[k] -= dp;

     cerr << " val " << val[0][0] << " " << val[0][1]
              << " " << val[1][0] << " " << val[1][1] << endl;

     cerr << " alpha[" << j << " " << k << "] = "
                   << (val[1][1] + val[0][0] - val[0][1] - val[1][0])/dp/dp
          << "   " << alpha[j][k] << endl;
   }
}

/*****************************************************************************/
int main(int arg, char *arc[])
{
  ClusterGenerator theGenerator;

  int N = 1e+4;

  if(strcmp(arc[1],"-process") == 0)
  {
    N = atoi(arc[2]);
  }

#ifndef Debug
  ifstream fileIn;

  if(strcmp(arc[3],"-pixel") == 0)
    fileIn.open("../out/pixels.dat");
  else
    fileIn.open("../out/strips.dat");
#endif

  int event = 0;

  ofstream fileFit;

  if(strcmp(arc[3],"-pixel") == 0)
    fileFit.open("../out/result_pixels.dat");
  else
    fileFit.open("../out/result_strips.dat");

  for(int ii = 0; ii < N; ii++)
  {
    ++event;

    if(event % 1000 == 0)
       cerr << " Event " << event << endl;

  // Read hit
#ifdef Debug
  if(strcmp(arc[3],"-pixel") == 0)
  {
    theGenerator.run(1, 0);
    ifstream fileIn("../out/pixels.dat");
  }
  else
  {
    theGenerator.run(1, 1);
    ifstream fileIn("../out/strips.dat");
  }
#endif

  Hit hit = readHit(fileIn);

#ifdef Debug
  fileIn.close();
#endif

  // Undo coupling, strips only
  if(strcmp(arc[3],"-pixel") == 0) hit.coupling = 0.;
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

    if(0)
    {
      // Linear programming
      LinearProgramming lp;
      int status;
  
      // FIXME
      vector<double> meas;
  
      for(vector<Pixel>::const_iterator m = hit.filledPixels.begin(); 
                                        m!= hit.filledPixels.end(); m++)
        meas.push_back(m->meas.y);
  
      TVectorD result = lp.solve(hit.coupling, meas, 1., status);
  
      solvedPixels = hit.filledPixels;
      for(int i = 0; i < n; i++)
        if(status == 0)
          solvedPixels[i].meas.y = result[i+1];
        else
          solvedPixels[i].meas.y = ups[i];
    }
    else
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

    // FIXME add first and last
/*
    if(result[0] > 1.) // 1 keV
    {
    Pixel pixel;
    Pixel::Meas meas;

    meas.x[0] = hit.filledPixels[0].meas.x[0] - 1; // position [pitch units]
    meas.x[1] = 50;
    meas.y    = result[0];

    pixel.meas = meas;

    solvedPixels.push_back(pixel);
    }

    if(result[n+1] > 1.) // 0.1 keV
    {
    Pixel pixel;
    Pixel::Meas meas;

    meas.x[0] = hit.filledPixels[n-1].meas.x[0] + 1; // position [pitch units]
    meas.x[1] = 50;
    meas.y    = result[n+1];

    pixel.meas = meas;

    solvedPixels.push_back(pixel);
    }
*/
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
//    hit.overflow  =  1e+9; // FIXME no limits
//    hit.threshold = -1e+9; // FIXME
  }

  ElossModel theModel(hit);

  HepMatrix errs(3,3);
  int nstep;
  
  vector<bool> isFix(3, false);
  if(strcmp(arc[3],"-strip") == 0)
    isFix[1] = true;

  NewtonMinimizer theMinimizer;
  theMinimizer.minimize(theModel, isFix, pars, errs, nstep);

  for(int k = 0; k < 2; k++)
    prec[2][k] = pars[k];

  // Write out
  fileFit << " " << n;

  // Go though the three methods
  for(int j = 0; j < 3; j++)
  {
    // Cartesian residuals
    for(int k = 0; k < 2; k++)
      fileFit << " " << (prec[j][k] - psim[k]) * hit.unit.pitch[k] * 1e+4;

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

    for(int k = 0; k < 2; k++)
      fileFit << " " << p[k] * 1e+4; 
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

//cerr << " cccc " << csim << " " << cdef << " " << crec[0] << " " << crec[1] << endl;
  }

  for(int k = 0; k < 2; k++)
    fileFit << " " << crec[k] - csim;
  fileFit << " " << cdef - csim;

  int ierr;
  errs.invert(ierr);

  // in um
  for(int k = 0; k < 3; k++)
    fileFit << " " << sqrt(errs[k][k]) * hit.unit.pitch[k] * 1e+4;

  fileFit << " " << sqrt(errs[2][2]);

  fileFit << " " << betaGamma
          << " " << csim    / hit.lambda
          << " " << crec[1] / hit.lambda << endl; 
  }

#ifndef Debug
  fileIn.close();
#endif

  fileFit.close();

  return 0;
}
