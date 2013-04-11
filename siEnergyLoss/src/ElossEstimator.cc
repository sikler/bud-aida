#include "../interface/Levels.h"
#include "../interface/ElossEstimator.h"

#include "../../DataFormats/interface/TTrack.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

#define sqr(x) ((x) * (x))

/*****************************************************************************/
ElossEstimator::ElossEstimator()
{
}

/*****************************************************************************/
ElossEstimator::~ElossEstimator()
{
}

/*****************************************************************************/
inline double ElossEstimator::sigmaD(double y)
{
  return Sigma0 + b * y;
}

/*****************************************************************************/
double ElossEstimator::logLikelihood(double epsilon, double y,
                                     double l, bool over)
{
  // calculate effective path
  double ls = l * (1 + a *log(l/l0));
  double Delta = epsilon * ls;

  double sD  = sigmaD(y);
  double s   = sqrt(sD *sD  + Noise*Noise);

  double chi2;

  if(!over)
  {
    if(Delta < y - Nu * s) chi2 = - 2*Nu * (Delta - y)/s - Nu*Nu;
                      else chi2 =      sqr((Delta - y))/ (s*s);
  }
  else
  {
    if(Delta < y + Nu * s) chi2 = - (Delta - y)/s + 1;
                      else chi2 = 0.;
  }

  return -2*log( exp(-chi2/2) / s);
}

/****************************************************************************/
double ElossEstimator::getChi2Gain
  (double gain, const vector<TSlimMeasurement> & hits)
{
  double chi2 = 0.;

  for(vector<TSlimMeasurement>::const_iterator hit = hits.begin();
                                               hit!= hits.end(); hit++)
  if(!hit->hasOverflow && hit->epsilon < 10 && hit->path < 2) // restrict
    chi2 += logLikelihood(hit->epsilon, gain * hit->energy,
                          hit->path, hit->hasOverflow);

  return chi2;
}

/****************************************************************************/
double ElossEstimator::getChi2Epsilon
  (double epsilon, const vector<TSlimMeasurement> & hits)
{
  double chi2 = 0.;

  for(vector<TSlimMeasurement>::const_iterator hit = hits.begin();
                                               hit!= hits.end(); hit++)
    chi2 += logLikelihood(epsilon, hit->energy, hit->path, hit->hasOverflow);

  return chi2;
}

/****************************************************************************/
void ElossEstimator::getValuesByVaryingGain
  (const vector<TSlimMeasurement> & hits,
   double gain, vector<double> & val)
{
  const double dg = 1e-4;

  // Values
  val[0] =  getChi2Gain(gain - dg, hits);
  val[1] =  getChi2Gain(gain     , hits);
  val[2] =  getChi2Gain(gain + dg, hits);

  // Derivatives
  double sec = (val[2] + val[0] - 2*val[1]) / (dg * dg);
  double fir = (val[2] - val[0]) / (2 * dg);

  val[0] = val[1]; // central value
  val[1] = fir;    // first derivative
  val[2] = sec;    // second derivative
}

/****************************************************************************/
void ElossEstimator::getValuesByVaryingEpsilon
  (const vector<TSlimMeasurement> & hits,
   double epsilon, vector<double> & val)
{
  const double de = 1e-3;

  // Values
  val[0] =  getChi2Epsilon(epsilon - de, hits);
  val[1] =  getChi2Epsilon(epsilon     , hits);
  val[2] =  getChi2Epsilon(epsilon + de, hits);

  // Derivatives
  double sec = (val[2] + val[0] - 2*val[1]) / (de * de);
  double fir = (val[2] - val[0]) / (2 * de);

  val[0] = val[1]; // central value
  val[1] = fir;    // first derivative
  val[2] = sec;    // second derivative
}

/*****************************************************************************/
void ElossEstimator::functionGain(double g, vector<double> & val)
{
  getValuesByVaryingGain(slim, g, val);
} 

/*****************************************************************************/
double ElossEstimator::functionGain(double gain)
{
  vector<double> val(3);

  functionGain(gain, val);

  return val[0];
}

/****************************************************************************/
void ElossEstimator::shft2(double &a, double &b, const double c)
{ a=b; b=c; }

/*****************************************************************************/
void ElossEstimator::shft3(double &a, double &b, double &c, const double d)
{ a=b; b=c; c=d; }

/*****************************************************************************/
double ElossEstimator::goldenSearch(double ax, double bx, double cx)
{
  double xmin; //,fmin;
  
  const double R = 0.61803399; // (sqrt(5.) - 1) / 2
  const double C = 1 - R;
  
  double x1,x2;
  double x0=ax;
  double x3=cx;

  if(fabs(cx-bx) > fabs(bx-ax))
  { x1 = bx; x2 = bx+C*(cx-bx); }
  else
  { x2 = bx; x1 = bx-C*(bx-ax); }

  double f1 = functionGain(x1);
  double f2 = functionGain(x2);

  while(fabs(x3-x0) > 1e-2)
  {
    if(f2 < f1)
    {
      shft3(x0,x1,x2, R*x2+C*x3);
      shft2(   f1,f2, functionGain(x2));
    }
    else
    {
      shft3(x3,x2,x1, R*x1+C*x0);
      shft2(   f2,f1, functionGain(x1));
    }
  }

  if(f1 < f2)
  { xmin=x1; } // fmin=f1; }
  else
  { xmin=x2; } // fmin=f2; }

  return xmin;
}

/*****************************************************************************/
double ElossEstimator::newtonMethodGain(pair<double,double> & value)
{
  int nStep = 0;

  double par = value.first; // input MeV/cm
  double dpar;

  vector<double> val(3);

  do
  { 
    functionGain(par, val);

    if(val[2] != 0.) dpar = - val[1]/fabs(val[2]);
                else dpar = 1.;                     // step up, for epsilon

    if(par + dpar > 0) par += dpar; // ok
                  else par /= 2;    // half

    nStep++;
  }
  while(fabs(dpar) > 1e-3 && nStep < 50);
  
  value.first  = par;
  value.second = 2/val[2]; // sigma2 was 2
  
  return val[0];
}

/****************************************************************************/
void ElossEstimator::loadGains(ifstream & fileGain)
{
  int d;
  float f;

  while(!fileGain.eof()) 
  {
    ChipId key;

    fileGain >> key.first
             >> key.second.first
             >> key.second.second
             >> d;

    if(!fileGain.eof())
      fileGain >> gains[key]
               >> f;
  }
}

/****************************************************************************/
void ElossEstimator::correctDeposits(TTrack * track)
{
  for(vector<Hit>::iterator hit = track->hits.begin();
                            hit!= track->hits.end(); hit++)
    for(vector<Pixel>::iterator pixel = hit->allPixels.begin();
                                pixel!= hit->allPixels.end(); pixel++)
      pixel->meas.y *= gains[hit->chipId.code];
}

/****************************************************************************/
void ElossEstimator::calibrateGains
  (map<ChipId, vector<TSlimMeasurement> > & hits, ofstream & fileGain)
{
  cerr << " calibrate gain with golden section search + Newton method" << endl;
  cerr << " number of keys : " << hits.size() << endl;

  for(map<ChipId, vector<TSlimMeasurement> >::iterator
      key = hits.begin(); key!= hits.end(); key++)
  if(key->second.size() > 1)
  {
    slim = key->second;

    pair<double,double> value;

    // Golden section search [0.1,10.] with center 1.0
    value.first =  goldenSearch(0.1, 1.0, 10.);

    // Newton method around value.first, get sigma2 = value.second
    newtonMethodGain(value);

    // If already pre-calibrated, multiply by previous
    if(gains.count(key->first) > 0)
    {
      float g = gains[key->first];

      value.first  *= g;
      value.second *= g*g;
    }

    fileGain << " " << key->first.first
             << " " << key->first.second.first
             << " " << key->first.second.second
             << " " << key->second.size()
             << " " << value.first
             << " " << sqrt(value.second) << endl;

    key->second.clear();
  }
}

/****************************************************************************/
// from Fitters::newtonMethodEpsilon
pair<double,double> ElossEstimator::getEstimate
  (const vector<TSlimMeasurement> & hits)
{
  bool allSaturated = false; // FIXME

  double epsilon = 3.0;
  double sig;

  if(!allSaturated)
  {
    double deps;
    vector<double> val(3);
    int nStep = 0;

    do
    {
      getValuesByVaryingEpsilon(hits, epsilon, val);

      deps = - val[1]/fabs(val[2]);

      // renorm
      if(val[2] == 0. || fabs(deps) > 1) deps = 1.;
  
      if(epsilon + deps > 0) epsilon += deps; // ok
                        else epsilon /= 2;    // half

      nStep++;
    }
    while(fabs(deps) > 1e-3 && nStep < 50);

    getValuesByVaryingEpsilon(hits, epsilon, val);

    sig = sqrt(2/val[2]);
  }
  else
  {
    int nStep = 0;

    vector<double> val(3);

    do
    {
      getValuesByVaryingEpsilon(hits, epsilon, val);

      if(val[1] != 0)
        epsilon += - val[0]/val[1];

      nStep++;
    }
    while(val[0] > 1e-3 && val[1] != 0 && nStep < 10);

    sig = sqr(epsilon*0.1);

    epsilon *= 1.1;
  }
 
  return pair<double,double>(epsilon, sig);
}

/****************************************************************************/
double ElossEstimator::getTruncMean(vector<double> & dedx)
{
  sort(dedx.begin(), dedx.end());

  vector<double> s(2,0.);

  int j = 0;
  for(vector<double>::const_iterator x = dedx.begin();
                                     x!= dedx.end(); x++)
  if(j++ < int(dedx.size())/2)
  {
    s[1] += *x;
    s[0] += 1.;
  }
  
  s[1] /= s[0];

  return s[1];
}

/****************************************************************************/
double ElossEstimator::getPowerMean(const vector<double> & dedx, double power)
{
  vector<double> s(2,0.);

  for(vector<double>::const_iterator x = dedx.begin();
                                     x!= dedx.end(); x++)
  {
    s[1] += pow(*x, power);
    s[0] += 1.;
  }

  s[1] /= s[0];

  return pow(s[1], 1/power);
}

/****************************************************************************/
pair<double,double> ElossEstimator::estimate
  (TTrack & recTrack, ofstream & fileEstimate)
{
  pair<double,double> estimate(0.,0.);

  if(recTrack.hits.size() > 0)
  {
    vector<TSlimMeasurement> track;
  
    bool hasOverflow = false;
  
    for(vector<Hit>::iterator hit = recTrack.hits.begin();
                              hit!= recTrack.hits.end(); hit++)
    {
      if(hit->hasOverflow) hasOverflow = true;

      TSlimMeasurement slim;
      slim.energy      = hit->charge;
      slim.path        = hit->length;
      slim.hasOverflow = hit->hasOverflow;
  
      track.push_back(slim);
    }
  
    // Fitter
    estimate = getEstimate(track);
    double p = recTrack.pt * cosh(recTrack.eta);
    fileEstimate <<  " estimates " << p
                 << " " << estimate.first << " " << estimate.second;
  
    vector<double> dedx;
    for(vector<TSlimMeasurement>::const_iterator deposit = track.begin();
                                                 deposit!= track.end();
                                                 deposit++)
      dedx.push_back(deposit->energy / deposit->path);

    fileEstimate << " " << getTruncMean(dedx)
                 << " " << getPowerMean(dedx, -2.)
                 << " " << getPowerMean(dedx, -1.)
                 << " " << getPowerMean(dedx,  1.)
                 << " " << hasOverflow
                 << endl;
  }
   
  return estimate;
}

