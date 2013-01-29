#include "../interface/ElossEstimator.h"

#include "../interface/ModelBichsel.h"
#include "../interface/MostProbable.h"

#include "../../DataFormats/interface/TTrack.h"

#include <cmath>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

#define nTracks 100000

#define Sqr(x) ((x) * (x))

struct Layer_t {
  double depth; // cm
  double gain;  // number
};

vector<Layer_t> layers;
int nLayers;

vector<vector<double> > tracks;

double mass[2] = {0.139, 0.493};
double p = 0.8;
double eps;

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
  double sigma0 = 2.0e-3;
  double b = 0.095;

  return sigma0 + b * y;
}

/*****************************************************************************/
double ElossEstimator::logLikelihood(double epsilon, double y, double l)
{
  double nu = 0.65;

  // noise
  double sn = 3e-3; // MeV

  // calculate effective path
  double a  = 0.07;   //
  double l0 = 300e-4; // cm
  double ls = l * (1 + a *log(l/l0));

  double Delta = epsilon * ls;

  double sD  = sigmaD(y);
  //double sDD = sigmaD(Delta);

  double s   = sqrt(sD *sD  + sn*sn);
  //double s0  = sqrt(sDD*sDD + sn*sn);

  double chi2;
  if(Delta < y - nu * s) chi2 = - 2*nu * (Delta - y)/s - nu*nu;
                    else chi2 =      Sqr((Delta - y))/ (s*s);

//  return -2*log( exp(-chi2/2) / (s/s0) );
  return -2*log( exp(-chi2/2) / s);
}

/****************************************************************************/
double ElossEstimator::getChi2(double gain, const vector<double> & hits, double l)
{
  double chi2 = 0.;

  for(vector<double>::const_iterator hit = hits.begin();
                                     hit!= hits.end(); hit++)
    chi2 += logLikelihood(eps, gain * (*hit), l);

  return chi2;
}

/****************************************************************************/
double ElossEstimator::getChi2(double epsilon, const vector<pair<double,double> > & hits)
{
  double chi2 = 0.;

  for(vector<pair<double,double> >::const_iterator hit = hits.begin();
                                                   hit!= hits.end(); hit++)
    chi2 += logLikelihood(epsilon, hit->first, hit->second);

  return chi2;
}

/****************************************************************************/
void ElossEstimator::calibrateGains(int pass)
{
  char name[256];
  sprintf(name,"../out/gains_%d.dat",pass);

  ofstream file(name);

  for(int j = 0; j < nLayers; j++)
  {
    cerr << "  layer[" << j << "] ";

    vector<double> hits;
/*
    for(vector<vector<double> >::const_iterator track = tracks.begin();
                                                track!= tracks.end(); track++)
*/
    for(int i = 0; i < nTracks; i++) // FIXME
      hits.push_back(tracks[i][j]); 

    double gain = 0.; 
    double chi2min = 1e+30;
    double dg = 1e-3;

    for(double g = 0.5; g < 1.5; g += dg)
    {
      double chi2 = getChi2(g, hits, layers[j].depth);

      if(chi2 < chi2min)
      { chi2min = chi2 ; gain = g; }
    }

    double g,val[3];
    g = gain - dg; val[0] =  getChi2(g, hits, layers[j].depth);
    g = gain     ; val[1] =  getChi2(g, hits, layers[j].depth);
    g = gain + dg; val[2] =  getChi2(g, hits, layers[j].depth);

    double secder = (val[2] + val[0] - 2*val[1]) / dg / dg;
    double sig = sqrt(2/secder);

    cerr << " \t" << 1/layers[j].gain << " " << gain << " " << sig << endl;
    file <<   " " << 1/layers[j].gain << " " << gain << " " << sig << endl;

    // correct
    for(vector<vector<double> >::iterator track = tracks.begin();
                                          track!= tracks.end(); track++)
//      (*track)[j] *= gain;
      (*track)[j] *= 1/layers[j].gain;
  }

  file.close();
}

/****************************************************************************/
pair<double,double> ElossEstimator::getEstimate
  (const vector<pair<double,double> > & track)
{
  vector<pair<double,double> > hits; // FIXME miert kell?
  for(vector<pair<double,double> >::const_iterator hit = track.begin();
                                                   hit!= track.end(); hit++)
    hits.push_back(*hit);

  double epsilon = 0.;
  double chi2min = 1e+30;
  double de = 1e-3;

  for(double e = 1.; e < 100.; e += de) // FIXME ?? Ez jobb ???
  {
    double chi2 = getChi2(e, hits);

    if(chi2 < chi2min)
    { chi2min = chi2 ; epsilon = e; }
  }

  double e,val[3];
  e = epsilon - de; val[0] =  getChi2(e, hits);
  e = epsilon     ; val[1] =  getChi2(e, hits);
  e = epsilon + de; val[2] =  getChi2(e, hits);

  double secder = (val[2] + val[0] - 2*val[1]) / de / de;
  double sig = sqrt(2/secder);
 
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
void ElossEstimator::estimate(TTrack & recTrack, ofstream & fileEstimate)
{
  if(recTrack.lhits.size() > 0)
  {
  vector<pair<double,double> > track;

  for(vector<Hit>::iterator hit = recTrack.lhits.begin();
                            hit!= recTrack.lhits.end(); hit++)
  {
    cerr << "  hit " << hit->charge << " " << hit->length << endl;

    pair<double,double> x(hit->charge*1e-3, hit->length); // MeV, cm

    track.push_back(x);
  }

  // Fitter
  pair<double,double> estimate = getEstimate(track);
  double p = recTrack.pt * cosh(recTrack.eta);
  fileEstimate <<  " estimates " << p
               << " " << estimate.first << " " << estimate.second;

  vector<double> dedx;
  for(vector<pair<double,double> >::const_iterator deposit = track.begin();
                                                   deposit!= track.end();
                                                   deposit++)
      dedx.push_back(deposit->first / deposit->second);

    fileEstimate << " " << getTruncMean(dedx)
         << " " << getPowerMean(dedx, -2.)
         << " " << getPowerMean(dedx, -1.)
         << " " << getPowerMean(dedx,  1.)
         << endl;

  cerr << " --------------------------------------------" << endl;
  }

/*
  char name[256];
  sprintf(name,"../out/estimates_%d.dat",pass);

  ofstream file(name);

  cerr << "  ";
*/

/*
  // This model, max likelihood
  int j = 0;
  for(vector<vector<double> >::const_iterator track = tracks.begin();
                                              track!= tracks.end(); track++)
  {
//    if((j+1) % (tracks.size()/40) == 0) cerr << ".";
//    file << " " << (j++ < nTracks ? 0 : 1);

    pair<double,double> estimate = getEstimate(*track);
//    file <<  " " << estimate.first << " " << estimate.second;

    vector<double> dedx;
    int j = 0;
    for(vector<double>::const_iterator hit = track->begin();
                                       hit!= track->end(); hit++)
      dedx.push_back(*hit / layers[j++].depth);
*
    file << " " << getTruncMean(dedx);
    file << " " << getPowerMean(dedx, -2.);
    file << " " << getPowerMean(dedx, -1.);
    file << " " << getPowerMean(dedx,  1.);

    file << endl;
*/
//  }

//  file.close();
}

/****************************************************************************/
//int main(int arg, char **arc)
//{
/*
  MostProbable mostProbable;
  eps = mostProbable.value(p / mass[0]);
  cerr << " epsilon = " << eps << endl;

  cerr << " calibrating gains [with pions only].." << endl;
  calibrateGains(pass);
  cerr << " [done]" << endl;

  cerr << " estimating energy loss..." << endl;
  estimateEloss(pass); 
  cerr << " [done]" << endl;
*/
//}
