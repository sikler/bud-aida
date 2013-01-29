#include "../interface/ModelBichsel.h"
#include "../interface/MostProbable.h"

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
double getFlatRandom()
{
  return drand48();
}

/****************************************************************************/
double getGaussRandom()
{
  return sqrt(-2*log(getFlatRandom())) * cos(M_PI*getFlatRandom());
}

/****************************************************************************/
void initializeLayers(int pass)
{
  int depth1[16] =
        {480,483,486,494, 501,514,525,547, 565,600,630,691, 749,884,1043,1714};
  int depth2[16] =
        {400,402,403,408, 412,419,425,436, 445,461,475,500, 520,560,594,666};
  int depth3[16] =
        {343,344,345,348, 350,354,358,364, 370,379,386,399, 410,428,443,471};
  int depth4[16] =
        {300,300,301,303, 304,307,310,314, 317,323,328,335, 341,352,360,375};

  for(int k = 0; k < 16; k ++) 
  {
    Layer_t layer;

    if(pass == 0) layer.depth = depth1[k] * 1e-4;
    if(pass == 1) layer.depth = depth2[k] * 1e-4;
    if(pass == 2) layer.depth = depth3[k] * 1e-4;
    if(pass == 3) layer.depth = depth4[k] * 1e-4;

    layer.gain  = 1 + 0.2 * (2*getFlatRandom() - 1);

    layers.push_back(layer);
  }

  // 300 300 300  300
  // 300 400 500  600
  // 300 500 700  900
  // 300 600 900 1200

/*
  double dlow=0, dmax=0, dd=0;

  if(pass == 0) { dlow = 300e-4; dmax =  301e-4; dd = (dmax - dlow) / 4; }
  if(pass == 1) { dlow = 300e-4; dmax =  600e-4; dd = (dmax - dlow) / 4; }
  if(pass == 2) { dlow = 300e-4; dmax =  900e-4; dd = (dmax - dlow) / 4; }
  if(pass == 3) { dlow = 300e-4; dmax = 1200e-4; dd = (dmax - dlow) / 4; }

  for(double d = dlow; d < dmax + 1e-6; d += dd)
  for(int k = 0; k < 4; k ++) 
  {
     Layer_t layer;

     layer.depth = d;
     layer.gain  = 1 + 0.2 * (2*getFlatRandom() - 1);

     layers.push_back(layer);
   } 
*/

  nLayers = layers.size();
}

/****************************************************************************/
void generateTracks(int n, int pid)
{
  double betaGamma = p / mass[pid];

  ModelBichsel theModel("Si");
  theModel.prepare(betaGamma);

  cerr << "  ";

  for(int i = 0; i < n; i++)
  {
    if((i+1) % int(n/40.) == 0) cerr << ".";

    vector<double> hits;
   
    for(int j = 0; j < nLayers; j++)
    {
      // get deposit, add noise
      double y = theModel.generate(layers[j].depth * 1e+4)
               + 3. * getGaussRandom(); // um, keV

      hits.push_back(y * 1e-3 * layers[j].gain);
    }

    tracks.push_back(hits);
  }
}

/*****************************************************************************/
inline double sigmaD(double y)
{
  double sigma0 = 2.0e-3;
  double b = 0.095;

  return sigma0 + b * y;
}

/*****************************************************************************/
double logLikelihood(double epsilon, double y, double l)
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
double getChi2(double gain, const vector<double> & hits, double l)
{
  double chi2 = 0.;

  for(vector<double>::const_iterator hit = hits.begin();
                                     hit!= hits.end(); hit++)
    chi2 += logLikelihood(eps, gain * (*hit), l);

  return chi2;
}

/****************************************************************************/
double getChi2(double epsilon, const vector<double> & hits)
{
  double chi2 = 0.;
  int j = 0;

  for(vector<double>::const_iterator hit = hits.begin();
                                     hit!= hits.end(); hit++)
    chi2 += logLikelihood(epsilon, (*hit), layers[j++].depth);

  return chi2;
}

/****************************************************************************/
void calibrateGains(int pass)
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
pair<double,double> getEstimate(const vector<double> & track)
{
  vector<double> hits;
  for(vector<double>::const_iterator hit = track.begin();
                                     hit!= track.end(); hit++)
    hits.push_back(*hit);

  double epsilon = 0.;
  double chi2min = 1e+30;
  double de = 1e-3;

  for(double e = 1.; e < 5.; e += de)
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
double getTruncMean(vector<double> & dedx)
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
double getPowerMean(const vector<double> & dedx, double power)
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
void estimateEloss(int pass)
{
  char name[256];
  sprintf(name,"../out/estimates_%d.dat",pass);

  ofstream file(name);

  cerr << "  ";

  // This model, max likelihood
  int j = 0;
  for(vector<vector<double> >::const_iterator track = tracks.begin();
                                              track!= tracks.end(); track++)
  {
    if((j+1) % (tracks.size()/40) == 0) cerr << ".";
    file << " " << (j++ < nTracks ? 0 : 1);

    pair<double,double> estimate = getEstimate(*track);
    file <<  " " << estimate.first << " " << estimate.second;

    vector<double> dedx;
    int j = 0;
    for(vector<double>::const_iterator hit = track->begin();
                                       hit!= track->end(); hit++)
      dedx.push_back(*hit / layers[j++].depth);

    file << " " << getTruncMean(dedx);
    file << " " << getPowerMean(dedx, -2.);
    file << " " << getPowerMean(dedx, -1.);
    file << " " << getPowerMean(dedx,  1.);

    file << endl;
  }

  file.close();
}

/****************************************************************************/
int main(int arg, char **arc)
{
  int pass = atoi(arc[1]);
  cerr << " pass = " << pass << endl;

  cerr << " initalizing layers..";
  initializeLayers(pass);
  cerr << " [" << nLayers << "]" << endl;

  cerr << " generating tracks..." << endl;
  generateTracks( nTracks, 0);
  cerr << " [" << nTracks << " done]" << endl;

  cerr << " generating tracks..." << endl;
  generateTracks( nTracks/3, 1);
  cerr << " [" << nTracks/3 << " done]" << endl;

  cerr << " altogether " << tracks.size() << endl;

  MostProbable mostProbable;
  eps = mostProbable.value(p / mass[0]);
  cerr << " epsilon = " << eps << endl;

  cerr << " calibrating gains [with pions only].." << endl;
  calibrateGains(pass);
  cerr << " [done]" << endl;

  cerr << " estimating energy loss..." << endl;
  estimateEloss(pass); 
  cerr << " [done]" << endl;
}
