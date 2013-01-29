#include "../interface/ElossModel.h"

#include "../interface/TouchedChannels.h"
#include "../interface/Levels.h"

#define sqr(x) ((x) * (x))

#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace CLHEP;

/*****************************************************************************/
ElossModel::ElossModel
  (Hit hit_) : hit(hit_)
{
}

/*****************************************************************************/
ElossModel::~ElossModel()
{
}

/*****************************************************************************/
void ElossModel::preparePixels
  (const vector<double>& pars)
{
  // Calculate endpoints
  double endpoint[2][2];
  for(int i=0; i<2; i++)
  for(int k=0; k<2; k++)
    endpoint[i][k] = pars[k] + (2*i-1) * 0.5 * hit.dpos[k];

  // Prepare pixels
  TouchedChannels theTouchedChannels(hit);
  hit.allPixels = theTouchedChannels.findChannels(endpoint);
}

/*****************************************************************************/
Hit ElossModel::getHit(const vector<double>& pars)
{
  // prepare pixels
  preparePixels(pars);

  getValue(pars); // also, fill theo

  return hit;
}

/*****************************************************************************/
double ElossModel::getSigma(double y)
{
  return Sigma0 + b * y;
}

/*****************************************************************************/
double ElossModel::getChi2(double y, double Delta,
                           double thr, double ovf,
                           vector<double> & der)
{
  double sigma, limit;

  if(y >= thr && y <= ovf)
  { // normal
    sigma = getSigma(y);
    limit = y - Nu*sigma;

    if(Delta < limit)
    {
      der[0] = -2*Nu*(Delta-y)/sigma - Nu*Nu;
      der[1] = -2*Nu/sigma;
    }
    else
    {
      der[0] = sqr((Delta - y)/sigma);
      der[1] =   2*(Delta - y)/sqr(sigma);
      der[2] =   2            /sqr(sigma);
    }
  }
  else
  {
    if(y < thr)
    { // left truncation, threshold
      sigma = getSigma(thr);
      limit = thr - sigma;

      if(Delta < limit)
      { 
      }
      else
      {  
        der[0] = sqr((Delta - thr)/sigma + 1);
        der[1] =  2*((Delta - thr)/sigma + 1)/sigma;
        der[2] =  2                      /sqr(sigma);
      }
    }
    else
    { // right censoring, saturation
      sigma = getSigma(ovf);
      limit = ovf + sigma;

      if(Delta < limit)
      {
        der[0] = -(Delta - ovf)/sigma + 1;
        der[1] = -1            /sigma; 
      }
      else
      {
      }
    }
  }

  return der[0];
}

/*****************************************************************************/
void ElossModel::calculate(const vector<double>& pars, double & val,
                           HepVector& beta, HepMatrix& alpha)
{
  val = 0.;

  for(int k = 0; k < 3; k++)
    beta[k] = 0.;

  for(int j = 0; j < 3; j++)
  for(int k = 0; k < 3; k++)
    alpha[j][k] = 0.;

  // prepare pixels
  preparePixels(pars);

  // take all pixels
  for(vector<Pixel>::iterator pixel = hit.allPixels.begin();
                              pixel!= hit.allPixels.end(); pixel++)
  {
    Pixel::Calc c = pixel->calc;

    double y       = pixel->meas.y;

    double epsilon = exp(pars[2]);
    double l       = c.l;
    double Delta   = epsilon * l; 
 
    vector<double> der(3, 0.);
    getChi2(y,Delta, hit.threshold, hit.overflow, der);

    // Value
    val += der[0];

    // Beta
    for(int k = 0; k < 2; k++)
      beta[k] += - der[1] * epsilon * c.dl_dP[k];

    beta[2] += - der[1] * Delta;
    
    // Alpha
    for(int j=0; j<2; j++)
    for(int k=0; k<2; k++)     
      alpha[j][k] += der[2] * sqr(epsilon) * c.dl_dP[j] * c.dl_dP[k];

    for(int k=0; k<2; k++)     
      alpha[k][2] += (der[2] * Delta + der[1]) * c.dl_dP[k] * epsilon;

    alpha[2][2] += der[2] * sqr(Delta) + der[1] * Delta;
  }

  // Mirror
  for(int k=0; k<2; k++)     
    alpha[2][k] = alpha[k][2];
}

/*****************************************************************************/
double ElossModel::getDerivatives
  (const vector<double>& pars, HepVector & beta, HepMatrix & alpha)
{
  double val;

  calculate(pars, val, beta, alpha);

  return val;
}

/*****************************************************************************/
double ElossModel::getValue(const vector<double>& pars)
{
  double val;
  int dim = pars.size();
  HepMatrix alpha(dim,dim);
  HepVector beta(dim,0.);

  calculate(pars, val, beta, alpha);

  return val;
}
