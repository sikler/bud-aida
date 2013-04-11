#include "../interface/FitStripCluster.h"

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

#include "TDecompLU.h"

#define sqr(x) ((x) * (x))

using namespace std;

enum { isNormal, isBelow, isOver };

/*****************************************************************************/
FitStripCluster::FitStripCluster
  (const vector<pair<double,int> > & b_,
   const double alpha, const double sigma) :
   b(b_), Coupling(alpha), Sigma(sigma)
{
  iSigma  = 1./Sigma;

  a00 =       sqr(      Coupling) * sqr(iSigma);
  a11 =       sqr(1 - 2*Coupling) * sqr(iSigma);
  a01 = Coupling*(1 - 2*Coupling) * sqr(iSigma);

  a0 =        Coupling  * iSigma;
  a1 = (1 - 2*Coupling) * iSigma;

  npar = b.size();
}

/*****************************************************************************/
FitStripCluster::~FitStripCluster()
{
}

/*****************************************************************************/
double FitStripCluster::getChi2(const vector<double> & x)
{
  vector<double> y(npar, 0.);

  for(int i = 0; i < npar; i++)
  {
    double dx = Coupling * x[i];

    if(i-1 >= 0)   y[i-1] += dx;
                   y[i  ] += x[i] - 2*dx;
    if(i+1 < npar) y[i+1] += dx;
  }

  double chi2 = 0.;

  for(int i = 0; i < npar; i++)
  {
    double q = (b[i].first - y[i]) * iSigma;

    if(b[i].second == isNormal) {                   chi2 +=     q*q; }
    if(b[i].second == isBelow ) { q -= 2; if(q < 0) chi2 += 0.5*q*q; }
    if(b[i].second == isOver  ) { q += 2; if(q > 0) chi2 += 0.5*q*q; }

    // Penalty for negatives
    if(x[i] < 0)
    {
      q = x[i] * iSigma;
      chi2 += q*q;
    }
  }

  return chi2;
}

/*****************************************************************************/
void FitStripCluster::getAlphaBeta(const vector<double> & x,
                              TMatrixD & alpha,
                              TVectorD & beta,
                             const vector<bool> & isFix)
{
  vector<double> y(npar, 0.);

  for(int i = 0; i < npar; i++)
  if(!isFix[i])
  {
    double dx = Coupling * x[i];

    if(i-1 >= 0)   y[i-1] += dx;
                   y[i  ] += x[i] - 2*dx;
    if(i+1 < npar) y[i+1] += dx;
  }

  for(int i = 0; i < npar; i++)
  {
    double q = (y[i] - b[i].first) * iSigma;

    int f = 0;
    if(b[i].second == isNormal) {                   f = 2; }
    if(b[i].second == isBelow ) { q += 2; if(q > 0) f = 1; }
    if(b[i].second == isOver  ) { q -= 2; if(q < 0) f = 1; }

    if(f > 0)
    {
      if(i-1 >= 0)
      if(!isFix[i-1])
      {
        alpha[i-1][i-1] += f * a00;

        if(!isFix[i])
        {
          alpha[i-1][i  ] += f * a01;
          alpha[i  ][i-1] += f * a01;
        }

        beta[i-1] += f*q*a0;
      }

      if(!isFix[i  ])
      {
        alpha[i  ][i  ] += f * a11;

        beta[i  ] += f*q*a1;
      }

      if(i+1 < npar) 
      if(!isFix[i+1])
      {
        alpha[i+1][i+1] += f * a00;

        if(!isFix[i])
        {
          alpha[i+1][i  ] += f * a01;
          alpha[i  ][i+1] += f * a01;
        }

        beta[i+1] += f*q*a0;
      }
    }

    // Penalty for negatives
    if(!isFix[i  ])
    if(x[i] < 0)
    {
      alpha[i][i] += 2*sqr(iSigma);

      q = x[i] * iSigma;
      beta[i] += 2 * q*iSigma;
    }
  }

  // Calculate average of non-zero matrix elements (with fabs)
  double s=0., n=0.;

  for(int i = 0; i < npar; i++) 
  for(int j = 0; j < npar; j++)
    if(alpha[i][j] > 0)
    {
      s += fabs(alpha[i][j]);
      n++;
    }

  s /= n;

  for(int i = 0; i < npar; i++)
  if(isFix[i] || alpha[i][i] == 0.)
  {
    alpha[i][i] = s;

    beta[i] = 0.;
  }
}

/*****************************************************************************/
int FitStripCluster::run(double & chi2, pair<double,double> & result,
                         vector<double> & x)
{
  for(vector<pair<double,int> >::const_iterator ib = b.begin();
                                                ib!= b.end(); ib++) 
  {
    x.push_back(ib->first);

//    if(ib->second == isOver) over = true;
  }

  vector<bool> isFix(npar, false);

  bool ok = false;

  TMatrixD hessian(npar,npar);

  int iter = 0;
  
  do
  {
    double old = getChi2(x);
    double diff;

    do
    {
      TMatrixD alpha(npar,npar); alpha.Zero();
      TVectorD beta(npar); beta.Zero();

      getAlphaBeta(x, alpha,beta, isFix);

      // Fix alpha if singular
      if(alpha.Determinant() == 0.) 
      {
        for(int i = 0; i < npar; i++)
        for(int j = 0; j < npar; j++)
          if(i != j) alpha[i][j] = 0.;
      }

      TDecompLU lu(alpha);

      bool ok; TVectorD delta = lu.Solve(beta, ok);
  
      vector<double> x_(npar);
      double lambda = 1., next;

      for(int i = 0; i < npar; i++)
         x_[i] = x[i] - lambda * delta[i];

      next = getChi2(x_);

      diff = old - next;

      if(diff > 0) x = x_; // FIXME

      old = next;

      hessian = alpha;

      iter++;
    }
    while((diff > 1e-6) && (iter < 100));

    chi2 = old;

    // Check if we have negatives
    double min = 1e+9; int mi = -1;
    for(int i = 0; i < npar; i++)
      if(x[i] < min)
      { min = x[i] ; mi = i; }

    if(min < 0)
    {
      x[mi] = 0.; isFix[mi] = true;
    }
    else
      ok = true;
  }
  while(!ok && iter < 100);

  if(iter >= 100) cerr << "[" << iter << "]" << endl;

  hessian *= 1./2;
  hessian.Invert();

  double var2 = 0.;
  for(int i = 0; i < npar; i++)
  if(!isFix[i])
  for(int j = 0; j < npar; j++)
  if(!isFix[j])
    var2 += hessian[i][j];

  var2 = fabs(var2);

  double sum = 0.; int num = 0;
  for(int i = 0; i < npar; i++)
  if(!isFix[i])
  { sum += x[i]; if(x[i] > 0) num++; }

  result.first  = sum;
  result.second = var2;

  return num;
}
