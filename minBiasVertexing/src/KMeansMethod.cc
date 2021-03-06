#include "../interface/KMeansMethod.h"

#include <cmath>

#define sqr(x) ((x) * (x))

using namespace std;

/*****************************************************************************/
double KMeansMethod::logProb(const pair<double,double> & p,
                             double mu, double P)
{
  return - 0.5 * sqr(p.first - mu) / p.second + log(P);
}

/*****************************************************************************/
void KMeansMethod::estimateAverages
  (const vector<pair<double,double> > & points,
   const TMatrixD & p,
   TVectorD & mu, TVectorD & P)
{
  int K = mu.GetNrows();

  // Recalculate averages
  for(int k = 0; k < K; k++)
  {
    double sum[3] = {0, 0, 0};

    for(unsigned int n = 0; n < points.size(); n++)
    {
      sum[2] += p(n,k) * points[n].first / points[n].second;
      sum[1] += p(n,k)                   / points[n].second;
      sum[0] += p(n,k);
    }

    mu(k) = sum[2]/sum[1];
     P(k) = sum[0]/points.size();
  }
}

/*****************************************************************************/
double KMeansMethod::estimateChiSquare
  (const vector<pair<double,double> > & points,
   const TVectorD & mu, const TVectorD & P,
   TMatrixD & p, bool estimateResponsibility)
{
  int K = mu.GetNrows();

  double chi2 = 0;
  for(unsigned int n = 0; n < points.size(); n++)
  {
    vector<double> logP(K);

    bool isFirst = true;
    double maxLogP = 0;
    int maxK = -1;

    for(int k = 0; k < K; k++)
    if(P(k) > 0)
    {
      logP[k] = logProb(points[n],mu(k),P(k));

      if(logP[k] > maxLogP || isFirst)
      { maxLogP = logP[k]; maxK = k; isFirst = false; }
    }

    if(estimateResponsibility)
      for(int k = 0; k < K; k++)
        if(k == maxK) p(n,k) = 1;
                 else p(n,k) = 0;
    else
      chi2 += -2 * maxLogP;
  }

  return chi2;
}

/*****************************************************************************/
double KMeansMethod::getChiSquare
  (const vector<pair<double,double> > & points,
   const TVectorD & mu, const TVectorD & P)
{
  TMatrixD p;

  return estimateChiSquare(points,mu,P, p,false);
}

/*****************************************************************************/
void KMeansMethod::estimateResponsibility
  (const vector<pair<double,double> > & points,
   const TVectorD & mu, const TVectorD & P,
   TMatrixD & p)
{
  estimateChiSquare(points,mu,P, p,true);
}

/*****************************************************************************/
double KMeansMethod::run
  (int K,
   const vector<pair<double, double> > & points,
   TVectorD mu, TVectorD P,
   VertexCollection & vertices)
{
  // Clear vertices
  vertices.clear();

  int N = points.size();

  TMatrixD p(N,K);

  double chi2 = 0.;
  double old_chi2;

  int iter = 0;

  do
  {
    estimateResponsibility(points,mu,P, p);
    estimateAverages      (points,p, mu,P);

    old_chi2 = chi2;
    chi2 = getChiSquare(points,mu,P);
  }
  while(fabs(chi2 - old_chi2) > 1e-3 && ++iter < 100);

  // Fill
  for(int k = 0; k < K; k++)
  {
    Vertex vertex;

    vertex.first = mu(k);

    for(int i = 0; i < N; i++)
    {
      double pmax = 0.;
      int kmax = -1;

      for(int j = 0; j < K; j++)
        if(p(i,j) > pmax)
        { pmax = p(i,j); kmax = j; }

      if(k == kmax)
        vertex.second.push_back(i);
    }

    vertices.push_back(vertex);
  }

  return chi2;
}

