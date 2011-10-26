#include "../interface/PairGroupMethod.h"

#include "TMatrixD.h"

#define sqr(x) ((x) * (x))

#include <iostream>
#include <cmath>
using namespace std;

/*****************************************************************************/
void PairGroupMethod::run(const vector<pair<double,double> > & points,
                                vector<pair<TVectorD, TVectorD> > & clusters,
                                unsigned int maxClusters)
{
  vector<pair<double,double> > p(points);
  vector<int> n(p.size(), 1);

  unsigned int nUse = p.size() - 1;
  vector<bool> use(p.size(), true);

  // Initiate distance matrix
  TMatrixD dist(p.size(),p.size());
  for(unsigned int i = 0  ; i < p.size() - 1; i++)
  for(unsigned int j = i+1; j < p.size()    ; j++)
    dist(i,j) =  sqr(p[i].first  - p[j].first ) /
                    (p[i].second + p[j].second);

  while(nUse > 0)
  {
    unsigned int imin=1e+3, jmin=1e+3;
    bool isFirst = true;
    double dmin = 0.;

    for(unsigned int i = 0  ; i < p.size() - 1; i++)
    if(use[i])
    for(unsigned int j = i+1; j < p.size()    ; j++)
      if(use[j])
      if(dist(i,j) < dmin || isFirst)
      {
        dmin = dist(i,j);
        imin = i; jmin = j;

        isFirst = false;
      }

cerr << " nUse = " << nUse << " dmin = " << sqrt(dmin) << endl;

    // Join imin and jmin
    double x  = (p[imin].first/p[imin].second + p[jmin].first/p[jmin].second) /
                (           1./p[imin].second +            1./p[jmin].second);
    double s2 = 1 /
                (           1./p[imin].second  +           1./p[jmin].second);

    // Update imin
    p[imin] = pair<double,double>(x,s2);
    n[imin] = n[imin] + n[jmin];

    // Update distances
    for(unsigned int i = 0; i < imin; i++)
    if(use[i])
      dist(i,imin) =  sqr(p[i].first  - p[imin].first ) /
                         (p[i].second + p[imin].second);


    for(unsigned int j = imin+1; j < p.size(); j++)
    if(use[j])
      dist(imin,j) =  sqr(p[imin].first  - p[j].first ) /
                         (p[imin].second + p[j].second);

    // Erase jmin
    use[jmin] = false;

    // Save to clusters
    if(nUse <= maxClusters)
    {
      int k = 0;

      TVectorD mu(nUse);
      TVectorD P(nUse);

      for(unsigned int i = 0  ; i < p.size(); i++)
      if(use[i])
      {
        // average position
        mu(k) = p[i].first;

        // probability, with weight
        P(k) = float(n[i])/nUse;

        k++;
      }

      clusters[nUse] = pair<TVectorD,TVectorD>(mu,P);
    }
    nUse--;
  }
}

