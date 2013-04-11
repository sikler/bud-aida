#include "../interface/VertexFinder.h"

#include "../interface/PairGroupMethod.h"
#include "../interface/FastPairwiseNearestNeighbor.h"
#include "../interface/GaussianMixture.h"
#include "../interface/KMeansMethod.h"
#include "../interface/OptimalTree.h"

#include "TMath.h"

#include <iostream>
#include <cmath>

using namespace std;

enum Generator { fix, single, poissonian };
enum Classification { divisive, kMeans, gaussianMixture, optimalTree };

#define nMethods 4
const char *methodName[nMethods] = { "div","kme","gau","opt" };

// FIXME
#define nAssoc 3
#define nMulti 500
#define maxVertex 100
#define nMin 2

/*****************************************************************************/
VertexFinder::VertexFinder()
{
  vector<int> v(nAssoc,0);

  for(int i = 0; i < nMethods; i++) // FIXME
  {
    vector<vector<int> > a;
    for(int j = 0; j < maxVertex; j++)
      a.push_back(v);

    asim.push_back(a);
    arec.push_back(a);
  }

  for(int i = 0; i < nMulti; i++) // FIXME
  {
    msim0.push_back(v); // no need!!! FIXME
    msim1.push_back(v);
    msim2.push_back(v);
    msim3.push_back(v);
  }


  fileMacro.open("../out/macro.dat");
}

/*****************************************************************************/
VertexFinder::~VertexFinder()
{
  ofstream fileRes("../out/result.dat");

  for(int method = kMeans; method <= optimalTree; method++)
  for(int j = 0; j < maxVertex; j++)
  {
    double nsim=0., nrec=0.;
    for(int i = 0; i < nAssoc; i++)
    {
      nsim += asim[method][j][i];
      nrec += arec[method][j][i];
    }

    for(int i = 0; i < nAssoc; i++)
      fileRes << " " << method // methodName[method]
              << " " << j // Mu
              << " " << i
              << " " << asim[method][j][i] / nsim
              << " " << arec[method][j][i] / nrec << endl;

    fileRes << endl << endl;
  }

  fileRes.close();

  fileRes.open("../out/multi.dat");
  for(int k = 0; k < nMulti; k++)
  {
    fileRes << " " << k
            << " " << float(msim0[k][1] + msim0[k][2]) /
                           (msim0[k][0] + msim0[k][1] + msim0[k][2])
            << " " << float(msim1[k][1] + msim1[k][2]) /
                           (msim1[k][0] + msim1[k][1] + msim1[k][2])
            << " " << float(msim2[k][1] + msim2[k][2]) /
                           (msim2[k][0] + msim2[k][1] + msim2[k][2])
            << " " << float(msim0[k][2]) /
                           (msim0[k][0] + msim0[k][1] + msim0[k][2])
            << " " << float(msim1[k][2]) /
                           (msim1[k][0] + msim1[k][1] + msim1[k][2])
            << " " << float(msim2[k][2]) /
                           (msim2[k][0] + msim2[k][1] + msim2[k][2])
            << " " << float(msim3[k][1] + msim3[k][2]) /
                           (msim3[k][0] + msim3[k][1] + msim3[k][2])
            << " " << float(msim3[k][2]) /
                           (msim3[k][0] + msim3[k][1] + msim3[k][2])
            << endl;
  }
}

/*****************************************************************************/
void VertexFinder::evaluatePerformance(const VertexCollection & sim,
                                       const VertexCollection & rec,
                                       vector<int> & asim,
                                       vector<int> & arec,
                                       vector<vector<int> > & msim, int method)
{
  // Look at simvertices
  for(VertexCollection::const_iterator s = sim.begin();
                                       s!= sim.end(); s++)
  {
    double dz = 1e+99;
    double fraction = -1;
    int loss = -1;

    int nassoc = 0;
    int nallsim = s->second.size();

    for(VertexCollection::const_iterator r = rec.begin();
                                         r!= rec.end(); r++)
    {
      int nallrec = r->second.size();
      int nshared = 0;

      for(vector<int>::const_iterator ir = r->second.begin();
                                      ir!= r->second.end(); ir++)
      for(vector<int>::const_iterator is = s->second.begin();
                                      is!= s->second.end(); is++)
        if(*ir == *is)
          nshared++;

      if(nshared > nallrec / 2)
      {
        nassoc++;
        dz = s->first - r->first;
        fraction = float(nshared) / nallsim;
        loss = nallsim - nshared;
      }
    }

    if(nallsim >= nMin) // accepted
    {
      if(nassoc < nAssoc)
        asim[nassoc]++;
    }

//    if(sim.size() == 1) // FIXME
      if(nallsim < nMulti && nassoc < nAssoc)
        msim[nallsim][nassoc]++;

    if(nassoc == 1)
    {
       fileMacro << methodName[method] << " " << sim.size() // FIXME Mu
                 << " " << dz << " " << fraction << " " << loss << endl;
    }
  }
//while(getchar() == 0);

  // Look at recvertices
  for(VertexCollection::const_iterator r = rec.begin();
                                       r!= rec.end(); r++)
  {
    int nallrec = r->second.size();

    int nassoc = 0;

    for(VertexCollection::const_iterator s = sim.begin();
                                         s!= sim.end(); s++)
    {
      int nallsim = s->second.size();
      int nshared = 0;

      for(vector<int>::const_iterator ir = r->second.begin();
                                      ir!= r->second.end(); ir++)
      for(vector<int>::const_iterator is = s->second.begin();
                                      is!= s->second.end(); is++)
        if(*ir == *is)
          nshared++;

      if(nshared > nallsim / 2 && nallsim >= nMin) // accepted
        nassoc++;
    }

    if(nallrec >= nMin) // reasonable
      if(nassoc < nAssoc)
        arec[nassoc]++;
  }
}

/*****************************************************************************/
double VertexFinder::getLambda(const VertexCollection & vertices)
{
  int n = 0;
  for(VertexCollection::const_iterator vertex = vertices.begin();
                                       vertex!= vertices.end(); vertex++)
    n += vertex->second.size();

  double chi2 = 0.;
  for(VertexCollection::const_iterator vertex = vertices.begin();
                                       vertex!= vertices.end(); vertex++)
  {
    int k = vertex->second.size();

    if(k > 0)
      chi2 += -2 * k * log(float(k)/n);
  }

  return chi2;
}

/*****************************************************************************/
void VertexFinder::numberOfVertices(const vector<pair<double,double> > & points,
                      const vector<pair<TVectorD,TVectorD> > & clusters,
                      VertexCollection & rec,
                      int classification)
{
  rec.clear();

  for(int K = 1; K <= int(points.size()); K++)
  {
    // Get results from hierarchical clustering
    TVectorD mu(K); mu = clusters[K].first;
    TVectorD P(K) ; P  = clusters[K].second;

    double chi2 = 0.;

    if(classification == kMeans)
    {
      KMeansMethod theKMeansMethod;
      chi2 = theKMeansMethod.run(K,points, mu,P, rec);
    }

    if(classification == gaussianMixture)
    {
      GaussianMixture theGaussianMixture;
      chi2 = theGaussianMixture.run(K,points, mu,P, rec);
    }

    // Stopping condition
    if(chi2 < getLambda(rec) ||
       TMath::Prob(chi2 - getLambda(rec), points.size()) > 1e-3)
    {
      for(int k = 0; k < K ; k++)
        cerr << "  vertex at " << mu(k) << endl;

      cerr << "  found vertices = " << K
           << " | chi2 = " << chi2 << endl;

      break;
    }
  }
}

/********************* ********************************************************/
VertexCollection VertexFinder::findVertices
  (const vector<pair<double, double> > & points,
   const VertexCollection & sim) // FIXME
{
  cerr << "\033[22;32m"
       << "  Vertex finding"
       << "\033[22;0m" << endl;

  // Initialize clusters for advanced methods
  unsigned int maxVertices = (unsigned int)(100); // ???
  vector<pair<TVectorD,TVectorD> > clusters;
  for(unsigned int i = 0; i <= maxVertices; i++)
  {
    TVectorD mu(i);
    TVectorD  P(i);

    clusters.push_back(pair<TVectorD,TVectorD>(mu,P));
  }

  cerr << "  ---- hierarchical clustering ----" << endl;
  unsigned int nOptimal;
  vector<vector<int> > lists;

  int dMax = 8; // ??
  FastPairwiseNearestNeighbor thePairGroupMethod(dMax);
  thePairGroupMethod.run(points, clusters, nOptimal, lists, maxVertices);

  cerr << "  ---- classification ---- "
       << " " << points.size() << " " << clusters.size()
       << endl;
  VertexCollection rec;

  // kMeans
  numberOfVertices(points, clusters, rec, kMeans);
  cerr << "  kMeans " << rec.size() << endl;
  evaluatePerformance(sim,rec, asim[kMeans][sim.size()],
                               arec[kMeans][sim.size()], msim1, kMeans);

  // gaussianMixture
  numberOfVertices(points, clusters, rec, gaussianMixture);
  cerr << "  gaussianMixture " << rec.size() << endl;
  evaluatePerformance(sim,rec, asim[gaussianMixture][sim.size()],
                               arec[gaussianMixture][sim.size()], msim2,
                               gaussianMixture);

  // optimalTree
  int K = nOptimal;
  OptimalTree theOptimalTree;
  TVectorD mu(K); mu = clusters[K].first;
  TVectorD P(K) ; P  = clusters[K].second;

  theOptimalTree.run(K,points, mu,P, lists, rec);
  cerr << "  optimalTree " << rec.size() << endl;
  evaluatePerformance(sim,rec, asim[optimalTree][sim.size()],
                               arec[optimalTree][sim.size()], msim3, optimalTree);

  return rec;
}


