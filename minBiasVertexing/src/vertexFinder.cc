#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include "TVectorD.h"
#include "TMath.h"

#include "../interface/PairGroupMethod.h"
#include "../interface/FastPairwiseNearestNeighbor.h"
#include "../interface/GaussianMixture.h"
#include "../interface/KMeansMethod.h"
#include "../interface/OptimalTree.h"

#undef Debug

#define sqr(x) ((x) * (x))

using namespace std;

#define SigmaZ 5.0

double Mu;
unsigned int Multi;

// z position, sequential ids of tracks
typedef pair<float, vector<int> > Vertex;

typedef vector<Vertex> VertexCollection;

enum Generator { fix, single, poissonian };
Generator generator;

enum Classification { divisive, kMeans, gaussianMixture, optimalTree };

#define nMethods 4
const char *methodName[nMethods] = { "div","kme","gau","opt" };

ifstream filePythia;
ofstream fileMicro;
ofstream fileMacro;

#define nAssoc 3
#define nMulti 500

bool optimize;
bool tree;
float zOffset;
int nTrkMin;
float zSeparation;

int dMax;

bool sensitivity;
char resultName[128];
float background, sysShift, ranShift;
int nMin;

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

/*****************************************************************************/
double poisson(double mu, double n)
{
  return pow(mu,n)*exp(-mu - lgamma(n+1));
}

/*****************************************************************************/
int getPoissonRandom(double mu)
{
  int n;
  double p;

  do
  {
    n = int(getFlatRandom() * 10*mu) + 1;
    p = poisson(mu,n);
  }
  while(getFlatRandom() > p);

  return n;
}

/*****************************************************************************/
vector<pair<double,double> > readEvent()
{
  int d, processId, n;

  filePythia >> d;
  filePythia >> processId;
  filePythia >> d;
  filePythia >> n;

  vector<pair<double,double> > tracks;

  for(int i = 0; i < n; i++)
  {
    int pid;
    double eta,pt;

    filePythia >> pid;
    filePythia >> eta;
    filePythia >> pt;

    if(fabs(eta) < 2.5 && pt > 0.1)
    {
      pair<double,double> track(eta,pt);
      tracks.push_back(track);
    }  
  }

  return tracks;
}

/*****************************************************************************/
void generatePoints(int K, vector<pair<double,double> > & points,
                    VertexCollection & vertices)
{
  vertices.clear();

#ifdef Debug
  ofstream file ("../out/z.dat");
  ofstream file0("../out/z0.gnu");
  cerr << " generate K = " << K << endl;
#endif

  // For each vertex
  for(int k = 0; k < K; k++)
  {
    vector<pair<double,double> > tracks;

    // Read Pythia event
    if(generator != single)
      tracks = readEvent();
    else
    { 
      do
        tracks = readEvent();
      while( tracks.size() <= Multi-5 || tracks.size() >  Multi+5 );
    }

    // Multiplicity
    int n = tracks.size();

    // Generate vertex
    double z0 = getGaussRandom() * SigmaZ;

#ifdef Debug
    if(n > 0)
    file0 << " set arrow from " << z0 << ", graph 0 to "
                                << z0 << ", graph 1 nohead lt 3" << endl;

    cerr << " z0 = " << z0 << " | n = " << n << endl;
#endif

    Vertex vertex;
    vertex.first = z0;

    // For each track
    for(int i = 0; i < n; i++)
    {
      // Eta
      double eta = tracks[i].first;

      // Pt
      double pt = tracks[i].second;

      // Sigma of reconstructed z
      double sigma = sqrt(sqr(50e-4) + sqr(100e-4/pt) * pow(cosh(eta),3.));

      // Generate reconstructed z0
      double z_vtx = z0 + sigma * getGaussRandom();

      // Modify assumed sigma (for sensitivity study)
      sigma = fabs(sigma * (1 + sysShift + ranShift * getGaussRandom()));

      // Put in
      vertex.second.push_back(points.size());
      points.push_back(pair<double,double>(z_vtx,sigma*sigma));

#ifdef Debug
      file << " " << z_vtx << " " << sigma << endl;
#endif

      // Add background
      if(getFlatRandom() < background)
      {
//        vertex.second.push_back(points.size());

        z_vtx = getGaussRandom() * SigmaZ;
        points.push_back(pair<double,double>(z_vtx,sigma*sigma));

#ifdef Debug
      file << " " << z_vtx << " " << sigma << endl;
#endif
      }
    } 

    vertices.push_back(vertex);
  }

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
      chi2 += k*(1 - 2 * log(float(k)/n));
  }

#ifdef Debug
  cerr << "\033[22;31m"
       << " chi2 --> " << chi2 << " +- " << sqrt(2*n)
       << "\033[22;0m" 
       << endl;

  file0.close();
  file.close();
#endif
}

/*****************************************************************************/
double getChi2(const VertexCollection & vertices)
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
      chi2 += k*(1 - 2 * log(float(k)/n));
  }

  return chi2;
}

/*****************************************************************************/
double getLambda(const VertexCollection & vertices)
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
void evaluatePerformance(const VertexCollection & sim,
                         const VertexCollection & rec,
                         vector<int> & asim,
                         vector<int> & arec,
                         vector<vector<int> > & msim, int method)
{
#ifdef Debug
  cerr << endl
       << "  sim = " << sim.size()
        << " rec = " << rec.size() << endl;
#endif

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
 
    if(sim.size() == 1)
      if(nallsim < nMulti && nassoc < nAssoc)
        msim[nallsim][nassoc]++; 

    if(nassoc == 1)
    {
      if(generator == single)
        fileMicro << methodName[method] << " " << Multi
                  << " " << dz << " " << fraction << " " << loss << endl; 

      if(generator == fix)
        fileMacro << methodName[method] << " " << Mu
                  << " " << dz << " " << fraction << " " << loss << endl; 
    } 
  }

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
void numberOfVertices(const vector<pair<double,double> > & points,
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

#ifdef Debug
    cerr << " chi2 --> " << chi2 << " / " <<   chi2 / getChi2(rec)
         << " " << chi2 /  (3*sqrt(2*points.size()) + getChi2(rec))
         << " p = " << (chi2 < getLambda(rec) ? 1 :
                        TMath::Prob(chi2 - getLambda(rec), points.size()))
         << endl;
#endif

    // Stopping condition
    if(chi2 < getLambda(rec) ||
       TMath::Prob(chi2 - getLambda(rec), points.size()) > 1e-3)
    {
#ifdef Debug
      ofstream file("../out/z.gnu");
      for(int k = 0; k < K ; k++)
        file << "set arrow from " << mu(k) << ",graph -0.03 to "
                                  << mu(k) << ",graph -0.01 lt 2" << endl;
      file.close();
 
      cerr << " vertices = " << K << endl
           << "\033[22;32m" << " chi2 ==> " << chi2 << "\033[22;0m" << endl;
#endif
      break;
    }
  }
}

/*****************************************************************************/
void options(int arg, char **arc)
{
  int i = 1;

  optimize = false;
  tree     = false;

  zOffset     = 3.0;
  nTrkMin     = 3;
  zSeparation = 0.3;

  dMax = 8;

  sprintf(resultName,"../out/result.dat");
  sensitivity = false;
  background = 0.;
  sysShift   = 0.;
  ranShift   = 0.;

  //
  nMin = 2;

  do
  {
    if(strcmp(arc[i],"-fix")     == 0) 
    {
      generator = fix;
      Mu = atof(arc[++i]);
    }

    if(strcmp(arc[i],"-poissonian") == 0)
    {
      generator = poissonian;
      Mu = atof(arc[++i]);
    }

    if(strcmp(arc[i],"-single")     == 0)
    {
      generator = single;
      Mu = 1.;
      Multi = atoi(arc[++i]);
    }

    if(strcmp(arc[i],"-file") == 0)
    {
      sensitivity = true;
      sprintf(resultName,"%s",arc[++i]);
    }

    if(strcmp(arc[i],"-bck") == 0)
      background = atof(arc[++i]);

    if(strcmp(arc[i],"-sys") == 0)
      sysShift = atof(arc[++i]);

    if(strcmp(arc[i],"-ran") == 0)
      ranShift     = atof(arc[++i]);

    if(strcmp(arc[i],"-nmin") == 0)
      nMin         = atoi(arc[++i]);

    // float zOffset, int nTrkMin, float zSeparation
    if(strcmp(arc[i],"-pars") == 0)
    {
      optimize = true;
      zOffset     = atof(arc[++i]);
      nTrkMin     = atoi(arc[++i]);
      zSeparation = atof(arc[++i]);
    }
  
    if(strcmp(arc[i],"-tree") == 0)
    {
      tree = true;
      dMax = atoi(arc[++i]);
    }

    i++;
  }
  while(i < arg);
}

/*****************************************************************************/
int main(int arg, char **arc)
{
  options(arg,arc);

  cerr << fixed << setprecision(3);

  filePythia.open("../data/minBias.dat");

  vector<vector<int> > asim(nMethods,vector<int>(nAssoc,0));
  vector<vector<int> > arec(nMethods,vector<int>(nAssoc,0));

  vector<vector<int> > msim0(nMulti,vector<int>(nAssoc,0));
  vector<vector<int> > msim1(nMulti,vector<int>(nAssoc,0));
  vector<vector<int> > msim2(nMulti,vector<int>(nAssoc,0));
  vector<vector<int> > msim3(nMulti,vector<int>(nAssoc,0));

  cerr << " Mu = " << Mu << " ";
  if(generator == single) cerr << " single = " << Multi << " ";

  int events;
  if(generator != single) events = 1e+4;
                     else events = 1e+5 / Multi;
  if(optimize || tree)    events = 1e+4;
  if(sensitivity) events = 1e+4;

  if(!optimize && !tree && !sensitivity)
  {
    if(generator == single)
      fileMicro.open("../out/micro.dat"); //,ios::app);

    if(generator == fix)
      fileMacro.open("../out/macro.dat"); //,ios::app);
  }

  for(int ii = 0; ii < events / Mu; ii++)
  {
    if(ii % 200  == 0) cerr << ".";

    // Generate
    vector<pair<double, double> > points;

    // Vertex collections
    VertexCollection sim;

    // Generate vertices and tracks
#ifdef Debug
    cerr << "--------------------------------------------------------" << endl;
#endif

    do
    {
      if(generator == poissonian)
        generatePoints(getPoissonRandom(Mu), points, sim);

      if(generator == fix)
        generatePoints(int(Mu + 0.5), points, sim);

      if(generator == single)
        generatePoints(int(Mu + 0.5), points, sim);
    }
    while(points.size() == 0);

    // Initialize clusters for advanced methods
    unsigned int maxVertices = (unsigned int)(10*Mu + 10);
    vector<pair<TVectorD,TVectorD> > clusters;
    for(unsigned int i = 0; i <= maxVertices; i++)
    {
      TVectorD mu(i);
      TVectorD  P(i);

      clusters.push_back(pair<TVectorD,TVectorD>(mu,P));
    }

    // Advanced methods
    if(!optimize)
    {
#ifdef Debug
      cerr << " ---- hierarchical clustering ----" << endl;
#endif
      unsigned int nOptimal;
      vector<vector<int> > lists;

      FastPairwiseNearestNeighbor thePairGroupMethod(dMax);
      thePairGroupMethod.run(points, clusters, nOptimal, lists, maxVertices);

#ifdef Debug
      cerr << " ---- classification ---- "
           << " " << points.size() << " " << clusters.size()
           << endl;
#endif
      VertexCollection rec;

      // kMeans
      numberOfVertices(points, clusters, rec, kMeans);
#ifdef Debug
      cerr << " kMeans";
#endif
      evaluatePerformance(sim,rec, asim[kMeans],
                                   arec[kMeans], msim1, kMeans);

      // gaussianMixture
      numberOfVertices(points, clusters, rec, gaussianMixture);
#ifdef Debug
      cerr << " gaussianMixture";
#endif
      evaluatePerformance(sim,rec, asim[gaussianMixture],
                                   arec[gaussianMixture], msim2,
                                   gaussianMixture);

      // optimalTree
      int K = nOptimal;
      OptimalTree theOptimalTree;
      TVectorD mu(K); mu = clusters[K].first;
      TVectorD P(K) ; P  = clusters[K].second;

#ifdef Debug
      cerr << " optimalTree";
#endif
      theOptimalTree.run(K,points, mu,P, lists, rec);
      evaluatePerformance(sim,rec, asim[optimalTree],
                                   arec[optimalTree], msim3, optimalTree);
    }

#ifdef Debug
    int ret = system("cd ../gnu ; gnuplot event.gnu"); ret++;
    cerr << " Look at ../gnu/event.eps..";
    while(getchar() == 0);
#endif
  }

  cerr << endl;

  filePythia.close();

  if(generator == single && !optimize && !tree && !sensitivity)
    fileMicro.close();

  if(generator == fix    && !optimize && !tree && !sensitivity)
    fileMacro.close();

  // Optimize
  if(optimize)
  {
    ofstream fileOpt("../out/optimize.dat"); //,ios::app);
    fileOpt
         << " " << int(Mu + 0.5)
         << " " << zOffset << " " << nTrkMin << " " << zSeparation
         << " " << float(asim[divisive][0] + asim[divisive][2]) /
    (asim[divisive][0] + asim[divisive][1] + asim[divisive][2])
         << " " << float(arec[divisive][0] + arec[divisive][2]) /
    (arec[divisive][0] + arec[divisive][1] + arec[divisive][2])
         << endl;
    fileOpt.close();
  }
  
  // Tree
  if(tree)
  {
    ofstream fileOpt("../out/tree.dat"); //,ios::app);
    fileOpt
         << " " << int(Mu + 0.5)
         << " " << dMax
         << " " << float(   asim[optimalTree][0] + asim[optimalTree][2]) /
    (asim[optimalTree][0] + asim[optimalTree][1] + asim[optimalTree][2])
         << " " << float(   arec[optimalTree][0] + arec[optimalTree][2]) /
    (arec[optimalTree][0] + arec[optimalTree][1] + arec[optimalTree][2])
         << endl;
    fileOpt.close();
  }

  // Write out
  if(generator != single && !optimize && !tree)
  {
  ofstream fileRes(resultName); //,ios::app);

  for(int method = kMeans; method <= optimalTree; method++)
  {
    double nsim=0., nrec=0.;
    for(int i = 0; i < nAssoc; i++)
    {
      nsim += asim[method][i];
      nrec += arec[method][i];
    }

    for(int i = 0; i < nAssoc; i++)
      fileRes << " " << methodName[method]
              << " " << Mu 
              << " " << i
              << " " << asim[method][i] / nsim
              << " " << arec[method][i] / nrec << endl;

    fileRes << endl << endl;
  }

  fileRes.close();
  }

  if(Mu < 1.5 && generator == fix && !optimize && !tree && !sensitivity)
  {
  ofstream fileRes("../out/multi.dat");
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

  return 0;
}
