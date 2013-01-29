#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;

#include "../../DataFormats/interface/TBunchCrossing.h"
#include "../../DataFormats/interface/TVertex.h"

#include "../interface/KalmanTracking.h"

#include "TVectorD.h"
#include "TMath.h"

#include "../../siEnergyLoss/interface/RecoClusters.h"
#include "../../siEnergyLoss/interface/ElossEstimator.h"

#include "../../minBiasVertexing/interface/PairGroupMethod.h"
#include "../../minBiasVertexing/interface/FastPairwiseNearestNeighbor.h"
#include "../../minBiasVertexing/interface/GaussianMixture.h"
#include "../../minBiasVertexing/interface/KMeansMethod.h"
#include "../../minBiasVertexing/interface/OptimalTree.h"


enum Classification { divisive, kMeans, gaussianMixture, optimalTree };

//#include "../interface/Hit.h"

#include "TFile.h"
#include "TTree.h"

const double sigmaZ = 5.;

/*****************************************************************************/
double getFlatRandom()
{
  return drand48();
}

/*****************************************************************************/
double getGaussRandom()
{
  return sqrt(-2*log(getFlatRandom())) * cos(M_PI*getFlatRandom());
}

/*****************************************************************************/
TVertex readGeneratedEvent(ifstream & fileEvents)
{
  int d, processId, n;

  fileEvents >> d;          // event number
  fileEvents >> processId;  // process id, Pythia
  fileEvents >> d;          // total number of particles
  fileEvents >> n;          // number of particles to follow

  // Initialize verrtex
  TVertex vertex;

  vertex.z = sigmaZ * getGaussRandom();
  vertex.processId = processId;

  // Add tracks
  for(int i = 0; i < n; i++)
  {
    int pid;
    double eta,pt;

    fileEvents >> pid; // particle id
    fileEvents >> eta; // eta
    fileEvents >> pt;  // transverse momentum

    if(fabs(eta) < 2.5 && pt > 0.1) // FIXME
    {
      TTrack track;

      track.pdgId  =  pid;
      track.charge = (pid > 0 ? 1 : -1);

      track.eta = eta;
      track.pt  = pt ;
      track.phi = drand48() * 2 * M_PI;

      track.z   = vertex.z; // FIMXE

      vertex.tracks.push_back(track);
    }
  }

  return vertex;
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

    // Stopping condition
    if(chi2 < getLambda(rec) ||
       TMath::Prob(chi2 - getLambda(rec), points.size()) > 1e-3)
    {
      for(int k = 0; k < K ; k++)
        cerr << "  vertex at " << mu(k) << endl;

      cerr << " vertices = " << K << endl
           << "\033[22;32m" << " chi2 ==> " << chi2 << "\033[22;0m" << endl;
 
      break;
    }
  }
}

/*****************************************************************************/
VertexCollection findVertices(const vector<pair<double, double> > & points)
{
  // Initialize clusters for advanced methods
  unsigned int maxVertices = (unsigned int)(100); // ???
  vector<pair<TVectorD,TVectorD> > clusters;
  for(unsigned int i = 0; i <= maxVertices; i++)
  {
    TVectorD mu(i);
    TVectorD  P(i);

    clusters.push_back(pair<TVectorD,TVectorD>(mu,P));
  }

cerr << " clusters " << clusters.size() << " " << points.size() << endl;

  cerr << " ---- hierarchical clustering ----" << endl;
  unsigned int nOptimal;
  vector<vector<int> > lists;

  int dMax = 8; // ??
  FastPairwiseNearestNeighbor thePairGroupMethod(dMax);
  thePairGroupMethod.run(points, clusters, nOptimal, lists, maxVertices);

  cerr << " ---- classification ---- "
       << " " << points.size() << " " << clusters.size()
       << endl;
  VertexCollection rec;

  // kMeans
  numberOfVertices(points, clusters, rec, kMeans);
  cerr << " kMeans " << rec.size() << endl;;

  // gaussianMixture
  numberOfVertices(points, clusters, rec, gaussianMixture);
  cerr << " gaussianMixture " << rec.size() << endl;

  // optimalTree
  int K = nOptimal;
  OptimalTree theOptimalTree;
  TVectorD mu(K); mu = clusters[K].first;
  TVectorD P(K) ; P  = clusters[K].second;

  theOptimalTree.run(K,points, mu,P, lists, rec);
  cerr << " optimalTree " << rec.size() << endl;

  return rec;
}

/*****************************************************************************/
int main(int arg, char **arc)
{
  TFile * file = new TFile("../out/bunchCrossings.root","recreate");
  TTree * tree = new TTree("trackTree","trackTree");

  TBunchCrossing * bunchCrossing = new TBunchCrossing();

  tree->Branch("bunchCrossing", "TBunchCrossing", &bunchCrossing, 16000, 2);

  // Take particles from Monte Carlo event generator
  ifstream fileEvents("../data/minBias.dat");
 
  // Initialize tracking
  KalmanTracking kalmanTracking;

  cerr << " simulating events..";

  RecoClusters recoClusters;
  ElossEstimator elossEstimator;

  ofstream fileEstimate("../out/estimates.dat");

//  do
  for(int i = 0; i < 1; i++) // number of bx!
  { 
    bunchCrossing->runNumber = 1;
    bunchCrossing->bxNumber  = i + 1;

    cerr << "\033[22;31m"
         << " Bunch crossing " << bunchCrossing->bxNumber
         << "\033[22;0m" << endl;

    vector<pair<double, double> > trackZ;

    for(int j = 0; j < 4; j++) // number of collisions per bx
    {
      cerr << "Vertex " << j + 1 << endl;
      // Read generted event
      TVertex simVertex = readGeneratedEvent(fileEvents);
     
      if(simVertex.tracks.size() > 0) // FIXME
      { 
      cerr << " b1" << endl;
    
      TVertex recVertex;

      cerr << " b2 " << simVertex.tracks.size() << endl;
      for(vector<TTrack>::iterator track = simVertex.tracks.begin();
                                   track!= simVertex.tracks.end(); track++)
      {
        TTrack recTrack;
        cerr << " Track " << int(track - simVertex.tracks.begin()) + 1 << endl;

//        if(kalmanTracking.process(*track, recTrack)) ??
        kalmanTracking.process(*track, recTrack); // FIXME empty track??
cerr << "  size " << recTrack.lhits.size() << endl;

        // for vertex finding
        trackZ.push_back(pair<double,double>(recTrack.z,0.1)); // FIXME sigma_z

        recVertex.tracks.push_back(recTrack);
      }

      // Do vertex finding! FIXME
      recVertex.z = simVertex.z;

      // Reconstruct clusters
      for(vector<TTrack>::iterator track = recVertex.tracks.begin();
                                   track!= recVertex.tracks.end(); track++)
      {
         cerr << " Track " << int(track - recVertex.tracks.begin()) + 1 << endl; 
//         RecoClusters recoClusters;
         recoClusters.run(*track);
//        siEnergyLoss.reconstruct(track); // or per hit?

         elossEstimator.estimate(*track, fileEstimate);
      }

      // Store 
      bunchCrossing->simVertices.push_back(simVertex);
      bunchCrossing->recVertices.push_back(recVertex);
      }
    }

    // Vertex finding
    VertexCollection recVertices;
    if(trackZ.size() > 0) // FIXME
      recVertices = findVertices(trackZ); 

    // update vertex positions FIXME
    for(vector<TVertex>::iterator 
        vertex = bunchCrossing->recVertices.begin();
        vertex!= bunchCrossing->recVertices.end(); vertex++)
    {
      float dzmin = 1e+99;
      vector<Vertex>::const_iterator bestVertex;
 
      for(vector<Vertex>::const_iterator recVertex = recVertices.begin();
                                         recVertex!= recVertices.end();
                                         recVertex++)
      if(fabs(vertex->z - recVertex->first) < dzmin)
      {
        dzmin = fabs(vertex->z - recVertex->first);
        bestVertex = recVertex;
      }

      cerr << " setting vertex z " << vertex->z
                                   << " -> " << bestVertex->first << endl;
      vertex->z = bestVertex->first;
    }

    // Fill tree
    tree->Fill();
    bunchCrossing->Clear();
  }
//  while(!fileEvents.eof());

  fileEstimate.close();

  fileEvents.close(); 

  file->Write();
  file->Close();

  cerr << " [done]" << endl;
}
