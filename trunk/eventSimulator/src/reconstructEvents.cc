#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;

#include "../interface/KalmanTracking.h"

#include "../../DataFormats/interface/Particles.h"
#include "../../DataFormats/interface/TBunchCrossing.h"
#include "../../DataFormats/interface/TVertex.h"
#include "../../DataFormats/interface/TSlimMeasurement.h"

#include "../../siEnergyLoss/interface/ClusterReco.h"
#include "../../siEnergyLoss/interface/ElossEstimator.h"
#include "../../siEnergyLoss/interface/MostProbable.h"

#include "../../minBiasVertexing/interface/VertexFinder.h"

#include "TFile.h"
#include "TTree.h"

char inputName[256], outputName[256];

enum { asPion, ifIdentified, notNeeded };
int collectHits;

/*****************************************************************************/
void printPropeller(int i)
{
  char c[4] = {'|', '/','-','\\'};

  if(i < 40) cerr << ".";
        else cerr << c[(i-40) % 4] << "\b";
}

/*****************************************************************************/
void options(int arg, char **arc)
{
  int i = 1;

  // default
  collectHits = notNeeded;

  do
  {
    if(strcmp(arc[i],"-i") == 0) sprintf(inputName ,"%s",arc[++i]);
    if(strcmp(arc[i],"-o") == 0) sprintf(outputName,"%s",arc[++i]);
    if(strcmp(arc[i],"-collectHits") == 0)
    {
      i++;

      if(strcmp(arc[i],"asPion")       == 0) collectHits = asPion;
      if(strcmp(arc[i],"ifIdentified") == 0) collectHits = ifIdentified;
    }

    i++;
  }
  while(i < arg);
}

/*****************************************************************************/
void reconstructEvents()
{
  // Input simu data
  TFile   * fileIn   = new TFile(inputName,  "read");
  TTree   * treeIn   = (TTree *) fileIn->Get("trackTree");
  TBranch * branchIn = treeIn->GetBranch("bunchCrossing");

  TBunchCrossing * bunchCrossing = new TBunchCrossing();
  branchIn->SetAddress(&bunchCrossing);

  // Input materials
  TTree   * treeMat   = (TTree *) fileIn->Get("materialTree");
  TBranch * branchMat = treeMat->GetBranch("materials");

  TBranch * branchB   = treeMat->GetBranch("magneticField");

  vector<TLayer> * materials = new vector<TLayer>;
  branchMat->SetAddress(&materials);
  branchMat->GetEntry();

  // Output reco data
  TFile * fileOut    = new TFile(outputName, "recreate");

  TTree * treeOut    = new TTree("recTrackTree","recTrackTree");
  treeOut->Branch("bunchCrossing", "TBunchCrossing", &bunchCrossing, 16000, 2);

  TTree * treeMaterialsOut = new TTree("materialTree","materialTree");
  treeMaterialsOut->Branch("materials", &materials);
  treeMaterialsOut->Fill();

  KalmanTracking kalmanTracking(*materials, 0., false);
  ClusterReco    clusterReco   (*materials); // pass

  ElossEstimator elossEstimator;

  if(collectHits != asPion)
  {
    ifstream fileGain;

    if(collectHits == ifIdentified)
      fileGain.open("../out/gains_asPion.dat");
    else // notNeeded
      fileGain.open("../out/gains_ifIdentified.dat");

    cerr << "  loading gains...";
    elossEstimator.loadGains(fileGain);
    cerr << " [done]" << endl;

    fileGain.close();
  }

  MostProbable   mostProbable;

  // Vertex finder
  VertexFinder vertexFinder;

  ofstream fileEstimate("../out/estimates.dat");

  map<ChipId, vector<TSlimMeasurement> > hits;

  for(int ievent = 0; ievent < treeIn->GetEntries(); ievent++) 
  {
    branchIn->GetEntry(ievent);

    cerr << "\033[22;31m"
         << " Bunch crossing " << bunchCrossing->bxNumber
         << "\033[22;0m" << endl;

    // Reconstruct clusters
    cerr << "  - reco clusters + energy loss ";

    int i = 0;
    for(vector<TVertex>::iterator
          vertex = bunchCrossing->recVertices.begin();
          vertex!= bunchCrossing->recVertices.end(); vertex++)
      for(vector<TTrack>::iterator track = vertex->tracks.begin();
                                   track!= vertex->tracks.end(); track++)
      {
        printPropeller(++i);

        // If gain correction is available, correct measured deposit values
        if(collectHits != asPion)
          elossEstimator.correctDeposits(&(*track)); 
      
        // Calculate cluster properties
        clusterReco.run(*track);

        // Estimate the most probable energy loss rate
        pair<double,double> estimate = 
           elossEstimator.estimate(*track, fileEstimate);

        double p = track->pt * cosh(track->eta);

        // Collect deposits for gain calibration
        if(collectHits != notNeeded)
        {
          int pid;

          if(collectHits == asPion)
            pid = pion;
          else // ifIdentified
          {
//            pid = mostProbable.surePid(p, estimate.first);
            pid = mostProbable.guessPid(p, estimate.first);

            // Should not be too off
            if(pid != unknown)
              if(fabs(estimate.first / mostProbable.getValue(p, pid) - 1) > 0.1)
                pid = unknown;
          }

          if(pid != unknown)
          for(vector<Hit>::iterator hit = track->hits.begin();
                                    hit!= track->hits.end(); hit++)
          {
            TSlimMeasurement slim;
  
            slim.energy      = hit->charge; // MeV
            slim.path        = hit->length; // cm
            slim.epsilon     = mostProbable.getValue(p, pid);
  
            slim.hasOverflow = hit->hasOverflow; // bool
    
            TChipId chipId = hit->chipId;
  
            if(hits.count(chipId.code) == 0)
               hits[chipId.code].reserve(100);
   
            hits[chipId.code].push_back(slim);
          }
        }
      }

    cerr << " [done]" << endl;

    // Collect tracks for vertex finding
    vector<pair<double, double> > trackZ;

    for(vector<TVertex>::const_iterator 
          vertex = bunchCrossing->recVertices.begin();
          vertex!= bunchCrossing->recVertices.end(); vertex++)
      for(vector<TTrack>::const_iterator track = vertex->tracks.begin();
                                         track!= vertex->tracks.end(); track++)
        trackZ.push_back(pair<double,double>(track->z, track->sigma_z));

    // bunchCrossing.simVertices -> sim // FIXME
    VertexCollection sim;
    int it = 0;
    for(vector<TVertex>::const_iterator
          vertex = bunchCrossing->simVertices.begin();
          vertex!= bunchCrossing->simVertices.end(); vertex++)
    {
      Vertex v;

      v.first = vertex->z; 

      for(vector<TTrack>::const_iterator track = vertex->tracks.begin();
                                         track!= vertex->tracks.end(); track++)
        v.second.push_back(it++); // FIXME

//      Vertex v(vertex->z, t);

      sim.push_back(v);
    }

    // Vertex finding
    VertexCollection recVertices;

    if(trackZ.size() > 0) // FIXME
      recVertices = vertexFinder.findVertices(trackZ, sim);

    // Update vertex positions FIXME
    cerr << " updating vertex positions " << endl;

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

      cerr << "  setting vertex z " << vertex->z
                                    << " -> " << bestVertex->first << endl;
      vertex->z = bestVertex->first;
    }

    // Fill tree
    treeOut->Fill();

    bunchCrossing->Clear();
  } 

  fileEstimate.close();

  fileIn->Close();

  fileOut->Write();
  fileOut->Close();

  if(collectHits != notNeeded)
  {
    cerr << " writing hits for gain calibration..";

    ofstream outFile("../out/hits.bin", ios::binary);

    cerr << " " << hits.size() << " chips";

    // hits
    {
      int size = hits.size();
      outFile.write((char *)&size, sizeof(size));

      for(map<ChipId, vector<TSlimMeasurement> >::const_iterator
          hit = hits.begin(); hit != hits.end(); ++hit)
      {
        outFile.write((char *)&hit->first , sizeof(hit->first ));

        int n = hit->second.size();
        outFile.write((char *)&n, sizeof(n));

        for(vector<TSlimMeasurement>::const_iterator
          slim = hit->second.begin(); slim!= hit->second.end(); slim++)
        outFile.write((char *)&(*slim), sizeof(*slim));
      }
    }

    outFile.close();
  }

  cerr << " [done]" << endl;
}

/*****************************************************************************/
int main(int arg, char **arc)
{
  gErrorIgnoreLevel = kSysError;

  options(arg, arc);

  reconstructEvents();
}

