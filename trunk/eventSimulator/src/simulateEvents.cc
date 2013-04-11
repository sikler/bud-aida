#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "../interface/gzstream.h"
#include "../interface/KalmanTracking.h"

using namespace std;

#include "../../DataFormats/interface/TBunchCrossing.h"
#include "../../DataFormats/interface/TVertex.h"

#include "TFile.h"
#include "TTree.h"

// FIXME
char inputName[256], outputName[256];

int nBunchCrossings, nCollisions;
double Mu, sigmaZ, ptMin;
float gainHalfWidth;
bool fix;

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
TVertex readGeneratedEvent(igzstream & fileEvents)
{
  int d, processId, n;

  fileEvents >> d;          // event number
  fileEvents >> processId;  // process id, Pythia
  fileEvents >> n;          // number of particles to follow
  fileEvents >> d;          // total number of particles

  // Initialize vertex
  TVertex vertex;

  vertex.z = sigmaZ * getGaussRandom();
  vertex.processId = processId;

  // Add tracks
  for(int i = 0; i < n; i++)
  {
    int pid;
    double eta,pt,phi;

    fileEvents >> pid; // particle id
    fileEvents >> eta; // eta
    fileEvents >> pt;  // transverse momentum
    fileEvents >> phi; // 

    if(pt > ptMin)
    {
      TTrack track;

      track.pdgId  =  pid;
      track.charge = (pid > 0 ? 1 : -1);

      track.eta = eta;
      track.pt  = pt ;
      track.phi = phi;

      track.z   = vertex.z; // FIMXE

      vertex.tracks.push_back(track);
    }
  }

  return vertex;
}

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

  // Set default values
  nBunchCrossings = 1; nCollisions = 1; Mu = 1.; sigmaZ = 1.;
  gainHalfWidth = 0.;
  ptMin = 0.050; // GeV/c
  fix = false;

  do
  {
    if(strcmp(arc[i],"-numberOfBunchCrossings") == 0 ||
       strcmp(arc[i],"-nbx") == 0) 
         nBunchCrossings = atoi(arc[++i]);

    if(strcmp(arc[i],"-numberOfCollisionsFix ") == 0 ||
       strcmp(arc[i],"-ncfix ") == 0)
    {
      fix = true;
      Mu = atoi(arc[++i]);
    }

    if(strcmp(arc[i],"-numberOfCollisionsPoissonMean") == 0 ||
       strcmp(arc[i],"-ncpoi ") == 0)
    {
      fix = false;
      Mu = atof(arc[++i]);
    }

    if(strcmp(arc[i],"-interactionRegionSigmaZ") == 0 ||
       strcmp(arc[i],"-irsig") == 0)
         sigmaZ          = atof(arc[++i]);

    if(strcmp(arc[i],"-randomGainHalfWidth") == 0 ||
       strcmp(arc[i],"-rghw") == 0)
         gainHalfWidth   = atof(arc[++i]);

    if(strcmp(arc[i],"-i") == 0) sprintf(inputName ,"%s",arc[++i]);
    if(strcmp(arc[i],"-o") == 0) sprintf(outputName,"%s",arc[++i]);

    i++;
  }
  while(i < arg);
}

/*****************************************************************************/
void simulateEvents()
{
  // Output file
  TFile * fileOut = new TFile(outputName, "recreate");

  // Tracks
  TTree * tree    = new TTree("trackTree","trackTree");

  TBunchCrossing * bunchCrossing = new TBunchCrossing();
  tree->Branch("bunchCrossing", "TBunchCrossing", &bunchCrossing, 16000, 2);

  // Materials
  TTree * treeMaterials = new TTree("materialTree","materialTree");

  vector<TLayer> materials;
  treeMaterials->Branch("materials", &materials);

  KalmanTracking kalmanTracking(materials, gainHalfWidth, true);

  treeMaterials->Fill();

  int nEvent = 0;

  // Input file
  igzstream fileEvents(inputName);

  for(int ibunx = 0; ibunx < nBunchCrossings; ibunx++)
  { 
    bunchCrossing->runNumber = 1;
    bunchCrossing->bxNumber  = ++nEvent;

    cerr << "\033[22;31m"
         << " Bunch crossing " << bunchCrossing->bxNumber
         << "\033[22;0m" << endl;

    if(fix) nCollisions = int(Mu+0.5);
       else nCollisions = getPoissonRandom(Mu);

    for(int j = 0; j < nCollisions; j++)
    {
      // Read generated event
      TVertex simVertex = readGeneratedEvent(fileEvents);
     
      if(simVertex.tracks.size() > 0)
      { 
        TVertex recVertex;

        cerr << "   - simu clusters + tracking    ";
  
        for(vector<TTrack>::iterator track = simVertex.tracks.begin();
                                     track!= simVertex.tracks.end(); )
        {
          TTrack recTrack;
  
          // simulate + reconstruct, true if at least 4 hits
          if(kalmanTracking.process(*track, recTrack))
          {
            printPropeller(int(track - simVertex.tracks.begin() + 1));
            recTrack.pdgId = track->pdgId; // FIXME
            recVertex.tracks.push_back(recTrack);

            track++;
          }
          else
          {
            simVertex.tracks.erase(track);
          }
        }
        cerr << " [done]" << endl;
  
        // Do vertex finding! FIXME
        recVertex.z = simVertex.z;
  
        // Store 
        if(simVertex.tracks.size() > 0 &&
           recVertex.tracks.size() > 0)
        {
          bunchCrossing->simVertices.push_back(simVertex);
          bunchCrossing->recVertices.push_back(recVertex);
        }
      }
    }

    // Fill tree
    tree->Fill();
    bunchCrossing->Clear();
  }

  fileEvents.close(); 

  fileOut->Write();
  fileOut->Close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
int main(int arg, char **arc)
{
  options(arg, arc);

  simulateEvents();
}

