#ifndef _TTrack_h_
#define _TTrack_h_

#include <utility>
#include "TObject.h"

#include "TPixelHit.h"
#include "TStripHit.h"

#include "TStripHit.h"

#include "Coord.h"

#include "Hit.h" 

class TTrack : public TObject
{
 public:
  TTrack();
  virtual ~TTrack();

  short int charge;
  float eta;
  float pt;
  float phi;

  float chi2;
  int    ndf;

  float z;
  float d0;

  std::pair<float,float> epsilon;

  std::vector<TPixelHit> pixelHits;
  std::vector<TStripHit> stripHits;

  std::vector<Coord> hits;

  std::vector<Hit> lhits;

  short int pdgId;

  ClassDef(TTrack,1)
};

#endif
