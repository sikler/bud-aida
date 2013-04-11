#ifndef _TTrack_h_
#define _TTrack_h_

#include <utility>

#include "TObject.h"

#include "Point.h"
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

  float z, sigma_z;
  float d0;

  std::pair<float,float> epsilon;

  std::vector<Point> points;
  std::vector<Hit>   hits;

  short int pdgId;

  ClassDef(TTrack,1)
};

#endif
