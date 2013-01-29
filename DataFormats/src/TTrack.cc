#include "../interface/TTrack.h"

ClassImp(TTrack)

TTrack::TTrack()
{
  pixelHits.clear();
  stripHits.clear();

  hits.clear();
  lhits.clear();
}

TTrack::~TTrack()
{
}
