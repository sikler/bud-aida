#ifndef _TPixelHit_h_
#define _TPixelHit_h_

#include <stdint.h>
#include <utility>
#include "TObject.h"

class TPixelHit : public TObject
{
 public:
  TPixelHit();
  virtual ~TPixelHit();

  std::pair<short int, short int> meas;
  std::pair<float,float>          pred;

  bool forCalib; // isStraight && isComplete && hasBigPixelsOnlyInside
  bool forEloss; // isComplete

  float x;          // pathlength [cm]
  float thickness;  //            [cm] 

  float Delta; // deposited energy [MeV]
  bool isSaturated; // any pixel with adc > 200?

  uint32_t detId;
  unsigned char chip;

  float r; // global r coordinate of the hit
  float p; // global momentum, updated state, at the hit

  unsigned char nChannels;

  ClassDef(TPixelHit,1)
};

#endif
