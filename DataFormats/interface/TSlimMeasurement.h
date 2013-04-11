#ifndef _TSlimMeasurement_h_
#define _TSlimMeasurement_h_

#include "TObject.h"

class TSlimMeasurement : public TObject
{
 public:
  TSlimMeasurement();
  virtual ~TSlimMeasurement();

  float energy;  // MeV, measured
  float path;    // cm,  measured
  float epsilon; // MeV/cm, predicted

  bool hasOverflow;

  ClassDef(TSlimMeasurement,1)
};

#endif

