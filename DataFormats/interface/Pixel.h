#ifndef _Pixel_h_
#define _Pixel_h_

#include <vector>

#include "TObject.h"
#include "Crossing.h"

/*****************************************************************************/
class Pixel : public TObject
{
 public:
  Pixel();
  virtual ~Pixel();

  class Meas // measured
  {
   public:
    int x[2]; // coordinate of lower left corner [pitch]
    double y; // deposited charge [MeV]

    Crossing c[2];
  };

  class Calc // calculated
  {
   public:
    bool isTouched;

    double l;        // true path [cm]
    double dl_dP[2]; //
  };

  Meas meas;
  Calc calc;

  ClassDef(Pixel,1)
};

#endif
