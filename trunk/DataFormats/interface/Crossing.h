#ifndef _Crossing_h_
#define _Crossing_h_

#include <vector>

#include "TObject.h"

/*****************************************************************************/
class Crossing : public TObject
{
 public:
  Crossing();
  virtual ~Crossing();

  double pos[2];
  int dir, typ;

  ClassDef(Crossing,1)
};

#endif
