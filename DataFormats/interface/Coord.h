#ifndef _Coord_h_
#define _Coord_h_

#include <utility>
#include "TObject.h"

#define nDims 3

class Coord : public TObject
{
 public:
  Coord();
  virtual ~Coord();

  int q;
  double x[nDims];
  double pt; // auxiliary, sometimes needed
  double pz;
  double eta;
  double p_; // total momentum
  double p[nDims];

  void print();

  ClassDef(Coord,1)  
};

#endif
