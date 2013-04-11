#ifndef _Point_h_
#define _Point_h_

#include <utility>
#include "TObject.h"

#define nDims 3

class Point : public TObject
{
 public:
  Point();
  virtual ~Point();

  int q;
  double x[nDims];
  double pt; // auxiliary, sometimes needed
  double pz;
  double eta;
  double p_; // total momentum
  double p[nDims];

  bool isPixel; // needed for plots

  void print();

  ClassDef(Point,1)  
};

#endif
