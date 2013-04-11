#ifndef _TLayer_h_
#define _TLayer_h_

#include <utility>
#include <iostream>
#include <map>

#include "TObject.h"

typedef std::pair<int, std::pair<int,int> > ChipId;

class TLayer : public TObject
{
 public:
  TLayer();
  virtual ~TLayer();

  int ilayer;

  double radius;     // cm
  double halfLength; // cm

  double thickness;  // cm
  double pitch[3];   // cm
  double diffSigma;  // cm

  int nrows;
  int ncolumns;

  double sigma_rphi;
  double sigma_z;

  double noise;      // MeV
  double threshold;  // MeV
  double overflow;   // MeV

  double coupling;

  double dz;
  double dphi;

  std::map<ChipId,float> gain;

  double B; // T

  bool isPixel;

  void print() const
  {
    std::cerr << " i="   << ilayer
              << " r="   << radius
              << " l/2=" << halfLength
              << " sigma_rphi=" << sigma_rphi
              << " sigma_z="    << sigma_z
              << " thickness="  << thickness
              << std::endl;
  }

  ClassDef(TLayer,1)
};

#endif

