#ifndef _Hit_h_
#define _Hit_h_

#include <iostream>
#include <vector>

#include "TObject.h"

#include "Pixel.h"
#include "TChipId.h"

class TChipId;

class Hit : public TObject
{
 public:
  Hit();
  virtual ~Hit();

  double length;                         // 3D [cm]
  double lambda, ulambda[2];             // 2D [cm]

  double dpos[2], pos_orig[2], pos[2];   // 2D [pitch]
  double shift[2];                       // shift due to drift [pitch]
  double charge_orig, charge;            // [MeV]

  double error[3][3];

  int ilayer;
  TChipId chipId;
  bool hasOverflow;

  std::vector<Pixel> filledPixels; // only those which have adc value
  std::vector<Pixel> allPixels;    // all: add missing touched

  void print()
  {
    std::cerr << " hit : l=" << length
                    << " lam="    << lambda
                    << " ulam[0]=" << ulambda[0]
                    << " ulam[1]=" << ulambda[1] << std::endl;
    std::cerr << " hit : dpos[0]=" << dpos[0]
                    << " dpos[1]=" << dpos[1] << std::endl;
    std::cerr << " hit : pos[0]=" << pos[0]
                    << " pos[1]=" << pos[1] << std::endl;
    std::cerr << " hit : ncha=" << allPixels.size()
                  << " charge=" << charge
                  << std::endl;
    for(std::vector<Pixel>::const_iterator pixel = allPixels.begin();
                                           pixel!= allPixels.end(); pixel++)
      std::cerr << "  cha : x=" << pixel->meas.x[0]
                       << ","   << pixel->meas.x[1]
                       << " y=" << pixel->meas.y
                       << " isT=" << (pixel->calc.isTouched ? "1" : "0")
                       << " l=" << pixel->calc.l
                       << std::endl;
  }

  ClassDef(Hit,1)
};

#endif
