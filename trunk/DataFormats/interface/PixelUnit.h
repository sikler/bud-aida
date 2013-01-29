#ifndef _PixelUnit_h_
#define _PixelUnit_h_

#include <iostream>
#include <vector>

#include "TObject.h"

/*****************************************************************************/
class PixelUnit : public TObject
{
 public:
  PixelUnit();
  virtual ~PixelUnit();

  int nrows;
  int ncolumns;

  double diffSigma; // cm, maximal Gaussian sigma for diffusion
  double pitch[3];  // cm

  void print()
  {
    std::cerr << " unit : nrows=" << nrows
                 << " ncolumnsm=" << ncolumns
                 << " diffSigma=" << diffSigma
                  << " pitch[0]=" << pitch[0]
                  << " pitch[1]=" << pitch[1]
                  << " pitch[2]=" << pitch[2] << std::endl;
  }

  ClassDef(PixelUnit,1)
};

#endif
