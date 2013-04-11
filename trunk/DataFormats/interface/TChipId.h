#ifndef _TChipId_h_
#define _TChipId_h_

#include <iostream>
#include <utility>

#include "TObject.h"

typedef std::pair<int, std::pair<int,int> > ChipId;

class TChipId : public TObject
{
 public:
  TChipId();
  TChipId(int ilayer, int iz, int iphi);
  virtual ~TChipId();

  int getiLayer() { return code.first;         }
  int getiZ()     { return code.second.first;  }
  int getiPhi()   { return code.second.second; }

  ChipId code;

  void print()
  {
    std::cerr << " TChipId : " << getiLayer() << " " << getiZ()
                        << " " << getiPhi() << std::endl;
  }

  ClassDef(TChipId,1)
};

#endif

