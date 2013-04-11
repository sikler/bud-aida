#include "../interface/TChipId.h"

ClassImp(TChipId)

TChipId::TChipId()
{
}

TChipId::TChipId(int ilayer, int iz, int iphi)
{
  code.first         = ilayer;
  code.second.first  = iz;
  code.second.second = iphi;
}

TChipId::~TChipId()
{
}
