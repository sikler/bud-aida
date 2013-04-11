#ifndef _TVertex_h_
#define _TVertex_h_

#include "TObject.h"
#include "TTrack.h"

#include <vector>

class TVertex : public TObject
{
 public:
  TVertex();
  virtual ~TVertex(); 

  std::vector<TTrack> tracks;

  float z, sigma_z;

  short int processId;

  ClassDef(TVertex,1)
};

#endif
