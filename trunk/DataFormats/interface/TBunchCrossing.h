#ifndef _TBunchCrossing_h_
#define _TBunchCrossing_h_

#include "TObject.h"
#include <vector>
#include <string>

class TVertex;

class TBunchCrossing : public TObject
{
 public:
  TBunchCrossing();
  virtual ~TBunchCrossing();

  int runNumber;
  int bxNumber;

  // vertices
  std::vector<TVertex> simVertices;
  std::vector<TVertex> recVertices;

  void Clear();

  ClassDef(TBunchCrossing,1)
};

#endif
