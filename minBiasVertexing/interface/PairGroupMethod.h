#ifndef _PairGroupMethod_
#define _PairGroupMethod_

#include <utility>
#include <vector>

#include "TVectorD.h"

class PairGroupMethod
{
 public:
  PairGroupMethod() {}

  void run(const std::vector<std::pair<double,double> > & points,
                 std::vector<std::pair<TVectorD, TVectorD> > & clusters,
                 unsigned int maxClusters); 
};

#endif
