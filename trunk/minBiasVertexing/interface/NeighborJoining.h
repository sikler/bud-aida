#ifndef _NeighborJoining_
#define _NeighborJoining_

#include <utility>
#include <vector>

#include "TVectorD.h"

class NeighborJoining
{
 public:
  NeighborJoining() {}

  void run(const std::vector<std::pair<double,double> > & points,
                 std::vector<std::pair<TVectorD, TVectorD> > & clusters,
                 unsigned int maxClusters); 
};

#endif
