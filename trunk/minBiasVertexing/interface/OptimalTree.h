#ifndef _OptimalTree_h_
#define _OptimalTree_h_

#include <utility>
#include <vector>

#include "TVectorD.h"
#include "TMatrixD.h"

typedef std::pair<float, std::vector<int> > Vertex;
typedef std::vector<Vertex> VertexCollection;

class OptimalTree
{
 public:
  OptimalTree() {}

  double run(int K, const std::vector<std::pair<double, double> > & points,
             TVectorD mu, TVectorD P,
             std::vector<std::vector<int> > lists,
             VertexCollection & vertices);

 private:
};

#endif
