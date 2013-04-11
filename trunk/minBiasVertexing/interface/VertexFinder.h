#ifndef _VertexFinder_h_
#define _VertexFinder_h_

#include <fstream>
#include <utility>
#include <vector>

#include "TVectorD.h"

typedef std::pair<float, std::vector<int> > Vertex;
typedef std::vector<Vertex> VertexCollection;

class VertexFinder
{
 public:
  VertexFinder();
  ~VertexFinder();

  VertexCollection findVertices
    (const std::vector<std::pair<double, double> > & points,
     const VertexCollection & sim);

 private:
  void evaluatePerformance(const VertexCollection & sim,
                           const VertexCollection & rec,
                           std::vector<int> & asim,
                           std::vector<int> & arec,
                           std::vector<std::vector<int> > & msim, int method);

  double getLambda(const VertexCollection & vertices);
 
  void numberOfVertices(const std::vector<std::pair<double,double> > & points,
                        const std::vector<std::pair<TVectorD,TVectorD> > & clusters,
                        VertexCollection & rec,
                        int classification);

  std::ifstream filePythia;
  std::ofstream fileMacro;

  std::vector<std::vector<std::vector<int> > > asim, arec;

  std::vector<std::vector<int> > msim0, msim1, msim2, msim3;
};

#endif
