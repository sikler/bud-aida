#ifndef _GaussianMixture_h_
#define _GaussianMixture_h_

#include <utility>
#include <vector>

#include "TVectorD.h"
#include "TMatrixD.h"

typedef std::pair<float, std::vector<int> > Vertex;
typedef std::vector<Vertex> VertexCollection;

class GaussianMixture
{
 public:
  GaussianMixture() {}

  double run(int K, const std::vector<std::pair<double, double> > & points,
             TVectorD mu, TVectorD P,
             VertexCollection & vertices);

 private:
  double logProb(const std::pair<double,double> & p, double mu, double P);
  void estimateAverages(const std::vector<std::pair<double,double> > & points,
                        const TMatrixD & p,
                        TVectorD & mu, TVectorD & P);
  double poisson(double mu, double n);
  double estimateChiSquare(const std::vector<std::pair<double,double> > & points,
                           const TVectorD & mu, const TVectorD & P,
                           TMatrixD & p, bool estimateResponsibility);
  double getChiSquare(const std::vector<std::pair<double,double> > & points,
                      const TVectorD & mu, const TVectorD & P);
  void estimateResponsibility(const std::vector<std::pair<double,double> > & points,
                              const TVectorD & mu, const TVectorD & P,
                              TMatrixD & p);
};

#endif
