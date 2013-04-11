#ifndef _FitStripCluster_h_
#define _FitStripCluster_h_

#include <vector>
#include <utility>

#include "TMatrixD.h"
#include "TVectorD.h"

class FitStripCluster
{
 public:
  FitStripCluster(const std::vector<std::pair<double,int> > & b_,
                  const double alpha, const double sigma);
  virtual ~FitStripCluster();

  int run(double & chi2, std::pair<double,double> & result,
          std::vector<double> & x);

 private:
  double getChi2   (const std::vector<double> & x);

  void getAlphaBeta(const std::vector<double> & x,
                    TMatrixD & alpha,
                    TVectorD & beta,
                    const std::vector<bool> & isFix);

  const std::vector<std::pair<double,int> > & b;

  double Coupling;
  double Sigma, iSigma;

  double a00,a11,a01, a0,a1;

  int npar;
};

#endif
