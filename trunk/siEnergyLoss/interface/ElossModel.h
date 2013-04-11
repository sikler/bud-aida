#ifndef _ElossModel_h_
#define _ElossModel_h_

#include "../../DataFormats/interface/Hit.h"

#include "TMatrixD.h"
#include "TVectorD.h"

#include <vector>

class TLayer;

class ElossModel
{
  public:
    ElossModel(Hit hit_, TLayer * unit_);
    ~ElossModel();

    Hit getHit(const std::vector<double> & pars);
    double getSigma(double y);
    double getDerivatives(const std::vector<double> & pars,
                          TVectorD & beta,
                          TMatrixD & alpha);
    double getValue(const std::vector<double> & pars);

  private:
    void preparePixels(const std::vector<double> & pars);
    double getChi2(double y, double Delta, double thr, double ovf,
                   std::vector<double> & der);

    void calculate(const std::vector<double> & pars, double & val,
                   TVectorD & beta,
                   TMatrixD & alpha);

    Hit hit;
    TLayer * unit;
};

#endif
