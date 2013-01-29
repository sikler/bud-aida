#ifndef _ElossModel_h_
#define _ElossModel_h_

#include "../../DataFormats/interface/Hit.h"

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

#include <vector>

class ElossModel
{
  public:
    ElossModel(Hit hit_);
    ~ElossModel();

    Hit getHit(const std::vector<double>& pars);
    double getDerivatives(const std::vector<double>& pars,
                         CLHEP::HepVector& beta,
                         CLHEP::HepMatrix& alpha);
    double getValue(const std::vector<double>& pars);

  private:
    void preparePixels(const std::vector<double>& pars);
    double getSigma(double y);
    double getChi2(double y, double Delta, double thr, double ovf,
                   std::vector<double> & der);

    void calculate(const std::vector<double>& pars, double & val,
                   CLHEP::HepVector & beta,
                   CLHEP::HepMatrix & alpha);

    Hit hit;
};

#endif
