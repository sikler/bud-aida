#ifndef _ElossEstimator_h_
#define _ElossEstimator_h_

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

#include <vector>
#include <utility>

#include <fstream>

class TTrack;

class ElossEstimator
{
 public:
  ElossEstimator();
  virtual ~ElossEstimator();

  void calibrateGains(int pass);
  void estimate(TTrack & recTrack, std::ofstream & fileEstimate);

 private:
  inline double sigmaD(double y);
  double logLikelihood(double epsilon, double y, double l);

  double getChi2(double gain,    const std::vector<double> & hits, double l);
  double getChi2(double epsilon, const std::vector<std::pair<double,double> > & hits);

  std::pair<double,double> getEstimate
    (const std::vector<std::pair<double,double> > & track);

  double getTruncMean(      std::vector<double> & dedx);
  double getPowerMean(const std::vector<double> & dedx, double power);
};

#endif
