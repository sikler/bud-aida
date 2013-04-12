#ifndef _ElossEstimator_h_
#define _ElossEstimator_h_

#include <fstream>
#include <vector>
#include <utility>
#include <map>

class TTrack;
class TLayer;

#include "../../DataFormats/interface/TSlimMeasurement.h"
#include "../../DataFormats/interface/TChipId.h" // FIXME

class ElossEstimator
{
 public:
  ElossEstimator();
  ElossEstimator(const std::vector<TLayer> & materials_);
  virtual ~ElossEstimator();

  void loadGains
    (std::ifstream & fileGain);

  void correctDeposits(TTrack * track);

  void calibrateGains
    (std::map<ChipId, std::vector<TSlimMeasurement> > & hits,
     std::ofstream & fileGain);

  std::pair<double,double>
    estimate(TTrack & recTrack, std::ofstream & fileEstimate);

 private:
  inline double sigmaD(double y);
  double logLikelihood(double epsilon, double y, double l, bool over);

  double getChi2Gain
    (double gain,    const std::vector<TSlimMeasurement> & hits);

  double getChi2Epsilon
    (double epsilon, const std::vector<TSlimMeasurement> & hits);

  void getValuesByVaryingGain
    (const std::vector<TSlimMeasurement> & hits,
     double gain, std::vector<double> & val);

  void getValuesByVaryingEpsilon
    (const std::vector<TSlimMeasurement> & hits,
     double epsilon, std::vector<double> & val);

  void   functionGain(double g, std::vector<double> & val);
  double functionGain(double gain);

  void shft2(double &a, double &b, const double c);
  void shft3(double &a, double &b, double &c, const double d);

  double goldenSearch(double ax, double bx, double cx);
  double newtonMethodGain(std::pair<double,double> & value);

  std::pair<double,double> getEstimate
    (const std::vector<TSlimMeasurement> & track);

  double getTruncMean(      std::vector<double> & dedx);
  double getPowerMean(const std::vector<double> & dedx, double power);

  std::vector<TSlimMeasurement> slim;

  std::map<ChipId, float> gains;
 
  std::vector<TLayer> materials;
};

#endif
