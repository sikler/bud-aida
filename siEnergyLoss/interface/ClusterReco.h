#ifndef _ClusterReco_h_
#define _ClusterReco_h_

#include <fstream>
#include <vector>

//#include "../../DataFormats/interface/TLayer.h"

class TTrack;
class Hit;
class TLayer;

class ClusterReco
{
 public:
  ClusterReco(const std::vector<TLayer> & materials_);
  virtual ~ClusterReco();

  void run(TTrack & recTrack);

 private:
  void fixHit(Hit & hit);

  void estimatePositionLF      (Hit & hit, double & epsilon);
  void estimatePositionWeighted(Hit & hit, double & epsilon);

  void estimateDeposit         (Hit & hit, TLayer * unit,
                                double epsilon, double * charge);

  std::pair<double,double> processHit(Hit & hit);

  std::ofstream filePixel;
  std::ofstream fileStrip;
  std::vector<std::ofstream *> fileFit;

  double cForOverflow;

  std::vector<TLayer> materials;
};

#endif
