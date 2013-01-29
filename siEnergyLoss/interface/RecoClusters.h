#ifndef _RecoClusters_h_
#define _RecoClusters_h_

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

#include <fstream>
#include <vector>
#include <utility>

class TTrack;
class Hit;
class ElossModel;

class RecoClusters
{
 public:
  RecoClusters();
  virtual ~RecoClusters();

  void run(TTrack & recTrack);

 private:
  void fixHit(Hit & hit);

  void estimatePositionLF      (Hit & hit, double & epsilon);
  void estimatePositionWeighted(Hit & hit, double & epsilon);

  void estimateDeposit         (Hit & hit, double epsilon, double * charge);

  std::pair<double,double> processHit(Hit & hit);

  std::ofstream fileFit;
};

#endif
