#ifndef _NewtonMinimizer_h_
#define _NewtonMinimizer_h_

#include "CLHEP/Matrix/Matrix.h"
#include <vector>

class ElossModel;

class NewtonMinimizer
{
 public:
   NewtonMinimizer();
   ~NewtonMinimizer();

   void minimize(ElossModel & theModel,
                 const std::vector<bool> & isFix,
                 std::vector<double> & par, CLHEP::HepMatrix & err,
                 int & nstep);

 private:
};

#endif
