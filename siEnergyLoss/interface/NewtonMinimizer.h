#ifndef _NewtonMinimizer_h_
#define _NewtonMinimizer_h_

#include "TMatrixD.h"

#include <vector>

class ElossModel;

class NewtonMinimizer
{
 public:
   NewtonMinimizer();
   ~NewtonMinimizer();

   void minimize(ElossModel & theModel,
                 const std::vector<bool> & isFix,
                 std::vector<double> & par, TMatrixD & err,
                 int & nstep);

 private:
};

#endif
