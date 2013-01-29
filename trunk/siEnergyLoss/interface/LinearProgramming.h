#ifndef _LinearProgramming_h_
#define _LinearProgramming_h_

#include "TVectorD.h"
#include <vector>

class LinearProgramming
{
 public:
   LinearProgramming();
   ~LinearProgramming();

   TVectorD solve(double a, std::vector<double> meas,
                  double c0, int & status);

 private:
};

#endif
