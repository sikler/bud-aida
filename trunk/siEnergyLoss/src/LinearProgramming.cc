#include "../interface/LinearProgramming.h"
#include "../interface/Levels.h"

#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TQpProbDens.h"
#include "TMehrotraSolver.h"

#include <iostream>

//  minimize    c^T x + ( 1/2 ) x^T Q x       
//  subject to                A x  = b  
//                    clo <=  C x <= cup
//                    xlo <=    x <= xup

using namespace std;

/*****************************************************************************/
LinearProgramming::LinearProgramming()
{
}

/*****************************************************************************/
LinearProgramming::~LinearProgramming()
{
}

/*****************************************************************************/
// a  = coupling
// c0 =  1.; minimize
// c0 = -1.; maximize
TVectorD LinearProgramming::solve
  (double a, vector<double> meas, double c0, int & status)
{
  const int nrVar  = meas.size() + 2; // number of all strips

  int nrEqual   = 0;
  int nrInEqual = 0;

  for(vector<double>::const_iterator strip = meas.begin();
                                     strip!= meas.end(); strip++)
  {
    if(*strip < Overflow) nrEqual++;
                     else nrInEqual++;
  }

  nrInEqual += 2; // first and last

  TVectorD    c(nrVar); c = c0;   // minimize or maximize
  TMatrixDSym Q(nrVar); Q.Zero(); // Q = 0

  TMatrixD A(nrEqual  , nrVar); A.Zero();
  TVectorD b(nrEqual);
  TMatrixD C(nrInEqual, nrVar); C.Zero();

  // All inequalities are for overflow, initially
  TVectorD clo (nrInEqual); TVectorD cup (nrInEqual);
  TVectorD iclo(nrInEqual); TVectorD icup(nrInEqual);

  // Initilalize again
  nrEqual = 0; nrInEqual = 0;

  // First strip (upper limit = Threshold)
  iclo(nrInEqual) = 0; clo(nrInEqual) = 0;
  icup(nrInEqual) = 1; cup(nrInEqual) = Threshold;

  C(nrInEqual, 0) = 1 - 2*a;
  C(nrInEqual, 1) = a;

  nrInEqual++;

  int i = 1;
  for(vector<double>::const_iterator strip = meas.begin();
                                     strip!= meas.end(); strip++)
  {
    if(*strip < Overflow)
    { // equality, =
      if(i+1 < nrVar) A(nrEqual,i+1) = a;
                      A(nrEqual,i  ) = 1 - 2*a;
      if(i-1 >=    0) A(nrEqual,i-1) = a;

      b(nrEqual) = *strip;

      nrEqual++;
    }
    else
    { // inequality, >=
      if(i+1 < nrVar) C(nrInEqual,i+1) = a;
                      C(nrInEqual,i  ) = 1 - 2*a;
      if(i-1 >=    0) C(nrInEqual,i-1) = a;

      iclo(nrInEqual) = 1; clo(nrInEqual) = Overflow;
      icup(nrInEqual) = 0; cup(nrInEqual) = 0;

      nrInEqual++;
    }

    i++;
  }

  // Last strip (upper limit = Threshold)
  C(nrInEqual, nrVar-1) = 1 - 2*a;
  C(nrInEqual, nrVar-2) = a;

  iclo(nrInEqual) = 0; clo(nrInEqual) = 0;
  icup(nrInEqual) = 1; cup(nrInEqual) = Threshold;

  nrInEqual++;

  // x_i >= 0
  TVectorD xlo (nrVar);  xlo = 0;
  TVectorD xup (nrVar);  xup = 0;
  TVectorD ixlo(nrVar); ixlo = 1;
  TVectorD ixup(nrVar); ixup = 0;

  // set up the problem
  TQpProbDens *qp = new TQpProbDens(nrVar,nrEqual,nrInEqual);
  TQpDataDens *prob = (TQpDataDens *)
                 qp->MakeData(c,Q,xlo,ixlo,xup,ixup,A,b,C,clo,iclo,cup,icup);

  // set up variables
  TQpVar      *vars  = qp->MakeVariables(prob);
  TQpResidual *resid = qp->MakeResiduals(prob);

  // solve
  TMehrotraSolver *s = new TMehrotraSolver(qp,prob);
  status = s->Solve(prob,vars,resid);

  const TVectorD amplitudes = vars->fX;

  delete qp;
  delete s;

  return amplitudes;
}
