#include "../interface/NewtonMinimizer.h"

#include "../interface/ElossModel.h"

//#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

#include <cmath>
#include <iostream>
#include <fstream>
//#include <iomanip>

using namespace std;
using namespace CLHEP;

/*****************************************************************************/
NewtonMinimizer::NewtonMinimizer()
{
}

/*****************************************************************************/
NewtonMinimizer::~NewtonMinimizer()
{
}

/*****************************************************************************/
void NewtonMinimizer::minimize(ElossModel & theModel, 
                               const vector<bool> & isFix,
                               vector<double> & par, HepMatrix & err,
                               int & nstep)
{
  nstep = 0;

  double old = theModel.getValue(par);

  HepVector g(3), h(3);

  double move;
  do
  {
    HepMatrix alpha(3,3, 0.);
    HepVector beta(3, 0.);

    // Get derivatives
    theModel.getDerivatives(par, beta,alpha);

    double sum[2] = {0.,0.};
    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      if(alpha[i][j] != 0.)
      { sum[1] += fabs(alpha[i][j]); sum[0] += 1.; }

    if(sum[0] != 0.) sum[1] /= sum[0];
                else sum[1] = 1.;

    for(int k = 0; k < 3; k++)
    if(isFix[k])
    {
      beta[k] = 0.;

      for(int j = 0; j < 3; j++)
      if(j != k)
      { alpha[j][k] = 0.; alpha[k][j] = 0.; }

      alpha[k][k] = sum[1];
    }

    // Conjugate gradient
    if(nstep++ == 0)
    { // initialize
      g = beta;
      h = g;
    }
    else
    { // update
      HepVector go = g;
      g = beta;

      double gamma = dot(g - go,g) / dot(go,go);
      h = g + gamma * h;
    }

    // Check if positive definite
    // FIXME
    bool isPosDef = (alpha[0][0] > 0 &&
                     alpha[0][0]*alpha[1][1] - alpha[0][1]*alpha[1][0] > 0 &&
                     alpha.determinant() > 0);

    HepVector dpar;
    double next = 0;
    vector<double> par_(3);

    if(isPosDef)
    {
      // Get Newton step
      dpar = solve(alpha, beta);

      for(int k = 0; k < 3; k++)
        par_[k] = par[k] + dpar[k]; 
 
      next = theModel.getValue(par_);
    }

    // Check 
    if(!isPosDef || next > old)
    {
      // Calculate step size;
      double lambda = dot(h,h) / fabs(dot(h, alpha * h)); // FIXME
      HepVector dpar = lambda * beta;

      double factor = 1.;

      do
      { 
        for(int k = 0; k < 3; k++)
          par_[k] = par[k] + factor * dpar[k]; 

        next = theModel.getValue(par_);

        factor /= 2;
      }
      while(next > old && factor > 1e-6);
    }

    if(next < old - 1e-3)
      par = par_;
    else
      next = old;

    move = next - old;

    old = next;

    err = alpha;
  }
  while(move < -1e-3 && nstep < 100);
}
