#include "../interface/NewtonMinimizer.h"
#include "../interface/ElossModel.h"

#include <cmath>

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompLU.h"

using namespace std;

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
                               vector<double> & par, TMatrixD & err,
                               int & nstep)
{
  nstep = 0;

  double old = theModel.getValue(par);

  TVectorD g(3), h(3);

  double move;
  do
  {
    TMatrixD alpha(3,3); alpha.Zero();
    TVectorD beta(3); beta.Zero();

    // Get derivatives
    theModel.getDerivatives(par, beta,alpha);

    if(beta[0] == 0. && beta[1] == 0. && beta[2] == 0.)
      break;

    // Guess diagonal elements by averaging nonzero elements of alpha
    double sum[2] = {0.,0.};
    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      if(alpha[i][j] != 0.)
      { sum[1] += fabs(alpha[i][j]); sum[0] += 1.; }

    if(sum[0] != 0.) sum[1] /= sum[0];
                else sum[1] = 1.;

    // Fix some elements (used for strips) or replace by a reasonable value
    for(int k = 0; k < 3; k++)
    if(isFix[k])
    {
      beta[k] = 0.;

      for(int j = 0; j < 3; j++)
      if(j != k)
      { alpha[j][k] = 0.; alpha[k][j] = 0.; }

      alpha[k][k] = sum[1];
    }

    // try to fix
    for(int k = 0; k < 3; k++)
    if(alpha[k][k] == 0.)
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
      TVectorD go = g;
      g = beta;

      double gamma = ((g - go)*g) / (go*go);
      h = g + gamma * h;
    }

    // Check if positive definite
    bool isPosDef = (alpha[0][0] > 0 &&
                     alpha[0][0]*alpha[1][1] - alpha[0][1]*alpha[1][0] > 0 &&
                     alpha.Determinant() > 0);

    double next = 0;
    vector<double> par_(3);

    if(isPosDef)
    {
      // Get Newton step
      bool ok; TDecompLU lu(alpha);
      TVectorD dpar = lu.Solve(beta, ok);

      if(fabs(dpar[0]) < 1 &&
         fabs(dpar[1]) < 1 &&
         fabs(dpar[2]) < 0.2)
        for(int k = 0; k < 3; k++)
          par_[k] = par[k] + dpar[k]; 
        else // only if the step is small
        {
          for(int k = 0; k < 3; k++)
            par_[k] = par[k]; 
        }

      next = theModel.getValue(par_);
    }

    // if not pos def or chi2 increased, try to step with conjugate gradient 
    if( !isPosDef || next > old )
    {
      // Calculate step size;
      double lambda = (h*h) / fabs(h * (alpha * h));
      TVectorD dpar = lambda * beta;

      double factor = 1.;

      do
      { 
        for(int k = 0; k < 3; k++)
        if(fabs(factor*dpar[0]) < 1 &&
           fabs(factor*dpar[1]) < 1 &&
           fabs(factor*dpar[2]) < 0.2)
          par_[k] = par[k] + factor * dpar[k]; 
        else
        {
          par_[k] = par[k];
        }

        next = theModel.getValue(par_);

        factor /= 2;
      }
      while(next > old && factor > 1e-6);
    }

    if(next < old - 1e-3) par  = par_;
                     else next = old;

    move = next - old;
    old = next;
    err = alpha;
  }
  while(move < -1e-3 && nstep < 100);
}

