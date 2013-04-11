#include "../interface/CubicSpline.h"

using namespace std;

/****************************************************************************/
CubicSpline::CubicSpline()
{
}

/****************************************************************************/
CubicSpline::~CubicSpline()
{
}

/****************************************************************************/

/*
Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi =
f(xi), with x1 < x2 < .. . < xN, and given values yp1 and ypn for the first
derivative of the interpolating function at points 1 and n, respectively,
this routine returns an array y2[1..n] that contains the second derivatives
of the interpolating function at the tabulated points xi. If yp1 and/or ypn
are equal to 1 x 1030 or larger, the routine is signaled to set the
corresponding boundary condition for a natural spline, with zero second
derivative on that boundary.
*/

void CubicSpline::prepare
  (double x[], double y[], int n, double yp1, double ypn, double y2[])
{ 
  int i,k;
  double p,qn,sig,un, u[n-1];
  
  if (yp1 > 0.99e30) // Lower boundary
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  
  for (i=1;i<n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  
  if (ypn > 0.99e30) // Upper boundary
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2;k>=0;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
}

/****************************************************************************/

/*
Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the
xais in order), and given the array y2a[1..n], which is the output from
spline above, and given a value of x, this routine returns a cubic-spline
interpolated value y.
*/

void CubicSpline::interpolate
  (double xa[], double ya[], double y2a[], int n, double x,
   double *y0, double *y1, double *y2)
{
  int klo,khi, k;
  double h,b,a;

  klo = 0; khi = n-1;

  if(x <= xa[klo] || x >= xa[khi])
  { 
    // x is outside of limits
    klo = 0; khi = 1;
  }
  else
    while (khi-klo > 1)
    {
      k = (khi+klo)/2;
      if (xa[k] > x) khi=k;
                else klo=k;
    }

  h=xa[khi]-xa[klo];
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;

  *y0 = a*ya[klo] +
        b*ya[khi] +
        ((a*a*a-a) * y2a[klo] +
         (b*b*b-b) * y2a[khi]) * (h*h)/6.0;

/*
  *y1 = (ya[khi] - ya[klo])/h -
        (3*a*a-1)/6 * h * y2a[klo] +
        (3*b*b-1)/6 * h * y2a[khi];

  *y2 = a*y2a[klo] + b*y2a[khi];
*/
}

