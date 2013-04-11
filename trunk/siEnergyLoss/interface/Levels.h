#ifndef _Levels_h_
#define _Levels_h_

// Silicon
const double Noise     = 1.5e-3; // MeV move to FIXME

// Energy loss parametrization
const double Nu     = 0.65;
const double a      = 0.07;
const double Sigma0 = 2e-3; // MeV
const double b      = 0.095;

// Reference path length
const double l0 = 450e-4;   // cm

// Bethe-Bloch, silicon
const double K = 0.307075e-3; // GeV cm^2 / g
const int    Z = 14;
const double A = 28.;

const double rho = 2.33;      // g/cm3

#endif
