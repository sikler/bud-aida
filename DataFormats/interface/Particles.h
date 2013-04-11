#ifndef _Particles_h_
#define _Particles_h_

enum Particles { unknown=0, elec=1, pion=2, kaon=3, prot=4, muon=5, nParts=6 };
static double mass[nParts] =
   { 0, 0.511e-3, 0.139570, 0.493677, 0.938272, 0.105658 }; // GeV
static int pdg[nParts] =
   { 0, 11      , 211     , 321     , 2212    , 13       };

#endif
