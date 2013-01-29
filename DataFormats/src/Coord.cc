#include "../interface/Coord.h"

#include <iostream>
using namespace std;

ClassImp(Coord)

Coord::Coord()
{
}

Coord::~Coord()
{
}

void Coord::print()
{
  cerr << " q="  << q
       << " p=" <<  p_;

  for(int i = 0; i < nDims; i++)
    if(i == 0) cerr << " p=(" << p[i];
          else cerr << ","    << p[i];

  cerr << ")" << endl;
}

