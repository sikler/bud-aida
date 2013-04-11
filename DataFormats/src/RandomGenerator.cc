#include "../interface/RandomGenerator.h"

#include <cmath>
#include <cstdlib>

/*****************************************************************************/
RandomGenerator::RandomGenerator()
{
}

/*****************************************************************************/
RandomGenerator::~RandomGenerator()
{
}

/*****************************************************************************/
double RandomGenerator::getFlatRandom()
{
  return drand48();
}

/****************************************************************************/
double RandomGenerator::getGaussRandom()
{
  return sqrt(-2*log(getFlatRandom())) * cos(M_PI*getFlatRandom());
}


