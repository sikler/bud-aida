#ifndef _RandomGenerator_h_
#define _RandomGenerator_h_

class RandomGenerator
{
  public:
    RandomGenerator();
    ~RandomGenerator();

    double getFlatRandom();
    double getGaussRandom();
};

#endif
