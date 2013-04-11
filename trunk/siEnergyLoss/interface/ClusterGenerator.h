#ifndef _ClusterGenerator_h_
#define _ClusterGenerator_h_

#include <fstream>
#include <vector>

class ModelBichsel;
class RandomGenerator;

class Pixel;
class Hit;
class TChipId;
class TLayer;

class ClusterGenerator
{
 public:
   ClusterGenerator();
   ~ClusterGenerator();

   Hit create(double betaGamma, double theta, double phi,
              TChipId & chipId, TLayer * mat, int ilayer);

 private:
   void addCoupling(Hit & hit, TLayer * unit);

   double getIntegral(double x, int x0, int x1, double sigma);

   void matchOrAdd(std::vector<Pixel> & pixels,
                   int x, int y, double Delta);

   ModelBichsel    * theModel;
   RandomGenerator * theRandom;
};

#endif
