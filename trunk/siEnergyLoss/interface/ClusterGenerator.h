#ifndef _ClusterGenerator_h_
#define _ClusterGenerator_h_

#include <fstream>
#include <vector>

class ModelBichsel;
class Pixel;
class Hit;

class ClusterGenerator
{
 public:
   ClusterGenerator();
   ~ClusterGenerator();

   Hit create(int type, double betaGamma, double theta, double phi);

//   void run(int n, int type);

 private:
   double getFlatRandom();
   double getGaussRandom();

   void addCoupling(Hit & hit);

   double getIntegral(double x, int x0, int x1, double sigma);

   void matchOrAdd(std::vector<Pixel> & pixels,
                   int x, int y, double Delta);

//   void writeHit(ModelBichsel & theModel, std::ofstream & file,
//                 int type, double & betaGamma);

   ModelBichsel * theModel;
};

#endif
