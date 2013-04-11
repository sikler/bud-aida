#ifndef _ModelBichsel_h_
#define _ModelBichsel_h_

#include <fstream>

class RandomGenerator;

class ModelBichsel
{
  public:
    ModelBichsel(const char * elem);
    ~ModelBichsel();

    void prepare(double betaGamma);

    double generate(double thickness);
    double generateDemo(double thickness, std::ofstream & file);

  private:
    void readMeanNumberOfCollisions(const char * elem);
    void readCumulativeCrossSection(const char * elem);

    #define mncRows 51
    double mnc[3][mncRows];

    #define ccsCols 10
    #define ccsRows 10001
    double cbg[ccsCols];
    double ccs[ccsCols][3][ccsRows];

    int ic;
    double xc, sigmaT;

    RandomGenerator * theRandom;
};

#endif
