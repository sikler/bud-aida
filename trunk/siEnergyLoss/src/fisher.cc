#include <iostream>
#include <fstream>
#include <cmath>

#include "../interface/Levels.h"

#include "../interface/MostProbable.h"
#include "../interface/ModelBichsel.h"

#define Sqr(x) ((x) * (x))

using namespace std;

/*****************************************************************************/
void fisher()
{
  // keV, um

  MostProbable theMostProbable;
  ModelBichsel theModel("Si");

  ofstream file("../out/fisher.dat");

  double bg[5] = {0.32, 0.56, 1.00, 3.16, 10.};
  for(int ii = 0; ii < 5; ii++)
  {
    double betaGamma = bg[ii];

    theModel.prepare(betaGamma);

    double epsilon = theMostProbable.value(betaGamma) * 1e-1; // keV/um

    for(double x = 50; x < 1000; x +=50)
    {
      double Delta = epsilon * x * (1 + a*log(x/l0)); // keV

      double s[2] = {0,0};

      for(int i = 0; i < 1000; i++)
      {
        double y = theModel.generate(x); 
        double sigmaD = 2.0 + 0.10 * y;
  
        if(Delta > y - Nu * sigmaD)     
          s[1] += 1/Sqr(sigmaD);   

        s[0] ++;
      } 

      s[1] /= s[0];

      file << " " << betaGamma << " " << x
           << " " << s[1] << " " << epsilon <<  endl;

      cerr << " " << betaGamma << " " << x
           << " " << s[1] << " " << epsilon << endl;
    }

    file << endl;
    cerr << endl;
  }

  file.close();
}

/*****************************************************************************/
int main(int arg, char **arc)
{
  fisher();

  return 0;
}
