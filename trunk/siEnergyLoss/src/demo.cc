#include <fstream>
#include <iostream>
#include <cmath>

#include "../interface/ModelBichsel.h"

using namespace std;

/*****************************************************************************/
int main()
{
  double betaGamma = 3.17;

  {
    ModelBichsel theModel("Si");
    theModel.prepare(betaGamma);

    ofstream file("/home/sikler/weightedMean/data/demo_Si.dat");

    for(int i = 0; i < 10; i++)
    {
      file << " # " << theModel.generateDemo(300., file) << endl;
      file << endl << endl;
    } 

    file.close();
  }

  {
    ModelBichsel theModel("Ne");
    theModel.prepare(betaGamma);


    ofstream file("/home/sikler/weightedMean/data/demo_Ne.dat");

    for(int i = 0; i < 10; i++)
    {
      file << " # " << theModel.generateDemo(1., file) << endl;
      file << endl << endl;
    }

    file.close();
  } 
}
