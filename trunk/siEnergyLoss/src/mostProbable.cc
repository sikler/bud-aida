#include <iostream>
#include <fstream>
#include <cmath>

#include "../interface/MostProbable.h"
#include "../interface/ModelBichsel.h"

#include "TH1F.h"
#include "TF1.h"

#define Sqr(x) ((x) * (x))

using namespace std;

/*****************************************************************************/
void process()
{
  // keV, um

  MostProbable theMostProbable;
  ModelBichsel theModel("Si");

  ofstream file("../out/mostProbable.dat");

  for(double betaGamma = 1.0; betaGamma < 10.; betaGamma *=1.1)
  {
    theModel.prepare(betaGamma);

    double epsilon = theMostProbable.value(betaGamma) * 1e-1; // keV/um

    double x = 300;
    {
      cerr << " bg = " << betaGamma << " " << x << " um" << endl;

      double Delta = epsilon * x * (1 + 0.07*log(x/300)); // keV

      TH1F * hisli = new TH1F("hisli","hisli", 100,    (Delta/2),    (2*Delta));
      TH1F * hislo = new TH1F("hislo","hislo", 100, log(Delta/2), log(2*Delta));

      for(int i = 0; i < 10000; i++)
      {
        double y = theModel.generate(x); 

        hisli->Fill(   (y));
        hislo->Fill(log(y));
      } 

      double mpli = hisli->GetBinCenter(hisli->GetMaximumBin());
      double mplo = hislo->GetBinCenter(hislo->GetMaximumBin());

      file << " " << betaGamma << " " << x
           << " " << mpli/30
           << " " << mplo - log(30)
           << " " << Delta << endl;

      cerr << " " << betaGamma << " " << x
           << " " << mpli/30
           << " " << mplo - log(30)
           << " " << Delta << endl;

      delete hisli;
      delete hislo;
    }

    //file << endl;
    //cerr << endl;
  }

  file.close();
}

/*****************************************************************************/
int main(int arg, char **arc)
{
  process();

  return 0;
}
