#include "../interface/ModelBichsel.h"
#include "../interface/Levels.h"

#include "TH1F.h"
#include "TH2F.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

static int N = 1e+5; // 1e+4

/*****************************************************************************/
void printToFile(TH1F & h1, const char* fileName)
{
  ofstream file(fileName);
  for(int i = 1; i <= h1.GetNbinsX(); i++) // avoid und and ovf
    file << " " << h1.GetBinCenter(i)
         << " " << h1.GetBinContent(i)
         << " " << h1.GetBinError(i)
         << endl;
  file.close();
}

/*****************************************************************************/
void printToFile(TH2F & h2, const char* fileName)
{
  ofstream file(fileName);

  for(int i = 1; i <= h2.GetNbinsX(); i++)
  {
    for(int j = 1; j <= h2.GetNbinsY(); j++)
      file << " " << h2.GetXaxis()->GetBinLowEdge(i)
           << " " << h2.GetYaxis()->GetBinLowEdge(j)
           << " " << h2.GetBinContent(i,j)
           << " " << h2.GetBinError(i,j)
           << " " << h2.GetXaxis()->GetBinCenter(i)
           << " " << h2.GetYaxis()->GetBinCenter(j)
           << " " << h2.GetXaxis()->GetBinWidth(i)
           << " " << h2.GetYaxis()->GetBinWidth(j)
           << endl;

      file << " " << h2.GetXaxis()->GetBinLowEdge(i)
           << " " << h2.GetYaxis()->GetXmax()
           << " 0 0 0. 0. 0. 0." << endl;
    file << endl;
  }

  for(int j = 1; j <= h2.GetNbinsY(); j++)
    file << " " << h2.GetXaxis()->GetXmax()
         << " " << h2.GetYaxis()->GetBinLowEdge(j)
         << " 0 0 0. 0. 0. 0." << endl;

    file << " " << h2.GetXaxis()->GetXmax()
         << " " << h2.GetYaxis()->GetXmax()
         << " 0 0 0. 0. 0. 0." << endl;

  file.close();
}

/*****************************************************************************/
int main()
{
  double betaGamma = 0.763;

  ModelBichsel theModel("Si");
  theModel.prepare(betaGamma);

  ofstream file("/home/sikler/dep.dat");

  for(int i = 0; i < 1e+4; i++)
    file << " " << theModel.generate(275.) << endl;

  file.close();

  double bg[5] = {0.32, 0.56, 1.00, 3.16, 10.};

  for(int ii = 0; ii < 5; ii++)
  {
    double betaGamma = bg[ii];

    cerr << "  betaGamma = " << betaGamma << endl;

    vector<double> deposit(N, 0.);

    ModelBichsel theModel("Si");
    theModel.prepare(betaGamma);

    TH2F histo("bichsel","bichsel",
                160,    0-5., 1600-5.,  // thickness [um]
                410,  -20-1.,  800-1.); // deposit   [keV]

    for(double x = 0; x < 1600-5; x +=10.)
    {
      cerr << "   depth = " << x << " um" << endl;
  
      for(int i = 0; i < N; i++)
      {
        // add last 10um, but only if x > 0
        if(x > 0)
          deposit[i] += theModel.generate(10.);

        // add noise
        double y = deposit[i] +
                   theModel.getGaussRandom() * Noise; 
        histo.Fill(x,y); // ke-
      }
    }

    char fileName[128];
    sprintf(fileName, "../data/bichsel_bg_%.2f.dat",betaGamma);
    printToFile(histo, fileName);
  }

  return 0;
}
