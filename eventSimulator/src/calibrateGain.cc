#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>

#include <map>

#include "../../DataFormats/interface/TSlimMeasurement.h"
#include "../../DataFormats/interface/TChipId.h"

#include "../../siEnergyLoss/interface/ElossEstimator.h"

using namespace std;

char inputName[256], outputName[256], previousGain[256];

map<ChipId, vector<TSlimMeasurement> > hits;

/*****************************************************************************/
void options(int arg, char **arc)
{
  int i = 1;

  sprintf(previousGain,"");

  do
  {
    if(strcmp(arc[i],"-i") == 0) sprintf(inputName ,"%s",arc[++i]);
    if(strcmp(arc[i],"-o") == 0) sprintf(outputName,"%s",arc[++i]);

    if(strcmp(arc[i],"-previousGain") == 0)
                               sprintf(previousGain,"%s",arc[++i]);

    i++;
  }
  while(i < arg);
}

/*****************************************************************************/
int main(int arg, char **arc)
{
  options(arg, arc);


  cerr << "\033[22;31m" << "Calibrating chip gains..."
       << "\033[22;0m"  << endl;

  // Read maps
  ifstream inFile(inputName, ios::binary);

  int size; inFile.read((char *)&size, sizeof(size));

  for(int i = 0; i < size; ++i)
  {
    ChipId key; inFile.read((char *)&key, sizeof(key));
    int n;      inFile.read((char *)&n  , sizeof(n  ));

    vector<TSlimMeasurement> value;
    for(int k = 0; k < n; k++)
    {
      TSlimMeasurement slim;
      inFile.read((char *)&slim, sizeof(slim));
      value.push_back(slim);
    }

    hits[key] = value;
  }

  inFile.close();

  ElossEstimator elossEstimator;

  if(strcmp(previousGain,"") != 0)
  {
    ifstream filePrev;

    filePrev.open(previousGain);

    cerr << "  loading gains...";
    elossEstimator.loadGains(filePrev);
    cerr << " [done]" << endl;

    filePrev.close();
  }

  ofstream fileGains(outputName);
  elossEstimator.calibrateGains(hits, fileGains);
  fileGains.close();

  cerr << " [done]" << endl;

  return 0;
}

