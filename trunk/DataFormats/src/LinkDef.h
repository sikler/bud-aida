#include "../interface/TBunchCrossing.h"
#include "../interface/TVertex.h"
#include "../interface/TTrack.h"

#include "../interface/TSlimMeasurement.h"

#include "../interface/Point.h"
#include "../interface/Hit.h"
#include "../interface/Pixel.h"
#include "../interface/Crossing.h"
#include "../interface/TLayer.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class TBunchCrossing+;
#pragma link C++ class TVertex+;
#pragma link C++ class TTrack+;

#pragma link C++ class TSlimMeasurement+;
#pragma link C++ class TChipId+;

#pragma link C++ class Point+;

#pragma link C++ class TLayer+;
#pragma link C++ class vector<TLayer>+;

#pragma link C++ class Hit+;
#pragma link C++ class Pixel+;
#pragma link C++ class Pixel::Meas+;
#pragma link C++ class Pixel::Calc+;
#pragma link C++ class Crossing+;

//#pragma link C++ class pair<bool,bool>+;
//#pragma link C++ class pair<short int,short int>+;
#pragma link C++ class pair<int, pair<int, int> >+;

#endif
