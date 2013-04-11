#ifndef _KalmanTracking_h_
#define _KalmanTracking_h_

#include <fstream>
#include <cmath>

#include "../../DataFormats/interface/Point.h"
#include "../../DataFormats/interface/TLayer.h"

#include "TMatrixD.h"
#include "TVectorD.h"

// Number of track parameters
#define nPars 5

// Dimension of a measurement
#define nMeas 2

// Global spatial dimension
#define nDims 3

using namespace std;


/*****************************************************************************/
class State
{
 public:
  State()
  {
    TMatrixD F_(nPars,nPars); F.ResizeTo(F_); F = F_;

    TVectorD x_(nPars); x.ResizeTo(x_); x = x_;
    TVectorD r_(nMeas); r.ResizeTo(r_); r = r_;

    TMatrixD C_(nPars,nPars); C.ResizeTo(C_); C = C_;
    TMatrixD R_(nMeas,nMeas); R.ResizeTo(R_); R = R_;

    TLayer material_;
    material_.radius     = 0.;
    material_.sigma_rphi = 0.;
    material_.sigma_z    = 0.;
    material_.thickness  = 0.;
    material = material_;
  }

  double m;

  TMatrixD F;

  TVectorD x; // state
  TVectorD r; // residual

  TMatrixD C; // covariance of state
  TMatrixD R; // covariance of residuals

  double getKappa() { return x[0]; }
  double getTheta() { return x[1]; }
  double getPsi()   { return x[2]; }
  double getRPhi()  { return x[3]; }
  double getZ()     { return x[4]; }

  void setKappa(double value) { x[0] = value; }
  void setTheta(double value) { x[1] = value; }
  void setPsi  (double value) { x[2] = value; }
  void setRPhi (double value) { x[3] = value; }
  void setZ    (double value) { x[4] = value; }

  TLayer material; // FIXME REMOVE!!!!

  double xpos()
  {
    double phi = 0;
    if(material.radius > 0) phi = getRPhi() / material.radius;
    return material.radius * cos(phi);
  }
  double ypos()
  {
    double phi = 0;
    if(material.radius > 0) phi = getRPhi() / material.radius;
    return material.radius * sin(phi);
  }
  double zpos() { return getZ(); }

  double chi2;

  void print()
  {
    cerr << " p="      << fabs(1/getKappa())
         << " theta="  << getTheta()
         << " psi="    << getPsi()
         << " phi="    << getRPhi() / material.radius
         << " radius=" << material.radius
         << endl;
  }
};

/*****************************************************************************/
class Layer
{
 public:
  State simulated;
  State measured;

  State predicted;
  State updated;
  State smoothed;
};

/*****************************************************************************/
class TVertex;
class TTrack;

class Point;
class State;

class ClusterGenerator;
class RandomGenerator;

class KalmanTracking
{
 public:
  KalmanTracking(std::vector<TLayer> & materials, float gainHalfWidth_, bool simu);
  virtual ~KalmanTracking();
 
  bool process(TTrack & simTrack, TTrack & recTrack); // TO BE REMOVED

//  bool simulate   (TTrack & simTrack);
//  bool reconstruct(TTrack & simTrack, TTrack & recTrack); // ??

 private:
  void readMaterial(vector<TLayer> & material);

  void convertToParameters(Point & coord, State & state);
  void convertToCoordinates(State & state, Point & coord);

  double multiScatt(State & state, int flag);
  double energyLoss(State & state, int flag);

  State makeRechit(State & state);

  bool propagate(State & state, const TLayer & material);

  double getAngle(double a[], double b[], double c[]);

  State fitCircle(vector<Layer> & layers);

  TMatrixD calculateF(State & rechit, TLayer & material);

  TMatrixD outerProduct(const TVectorD & a, const TVectorD & b);

  TMatrixD calculateQ(State & state, const TMatrixD & F);

  double getChi2(const State & state);

  void fit   (vector<Layer> & layers);

  void set(TMatrixD & a, TMatrixD & b);

  void smooth(vector<Layer> & layers);

  State generate(TTrack & track);

  void simulate(TTrack & track, vector<Layer> & layers);
  vector<Point> reconstruct(vector<Layer> & layers, double & c, int pdgId);

  // Material: silicon
  // Units: cm, GeV, T
  double B;
  double X0;

  float gainHalfWidth;

  vector<TLayer> material;

  TMatrixD gR;

  std::ofstream fileChi;

  ClusterGenerator * clusterGenerator;
  RandomGenerator  * theRandom;
};

#endif
