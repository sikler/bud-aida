#ifndef _KalmanTracking_
#define _KalmanTracking_

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "../../DataFormats/interface/Coord.h"

// Number of track parameters
#define nPars 5

// Dimension of a measurement
#define nMeas 2

// Global spatial dimension
#define nDims 3

using namespace std;
using namespace CLHEP;

/*****************************************************************************/
class Material
{
 public:
  double radius;
  double halfLength; // FIXME new
  double sigma_rphi;
  double sigma_z;
  double thickness;

  void print() const
  {
    cerr << " r=" << radius
         << " l/2=" << halfLength
         << " sigma_rphi=" << sigma_rphi
         << " sigma_z=" << sigma_z
         << " thickness=" << thickness
         << endl;
  }
};

/*****************************************************************************/
/*
class Coord
{
 public:
  int q;
  double x[nDims];
  double pt; // auxiliary, sometimes needed
  double pz;
  double eta;
  double p_; // total momentum
  double p[nDims];

  void print()
  {
    cerr << " q="  << q
         << " p=" <<  p_;

    for(int i = 0; i < nDims; i++)
      if(i == 0) cerr << " p=(" << p[i];
            else cerr << ","    << p[i];

    cerr << ")" << endl;
  }
};
*/

/*****************************************************************************/
class State
{
 public:
  State()
  {
    HepMatrix F_(nPars,nPars); F = F_;

    HepVector x_(nPars); x = x_;
    HepVector r_(nMeas); r = r_;

    HepMatrix C_(nPars,nPars); C = C_;
    HepMatrix R_(nMeas,nMeas); R = R_;

    Material material_;
    material_.radius     = 0.;
    material_.sigma_rphi = 0.;
    material_.sigma_z    = 0.;
    material_.thickness  = 0.;
    material = material_;
  }

  double m;

  HepMatrix F;

  HepVector x; // state
  HepVector r; // residual

  HepMatrix C; // covariance of state
  HepMatrix R; // covariance of residuals

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

  Material material;

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

class Coord;
class State;

class ClusterGenerator;

class KalmanTracking
{
 public:
  KalmanTracking();
 
  bool process(TTrack & simTrack, TTrack & recTrack);

 private:
  void readMaterial(vector<Material> & material);

  void convertToParameters(Coord & coord, State & state);
  void convertToCoordinates(State & state, Coord & coord);

  double getFlatRandom();
  double getGaussRandom();

  double multiScatt(State & state, int flag);
  double energyLoss(State & state, int flag);

  State makeRechit(State & state);

  bool propagate(State & state, const Material & material);

  double getAngle(double a[], double b[], double c[]);

  State fitCircle(vector<Layer> & layers);

  CLHEP::HepMatrix calculateF(State & rechit, Material & material);
  CLHEP::HepMatrix calculateQ(State & state, const HepMatrix & F);

  double getChi2(const State & state);

  void fit   (vector<Layer> & layers);
  void smooth(vector<Layer> & layers);

  State generate(TTrack & track);

  void simulate(TTrack & track, vector<Layer> & layers);
  vector<Coord> reconstruct(vector<Layer> & layers, double & c);

  // Material: silicon
  // Units: cm, GeV, T
  double B;
  double rho; // g cm^-3
  double X0;

  vector<Material> material;

  HepMatrix gR;

  ClusterGenerator * clusterGenerator;
};

#endif
