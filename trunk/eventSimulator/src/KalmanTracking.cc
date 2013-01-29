#include "../interface/KalmanTracking.h"

#include "../../DataFormats/interface/TVertex.h"
#include "../../DataFormats/interface/TTrack.h"

//#include "../interface/Hit.h"
#include "../../siEnergyLoss/interface/ClusterGenerator.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace CLHEP;

//#define Debug

#define sqr(x) ((x) * (x))

enum Flag { Random, Mean, Sigma };

// FIXME 
const double recMass = 0.139;
const bool useOrigin = true;
const bool useLogLikelihood = false;
//const double maxChi = 500000.;
const double maxChi = 100.;

/*****************************************************************************/
KalmanTracking::KalmanTracking()
{
  // silicon
  rho = 2.329;        // g cm^-3
  X0  = 21.82 / rho;

  readMaterial(material);

  clusterGenerator = new ClusterGenerator();
}

/*****************************************************************************/
void KalmanTracking::readMaterial(vector<Material> & material)
{
  cerr << " reading material description..";

  ifstream file;
// FIXME
/*
  if(experiment == atlas)      file.open("../data/atlas.dat");
  if(experiment == cms  )      file.open("../data/cms.dat"  );
  if(experiment == alice)      file.open("../data/alice.dat");
  if(experiment == perfect)    file.open("../data/perfect.dat");
  if(experiment == atlas_full) file.open("../data/atlas_full.dat");
*/
  file.open("../data/cms.dat"  );

  // Read magnetic field
  file >> B;

  while(file.eof() == false)
  {
    float f;
    file >> f;

    if(file.eof() == false)
    {
      Material m;
      int d;

      m.radius = f;

      file >> f; m.halfLength = f;
 
      file >> d; m.sigma_rphi = d * 1e-4; // um -> cm
      file >> d; m.sigma_z    = d * 1e-4; // um -> cm
      file >> f; m.thickness  = (f * 1e-2) * X0;   // x/X0 -> x

//      m.print();

      material.push_back(m);
    }
  }

  file.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
void KalmanTracking::convertToParameters(Coord & coord, State & state)
{ // x,p ->
  double pz =                                              coord.p[2];
  double p  = sqrt(sqr(coord.p[0]) + sqr(coord.p[1]) + sqr(coord.p[2]));
  double r  = sqrt(sqr(coord.x[0]) + sqr(coord.x[1]));

  state.material.radius = r;
  state.setKappa(coord.q / p);
  state.setTheta(acos(pz / p));

  double phi;
  if(r > 0) phi = atan2(coord.x[1], coord.x[0]);
       else phi = 0.;
  
  state.setPsi  (atan2(coord.p[1], coord.p[0]));
  state.setRPhi (r * phi); 
  state.setZ    (coord.x[2]);
}

/*****************************************************************************/
void KalmanTracking::convertToCoordinates(State & state, Coord & coord)
{ // kappa, theta, psi, rphi, z ->
  double r   = state.material.radius; 

  coord.q    = (state.getKappa() > 0 ? 1 : -1);

  double phi;
  if(r > 0) phi = state.getRPhi()/r;
       else phi = 0.;

  coord.x[0] = r  * cos(phi);
  coord.x[1] = r  * sin(phi);
  coord.x[2] = state.getZ();

  coord.p_   = fabs(1. / state.getKappa());
  coord.pt   = coord.p_ * sin(state.getTheta());
  coord.pz   = coord.p_ * cos(state.getTheta());

  coord.eta = -log(tan(state.getTheta()/2));

  coord.p[0] = coord.pt * cos(state.getPsi());
  coord.p[1] = coord.pt * sin(state.getPsi());
  coord.p[2] = coord.p_ * cos(state.getTheta());
}

/*****************************************************************************/
double KalmanTracking::getFlatRandom()
{
  return drand48();
}

/****************************************************************************/
double KalmanTracking::getGaussRandom()
{
  return sqrt(-2*log(getFlatRandom())) * cos(M_PI*getFlatRandom());
}

/*****************************************************************************/
double KalmanTracking::multiScatt(State & state, int flag)
{
  Coord coord;
  convertToCoordinates(state, coord); // just to get p_

  double E = sqrt(sqr(coord.p_) + sqr(state.m));
  double beta = coord.p_ / E;

  double phi = 0.;
  if(state.material.radius > 0) phi = state.getRPhi() / state.material.radius;
  double d = state.material.thickness /
             fabs(cos(state.getPsi() - phi) / sin(state.getTheta()));
  
  double mean  = 0.;
  double sigma = 13.6e-3/(beta * coord.p_) * sqrt(d / X0);

  if(flag == Sigma) return sigma;

  double alpha, u[2];

  // Rotate (pz,pt)
  alpha = mean;
  if(flag == Random)
    alpha += sigma * getGaussRandom();

  u[0] =   coord.pz*cos(alpha) + coord.pt*sin(alpha); // pz'
  u[1] = - coord.pz*sin(alpha) + coord.pt*cos(alpha); // pt'

  coord.p[2] = u[0];

  coord.p[0] *= u[1] / coord.pt;
  coord.p[1] *= u[1] / coord.pt;

  // Rotate (px,py)
  alpha = mean;
  if(flag == Random)
    alpha += sigma * getGaussRandom();

  u[0] =   coord.p[0]*cos(alpha) + coord.p[1]*sin(alpha);
  u[1] = - coord.p[0]*sin(alpha) + coord.p[1]*cos(alpha);

  for(int k = 0 ; k < 2; k++)
    coord.p[k] = u[k];

  convertToParameters(coord, state);

  return 0.;
}

/*****************************************************************************/
double KalmanTracking::energyLoss(State & state, int flag)
{
  Coord coord;
  convertToCoordinates(state, coord); // just to get pt

  const double K = 0.307075e-3;
  const double Z = 14.;
  const double A = 28.;
  const double me = 0.511e-3;

  double E = sqrt(sqr(coord.p_) + sqr(state.m));
  double beta = coord.p_ / E;
  double gamm = 1/sqrt(1 - sqr(beta));

  double phi = 0.;
  if(state.material.radius > 0) phi = state.getRPhi() / state.material.radius;
  double d = state.material.thickness /
             fabs(cos(state.getPsi() - phi) / sin(state.getTheta()));

  double xi = K/2 * Z/A * rho * d / sqr(beta);
  double Tmax = 2*me*sqr(beta*gamm)/(1 + 2*gamm*me/state.m + sqr(me/state.m));
  double I = 16 * pow(Z,0.9) * 1e-9;

  // Most probable value
  double mean  = - xi * (log(2*me*sqr(beta*gamm*Tmax/I)) - 2*sqr(beta));
  double sigma = 4.018 / (2 * sqrt(2 * log(2))) * xi;

  if(flag == Sigma)
    return coord.q * (sigma / beta) / sqr(coord.p_);

  double dE  = mean;

  if(flag == Random)
    dE += sigma * getGaussRandom(); 

  double dp_ = dE / beta;

  // Slow down
  for(int k = 0 ; k < 2; k++)
    coord.p[k] *= (coord.p_ + dp_) / coord.p_;

  convertToParameters(coord, state);

  return 0.;
}

/*****************************************************************************/
State KalmanTracking::makeRechit(State & state)
{
  State rechit = state;

// MODIFICATION STARTS HERE
  float factor = 1.;
  if(getFlatRandom() < 1e-2) factor = 2.; // outliers, another FIXME

  rechit.setRPhi(state.getRPhi() +
                 state.material.sigma_rphi * factor * getGaussRandom());

  if(state.material.sigma_z < 5.) // smaller than 5 cm 
    rechit.setZ(state.getZ() +
                state.material.sigma_z     * factor * getGaussRandom());
  else
    rechit.setZ(state.getZ());
// MODIFICATION ENDS HERE

  return rechit;
}

/*****************************************************************************/
// Propagate to barrel, r1
bool KalmanTracking::propagate(State & state, const Material & material)
{
  double r1 = material.radius;

  Coord coord;
  convertToCoordinates(state, coord);

  // First look at two dimensional projection

  // radius of the circle, r2
  double r2 = coord.pt / (0.3e-2 * B);

  double n[2], C[2];
  n[0] =   coord.q * coord.p[1] / coord.pt;
  n[1] = - coord.q * coord.p[0] / coord.pt;

  // Center of the circle
  for(int k = 0 ; k < 2; k++)
    C[k] = coord.x[k] + n[k] * r2;

  // distance beamline - C
  double r12 = sqrt(sqr(C[0]) + sqr(C[1]));

  if(r12 > r1 + r2)
  { // disjoint
    return false;
  }
  if(r12 < fabs(r1 - r2)) 
  { // contains
    return false;
  }

  double P[3]; 
  double phi  = atan2(C[1],C[0]);
  double dphi = coord.q * acos((sqr(r1) + sqr(r12) - sqr(r2))/(2*r1*r12));
  P[0] = r1 * cos(phi + dphi);
  P[1] = r1 * sin(phi + dphi);

  double u[2];
  double alpha = asin( (n[0]*(P[1]-C[1]) - n[1]*(P[0]-C[0])) / r2 );
  u[0] =   coord.p[0]*cos(alpha) + coord.p[1]*sin(alpha);
  u[1] = - coord.p[0]*sin(alpha) + coord.p[1]*cos(alpha);

  P[2] =   coord.x[2] + coord.q * r2 * alpha / tan(state.getTheta());

  for(int k = 0 ; k < 3; k++)
    coord.x[k] = P[k];

  // p[2] is unchanged
  for(int k = 0 ; k < 2; k++) 
    coord.p[k] = u[k];

  if(fabs(coord.x[2]) < material.halfLength) // FIXME
  {
    convertToParameters(coord, state);

    state.material= material;

    return true;
  }
  else
  {
    // we are outside halfLengts at this radius
    return false;
  }
}

/*****************************************************************************/
double KalmanTracking::getAngle(double a[], double b[], double c[])
{
  double r2 = sqr(a[0] - b[0]) + sqr(a[1] - b[1]);  

  return asin(( (a[0]-b[0])*(c[1]-b[1]) -
                (a[1]-b[1])*(c[0]-b[0])) / r2);
}

/*****************************************************************************/
// Input  : layers[1-3].measured (first three measured points) 
// Output : layers[0].predicted  (vertex)
State KalmanTracking::fitCircle(vector<Layer> & layers)
{ 
  vector<State> rechits;
  rechits.push_back(layers[1].measured);
  rechits.push_back(layers[2].measured);
  rechits.push_back(layers[3].measured);

  // Use first three points
  double r, a[3],b[3],c[3];

  r = rechits[0].material.radius;
  a[0] = r * cos(rechits[0].getRPhi() / r);
  a[1] = r * sin(rechits[0].getRPhi() / r);
  a[2] = rechits[0].getZ();

  r = rechits[1].material.radius;
  b[0] = r * cos(rechits[1].getRPhi() / r);
  b[1] = r * sin(rechits[1].getRPhi() / r);
  b[2] = rechits[1].getZ();

  r = rechits[2].material.radius;
  c[0] = r * cos(rechits[2].getRPhi() / r);
  c[1] = r * sin(rechits[2].getRPhi() / r);
  c[2] = rechits[2].getZ();

  double ab[2], ac[2];
  for(int k = 0; k < 2; k++)
  {
    ab[k] = a[k] - b[k];
    ac[k] = a[k] - c[k];
  } 

  double z12 = sqr(a[0]) + sqr(a[1]);
  double z22 = sqr(b[0]) + sqr(b[1]);
  double z32 = sqr(c[0]) + sqr(c[1]);

  double area = -(ab[0]*ac[1] - ab[1]*ac[0])/2;

  // Center of the circle
  double C[2];
  C[0] = -(z12*(b[1]-c[1]) + z22*(c[1]-a[1]) + z32*(a[1]-b[1])) / (4 * area);
  C[1] =  (z12*(b[0]-c[0]) + z22*(c[0]-a[0]) + z32*(a[0]-b[0])) / (4 * area);

  // Radius
  double r2 = sqrt(sqr(C[0] - a[0]) + sqr(C[1] - a[1]));

  // p_T
  double pt = 0.3e-2 * B * r2;

  // Signed impact parameter
  double d0 = sqrt(sqr(C[0]) + sqr(C[1])) - r2;

  Coord coord;

  for(int k = 0; k < 2; k++)
    coord.x[k] = d0 * C[k] / sqrt(sqr(C[0]) + sqr(C[1]));

  // z posiiton of closest approaxh
  double aCx = getAngle(a,C,coord.x);
  double xCb = getAngle(coord.x,C,b);
  double aCc = getAngle(a,C,c);

  coord.x[2] = (aCx*b[2] + xCb*a[2]) / (aCx + xCb);

  coord.q = (area > 0 ? 1 : -1);

  double phi = atan2(C[1],C[0]) + coord.q * M_PI/2;

  coord.p[0] = pt * cos(phi);
  coord.p[1] = pt * sin(phi);

  coord.p[2] = pt *(c[2] - a[2]) / fabs(r2 * aCc);

  State vertex;
  vertex.m = recMass;
  convertToParameters(coord, vertex);

#ifdef Debug
  cerr << " pt = " << pt << "; eta = "
       << -log(tan(vertex.getTheta()/2)) << "; d0 = " << d0 << endl;
#endif

  if(useOrigin)
    vertex.material.radius = 0.;

  return vertex;
}

/*****************************************************************************/
HepMatrix KalmanTracking::calculateF(State & rechit, Material & material)
{
  HepMatrix F(nPars,nPars, 0);

  // Calculate derivatives
  State refer = rechit;

  propagate(refer, material);

  double dx = 1e-3;

  for(int j = 0; j < nPars; j++)
  {
    State state = rechit;
    state.x[j] += dx;

    propagate(state, material);

    for(int i = 0; i < nPars; i++)
      F[i][j] = (state.x[i] - refer.x[i]) / dx;
  }

  return F;
}

/*****************************************************************************/
HepMatrix KalmanTracking::calculateQ(State & state, const HepMatrix & F)
{
  HepMatrix Q(nPars,nPars);
  HepVector k(nPars), t(nPars), p(nPars);

  for(int i = 0; i < nPars; i ++)
  {
    k[i] = F[i][0]; // 0th column
    t[i] = F[i][1]; // 1st column
    p[i] = F[i][2]; // 2nd column
  }

  double sigma_kappa = energyLoss(state, Sigma);
  double sigma_theta = multiScatt(state, Sigma);
  double sigma_psi   = sigma_theta;

  Q = (k * k.T()) * sqr(sigma_kappa) +
      (t * t.T()) * sqr(sigma_theta) +
      (p * p.T()) * sqr(sigma_psi  );

  return Q;
}

/*****************************************************************************/
double KalmanTracking::getChi2(const State & state)
{
  int ierr;

  double chi2 = (state.r.T() * state.R.inverse(ierr) * state.r)[0];

  if(useLogLikelihood)
    chi2 += log(2 * M_PI * state.R.determinant());

  return chi2;
}

/*****************************************************************************/
void KalmanTracking::fit(vector<Layer> & layers)
{
  // Constant measurement and identity matrices
  HepMatrix H(nMeas,nPars,0.); H[0][3] = 1.; H[1][4] = 1.;
  HepMatrix Ip(nPars,nPars, 1);
  HepMatrix Im(nMeas,nMeas, 1);

  int ierr;

#ifdef Debug
  ofstream file("../out/fit.dat",ios::app);
#endif

  for(vector<Layer>::iterator layer = layers.begin();
                              layer!= layers.end(); layer++)
  {
    if(layer == layers.begin())
    {
      // Fit circle to first three hits (triplet)
      layer->updated = fitCircle(layers);

      // Initial guess on covariance
      HepMatrix C(nPars,nPars,0.);
      for(int i = 0; i < nPars; i++)
      for(int j = 0; j < nPars; j++)
      {
        if(i == j) C[i][j] = 1e+6;
              else C[i][j] = 1e+5;
      }
  
      layer->updated.C = C;
    }
    else
    {
      /************
      * Predicted *
      ************/
  
      // Target material
      Material material = layer->simulated.material;
  
      // Start with previous updated state
      State state = (layer-1)->updated;

      // Transient F, Q, V, G
      HepMatrix F = calculateF(state, material);

      HepMatrix Q(nPars,nPars, 0.);
      if(layer > layers.begin() + 1)
      Q = calculateQ(state, F);

      HepMatrix V(nMeas,nMeas, 0.);
      V[0][0] = sqr(layer->simulated.material.sigma_rphi);
      V[1][1] = sqr(layer->simulated.material.sigma_z   );

      HepMatrix G(nMeas,nMeas, 0.);
      G[0][0] = 1. / sqr(layer->simulated.material.sigma_rphi);
      G[1][1] = 1. / sqr(layer->simulated.material.sigma_z   );

      (layer-1)->updated.F = F;

      // Try to propagate
      if(!propagate(state, material)) 
      {
#ifdef Debug
        cerr << " could not propagate!" << endl;
#endif

        layers.erase(layer,layers.end());

        break;
      }

      // Energy loss
      energyLoss(state, Mean);

      // Extrapolation of the covariance
      state.C = F * state.C * F.T() + Q; 

      // Measured
      HepMatrix m = H * layer->measured.x;

      // Residuals of predictions
      state.r = m - H * state.x;

      // Covariance of predicted residuals
      state.R = V + H * state.C * H.T();

      // Chi2 increment
      state.chi2 = getChi2(state);

      if(state.chi2 > maxChi * maxChi)
      {
#ifdef Debug
        cerr << " chi too big = " << sqrt(state.chi2) << endl;
#endif

        layers.erase(layer,layers.end());
        break;
      }

      // Set predicted
      layer->predicted = state;

      /************
      * Filtering *
      ************/
  
      // Kalman gain matrix
      HepMatrix K =
        state.C * H.T() * (V + H * state.C * H.T()).inverse(ierr);
  
      // Update of the state vector
      state.x += K * (m - H * state.x);
  
      // Update of the covariance matrix 
      state.C = (Ip - K * H) * state.C;

      // Filtered redisuals
      state.r = m - H * state.x;
  
      // Covariance matrix of filtered residuals
      state.R = (Im - H * K) * V;
  
      // Chi2 increment
      state.chi2 = getChi2(state);

      if(state.chi2 > maxChi * maxChi)
      {
#ifdef Debug
        cerr << " chi too big = " << sqrt(state.chi2) << endl;
#endif

        layers.erase(layer,layers.end());
        break;
      }
  
#ifdef Debug
      Coord coord;
      convertToCoordinates(state, coord);
  
      cerr <<  " fitted pt = " << coord.pt
           << "; eta = " << coord.eta
           << "; chi2 = " << layer->predicted.chi2 << "/" 
           << state.chi2 << endl;
#endif
   
      // Set updated 
      layer->updated = state;
    }

#ifdef Debug
    file << " " << layer->updated.xpos()
         << " " << layer->updated.ypos()
         << " " << layer->updated.zpos()
         << endl;
#endif
  }

#ifdef Debug
  file << endl << endl;
  file.close();
#endif

//  while(getchar() == 0);
}

/*****************************************************************************/
void KalmanTracking::smooth(vector<Layer> & layers)
{
  HepMatrix H(nMeas,nPars,0.); H[0][3] = 1.; H[1][4] = 1.;

  // Global covariance
//  HepMatrix R(layers.size()-1, layers.size()-1);
  HepMatrix R(nMeas*(layers.size()-1),
              nMeas*(layers.size()-1));
  vector<HepMatrix> C(layers.size()-1);

#ifdef Debug
  ofstream file("../out/smoothing.dat",ios::app);
#endif

  for(vector<Layer>::iterator layer = layers.end()   - 1;
                              layer > layers.begin(); layer--) // FIXME was >
  {
    HepMatrix V(nMeas,nMeas, 0.);
    V[0][0] = sqr(layer->simulated.material.sigma_rphi);
    V[1][1] = sqr(layer->simulated.material.sigma_z   );


    if(layer == layers.end() - 1)
    {
      layer->smoothed = layer->updated;

      int i = layer - (layers.begin()+1);
      C[i] = layer->smoothed.C;

      HepMatrix r = V - (H * C[i] * H.T());

      for(int i1 = 0; i1 < nMeas; i1++)
      for(int j1 = 0; j1 < nMeas; j1++)
        R[nMeas*i + i1][nMeas*i + j1] = r[i1][j1];
    }
    else
    {
      // Smoother gain matrix
      int ierr;
      HepMatrix A = layer->updated.C * 
                           layer->updated.F.T() *
                          (layer+1)->predicted.C.inverse(ierr);

      State state = layer->updated;
  
      state.x += A * ( (layer+1)->smoothed.x - (layer+1)->predicted.x );
      state.C += A * ( (layer+1)->smoothed.C - (layer+1)->predicted.C ) * A.T();

      // Measured
      HepMatrix m = H * layer->measured.x;
  
      state.r = m - H * state.x;
  
/*
      state.R += - H * A *
                   ((layer+1)->smoothed.C - (layer+1)->predicted.C ) *
                   A.T() * H.T();
*/
      state.R = V - H * state.C * H.T();

      {
        int i = layer - (layers.begin()+1);
        int n = layers.size() - 1;

        C[i] = state.C; 

        for(int j = i+1; j < n; j++)
          C[j] = A * C[j];

        for(int j = i; j < n; j++)
        {
          HepMatrix r;

          if(i == j) r = V - (H * C[j] * H.T());
                else r =   - (H * C[j] * H.T());

          for(int i1 = 0; i1 < nMeas; i1++)
          for(int j1 = 0; j1 < nMeas; j1++)
          {
            R[nMeas*i + i1][nMeas*j + j1] = r[i1][j1];
            R[nMeas*j + i1][nMeas*i + j1] = r[i1][j1];
          }
        }
      }
  
      // Chi2 increment
      state.chi2 = getChi2(state);

      if(state.chi2 > maxChi * maxChi)
      {
// FIXME< new
//#ifdef Debug
        cerr << " chi too big = " << sqrt(state.chi2) << endl;
//#endif

        layers.erase(layer,layers.end());
        break;
      }


      layer->smoothed = state;

#ifdef Debug
    Coord coord;
    convertToCoordinates(state, coord);

    cerr <<  " smooth pt = " << coord.pt
         << "; eta = " << coord.eta
         << "; chi2 = " << state.chi2 << " " << state.material.radius << endl;
#endif
    }
#ifdef Debug
    file << " " << layer->updated.xpos()
         << " " << layer->updated.ypos()
         << " " << layer->updated.zpos()
         << endl;
#endif
  }
  
  gR = R;

#ifdef Debug
  file << endl << endl;
  file.close();
#endif
}

/*****************************************************************************/
State KalmanTracking::generate(TTrack & track)
{
  Coord coord;

  // Starts from origo
  coord.x[0] = 0.; // FIXME
  coord.x[1] = 0.;
  coord.x[2] = track.z;

  coord.pt = track.pt;

  coord.p[0] = coord.pt * cos(track.phi);
  coord.p[1] = coord.pt * sin(track.phi);
  coord.p[2] = coord.pt * sinh(track.eta);

  // Charge
  coord.q = track.charge;

  State state;
  convertToParameters(coord, state);

  state.material.radius = 0.;

//  state.m = genMass;
  state.m = 0.139; // FIXME // MASSS!

  if(abs(track.pdgId) ==  211) state.m = 0.139;
  if(abs(track.pdgId) ==  321) state.m = 0.493;
  if(abs(track.pdgId) == 2212) state.m = 0.938;

#ifdef Debug
  state.print();
#endif

  return state;
}

/*****************************************************************************/
void KalmanTracking::simulate(TTrack & track, vector<Layer> & layers)
{
#ifdef Debug
  ofstream file("../out/simulation.dat",ios::app);
#endif

  // Vertex
  Layer layer;
  layer.simulated = generate(track);

  // Layers
  layers.push_back(layer);

  bool ok = true;
  for(vector<Material>::const_iterator m = material.begin();
                                       m!= material.end() && ok; m++)
  {
    Layer layer = layers.back();

    if(propagate(layer.simulated, *m))
    {
      multiScatt(layer.simulated, Random);

      {
      State state = layer.simulated;

      double phi = 0.;
      if(state.material.radius > 0)
        phi = state.getRPhi() / state.material.radius - state.getPsi();

//      double theta = state.getTheta() - state.getPsi(); // FIXME
      double theta = state.getTheta() - M_PI/2; // FIXME
// force PIXEL = 0

double p = fabs(1/layer.simulated.getKappa());
double betaGamma = p / layer.simulated.m;
//cerr << " p = " << p  << " m = " << layer.simulated.m << endl;

      track.lhits.push_back(
           clusterGenerator->create(0, betaGamma, theta, phi)); // FIXME

cerr << " generated " << track.lhits.back().filledPixels.size()
               << " " << track.lhits.back().allPixels.size();
state.print();
/*
               << " theta,psi,rphi " << state.getTheta();
               << " " << state.getPsi()
               << " " << state.getRPhi()
               << endl;
*/
      }

      energyLoss(layer.simulated, Random);
//  FIXME    call ClusterGenerator here!!!! FIXME

      layer.measured = makeRechit(layer.simulated);

      layers.push_back(layer);

      // FIXME
      Coord coord;
      convertToCoordinates(layer.simulated, coord);

      track.hits.push_back(coord); // Push FIXME

#ifdef Debug
      file << " " << layer.simulated.xpos()
           << " " << layer.simulated.ypos()
           << " " << layer.simulated.zpos()
           << " " << layer.measured.xpos()
           << " " << layer.measured.ypos()
           << " " << layer.measured.zpos()
           << endl;

      cerr <<  " pt : " << coord.pt
           << "; eta : " << coord.eta << endl;
#endif
    }
    else
      ok = false;
  }

#ifdef Debug
  file << endl << endl;
  file.close();

//  while(getchar() == 0);
#endif

  cerr << " lhits = " << track.lhits.size() << endl;
}

/*****************************************************************************/
vector<Coord> KalmanTracking::reconstruct(vector<Layer> & layers, double & c)
{
  // Fitting
  fit(layers);

  // Smoothing
  smooth(layers);

  double chi2[4] = {0,0,0,0};

  int nLayers = 1;

  for(vector<Layer>::iterator layer = layers.begin() + 1;
                              layer!= layers.end(); layer++)
    if(layer->updated.chi2 < 50) // estimateCut
    {
      chi2[0] += layer->predicted.chi2;
      chi2[1] += layer->updated.chi2;
      chi2[2] += layer->smoothed.chi2;

      nLayers++;
    }

  HepVector gr(nMeas * (layers.size() - 1));
  for(vector<Layer>::iterator layer = layers.begin() + 1;
                              layer!= layers.end(); layer++)
    for(int i1 = 0; i1 < nMeas; i1++)
      gr[nMeas*(layer - (layers.begin()+1)) + i1] = (layer->smoothed.r)[i1];

  int ierr;
  chi2[3] = (gr.T() * gR.inverse(ierr) * gr)[0];

  Coord coord;
  convertToCoordinates(layers[layers.size()-1].smoothed, coord);

  c = chi2[0];
 
  // Copy to coord
  vector<Coord> hits;
  for(vector<Layer>::iterator layer = layers.begin() + 0; // FIXME!!!
                              layer!= layers.end(); layer++)
  {
    Coord coord;

    if(layer == layers.begin())
      convertToCoordinates(layer->updated,  coord);
    else
      convertToCoordinates(layer->smoothed, coord);

     hits.push_back(coord);
  }

  return hits;
}

/*****************************************************************************/
bool KalmanTracking::process(TTrack & simTrack, TTrack & recTrack)
{
  vector<Layer> layers;

  simulate(simTrack, layers);

  if(layers.size() >= 4)
  {
    double c;

    recTrack.hits = reconstruct(layers, c);

    // FIXME
    recTrack.d0 = sqrt(sqr(recTrack.hits[0].x[0])
                     + sqr(recTrack.hits[0].x[1]));
    recTrack.z  =          recTrack.hits[0].x[2];

    recTrack.ndf  = recTrack.hits.size() - 3 - 1; // FIXME
    recTrack.chi2 = c;

    // FIXME
    recTrack.lhits = simTrack.lhits; 

    // FIXME !!!!
    recTrack.eta = simTrack.eta;
    recTrack.pt  = simTrack.pt;

    return true;
  }
  else
    return false;
}

