#include "../interface/KalmanTracking.h"

#include "../../DataFormats/interface/TVertex.h"
#include "../../DataFormats/interface/TTrack.h"
#include "../../DataFormats/interface/TChipId.h"
#include "../../DataFormats/interface/RandomGenerator.h"

#include "../../siEnergyLoss/interface/Levels.h"
#include "../../siEnergyLoss/interface/ClusterGenerator.h"

#include "../../DataFormats/interface/Particles.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>

using namespace std;

//#define Debug

#define sqr(x) ((x) * (x))

enum Flag { Random, Mean, Sigma };

const bool useOrigin = true;
const bool useLogLikelihood = false;
const double maxChi = 100.;

/*****************************************************************************/
KalmanTracking::KalmanTracking
  (vector<TLayer> & mat, float gainHalfWidth_, bool simu)
{
  // silicon
  X0  = 21.82 / rho;

  if(simu)
  {
    gainHalfWidth = gainHalfWidth_;

    readMaterial(material);
    mat = material;

    clusterGenerator = new ClusterGenerator();

    fileChi.open("../out/chi.dat");
  }
  else
  {
    material = mat;
    B = material[0].B; // FIXME

    cerr << "  read " << material.size() << " layers" << endl;
    cerr << "  magnetic field B = " << B << endl;

    fileChi.open("../out/chi_.dat");
  }

  theRandom = new RandomGenerator();
}

/*****************************************************************************/
KalmanTracking::~KalmanTracking()
{
  delete theRandom;
}

/*****************************************************************************/
void KalmanTracking::readMaterial(vector<TLayer> & material)
{
  cerr << " reading material description..";

  ifstream file;
  file.open("../data/cms.dat");

  ofstream fileGain("../out/gains.orig");

  // Read magnetic field
  file >> B;

  int ilayer = 0;

  bool stop;

  do
  {
    string s;
    file >> s;

    stop = (s.compare("stop") == 0);

    if(!stop)
    {
      TLayer m;
      int d;

      m.isPixel = (s.compare("p") == 0); // p (pixel) or s (strip)

      m.ilayer = ilayer++;

      m.nrows    = 100;
      m.ncolumns = 100;

      float f;
      file >> m.radius;
      file >> m.halfLength;

      for(int k = 0; k < 3; k++) 
        file >> m.pitch[k]; // cm
 
      file >> d; m.sigma_rphi =  d * 1e-4;       // um   -> cm
      file >> d; m.sigma_z    =  d * 1e-4;       // um   -> cm
      file >> f; m.thickness  = (f * 1e-2) * X0; // x/X0 -> cm
      file >> d; m.diffSigma  =  d * 1e-4;       // um   -> xm

      file >> f; m.noise     = f*1e-3; // keV -> MeV
      file >> f; m.threshold = f*1e-3; // keV -> MeV
      file >> f; m.overflow  = f*1e-3; // keV -> MeV

      m.coupling = 0.;
      if(m.isPixel) file >> s; // pixel should not have coupling
               else file >> m.coupling;

      file >> m.dz;                   // dz
      file >> d; m.dphi = (2*M_PI)/d; // dphi

      // Find chip, set gains
      for(double z   = -m.halfLength + m.dz/2; z < m.halfLength; z += m.dz)
      for(double phi = m.dphi/2; phi < 2*M_PI; phi += m.dphi)
      {
        int iz   = floor(z / m.dz);
        int iphi = floor(phi/ m.dphi);

        TChipId chipId(m.ilayer,iz,iphi);

        m.gain[chipId.code] =
           1 + gainHalfWidth * (2*theRandom->getFlatRandom()-1);

        fileGain << " " << chipId.getiLayer()
                 << " " << chipId.getiZ()
                 << " " << chipId.getiPhi() 
                 << " " << m.gain[chipId.code]
                 << endl;
      }

      m.B = B; // FIXME

      material.push_back(m);
    }
  }
  while(!stop);

  file.close();

  fileGain.close();

  cerr << " [done]" << endl;
}

/*****************************************************************************/
void KalmanTracking::convertToParameters(Point & coord, State & state)
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
void KalmanTracking::convertToCoordinates(State & state, Point & coord)
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

  coord.isPixel = state.material.isPixel;
}

/*****************************************************************************/
double KalmanTracking::multiScatt(State & state, int flag)
{
  Point coord;
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
    alpha += sigma * theRandom->getGaussRandom();

  u[0] =   coord.pz*cos(alpha) + coord.pt*sin(alpha); // pz'
  u[1] = - coord.pz*sin(alpha) + coord.pt*cos(alpha); // pt'

  coord.p[2] = u[0];

  coord.p[0] *= u[1] / coord.pt;
  coord.p[1] *= u[1] / coord.pt;

  // Rotate (px,py)
  alpha = mean;
  if(flag == Random)
    alpha += sigma * theRandom->getGaussRandom();

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
  Point coord;
  convertToCoordinates(state, coord); // just to get pt

  const double me = mass[elec]; 

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
    dE += sigma * theRandom->getGaussRandom(); 

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

  float factor = 1.;

  // sometimes add another component with factor 2 higher RMS
  if(theRandom->getFlatRandom() < 1e-2) factor = 2.;

  rechit.setRPhi(state.getRPhi() +
    state.material.sigma_rphi * factor * theRandom->getGaussRandom());

  if(state.material.sigma_z < 5.) // smaller than 5 cm 
    rechit.setZ(state.getZ() +
       state.material.sigma_z * factor * theRandom->getGaussRandom());
  else
    rechit.setZ(state.getZ());

  return rechit;
}

/*****************************************************************************/
// Propagate to barrel, r1
bool KalmanTracking::propagate(State & state, const TLayer & material)
{
  double r1 = material.radius;

  Point coord;
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

  if(fabs(coord.x[2]) < material.halfLength)
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

  Point coord;

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
  vertex.m = mass[pion]; // recMass
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
TMatrixD KalmanTracking::calculateF(State & rechit, TLayer & material)
{
  TMatrixD F(nPars,nPars); F.Zero();

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
TMatrixD KalmanTracking::outerProduct(const TVectorD & a, const TVectorD & b)
{
  int na = a.GetNoElements();
  int nb = b.GetNoElements();

  TMatrixD M(na,nb);

  for(int i = 0; i < na; i++)
  for(int j = 0; j < nb; j++)
    M[i][j] = a[i] * b[j];

  return M;
}

/*****************************************************************************/
TMatrixD KalmanTracking::calculateQ(State & state, const TMatrixD & F)
{
  TMatrixD Q(nPars,nPars);
  TVectorD k(nPars), t(nPars), p(nPars);

  for(int i = 0; i < nPars; i ++)
  {
    k[i] = F[i][0]; // 0th column
    t[i] = F[i][1]; // 1st column
    p[i] = F[i][2]; // 2nd column
  }

  double sigma_kappa = energyLoss(state, Sigma);
  double sigma_theta = multiScatt(state, Sigma);
  double sigma_psi   = sigma_theta;

  Q = outerProduct(k, k) * sqr(sigma_kappa) +
      outerProduct(t, t) * sqr(sigma_theta) +
      outerProduct(p, p) * sqr(sigma_psi  );

  return Q;
}

/*****************************************************************************/
double KalmanTracking::getChi2(const State & state)
{
  TMatrixD In(TMatrixD::kInverted, state.R);
  double chi2 = (state.r * (In * state.r));

  if(useLogLikelihood)
    chi2 += log(2 * M_PI * state.R.Determinant());

  return chi2;
}

/*****************************************************************************/
void KalmanTracking::fit(vector<Layer> & layers)
{
  // Constant measurement and identity matrices
  TMatrixD H(nMeas,nPars);
    H.Zero(); H[0][3] = 1.; H[1][4] = 1.;
  TMatrixD Ip(nPars,nPars);
   Ip.Zero(); for(int i = 0; i < nPars; i++) Ip[i][i] = 1.;
  TMatrixD Im(nMeas,nMeas);
   Im.Zero(); for(int i = 0; i < nMeas; i++) Im[i][i] = 1.;

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
      TMatrixD C(nPars,nPars); C.Zero();
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
      TLayer material = layer->simulated.material;
  
      // Start with previous updated state
      State state = (layer-1)->updated;

      // Transient F, Q, V, G
      TMatrixD F = calculateF(state, material);

      TMatrixD Q(nPars,nPars); Q.Zero();
      if(layer > layers.begin() + 1)
      Q = calculateQ(state, F);

      TMatrixD V(nMeas,nMeas); V.Zero();
      V[0][0] = sqr(layer->simulated.material.sigma_rphi);
      V[1][1] = sqr(layer->simulated.material.sigma_z   );

      TMatrixD G(nMeas,nMeas); G.Zero();
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
      TMatrixD FT(TMatrixD::kTransposed, F);
      state.C = F * state.C * FT + Q; 

      // Measured
      TVectorD m = H * layer->measured.x;

      // Residuals of predictions
      state.r = m - H * state.x;

      // Covariance of predicted residuals
      TMatrixD HT(TMatrixD::kTransposed, H);
      state.R = V + H * (state.C * HT);

      // Chi2 increment
      state.chi2 = getChi2(state);

      if(state.chi2 > maxChi * maxChi)
      {
#ifdef Debug
        cerr << " (b) chi too big = " << sqrt(state.chi2) << endl;
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
      TMatrixD K =
        state.C * HT * (V + H * state.C * HT).Invert();
  
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
        cerr << " (c) chi too big = " << sqrt(state.chi2) << endl;
#endif

        layers.erase(layer,layers.end());
        break;
      }
  
#ifdef Debug
      Point coord;
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
}

/*****************************************************************************/
void KalmanTracking::set(TMatrixD & a, TMatrixD & b)
{
  a.ResizeTo(b);
  a = b;
}

/*****************************************************************************/
void KalmanTracking::smooth(vector<Layer> & layers)
{
  TMatrixD H(nMeas,nPars); H.Zero(); H[0][3] = 1.; H[1][4] = 1.;

  // Global covariance
  TMatrixD R(nMeas*(layers.size()-1),
             nMeas*(layers.size()-1));
  vector<TMatrixD> C(layers.size()-1);

#ifdef Debug
  ofstream file("../out/smoothing.dat",ios::app);
#endif

  for(vector<Layer>::iterator layer = layers.end()   - 1;
                              layer > layers.begin(); layer--)
  {
    TMatrixD V(nMeas,nMeas); V.Zero();
    V[0][0] = sqr(layer->simulated.material.sigma_rphi);
    V[1][1] = sqr(layer->simulated.material.sigma_z   );

    if(layer == layers.end() - 1)
    {
      layer->smoothed = layer->updated;

      int i = layer - (layers.begin()+1);

      set(C[i], layer->smoothed.C);

      TMatrixD HT(TMatrixD::kTransposed, H);
      TMatrixD r(V - H * C[i] * HT);

      for(int i1 = 0; i1 < nMeas; i1++)
      for(int j1 = 0; j1 < nMeas; j1++)
        R[nMeas*i + i1][nMeas*i + j1] = r[i1][j1];
    }
    else
    {
      // Smoother gain matrix
      TMatrixD In(TMatrixD::kInverted, (layer+1)->predicted.C);

      TMatrixD FT(TMatrixD::kTransposed, layer->updated.F);
      TMatrixD A = layer->updated.C * FT * In;

      State state = layer->updated;
  
      state.x += A * ( (layer+1)->smoothed.x - (layer+1)->predicted.x );

      TMatrixD AT(TMatrixD::kTransposed, A);
      state.C += A * ( (layer+1)->smoothed.C - (layer+1)->predicted.C ) * AT;

      // Measured
      TVectorD m = H * layer->measured.x;
  
      state.r = m - H * state.x;
  
      TMatrixD HT(TMatrixD::kTransposed, H);
      state.R = V - H * state.C * HT;

      {
        int i = layer - (layers.begin()+1);
        int n = layers.size() - 1;

        set(C[i], state.C); 

        for(int j = i+1; j < n; j++)
          C[j] = A * C[j];

        for(int j = i; j < n; j++)
        {
          TMatrixD r(nMeas,nMeas); r.Zero();

          TMatrixD HT(TMatrixD::kTransposed, H);
          if(i == j) r = V - (H * C[j] * HT);
                else r = r - (H * C[j] * HT);

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
#ifdef Debug
        cerr << " (a) chi too big = " << sqrt(state.chi2) << endl;
#endif

        layers.erase(layer,layers.end());
        break;
      }

      layer->smoothed = state;

#ifdef Debug
    Point coord;
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
  
  // copy with care
  int n = nMeas*(layers.size() - 1);
  TMatrixD Rp(n,n);

  for(int i = 0; i < n; i++)
  for(int j = 0; j < n; j++)
   Rp[i][j] = R[i][j];

//  gR = R;
  set(gR, Rp);

#ifdef Debug
  file << endl << endl;
  file.close();
#endif
}

/*****************************************************************************/
State KalmanTracking::generate(TTrack & track)
{
  Point coord;

  // Starts from origo
  coord.x[0] = 0.;
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

  // Default
  state.m = mass[pion];

  // Ovewrite if there is a pdgId
  for(int i = 1; i < nParts; i++)
    if(abs(track.pdgId) == pdg[i]) state.m = mass[i];

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

  for(vector<TLayer>::iterator m = material.begin(); // iterator??
                               m!= material.end() && ok; m++)
  {
    Layer layer = layers.back();

    if(propagate(layer.simulated, *m))
    {
      multiScatt(layer.simulated, Random);

      State state = layer.simulated;

      double phi = 0.;
      if(state.material.radius > 0)
        phi = state.getRPhi() / state.material.radius - state.getPsi();

      double theta = state.getTheta() - M_PI/2;

      double p = fabs(1/layer.simulated.getKappa());
      double betaGamma = p / layer.simulated.m;

      // FIXME
      Point coord;
      convertToCoordinates(layer.simulated, coord);

      int ilayer, iz, iphi;
      {
        ilayer = int(m - material.begin());

        double phi = atan2(coord.x[1], coord.x[0]);
        double z   =       coord.x[2];

        iz   = floor(z / m->dz);
        iphi = floor( (phi < 0 ? phi + 2*M_PI : phi)/ m->dphi);
      }
      TChipId chipId(ilayer, iz, iphi);

      track.hits.push_back(
           clusterGenerator->create(betaGamma, theta, phi, chipId,
                                    &(*m), int(m - material.begin()) ) );

      energyLoss(layer.simulated, Random);

      layer.measured = makeRechit(layer.simulated);

      layers.push_back(layer);

      track.points.push_back(coord); // Push

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
}

/*****************************************************************************/
vector<Point> KalmanTracking::reconstruct(vector<Layer> & layers, double & c, int pdgId)
{
  // Fitting
  fit(layers);

  // Smoothing
  smooth(layers);

  vector<double> chi2(4,0.);

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

  TVectorD gr(nMeas * (layers.size() - 1));
  for(vector<Layer>::iterator layer = layers.begin() + 1;
                              layer!= layers.end(); layer++)
    for(int i1 = 0; i1 < nMeas; i1++)
      gr[nMeas*(layer - (layers.begin()+1)) + i1] = (layer->smoothed.r)[i1];

  TMatrixD In(TMatrixD::kInverted, gR);

  chi2[3] = (gr * (In * gr));

  Point coord;
  convertToCoordinates(layers[layers.size()-1].smoothed, coord);

  if(layers.size() > 1)
  {
    Point coord;
    convertToCoordinates(layers[layers.size()-1].smoothed, coord);

    fileChi << " " << layers[1].smoothed.getKappa()
            << " " << sqrt(chi2[0])
            << " " << sqrt(chi2[1])
            << " " << sqrt(chi2[2])
            << " " << sqrt(chi2[3])
            << " " << nLayers - 1
            << " " << coord.p_
            << " " << pdgId
            << endl;
  }

  c = chi2[0];
 
  // Copy to Point
  vector<Point> points;
  for(vector<Layer>::iterator layer = layers.begin(); // FIXME!!!
                              layer!= layers.end(); layer++)
  {
    Point coord;

    if(layer == layers.begin())
      convertToCoordinates(layer->updated,  coord);
    else
      convertToCoordinates(layer->smoothed, coord);

     points.push_back(coord);
  }

  return points;
}

/*****************************************************************************/
bool KalmanTracking::process(TTrack & simTrack, TTrack & recTrack)
{
  vector<Layer> layers;

  simulate(simTrack, layers);

  if(layers.size() >= 4)
  {
    double c;

    // FIXME all below should be in the reco step, simple copy for now
    recTrack.points = reconstruct(layers, c, simTrack.pdgId);

    // Take updated (eta,pt) at beam line (radius = 0)
    Point coord;
    convertToCoordinates(layers[0].updated, coord);

    recTrack.eta = coord.eta;
    recTrack.pt  = coord.pt;

    recTrack.d0 = sqrt(sqr(recTrack.points[0].x[0])
                     + sqr(recTrack.points[0].x[1]));
    recTrack.z  =          recTrack.points[0].x[2];

//    recTrack.sigma_z = 0.1;
//    recTrack.sigma_z = sqrt(layers[1].updated.C[4][4]); // FIXME
    recTrack.sigma_z = 100e-4/recTrack.pt * pow(cosh(recTrack.eta),1.5);

    recTrack.ndf  = recTrack.points.size() - 3 - 1;
    recTrack.chi2 = c;

    recTrack.hits = simTrack.hits; 

    return true;
  }
  else
    return false;
}

/*****************************************************************************/
/*
bool KalmanTracking::simulate(TTrack & simTrack)
{
  vector<Layer> layers;

  simulate(simTrack, layers);


  simTrack.layers = layers; 

  return true;
}
*/

/*****************************************************************************/
/*
bool KalmanTracking::reconstruct(TTrack & simTrack, TTrack & recTrack)
{
  vector<Layer> layers = simTrack.layers;

  if(layers.size() >= 4)
  {
    double c;

    recTrack.hits = reconstruct(layers, c, simTrack.pdgId);

    recTrack.d0 = sqrt(sqr(recTrack.hits[0].x[0])
                     + sqr(recTrack.hits[0].x[1]));
    recTrack.z  =          recTrack.hits[0].x[2];

    recTrack.ndf  = recTrack.hits.size() - 3 - 1;
    recTrack.chi2 = c;

    recTrack.hits = simTrack.hits;

    recTrack.eta = simTrack.eta;
    recTrack.pt  = simTrack.pt;

    return true;
  }
  else
    return false;
}
*/
