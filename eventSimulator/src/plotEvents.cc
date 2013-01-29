#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "../../DataFormats/interface/TBunchCrossing.h"
#include "../../DataFormats/interface/TVertex.h"
#include "../../DataFormats/interface/TTrack.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace CLHEP;

/*****************************************************************************/
void printHelix(const double * p1_, const double * p2_, const double * v2_,
                ofstream& outFile, int charge)
{
  Hep3Vector p1(p1_[0], p1_[1], p1_[2]);
  Hep3Vector p2(p2_[0], p2_[1], p2_[2]);
  Hep3Vector v2(v2_[0], v2_[1], v2_[2]);

  Hep3Vector dp = p2 - p1;
  Hep3Vector n2(-v2.y(),v2.x(),0.);
  n2 = n2.unit();

  double r = -0.5 * (dp.x()*dp.x() + dp.y()*dp.y()) /
                    (dp.x()*n2.x() + dp.y()*n2.y());
  Hep3Vector c = p2 + r * n2;

  double dphi = sqrt(2 * 0.1 / fabs(r)); // allowed deflection: 0.1 cm

  double phi = acos(( (p1-c).x()*(p2-c).x() +
                      (p1-c).y()*(p2-c).y() )/(r*r));

  if(dp.x()*v2.x() + dp.y()*v2.y() < 0) phi = 2*M_PI - phi;

  int nstep = (int)(phi/dphi) + 1;

  if(nstep > 1)
  {
    dphi = phi / nstep;
    double dz = (p2 - p1).z() / nstep;

    Hep3Vector P0 = p1;
    Hep3Vector P1;

    charge = ((p1 - c).x() * (p2 - c).y() - (p1 - c).y() * (p2 - c).x() > 0 ?
              -1 : 1);
    if(dp.x()*v2.x() + dp.y()*v2.y() < 0) charge = -charge;

    outFile << ", Line[{{"<<P0.x()<<","<<P0.y()<<","<<P0.z()<<"}" ;

    for(int i = 0; i < nstep; i++)
    {
      double a = -charge * (i+1)*dphi;
      double z = p1.z() + (i+1)*dz;

      double x = c.x() + cos(a)*(p1 - c).x() - sin(a)*(p1 - c).y();
      double y = c.y() + sin(a)*(p1 - c).x() + cos(a)*(p1 - c).y();

      P1 = Hep3Vector(x,y,z);

      outFile << ", {"<<P1.x()<<","<<P1.y()<<","<<P1.z()<<"}";
      P0 = P1;
    }
    outFile << "}]" << endl;
  }
  else
  {
    Hep3Vector P0 = p1;
    Hep3Vector P1 = p2;

    outFile << ", Line[{{"<<P0.x()<<","<<P0.y()<<","<<P0.z()<<"}"
                  << ", {"<<P1.x()<<","<<P1.y()<<","<<P1.z()<<"}}]" << endl;
  }
}

/*****************************************************************************/
void plotDetUnit(const Coord & hit, ofstream & file)
{
  double phi = atan2(hit.x[1],hit.x[0]);
  double r   = sqrt (hit.x[0]*hit.x[0] + hit.x[1]*hit.x[1]);
  double z   = hit.x[2];

  double dz  = 10;

  int iphi = floor(phi/ (2*M_PI/32));
  int iz   = floor(z / dz);

  double phi0 =  iphi    * 2*M_PI/32;
  double phi1 = (iphi+1) * 2*M_PI/32;

  double z0   =  iz    * dz;
  double z1   = (iz+1) * dz;

  Hep3Vector p00(r*cos(phi0),r*sin(phi0),z0);
  Hep3Vector p01(r*cos(phi0),r*sin(phi0),z1);
  Hep3Vector p11(r*cos(phi1),r*sin(phi1),z1);
  Hep3Vector p10(r*cos(phi1),r*sin(phi1),z0);

  file << ", If[sd, {RGBColor[0.4,0.4,0.4]"
       << ", Line[{{"<<p00.x()<<","<<p00.y()<<","<<p00.z()<<"}, "
                <<"{"<<p01.x()<<","<<p01.y()<<","<<p01.z()<<"}, "
                <<"{"<<p11.x()<<","<<p11.y()<<","<<p11.z()<<"}, "
                <<"{"<<p10.x()<<","<<p10.y()<<","<<p10.z()<<"}, "
                <<"{"<<p00.x()<<","<<p00.y()<<","<<p00.z()<<"}}]}]"
       << endl;
}

/*****************************************************************************/
void plotPixel(ofstream & fileC, const Pixel & pixel, const double * pitch,
               double s0, double s1)
{
  vector<double> x(2);
  vector<double> y(2);
  vector<double> z(2);

  x[0] = (pixel.meas.x[0]     - s0) * pitch[0];
  x[1] = (pixel.meas.x[0] + 1 - s0) * pitch[0];

  y[0] = (pixel.meas.x[1]     - s1) * pitch[1];
  y[1] = (pixel.meas.x[1] + 1 - s1) * pitch[1];

  z[0] = -pitch[2]/2;
  z[1] =  pitch[2]/2;

  for(int i = 0; i < 2; i++)
  fileC<< ", Line[{{"<<x[0]<<","<<y[0]<<","<<z[i]<<"},"
               << "{"<<x[0]<<","<<y[1]<<","<<z[i]<<"},"
               << "{"<<x[1]<<","<<y[1]<<","<<z[i]<<"},"
               << "{"<<x[1]<<","<<y[0]<<","<<z[i]<<"},"
               << "{"<<x[0]<<","<<y[0]<<","<<z[i]<<"}}]";

  for(int i = 0; i < 2; i++)
  fileC<< ", Line[{{"<<x[i]<<","<<y[0]<<","<<z[0]<<"},"
               << "{"<<x[i]<<","<<y[0]<<","<<z[1]<<"},"
               << "{"<<x[i]<<","<<y[1]<<","<<z[1]<<"},"
               << "{"<<x[i]<<","<<y[1]<<","<<z[0]<<"},"
               << "{"<<x[i]<<","<<y[0]<<","<<z[0]<<"}}]";
/*
  for(int i = 0; i < 2; i++)
  fileC<< ", Line[{{"<<x[0]<<","<<y[i]<<","<<z[0]<<"},"
               << "{"<<x[0]<<","<<y[i]<<","<<z[1]<<"},"
               << "{"<<x[1]<<","<<y[i]<<","<<z[1]<<"},"
               << "{"<<x[1]<<","<<y[i]<<","<<z[0]<<"},"
               << "{"<<x[0]<<","<<y[i]<<","<<z[0]<<"}}]";
*/
}

/*****************************************************************************/
void plotBunchCrossing(TBunchCrossing * bunx, ofstream & file, ofstream & fileC)
{
  // start graphics 
  file << "Graphics3D[";
  fileC<< "Graphics3D[";

  // start physics
  file << "{AbsolutePointSize[5]";
  fileC<< "{AbsolutePointSize[5]";

// FIXME
double dx,dy,dz;

  // simVertices
  file << ", If[sv, {RGBColor[0.5,0.5,0.5]";
  for(vector<TVertex>::const_iterator vertex = bunx->simVertices.begin();
                                      vertex!= bunx->simVertices.end();
                                      vertex++)
  {
    int iv = int(vertex - bunx->simVertices.begin()) + 1;

    file << ", If[(iv==0 || iv=="<<iv<<"), {AbsolutePointSize[6]";

    // Vertex
    file << ", RGBColor[0.5,0.5,0.5], Point[{" << 0
                                        << "," << 0
                                        << "," << vertex->z << "}]" << endl;
    file << ", Text[StyleForm[\"" << iv << "\", FontFamily->\"Helvetica\"";

    ostringstream info;
        info << setprecision(3) << "&z0=-"
             << "&d0=-"
             << "&chi2=-"
             << "&ndf=-"
             << "&eta=-"
             << "&pt=-";

    ostringstream url;
//    url << "file:///home/sikler/aida/bud-aida/eventSimulator/view/passClick.html";
    url << "http://localhost:8000/view/passClick.html";

    ostringstream vth;
    vth << "?iv=" << iv
        << "&it=-1"
        << "&ih=-1";// FIXME

    file << ", URL->\""
         << url.str() << vth.str() << info.str() << ",target=passClick\"]";

    file << ", {" << 0
           << "," << 0
           << "," << vertex->z << "}" << ", {2,1}]" << endl;
    file << "}]";
  }
  file << "}]";

  // simTracks
  file << ", If[st, {RGBColor[0.5,0.5,0.5]";
  fileC<< ", If[(1==1), {RGBColor[0,0,0]";
  for(vector<TVertex>::const_iterator vertex = bunx->simVertices.begin();
                                      vertex!= bunx->simVertices.end();
                                      vertex++)
  {
    int iv = int(vertex - bunx->simVertices.begin()) + 1;

    file << ", If[(iv==0 || iv=="<<iv<<"), {AbsolutePointSize[6]";

    // Vertex
/*
    file << ", RGBColor[0.5,0.5,0.5], Point[{" << 0
                                        << "," << 0
                                        << "," << vertex->z << "}]" << endl;
    file << ", Text[StyleForm[\"" << iv << "\", FontFamily->\"Helvetica\"";

    ostringstream info;
        info << setprecision(3) << "&z0=-"     
             << "&d0=-"    
             << "&chi2=-"   
             << "&ndf=-"    
             << "&eta=-"    
             << "&pt=-";

    ostringstream url;
    url << "file:///home/sikler/aida/bud-aida/eventSimulator/view/passClick.html";

    ostringstream vth;
    vth << "?iv=" << iv
        << "&it=-1"
        << "&ih=-1";// FIXME

    file << ", URL->\""
         << url.str() << vth.str() << info.str() << ",target=passClick\"]";    

    file << ", {" << 0
           << "," << 0
           << "," << vertex->z << "}" << ", {2,1}]" << endl;
*/
    //

    fileC<< ", If[(iv=="<<iv<<"), {AbsolutePointSize[6]";
    for(vector<TTrack>::const_iterator track = vertex->tracks.begin();
                                       track!= vertex->tracks.end(); track++)
    {
      int it = int(track - vertex->tracks.begin()) + 1;

      file << ", If[(it==0 || it=="<<it<<"), {AbsolutePointSize[6]";
      fileC<< ", If[(it=="<<it<<"), {AbsolutePointSize[6]";
      for(vector<Coord>::const_iterator hit = track->hits.begin();     
                                        hit!= track->hits.end(); hit++)
      {
        int ih = int(hit - track->hits.begin()) + 1;

        plotDetUnit(*hit, file);

        file << ", RGBColor[0.6,0.6,1.0]" << endl; // Point

        file << ", Point[{" << hit->x[0]
                     << "," << hit->x[1]
                     << "," << hit->x[2] << "}]" << endl;

        ostringstream info;
        info << setprecision(3) << "&z0="     << track->z
             << "&d0="     << track->d0
             << "&chi2="   << track->chi2
             << "&ndf="    << track->ndf
             << "&eta="    << hit->eta
             << "&pt="     << hit->pt;

        ostringstream url;
//        url << "file:///home/sikler/aida/bud-aida/eventSimulator/view/passClick.html";
        url << "http://localhost:8000/view/passClick.html";

        ostringstream vth;
        vth << "?iv=" << iv
            << "&it=" << it
            << "&ih=" << ih; // FIXME

        file << ", Text[StyleForm[\"" << it << "\", FontFamily->\"Helvetica\"";

        file << ", URL->\""
             << url.str() << vth.str() << info.str() << ",target=passClick\"]";

        file << ", {" << hit->x[0]
               << "," << hit->x[1]
               << "," << hit->x[2] << "}" << ", {2,1}]" << endl;

        file << ", RGBColor[1.0,0.5,0.5]" << endl; // Line

        // vertex to first
        if(hit == track->hits.begin())
        {
          double pos[3] = {0.,0.,vertex->z};
          printHelix(pos, (hit  )->x, (hit  )->p, file, track->charge);
        }

        // if not last, draw helix
        if(hit != track->hits.end() - 1)
          printHelix(hit->x, (hit+1)->x, (hit+1)->p, file, track->charge);
      }
      file << "}]";

      // lhits
      for(vector<Hit>::const_iterator hit = track->lhits.begin();     
                                      hit!= track->lhits.end(); hit++)
      {
        int ih = int(hit - track->lhits.begin()) + 1;

        fileC<< ", If[(ih=="<<ih<<"), {AbsolutePointSize[6]";
/*
        fileC<< ", Text[StyleForm[" 
             << hit->allPixels.size() << "], {0,0,1}, {0,0}]" << endl;
*/

        double s0=0,s1=0,n=0;
        for(vector<Pixel>::const_iterator pixel = hit->allPixels.begin();
                                          pixel!= hit->allPixels.end(); pixel++)
        {
          s0 += pixel->meas.x[0] + 0.5;
          s1 += pixel->meas.x[1] + 0.5;
          n++;
        }

        s0/=n; s1/=n;
//s0 = 0;
//s1 = 0;

        dx = hit->unit.pitch[0];
        dy = hit->unit.pitch[1];
        dz = hit->unit.pitch[2];

        for(vector<Pixel>::const_iterator pixel = hit->allPixels.begin();
                                          pixel!= hit->allPixels.end(); pixel++)
        {
          fileC<< ", RGBColor[0.0,0.0,0.0]" << endl;
          fileC<< ", Text[StyleForm[\""
          // Put energy deposit per pixel here
                << int(pixel->meas.y) << "\", FontFamily->\"Helvetica\"], {" 
                   << (pixel->meas.x[0] + 0.5 - s0)*dx
            << "," << (pixel->meas.x[1] + 0.5 - s1)*dy 
            << "," << "0}, {0,0}]" << endl;

          fileC<< ", RGBColor[0.4,0.4,0.4]" << endl;
          plotPixel(fileC, *pixel, hit->unit.pitch, s0, s1);
        }

        fileC<< ", RGBColor[0.5,0,0]" << endl;
        fileC<< ", Line[{{" << (hit->pos_orig[0] - s0 - hit->dpos[0]/2)*dx
                      <<"," << (hit->pos_orig[1] - s1 - hit->dpos[1]/2)*dy
                      <<"," <<  dz/2 << "}"
                   <<", {"  << (hit->pos_orig[0] - s0 + hit->dpos[0]/2)*dx
                      <<"," << (hit->pos_orig[1] - s1 + hit->dpos[1]/2)*dy
                       <<","<< -dz/2 << "}}]";

        fileC<< "}]";
      }
      fileC<< "}]";
    }
    file << "}]";
    fileC<< "}]";
  }
  file << "}]";
  fileC<< "}]";

  // recVertices
  file << ", If[rv, {AbsolutePointSize[6]";
  for(vector<TVertex>::const_iterator vertex = bunx->recVertices.begin();
                                      vertex!= bunx->recVertices.end();
                                      vertex++)
  {
    int iv = int(vertex - bunx->recVertices.begin()) + 1;

    file << ", If[(iv==0 || iv=="<<iv<<"), {AbsolutePointSize[6]";

    // Vertex
    file << ", RGBColor[0.8,0.2,0.2], Point[{" << 0       
                                        << "," << 0
                                        << "," << vertex->z << "}]" << endl;
    file << ", Text[StyleForm[\"" << iv << "\", FontFamily->\"Helvetica\"";

    ostringstream info;
        info << setprecision(3) << "&z0=-"
             << "&d0=-"
             << "&chi2=-"
             << "&ndf=-"
             << "&eta=-"
             << "&pt=-";

    ostringstream url;
//    url << "file:///home/sikler/aida/bud-aida/eventSimulator/view/passClick.html";
    url << "http://localhost:8000/view/passClick.html";

    ostringstream vth;
    vth << "?iv=" << iv
        << "&it=-1"
        << "&ih=-1";// FIXME

    file << ", URL->\""
         << url.str() << vth.str() << info.str() << ",target=passClick\"]";

    file << ", {" << 0
           << "," << 0
           << "," << vertex->z << "}" << ", {-1,1}]" << endl;
    file << "}]";
  }
  file << "}]";

  // recTracks
  file << ", If[rt, {AbsolutePointSize[6]";
  for(vector<TVertex>::const_iterator vertex = bunx->recVertices.begin();
                                      vertex!= bunx->recVertices.end();
                                      vertex++)
  {
    int iv = int(vertex - bunx->recVertices.begin()) + 1;

    file << ", If[(iv==0 || iv=="<<iv<<"), {AbsolutePointSize[6]";

/*
    // Vertex
    file << ", RGBColor[0.8,0.2,0.2], Point[{" << 0 
                                        << "," << 0
                                        << "," << vertex->z << "}]" << endl;
    file << ", Text[StyleForm[\"" << iv << "\", FontFamily->\"Helvetica\"";

    ostringstream info;
        info << setprecision(3) << "&z0=-"
             << "&d0=-"
             << "&chi2=-"
             << "&ndf=-"
             << "&eta=-"
             << "&pt=-";

    ostringstream url;
    url << "file:///home/sikler/aida/bud-aida/eventSimulator/view/passClick.html";

    ostringstream vth;
    vth << "?iv=" << iv
        << "&it=-1" 
        << "&ih=-1";// FIXME

    file << ", URL->\""
         << url.str() << vth.str() << info.str() << ",target=passClick\"]"; 

    file << ", {" << 0
           << "," << 0
           << "," << vertex->z << "}" << ", {-1,1}]" << endl;
*/

    for(vector<TTrack>::const_iterator track = vertex->tracks.begin();
                                       track!= vertex->tracks.end(); track++)
    {
      int it = int(track - vertex->tracks.begin()) + 1;

      file << ", If[(it==0 || it=="<<it<<"), {AbsolutePointSize[6]";
      for(vector<Coord>::const_iterator hit = track->hits.begin();
                                        hit!= track->hits.end(); hit++)
      {
        if(hit != track->hits.begin())
        {
        int ih = int(hit - track->hits.begin()) + 1;

        plotDetUnit(*hit, file);

        file << ", RGBColor[0.2,0.2,0.6]" << endl; // Point
        
        file << ", Point[{" << hit->x[0]
                     << "," << hit->x[1]
                     << "," << hit->x[2] << "}]" << endl;

        ostringstream info;
        info << setprecision(3) << "&z0="     << track->z
             << "&d0="     << track->d0
             << "&chi2="   << track->chi2
             << "&ndf="    << track->ndf
             << "&eta="    << hit->eta
             << "&pt="     << hit->pt;

        ostringstream url;
//        url << "file:///home/sikler/aida/bud-aida/eventSimulator/view/passClick.html";
        url << "http://localhost:8000/view/passClick.html";

        ostringstream vth;
        vth << "?iv=" << iv
            << "&it=" << it
            << "&ih=" << ih-1; // FIXME

        file << ", Text[StyleForm[\"" << it << "\", FontFamily->\"Helvetica\"";

	file << ", URL->\""
             << url.str() << vth.str() << info.str() << ",target=passClick\"]"; 

        file << ", {" << hit->x[0] 
               << "," << hit->x[1]
               << "," << hit->x[2] << "}" << ", {-1,1}]" << endl;
        }

        file << ", RGBColor[0.5,0,0]" << endl; // Line

        // vertex to first
        if(hit == track->hits.begin())
        {
          double pos[3] = {0.,0.,vertex->z};
          printHelix(pos, (hit  )->x, (hit  )->p, file, track->charge);
        }

        // if not last, draw helix
        if(hit != track->hits.end() - 1)
          printHelix(hit->x, (hit+1)->x, (hit+1)->p, file, track->charge);
      }
      file << "}]";
    }
    file << "}]";
  }
  file << "}]";

  // region (tracker + ecal) // FIXME
  int mz = 300;
  
  // beam line
  file << ", RGBColor[0.0,1.0,0.0]";
  
  for(int z = -mz; z < mz; z += mz/30)
    file << ", Line[{{0,0,"<<z<<"}, {0,0,"<<z+mz/30<<"}}]" << endl; 

  // stop physics
  file << "}";

  // options
  file << ", PlotRange->All";
  file << ", Boxed->False";
  file << "]";

  // stop cluster
  fileC<< "}";

  // options
//  fileC<<", ViewPoint -> {0, 0, 1000}";
//  fileC<<", ViewVertical -> {0, 1, 0}";
//  fileC<< ", PlotRange->All";
  fileC<< ", PlotRange -> {{"<<-2*dx<<","<<2*dx<<"}, {"
                             <<-4*dy<<","<<4*dy<<"}, {"
                             <<-dz/2<<","<<dz/2<<"}}";
//  fileC<< ", Boxed->True";
  fileC<< ", Boxed->False";
  fileC<< "]";
}

/*****************************************************************************/
int main()
{
  // Open file
  TFile * fileData = new TFile("../out/bunchCrossings.root", "read");

  // Get tree
  TTree * tree = (TTree *) fileData->Get("trackTree");
  TBranch * branch = tree->GetBranch("bunchCrossing");

  // Set bunch crossing
  TBunchCrossing * bunx = new TBunchCrossing();
  branch->SetAddress(&bunx);

  cerr << " plotting bunch crossing..";

  branch->GetEntry(0);

  // Open output
  ofstream file("../out/bunchCrossing.m");
  ofstream fileC("../out/clusters.m");

  file << fixed << setprecision(3);

  plotBunchCrossing(bunx, file, fileC);

  file.close();
  fileC.close();

  cerr << " [done]" << endl;

  bunx->Clear();
}

