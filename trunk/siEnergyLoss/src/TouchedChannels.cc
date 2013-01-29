#include "../interface/TouchedChannels.h"

#include <iostream>
#include <cmath>

#define sqr(x) ((x) * (x))

using namespace std;

enum outlier { Vert = 0, Horiz = 1, Endpoint = 2};

/*****************************************************************************/
TouchedChannels::TouchedChannels (const Hit & hit_) : hit(hit_)
{
}

/*****************************************************************************/
TouchedChannels::~TouchedChannels()
{
}

/*****************************************************************************/
bool TouchedChannels::isInside(double point[])
{
  return (point[0] > 0 && point[0] < hit.unit.nrows && 
          point[1] > 0 && point[1] < hit.unit.ncolumns);
}

/****************************************************************************/
void TouchedChannels::addEndpoint(double point[], vector<Crossing>& crossing)
{
  Crossing c;
  for(int k=0; k<2; k++)
    c.pos[k] = point[k];
  c.dir = -1;

  crossing.push_back(c); 
}

/****************************************************************************/
void TouchedChannels::lookForCrossing
  (int dir, double point[][2],
   const vector<double>& divisionLines, vector<Crossing>& crossing)
{
  double dpos[2];
  for(int k=0; k<2; k++)
    dpos[k] = (point[1][k] - point[0][k]);
     
  for(vector<double>::const_iterator div = divisionLines.begin();
                                     div!= divisionLines.end(); div++)      
  {
    double position[2];
      
    position[  dir] = *div;
    position[1-dir] = point[0][1-dir] +
                          dpos[1-dir]*(*div - point[0][dir])/dpos[dir];
      
    if(isInside(position))
    {
      Crossing c;
      for(int k=0; k<2; k++)
        c.pos[k] = position[k];
      c.dir = dir;

      crossing.push_back(c);
    } 
  }
}

/****************************************************************************/
void TouchedChannels::sortCrossings(vector<Crossing>& crossing)
{
  bool change;

  do
  {
    change = false;
   
    for(vector<Crossing>::iterator ic = crossing.begin();
                                   ic!= crossing.end() - 1; ic++)
    if((*ic).pos[dir] > (*(ic+1)).pos[dir])
    {
      swap(*ic, *(ic+1));
      change = true;
    }
  }
  while(change);
}

/****************************************************************************/
void TouchedChannels::matchOrAddPixel
  (double dpos[2], Crossing c[2], vector<Pixel>& pixels)
{
  // Calculate position
  pair<float,float> pix( (c[1].pos[0]  + c[0].pos[0])/2,
                         (c[1].pos[1]  + c[0].pos[1])/2 );

  int x = (int) pix.first;
  int y = (int) pix.second;

  // Calculate partial path
  Pixel::Calc calc;
  for(int k = 0; k < 2; k++)
    calc.dl_dP[k] = 0.;

  calc.l = fabs(c[1].pos[dir] - c[0].pos[dir])/dpos[dir] * hit.length; 

  // Calculate normal vectors
  for(int i=0; i<2; i++)
  {
    int n = c[i].dir;

    if(n >= 0)
    {
      double normal[2];

      normal[    n] = (c[1-i].pos[n] > c[i].pos[n] ? 1 : -1);
      normal[1 - n] = 0;

      double scalar = 0.; 
      for(int k = 0; k < 2; k++)
        scalar += normal[k] * hit.ulambda[k];
      scalar = fabs(scalar);

      for(int k = 0; k < 2; k++)
        calc.dl_dP[k] += normal[k] / scalar *
                         hit.unit.pitch[k] * hit.length / hit.lambda;
    }
  }
 
  calc.isTouched = true;

  bool isFound = false;
  for(vector<Pixel>::iterator pixel = pixels.begin();
                              pixel!= pixels.end(); pixel++)
    if(pixel->meas.x[0] == x &&
       pixel->meas.x[1] == y)
    {
      pixel->calc = calc;
      isFound = true;

      break;
    }

  if(!isFound)
  {
    Pixel pixel;
    pixel.meas.x[0] = x;
    pixel.meas.x[1] = y;
    pixel.meas.y = 0;

    pixel.meas.c[0] = c[0];
    pixel.meas.c[1] = c[1];

    pixel.calc = calc;

    pixels.push_back(pixel);
  }
}

/****************************************************************************/
vector<double> TouchedChannels::getDivisionLines(int dir, double point[2][2])
{
  int e[2];
  e[0] = int(floor(min(point[0][dir],point[1][dir])) + 1);
  e[1] = int(floor(max(point[0][dir],point[1][dir]))    );

  int high;
  if(dir == 0) high = hit.unit.nrows;
          else high = hit.unit.ncolumns;

  for(int i=0; i<2; i++)
  {
    if(e[i] < 0)    e[i] = 0;
    if(e[i] > high) e[i] = high;
  }

  vector<double> divisionLines;

  for(int i=e[0]; i<=e[1]; i++)
  {
    double j = i;
    if(i == 0   ) j += 1e-3; // FIXME
    if(i == high) j -= 1e-3; // FIXME

    divisionLines.push_back(j); 
  }

  return divisionLines;
}

/****************************************************************************/
// Outlier start
void TouchedChannels::selectEndpoints(int line[2], int point[2], int index[2])
{
  if( (line[0] <= 1 && line[1] <= 1) ||
      (line[0] >= 2 && line[1] >= 2) )
  { // crosses pixel outside (vvhh, hhvv)
    if( (point[0] <= 1 && point[1] >= 3) ||
        (point[0] == 2 || point[1] == 2) )
    { index[0] = -1; index[1] = -2; }
    else
    {
      if(point[0] <= 1)
      { index[0] = -2; index[1] = 1; }
      else
      { index[0] = -1; index[1] = 0; }
    }
  }
  else
  { // crosses pixel inside (vhvh, hvhv, vhhv, hvvh)
    if(point[0] <= 1)
    { index[0] = -1; index[1] = 1; } // take higher point
    else
    { index[0] = -2; index[1] = 0; } // take lower point
  }
}

/*****************************************************************************/
void TouchedChannels::fillLinesPoints(const vector<Crossing> & cross,
 int  line[2], int  point[2],
 int iline[4], int ipoint[2])
{
  int il=0,ik=0, k=0,l=0, j=0;

  for(vector<Crossing>::const_iterator c = cross.begin();
                                       c!= cross.end(); c++)
  {
    // line : sequential number of vertical lines among all lines
    if(c->typ == Vert) line[il++] = k;
 
    if(c->typ != Endpoint)
    {
      // iline: kth line is cross number j
      iline[k++]  = j;
    }
    else
    {
      //  point: sequential number of points among lines
       point[ik++] = k;

      // ipoint: lth point is cross number j
      ipoint[l++]  = j;
    }
 
    j++;
  }
}

/*****************************************************************************/
void TouchedChannels::calculateOutlierDistance
  (double endpoint[2][2], Pixel& pixel) 
{
  double dist = 1e+99;
  vector<double> vdist(2, 0.);

  // Look at both endpoints, find closest part (vertex or edge) of the pixel
  for(int i = 0; i < 2; i++)
  {
    vector<double> a(2, 0.);

    // Check vertices
    for(int dx = 0; dx < 2; dx++)
    for(int dy = 0; dy < 2; dy++)
    {
      vector<double> a(2, 0.);
      a[0] = ((pixel.meas.x[0] + dx) - endpoint[i][0]) * hit.unit.pitch[0];
      a[1] = ((pixel.meas.x[1] + dy) - endpoint[i][1]) * hit.unit.pitch[1];

      double d = sqrt(sqr(a[0]) + sqr(a[1])); 

      if(d < dist) 
      { dist = d; vdist = a; }
    } 

    // Check horizontal edges
    for(int dy = 0; dy < 2; dy++)
    {
      for(int dx = 0; dx < 2; dx++) // projection on side?
      if(endpoint[i][0] > pixel.meas.x[0] &&
         endpoint[i][0] < pixel.meas.x[0] + dx)
      {
        vector<double> a(2, 0.);
        a[1] = ((pixel.meas.x[1] + dy) - endpoint[i][1]) * hit.unit.pitch[1];

        double d = fabs(a[1]);

        if(d < dist) 
        { dist = d; vdist = a; }
      }
    }

    // Check vertical edges
    for(int dx = 0; dx < 2; dx++)
    {
      for(int dy = 0; dy < 2; dy++) // projection on side?
      if(endpoint[i][1] > pixel.meas.x[1] &&
         endpoint[i][1] < pixel.meas.x[1] + dy)
      {
        vector<double> a(2, 0.);
        a[0] = ((pixel.meas.x[0] + dx) - endpoint[i][0]) * hit.unit.pitch[0];

        double d = fabs(a[0]);

        if(d < dist) 
        { dist = d; vdist = a; }
      }
    }
  }

  // Look for projections of pixel vertices
  for(int dx = 0; dx < 2; dx++)
  for(int dy = 0; dy < 2; dy++)
  {
    double ca[2], ba[2];
  
    ca[0] = pixel.meas.x[0] + dx - endpoint[0][0];
    ca[1] = pixel.meas.x[1] + dy - endpoint[0][1];

    for(int k = 0; k < 2; k++)
      ba[k] = endpoint[1][k] - endpoint[0][k];

    double lam = (ca[0]*ba[0] + ca[1]*ba[1]) /
                 (ba[0]*ba[0] + ba[1]*ba[1]);

    if(lam > 0 && lam < 1) // is inside
    {
      double p[2];
      for(int k = 0; k < 2; k++)
        p[k] = endpoint[0][k] + lam*(endpoint[1][k] - endpoint[0][k]);

      vector<double> a(2, 0.);
      a[0] = ((pixel.meas.x[0] + dx) - p[0]) * hit.unit.pitch[0];
      a[1] = ((pixel.meas.x[1] + dy) - p[1]) * hit.unit.pitch[1];

      double d = sqrt(sqr(a[0]) + sqr(a[1]));

      if(d < dist)
      { dist = d; vdist = a; }
    }
  }

  pixel.calc.l = - dist;

  for(int k = 0; k < 2; k++)
    pixel.calc.dl_dP[k] = vdist[k] / dist * hit.unit.pitch[k];
}

/****************************************************************************/
vector<Pixel> TouchedChannels::findChannels(double point[2][2])
{
  // Collect crossing points
  vector<Crossing> crossing;

  // Starting points
  for(int i=0; i<2; i++)
    if(isInside(point[i]))
      addEndpoint(point[i], crossing);  

  // Crossing points in x and y directions 
  for(int direction = 0; direction < 2; direction++)
  {
    vector<double> divisionLines = getDivisionLines(direction, point);
    lookForCrossing(direction, point, divisionLines, crossing);
  }

  // Choose dominant direction
  double dpos[2];
  for(int k=0; k<2; k++)
    dpos[k] = (point[1][k] - point[0][k]);
  
  if(fabs(dpos[0]) > fabs(dpos[1])) dir = 0;
                               else dir = 1;

  // Copy
  vector<Pixel> allPixels;
  for(vector<Pixel>::const_iterator pixel = hit.filledPixels.begin();
                                    pixel!= hit.filledPixels.end(); pixel++)
    allPixels.push_back(*pixel);

  // Clean allPixels
  Pixel::Calc calc;
  calc.isTouched = false;
  calc.l = 0.;
  for(int k = 0; k < 2; k++)
    calc.dl_dP[k] = 0.;


  for(vector<Pixel>::iterator pixel = allPixels.begin();
                              pixel!= allPixels.end(); pixel++)
    pixel->calc = calc;

  // Sort crossings along dir
  if(crossing.size() > 0)
  {
    sortCrossings(crossing);

    // Process touched channels
    for(vector<Crossing>::iterator ic = crossing.begin();
                                   ic!= crossing.end() - 1;
                                   ic++)
      matchOrAddPixel(dpos, &(*ic), allPixels); // fill, add missing touched
  }

  // Take outliers
  for(vector<Pixel>::iterator pixel = allPixels.begin();
                              pixel!= allPixels.end(); pixel++)
    if(pixel->calc.isTouched == false)
      calculateOutlierDistance(point, *pixel);

  return allPixels;
}
