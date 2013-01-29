#ifndef _TouchedChannels_h_
#define _TouchedChannels_h_

#include "../../DataFormats/interface/Hit.h"
#include <vector>

class TouchedChannels
{
  public:
    TouchedChannels(const Hit & hit_);
    ~TouchedChannels();

    std::vector<Pixel> findChannels(double endpoint[2][2]);

  private:
    bool isInside(double point[]);
    void addEndpoint(double point[], std::vector<Crossing>& crossing);
    void lookForCrossing
      (int dir, double point[][2], const std::vector<double>& divisionLines,
       std::vector<Crossing>& crossing);
    void sortCrossings(std::vector<Crossing>& crossing);
    void matchOrAddPixel(double dpos[2], Crossing c[2],
                         std::vector<Pixel>& pixels);
    std::vector<double> getDivisionLines (int dir, double point[2][2]);

    // outliers
    void selectEndpoints(int line[2], int point[2], int index[2]);
    void fillLinesPoints(const std::vector<Crossing> & cross,
                        int  line[2], int  point[2],
                        int iline[4], int ipoint[2]);

    void calculateOutlierDistance(double endpoint[2][2], Pixel& pixel);

    const Hit hit;

    int dir; // dominant direction
};

#endif
