#ifndef _CubicSpline_h_
#define _CubicSpline_h_

class CubicSpline
{
  public:
    CubicSpline();
    ~CubicSpline();
    void prepare
      (double x[], double y[], int n, double yp1, double ypn, double y2[]);
    void interpolate
      (double xa[], double ya[], double y2a[], int n, double x,
       double *y0, double *y1, double *y2);
  private:
};

#endif
