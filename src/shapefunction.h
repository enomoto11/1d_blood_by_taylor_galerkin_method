#ifndef _SHAPEFUNCTION_H_
#define _SHAPEFUNCTION_H_

#include <iostream>
#include <vector>

using namespace std;

class ShapeFunction1D
{
public:
  static void P2_N(vector<double> &N, const double &g1)
  {
    N.at(0) = 5e-1 * (1e0 - g1);
    N.at(1) = 5e-1 * (1e0 + g1);
  }
  static void P2_dNdr(vector<double> &dNdr, const double &g1)
  {
    dNdr.at(0) = -5e-1;
    dNdr.at(1) = 5e-1;
  }

  static void P2_dNinvdr(vector<double> &dNinvdr, const double &g1)
  {
    dNinvdr.at(0) = 5e-1 * (1e0 - g1) * (1e0 - g1);
    dNinvdr.at(1) = 5e-1 * (1e0 + g1) * (1e0 + g1);
  }
};

#endif