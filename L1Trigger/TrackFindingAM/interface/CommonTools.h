#ifndef _COMMONTOOLS_H_
#define _COMMONTOOLS_H_

#include <cmath>
#include <iostream>

using namespace std;

enum HW_SIGN_TYPE {UNSIGNED, SIGNED};


class CommonTools{

 public:
  static bool hardwareSimulation;
  static double binning(double fNumber, int nMSBpowOfTwo, int nBits, HW_SIGN_TYPE signType);
  static void binCordic(double X, double Y, double &result_R, double &result_PHI);
};

#endif
