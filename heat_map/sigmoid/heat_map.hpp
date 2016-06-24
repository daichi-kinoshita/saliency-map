#ifndef __HEAT_MAP_HPP
#defien __HEAT_MAP_HPP

#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

class Map
{
public :
  Map(char*);
  ~Map(void);

  void writeHeatMap(char*, int, double, );
  double**;
  
private :
  int** data;
  int rows, cols;
  int max;

  int linearToneCurve_R(int);
  int linearToneCurve_G(int);
  int linearToneCurve_B(int);

  int sigmoidToneCurve_R(int, double);
  int sigmoidToneCurve_G(int, double);

  static const int HEAT_MAP_LINEAR = 1;
  static const int HEAT_MAP_SIGMOID = 2;

  static const int I_MAX = 255; //max intense
  static const double SIGMOID_GAIN = 0.03; //sigmoid gain
  static const double SIGMOID_TH_R = 0.5;  //sigmoid threshold red
  static const double SIGMOID_TH_G = 0.5;  //sigmoid threshold green
};


#endif
