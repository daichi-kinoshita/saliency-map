#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

using namespace std;

static const int I_MAX = 255;

int toneCurve_R(int val)
{
   if (val < 128) {
    return 0;
  }
  else if (val < 192){
    return -2.0*I_MAX + 4.0*val;
  }
  else {
    return I_MAX;
  }
}

int toneCurve_G(int val)
{
  if (val < 64) {
    return 4.0*val;
  }
  else if (val < 192) {
    return I_MAX;
  }
  else {
    return 4.0*(I_MAX - val);
  }
}

int toneCurve_B(int val)
{
   if (val < 64) {
    return I_MAX;
  }
  else if (val < 128) {
    return 2*I_MAX - 4.0*val;
  }
  else {
    return 0;
  }   
}

int main(int argc, char *argv[])
{
  ofstream output(argv[1]);

  //write data
  int r, g, b;
  for (int y=0;y<=255;++y) {
    r = toneCurve_R(y);
    g = toneCurve_G(y);
    b = toneCurve_B(y);
    
    output << y << " " << r << " " << g << " " << b << endl;
  }

  return 0;
}
