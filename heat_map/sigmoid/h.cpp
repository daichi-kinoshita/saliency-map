#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <cstdlib>

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

int main(int argc, char* argv[])
{
  ofstream output(argv[1]);
  
  string magic_num = "P3";
  string comment = "# ascii color scale";
  int rows = 255;
  int cols = 50;
  int max = 255;

  output << magic_num << endl;
  output << comment << endl;
  output << cols << " " << rows << endl;
  output << max << endl;

  int r, g, b;
  for (int y=rows;y>=0;--y) {
    for (int x=0;x<cols;++x) {
      if (y < 48) {
	output << 255 << " " << 255 << " " << 255 << endl;
      }
      else {
	r = toneCurve_R(y);
	g = toneCurve_G(y);
	b = toneCurve_B(y);

	output << r << " " << g << " " << b << endl;
      }
    }
  }

  return 0;
}
