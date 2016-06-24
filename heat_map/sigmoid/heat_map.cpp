#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

class Map
{
public :
  Map(char*);
  ~Map(void);

  void setSigmoidGain(double);
  void setSigmoidThreshold_R(double);
  void setSigmoidThreshold_G(double);

  void writeHeatMap(char*, int, int);
  
private :
  int** data;
  int rows, cols;
  int max;

  double a;
  double th_r;
  double th_g;

  int linearToneCurve_R(int);
  int linearToneCurve_G(int);
  int linearToneCurve_B(int);

  int sinToneCurve_R(int);
  int sinToneCurve_G(int);
  int sinToneCurve_B(int);

  int sigmoidToneCurve_R(int, double);
  int sigmoidToneCurve_G(int, double);

  static const int HEAT_MAP_LINEAR  = 1;
  static const int HEAT_MAP_SIN     = 2;
  static const int HEAT_MAP_SIGMOID = 3;

  static const int I_MAX = 255; //max intense

  static const double SIGMOID_GAIN = 0.03; //sigmoid gain
  static const double SIGMOID_TH_R = 0.5;  //sigmoid threshold red
  static const double SIGMOID_TH_G = 0.5;  //sigmoid threshold green
};

Map::Map(char* file_name)
{
  ifstream input(file_name);

  a = SIGMOID_GAIN;
  th_r = SIGMOID_TH_R;
  th_g = SIGMOID_TH_G;

  //read header
  string  magic_num;
  string comment;
  getline(input, magic_num);
  getline(input, comment);
  input >> cols >> rows;
  input >> max;

  //read data
  data = new int*[rows];
  for (int y=0;y<rows;++y) {
    data[y] = new int[cols];
  }

  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
      input >> data[y][x];
    }
  }

}

Map::~Map()
{
  for (int i=0;i<rows;++i) {
    delete[] data[i];
  }
  delete[] data;
}

void Map::setSigmoidGain(double a)
{
  this->a = a;
}

void Map::setSigmoidThreshold_R(double th_r)
{
  this->th_r = th_r;
}

void Map::setSigmoidThreshold_G(double th_g)
{
  this->th_g = th_g;
}

int Map::linearToneCurve_R(int val)
{
  if (val < (I_MAX / 2.0)) {
    return 0;
  }
  else if (val < 192){
    return 4.0*val - 2.0*I_MAX;
  }
  else {
    return I_MAX;
  }
}

int Map::linearToneCurve_G(int val)
{
  if (val < (I_MAX / 4.0)) {
    return 4.0*val;
  }
  else if (val < 192) {
    return I_MAX;
  }
  else {
    return -4.0*val + 4.0*I_MAX;
  }
}

int Map::linearToneCurve_B(int val)
{
  if (val < (I_MAX / 4.0)) {
    return I_MAX;
  }
  else if (val < (I_MAX / 2.0)) {
    return -4.0*val + 2.0*I_MAX;
  }
  else {
    return 0;
  }   
}

int Map::sigmoidToneCurve_R(int val, double a)
{
  return I_MAX * (1.0 / (1.0 + exp(-a * (val - SIGMOID_TH_R*I_MAX)) ));
}

int Map::sigmoidToneCurve_G(int val, double a)
{
  return I_MAX - (I_MAX * (1.0 / (1.0 + exp(-a * (val - SIGMOID_TH_G*I_MAX)) )));
}

//th: threshold a: sigmoid gain
void Map::writeHeatMap(char* file_name, int th=0, int type=HEAT_MAP_LINEAR)
{

  ofstream output(file_name);

  //write header 
  string magic_num = "P3";
  string comment = "# ascii color scale";
  output << magic_num << endl;
  output << comment << endl;
  output << cols << " " << rows << endl;
  output << max << endl;

  //write data
  int r, g, b;
  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
 
      if (data[y][x] < th) {
	output << I_MAX << " " << I_MAX << " " << I_MAX << endl;
      }
      else {

	switch(type) {
	case HEAT_MAP_LINEAR:
	  r = linearToneCurve_R(data[y][x]);
	  g = linearToneCurve_G(data[y][x]);
	  b = linearToneCurve_B(data[y][x]);
	  break;

	case HEAT_MAP_SIGMOID:
	  r = sigmoidToneCurve_R(data[y][x], a);
	  g = sigmoidToneCurve_G(data[y][x], a);
	  b = 0;
	  break;

	default :
	  r = 0;
	  g = 0;
	  b = 0;
	}

	output << r << " " << g << " " << b << endl;
      }
    }
  }
}

// argv[1]: input file name(.pgm)
// argv[2]: output file name(.pnm)
int main(int argc, char *argv[])
{
  Map img(argv[1]);
  
  img.writeHeatMap(argv[2], 48);

  return 0;
}
