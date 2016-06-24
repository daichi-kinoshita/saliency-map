#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

using namespace std;

class Map
{
public :
  Map(char*);
  ~Map(void);

  void writeHeatMap(char*, int, double);

private :
  int** data;
  int rows, cols;
  int max;
  
  void normalizeMap(void);

  int toneCurve_R(int);
  int toneCurve_G(int);
  int toneCurve_B(int);

  static const int I_MAX = 255;
};

Map::Map(char* file_name)
{
  ifstream input(file_name);

  //read header
  string magic_num;
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

int Map::toneCurve_R(int val)
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

int Map::toneCurve_G(int val)
{
  if (val < 64) {
    return 4.0*val;
  }
  else if (val < 192) {
    return I_MAX;
  }
  else {
    return 4.0*I_MAX - 4.0*val;
  }
}

int Map::toneCurve_B(int val)
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

//th: threshold a: sigmoid gain
void Map::writeHeatMap(char* file_name, int th=48, double a=0.05)
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
	output << 255 << " " << 255 << " " << 255 << endl;
      }
      else {
	r = toneCurve_R(data[y][x]);
	g = toneCurve_G(data[y][x]);
	b = toneCurve_B(data[y][x]);

	output << r << " " << g << " " << b << endl;
      }
    }
  }
}

int main(int argc, char *argv[])
{
  Map img(argv[1]);
  
  img.writeHeatMap(argv[2], 30);

  return 0;
}
