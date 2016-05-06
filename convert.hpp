#ifndef __CONVERT_HPP
#define __CONVERT_HPP

#include <fstream>
#include "map.hpp"

using namespace std;

bool getImage(char*, Map2D&);
bool saveImage(char*, Map2D);
bool saveImage(char*, vector<Map2D>);

bool getImage(char* pic_name, Map2D& image)
{
  ifstream pic_stream(pic_name);
  if (pic_stream == NULL) {
    return false;
  }

  //header
  string magic_number;
  int rows;
  int cols;
  int max;
  pic_stream >> magic_number;
  pic_stream >> cols >> rows;
  pic_stream >> max;

  image.set(rows, cols);

  //data (0.0-1.0)
  double data;
  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
      pic_stream >> data;
      image.data[y][x] = data;
    }
  }
  
  return true;
}


bool saveImage(char* file_name, Map2D image)
{
  if (file_name == NULL) {
    return false;
  }
  if (image.empty()) {
    return false;
  }
  
  ofstream pic_stream(file_name, ios::out);
   
  image = Map2D::abs(image); 
  image.normalizeRange(255);

  int rows = image.rows;
  int cols = image.cols;
  int max = (int)image.max();

  //header
  pic_stream << "P2" << endl;
  pic_stream << cols << " " << rows << endl;
  pic_stream << max << endl;
  
  //data (0-255)
  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
      pic_stream << (int)image.data[y][x] << " ";
    }
    pic_stream << endl;
  }

  return true;
}

/* save image (.pnm) */
bool saveImage(char* file_name, vector<Map2D> image)
{
  if (file_name == NULL) {
    return false;
  }

  for (int c=0;c<image.size();++c) {
    if (image[c].empty()) {
      return false;
    }
  }
  
  ofstream pic_stream(file_name, ios::out);
   
  int rows = image[0].rows;
  int cols = image[0].cols;

  //header
  pic_stream << "P3" << endl;
  pic_stream << cols << " " << rows << endl;
  pic_stream << 255 << endl;
  
  //data (0-255)

  for (int c=0;c<image.size();++c) {
    image[c].normalizeRange(255);
  }

  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
      for (int c=0;c<image.size();++c) {
	pic_stream << (int)image[c].data[y][x] << " ";
      }
    }
    pic_stream << endl;
  }

  return true;
}

#endif
