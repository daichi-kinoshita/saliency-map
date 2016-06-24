#ifndef __IO_IMAGE_HPP
#define __IO_IMAGE_HPP

#include <fstream>
#include <cstdlib>
#include "map.hpp"

using namespace std;

namespace io
{

  class Image
  {
  public :
    static bool GetPGMImage(char*, Map2D&);
    static bool GetPNMImage(char*, vector<Map2D>&);
    static bool SavePGMImage(char*, Map2D);
    static bool SavePNMImage(char*, vector<Map2D>);

  private :
    Image(){};
    
  };

  bool io::Image::GetPGMImage(char* pic_name, Map2D& image)
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

    image.Resize(rows, cols);

    //data (0-255)
    double data;
    for (int y=0;y<rows;++y) {
      for (int x=0;x<cols;++x) {
	pic_stream >> data;
	image.data[y][x] = data;
      }
    }
  
    return true;
  }

  bool io::Image::GetPNMImage(char* pic_name, vector<Map2D>& image)
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

    for (int c=0;c<3;++c) {
      image.push_back(Map2D(rows, cols));
    }

    //data (0-255)
    double data;

    for (int y=0;y<rows;++y) {
      for (int x=0;x<cols;++x) {

	for (int c=0;c<3;++c) {
	  pic_stream >> data;	
	  image[c].data[y][x] = data;
	}

      }
    }
  
    return true;
  }

  bool io::Image::SavePGMImage(char* file_name, Map2D image)
  {
    if (file_name == NULL) {
      return false;
    }
    if (image.Empty()) {
      return false;
    }
  
    ofstream pic_stream(file_name, ios::out);
  
    image.ClampZero(0);
    image = Map2D::Abs(image);
    image.NormalizeRange(255);

    int rows = image.Rows();
    int cols = image.Cols();
    int max = (int)image.Max();

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

  bool io::Image::SavePNMImage(char* file_name, vector<Map2D> image)
  {
    if (file_name == NULL) {
      return false;
    }

    for (int c=0;c<image.size();++c) {
      if (image[c].Empty()) {
	return false;
      }
    }
  
    ofstream pic_stream(file_name, ios::out);
   
    int rows = image[0].Rows();
    int cols = image[0].Cols();

    //header
    pic_stream << "P3" << endl;
    pic_stream << cols << " " << rows << endl;
    pic_stream << 255 << endl;
  
    //data (0-255)

    for (int c=0;c<image.size();++c) {
      image[c].ClampZero(0);
      image[c].NormalizeRange(255);
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

  class Movie
  {

  };

}

#endif
