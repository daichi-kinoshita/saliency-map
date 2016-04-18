#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

#define CMD_LINE_NUM 3
#define COLOR_NUM 3
#define COLOR "BGR"

bool saveImage(char*, Mat&);
bool saveImage(char*, vector<vector<double> >&);

bool saveImage(char* file_name, Mat& image)
{
  if (file_name == NULL) {
    return false;
  }
  if (image.empty()) {
    return false;
  }
  
  ofstream pic_stream(file_name, ios::out);
  if (pic_stream == NULL) {
    return false;
  }

  int rows = image.rows;
  int cols = image.cols;
  int channels = image.channels();
  int step = image.step;
  int elem_size = image.elemSize();
  double max;
  minMaxLoc(image, NULL, &max);
  
  //header
  pic_stream << "P2" << endl;
  pic_stream << cols << " " << rows << endl;
  pic_stream << (int)max << endl;
  
  //data
  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
      pic_stream << (int)image.data[y*step + x*elem_size] << " ";
    }
    pic_stream << endl;
  }
  
  return true;
}

bool saveImage(char* file_name, vector<vector<double> >& image)
{
  if (file_name == NULL) {
    return false;
  }
  if (image.empty()) {
    return false;
  }
  
  ofstream pic_stream(file_name, ios::out);
  
  int rows = image.size();
  int cols = image.front().size();
  double max = 0;
  for (int y=0;y<rows;++y) {
    double tmp = *max_element(image[y].begin(), image[y].end());
    if (tmp > max) {
      max = tmp;
    }
  }
  
  //header
  pic_stream << "P2" << endl;
  pic_stream << cols << " " << rows << endl;
  pic_stream << (int)max << endl;
  
  //data
  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
      pic_stream << image[y][x];
    }
    pic_stream << endl;
  }
  
  return true;
}

int main(int argc, char* argv[])
{
  if (argc != CMD_LINE_NUM) {
    cerr << "ERROR: usage " << argv[0] << " <input file name> <output directory>" << endl;
    return -1;
  }
  
  Mat frame;
  int frame_cnt = 1;
  char file_name[64];
  
  /* init input movie */
  String input_path = argv[1];
  double fps = 30;
  int codec = VideoWriter::fourcc('X', 'V', 'I', 'D');
  VideoCapture input(input_path);

  if (!input.isOpened()) {
    cerr << "ERROR: " << argv[1] << " open failed" << endl;
    return -1;
  }

  for(;;) {

    input >> frame;

    if (frame.empty()) {
      break;
    }
    
    vector<Mat> color;
    split(frame, color);

    for (int i=0;i<COLOR_NUM;++i) {
      sprintf(file_name, "%s/%c/frame%05d.pgm", argv[2], COLOR[i], frame_cnt);
	if (!saveImage(file_name, color[i])) {
	  cerr << "ERROR: " << file_name << " open failed" << endl;
	  break;
	}
    }

    frame_cnt++;
  }

  return 0;
}
