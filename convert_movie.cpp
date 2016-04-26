#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <opencv2/opencv.hpp>
#include "convert.hpp"

using namespace std;
using namespace cv;

#define CMD_LINE_NUM 3

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

    char color_name[3] = {'R', 'G', 'B'};

    for (int i=0;i<3;++i) {
      sprintf(file_name, "input/%s/%c/frame%05d.pgm", argv[2], color_name[i], frame_cnt);
	if (!saveImage(file_name, color[i])) {
	  cerr << "ERROR: " << file_name << " open failed" << endl;
	  break;
	}
    }

    frame_cnt++;
  }

  return 0;
}
