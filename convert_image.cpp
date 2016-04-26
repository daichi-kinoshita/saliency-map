#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <opencv2/opencv.hpp>
#include "convert.hpp"

using namespace std;
using namespace cv;

#define CMD_LINE_NUM 4

int main(int argc, char* argv[])
{
  if (argc != CMD_LINE_NUM) {
    cerr << "ERROR: usage " << argv[0] << " <input directory> <output file name> <fps>" << endl;
    return -1;
  }

  Mat frame;
  int frame_cnt = 1;
  char file_name[64];

  /* input image rows,cols */
  sprintf(file_name, "%s/frame%05d.pgm", argv[1], frame_cnt);
  frame = imread(file_name);
  if (frame.empty()) {
    cerr << "ERROR: " << file_name << " open failed" << endl;
    return -1;
  }
  int rows = frame.rows;
  int cols = frame.cols;

  /* init output movie */
  String output_path = argv[2];
  double fps = atoi(argv[3]);
  int codec = VideoWriter::fourcc('X', 'V', 'I', 'D');
  VideoWriter output;
  output.open(output_path, codec, fps, Size(cols, rows));

  for(;;) {

    sprintf(file_name, "%s/frame%05d.pgm", argv[1], frame_cnt);
    frame = imread(file_name);
    
    if (frame.empty()) {
      break;
    }
    
    output << frame;
    
    frame_cnt++;
  }
  
  return 0;
}
