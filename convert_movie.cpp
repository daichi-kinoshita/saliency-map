#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

#define CMD_LINE_NUM 3

/* save image (.pnm) */
bool saveImage(char* file_name, Mat image)
{
  if (file_name == NULL) {
    return false;
  }
  if (image.empty()) {
    return false;
  }
  
  ofstream pic_stream(file_name, ios::out);
   
  int rows = image.rows;
  int cols = image.cols;
  int step = image.step;
  int elem_size = image.elemSize();
  int channels = image.channels();
  double max;

  minMaxLoc(image, NULL, &max);

  //header
  pic_stream << "P3" << endl;
  pic_stream << cols << " " << rows << endl;
  pic_stream << (int)max << endl;

  //data (0-255)
  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
      for (int c=0;c<channels;++c) {
	pic_stream << (int)image.data[y * step + x * elem_size + c] << " ";
      }
      pic_stream << endl;
    }
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


    sprintf(file_name, "%s/RGB/frame%05d.pnm", argv[2], frame_cnt);
    
    imwrite(file_name, frame);
    /*if (!saveImage(file_name, frame)) {
       cerr << "ERROR: " << file_name << " open failed" << endl;
	break;
	}*/

    frame_cnt++;
  }

  return 0;
}
