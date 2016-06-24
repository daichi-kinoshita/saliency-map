#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <opencv2/opencv.hpp>
#include "progress.hpp"

using namespace std;
using namespace cv;

static const int CMD_LINE_NUM = 6;

int main(int argc, char* argv[])
{
  if (argc != CMD_LINE_NUM) {
    cerr << "ERROR: usage " << argv[0] << " <input file> <output file> <fps> <sampling interval> <sampling frames>" << endl;
    return -1;
  }

  Mat frame;

  /* input */
  String input_path = argv[1];
  VideoCapture input(input_path);
  int input_fps = (int)(input.get(CV_CAP_PROP_FPS) + 0.5);
  int frame_num = input.get(CV_CAP_PROP_FRAME_COUNT);

  cout << frame_num << endl;

  /* output */
  String output_path = argv[2];
  int output_fps = atoi(argv[3]);
  int cols = input.get(CV_CAP_PROP_FRAME_WIDTH);
  int rows = input.get(CV_CAP_PROP_FRAME_HEIGHT);
  VideoWriter output;
  output.open(output_path, CV_FOURCC('P','I','M','1'), output_fps, Size(cols , rows));

  int sampling_interval = atoi(argv[4]);
  int sampling_frames = atoi(argv[5]);

  if (sampling_interval == 0) {
     cerr << "ERROR: sampling interval < 1" << endl;
    return -1;
  }

  cout << "[down-sampling]" << endl;
  cout << "input file : " << argv[1] << endl;
  cout << "output file: " << argv[2] << endl;

  int frame_cnt = 0;
  
  for(;;) {

    if (frame_cnt % sampling_interval == 0) {
      for (int i=0;i<sampling_frames;++i) {
	
	input >> frame;
	if (frame.empty()) {
	  break;
	}
	frame_cnt++;
	
	output << frame;
	
	DispProgressBar((double)frame_cnt / frame_num * 100.0);
      }
    }

    input >> frame; 
    if (frame.empty()) {
      break;
    }
    frame_cnt++;
      
  }

  cout << endl;

  return 0;

}
