#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

static const int CMD_LINE_NUM = 5;

void dispProgressBar(double);

void dispProgressBar(double percent)
{
  int max_bar_len = 25;
  string progress_bar;

  int bar_len = percent / (100 / max_bar_len);

  for (int i=0;i<max_bar_len;++i) {
    if (i < bar_len) {
      progress_bar += '=';
    }
    else if (i == bar_len) {
      progress_bar += '>';
    }
    else {
      progress_bar += ' ';
    }
  }

  printf("\r");
  cerr << "[" << progress_bar << "] " << (int)(percent + 0.5) << "%";
}

int main(int argc, char* argv[])
{
  if (argc != CMD_LINE_NUM) {
    cerr << "ERROR: usage " << argv[0] << " <input file> <output file> <fps> <sampling interval>" << endl;
    return -1;
  }

  Mat frame;

  /* input */
  String input_path = argv[1];
  VideoCapture input(input_path);
  int input_fps = (int)(input.get(CV_CAP_PROP_FPS) + 0.5);
  int frame_num = input.get(CV_CAP_PROP_FRAME_COUNT);

  /* output */
  String output_path = argv[2];
  int output_fps = atoi(argv[3]);
  int cols = input.get(CV_CAP_PROP_FRAME_WIDTH);
  int rows = input.get(CV_CAP_PROP_FRAME_HEIGHT);
  int codec = VideoWriter::fourcc('P', 'I', 'M', '1');
  VideoWriter output;
  output.open(output_path, codec, output_fps, Size(cols , rows));

  int sampling_interval = atoi(argv[4]);

  if (sampling_interval == 0) {
     cerr << "ERROR: sampling interval < 1" << endl;
    return -1;
  }

  cout << "[dawn-sampling]" << endl;
  cout << "input file : " << argv[1] << endl;
  cout << "output file: " << argv[2] << endl;

  int frame_cnt = 0;
  
  for(;;) {
      input >> frame;
      
      if (frame.empty()) {
	break;
      }
      frame_cnt++;

      if (frame_cnt % sampling_interval == 0) {
	dispProgressBar((double)frame_cnt / frame_num * 100.0);
	output << frame;
      }

  }

  cout << endl;

  return 0;

}
