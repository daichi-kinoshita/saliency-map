#include <iostream>
#include <fstream>
#include <vector>
#include <opencv2/opencv.hpp>
#include "sm.hpp"

using namespace std;
using namespace cv;
using namespace sm;

static const int CMD_LINE_NUM = 3;

int main(int argc, char* argv[])
{
  if (argc != CMD_LINE_NUM) {
    cerr << "ERROR: usage " << argv[0] << " <input> <output>" << endl;
    return -1;
  }

  /* init */
  String input_path = argv[1];
  String output_path = argv[2];
  double fps = 30.00;
  int codec = cv::VideoWriter::fourcc('X', 'V', 'I', 'D');
  VideoCapture input(input_path);
  VideoWriter output;
  output.open(output_path, codec, fps, Size(640, 480));

  if (!input.isOpened()) {
    cerr << "ERROR: Can't open " << argv[1] << endl;
  }

  vector<Mat> video;
  Saliency *current = NULL;
  Saliency *prev    = NULL;
  
  Mat frame;


  while(true) {
    /* input frame */
    input >> frame;
    
    if (frame.empty()) {
      delete current;
      break;
    }

    /* update frame*/
    current = new Saliency(frame, prev);
    delete prev;
    prev = current;

    Mat tmp = current->getSaliencyMap();

    video.push_back(convertTo256(tmp));
  }


  /* output frame */
  for(int i=0;i<video.size();i++) {
    output << video[i];
   }

 
  return 0;
}
