#include <iostream>
#include <fstream>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace std;
  
class SM
{
private :
/* fixed number */
  static const int SCALES    = 9;                   //pyramid scale
  static const int CS_NUM    = 6;                   //number of center-surroud map
  static const int ORIENTS   = 4;                   //number of orientation

  static const int CENTER[2] = {2, 4};              //center = sigma 2 or 3 or 4
  static const int DELTA[2]  = {3, 4};              //delta  = sigma +3 or +4;
  
  static const int SHIFT_X[4] = {1, 1, 0, -1};
  static const int SHIFT_Y[4] = {0, 1, 1,  1};

  //gaussian filter
  static const int GAUSSIAN_RADIUS = 8;                   //radius
  static const int GAUSSIAN_SIZE = 2*GAUSSIAN_RADIUS + 1; //size
  static const double GAUSSIAN_SIGMA = 1.3;               //parameter: sigma
 
  //gabor filter
  static const int GABOR_RADIUS = 8;                    //radius
  static const int GABOR_SIZE   = 2*GABOR_RADIUS + 1;   //size
  static const double GABOR_SIGMA  = 1.3;               //parameter: sigma
  static const double GABOR_LAMBDA = GABOR_RADIUS + 1;  //parameter: lambda    

  static const int STEP = 8;

  static const int SM_SIZE = 0;

  static const int ITERATION = 18;


public :
  
  static Mat CalcSaliencyMap(Mat);
  static vector<Mat> Color(Mat);
  static vector<Mat> Intense();
  static vector<Mat> Orientation(Mat);
  static vector<Mat> Fliker(Mat, Mat);
  static vector<Mat> Motion(Mat, Mat);

};

int main()
{

  return 0;
}
