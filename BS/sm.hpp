/*****************************
Saliency Map library.

sm
{
class Saliency
class M_Saliency : Saliency

NSS function
}

*******************************/

#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

typedef vector<Mat> vecMat;

namespace sm
{
  /* fixed number */
  const int SCALES    = 9;                   //pyramid scale
  const int CS_NUM    = 6;                   //number of center-surroud map
  const int ORIENTS   = 4;                   //number of orientation

  const int CENTER[2] = {2, 4};              //center = sigma 2 or 3 or 4
  const int DELTA[2]  = {3, 4};              //delta  = sigma +3 or +4;
  
  const int SHIFT_X[4] = {1, 1, 0, -1};
  const int SHIFT_Y[4] = {0, 1, 1,  1};

  //gaussian filter
  const int GAUSSIAN_RADIUS = 8;                   //radius
  const int GAUSSIAN_SIZE = 2*GAUSSIAN_RADIUS + 1; //size
  const double GAUSSIAN_SIGMA = 1.3;               //parameter: sigma
 
  //gabor filter
  const int GABOR_RADIUS = 8;                    //radius
  const int GABOR_SIZE   = 2*GABOR_RADIUS + 1;   //size
  const double GABOR_SIGMA  = 1.3;               //parameter: sigma
  const double GABOR_LAMBDA = GABOR_RADIUS + 1;  //parameter: lambda    

  const int STEP = 8;

  const int SM_SIZE = 0;

  const int ITERATION = 18;
  
  /*** Saliency class **********************************************/
  class Saliency
  {
  public:

    /*feature maps*/
    vecMat I;          //intense
    vecMat R;          //red
    vecMat G;          //green
    vecMat B;          //blue
    vecMat Y;          //yellow
    vecMat F;          //flicker
    vector<Mat_<float> > O[ORIENTS]; //orientation
    vector<Mat_<float> > M[ORIENTS]; //motion

    static const vecMat gabor_f;

    Size sigma_size;

    const vecMat  makeGaussianPyramid(const Mat&);
    vecMat makeCSPyrmd(const vecMat&);
    vecMat makeCSPyrmd(const vecMat&, const vecMat&);
    Mat centerSurround(const Mat&, const Mat&);
    void normalizePeak(Mat&);
    void normalizeRange(Mat&);

  public:
    Saliency();
    Saliency(const Mat&, const Saliency*);
    //Saliency(const Mat&, const vecMat&);

    Mat getSaliencyMap(void);

  };
  /*********************************************************************/
  
  Mat convertTo256(const Mat&);
  Mat HeatMap(const Mat&);
  const vecMat makeGaborFilter(void);
  const Mat makeGaussianFilter(void);
  
  Mat setGazePoint(const Mat&, const Mat&);

  //double NSS(const Mat&, vector<double>&, vector<double>&);

  
};
