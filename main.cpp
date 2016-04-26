#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include "convert.hpp"
#include "map.hpp"

using namespace std;

static const int CMD_LINE_NUM = 2;


/* Filter class */
class Filter
{
public :
  Filter(int);

  Map2D filter2D(Map2D);

protected :
  Map2D kernel;
  double sigma;   // filter sigma
  int k_size;     // filter kernel size
};

Filter::Filter(int k_size)
{
  sigma = 1.0;

  if (k_size % 2 == 0) {
    k_size++;
  }

  this->k_size = k_size;
  kernel = Map2D(k_size, k_size);
}

Map2D Filter::filter2D(Map2D image)
{
  Map2D _image = image;

  int rows = image.rows;
  int cols = image.cols;
  int h_k_size = (int)(k_size / 2);
  double sum;
  double value;

  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {

      sum = 0.0;
      for (int dy=-h_k_size;dy<=h_k_size;++dy) {
	for (int dx=-h_k_size;dx<=h_k_size;++dx) {
	  int yy = y + dy;
	  int xx = x + dx;

	  if (yy < 0) {
	    yy = 0;
	  }else if (yy >= rows) {
	    yy = rows - 1;
	  }

	  if (xx < 0) {
	    xx = 0;
	  }else if (xx >= cols) {
	    xx = cols - 1;
	  }

	  sum += image.data[yy][xx] * kernel.data[dy+h_k_size][dx+h_k_size];
	}
      }

      _image.data[y][x] = sum;
    }
  }

  _image = Map2D::abs(_image);

  return _image;
}

/* Gaussian class */
class Gaussian : public Filter
{
public :
  Gaussian(int);
  Gaussian(int, double);
  Gaussian(const Gaussian&);

private :
  void makeGaussianKernel(void);
};

Gaussian::Gaussian(int k_size) : Filter(k_size)
{
  sigma = 0.3 * (k_size / 2.0 - 1) + 0.5;

  makeGaussianKernel();
}

Gaussian::Gaussian(int k_size, double sigma) : Filter(k_size)
{
  this->sigma = sigma;

  makeGaussianKernel();
}

Gaussian::Gaussian(const Gaussian& obj) : Filter(obj){}

void Gaussian::makeGaussianKernel(void)
{
  int h_k_size = (int)(k_size / 2);

  double tmp1 = 1.0 / (2.0 * sigma * sigma);
  double tmp2 = tmp1 / M_PI;
  double value;

  for (int y=-h_k_size;y<=h_k_size;++y) {
    for (int x=-h_k_size;x<=h_k_size;++x) {
      value = tmp2 * exp(-1.0 * (x*x + y*y) * tmp1);
      kernel.data[y+h_k_size][x+h_k_size] = value;
    }
  } 

}

/* GaborFilter class */
class Gabor : public Filter
{
public :
  Gabor(int, double, double, double, double);
  //size, theta, sigma, gamma, psi, band_width
  Gabor(const Gabor&);

private :
  void makeGaborKernel(double, double, double, double);
};

Gabor::Gabor(int k_size, double theta, double gamma=1.0, double psi=0.0, double band_w=1.0) : Filter(k_size)
{
  double lambda = k_size / 2.0;

  sigma = sqrt(log(2.0) / 2.0) / M_PI * lambda;
  sigma *= (pow(2, band_w) + 1.0 ) / (pow(2, band_w) - 1.0);

  makeGaborKernel(theta, lambda, gamma, psi);
}

Gabor::Gabor(const Gabor& obj) : Filter(obj){}

void Gabor::makeGaborKernel(double theta, double lambda, double gamma, double psi)
{ 
  int h_k_size = (int)(k_size / 2.0);

  double tmp1 = -1.0 / (2.0 * sigma * sigma);
  double tmp2 = ((2.0 * M_PI) / lambda);
  double value;

  double x_theta;
  double y_theta;

  for (int y=-h_k_size;y<=h_k_size;++y) {
    for (int x=-h_k_size;x<=h_k_size;++x) {
      x_theta = x*cos(theta) + y*sin(theta);
      y_theta = -x*sin(theta) + y*cos(theta);

      value = exp((x_theta*x_theta + gamma*gamma*y_theta*y_theta) * tmp1); 
      value *= cos(tmp2 * x_theta + psi);
      kernel.data[y+h_k_size][x+h_k_size] = value;
    }
  }

}

/* GaborBank class */
class GaborBank
{
public :
  vector<Gabor> kernel;
  int size;
  
  GaborBank(int, int);
  GaborBank(const GaborBank&);
};

GaborBank::GaborBank(int k_size, int orientation)
{
  size = orientation;

  for (int i=0;i<orientation;++i) {
    Gabor tmp(k_size, (double)i / orientation * M_PI);
    kernel.push_back(tmp);
  }
}

GaborBank::GaborBank(const GaborBank& obj)
{
  this->size = obj.size;
  this->kernel = obj.kernel;
}

/* calculate SaliencyMap class */
class SM
{
public :
  static Map2D calcSaliency(Map2D, Map2D, Map2D);
  //static Map2D calcSaliency(Map2D, Map2D, Map2D, Map2D, Map2D, Map2D);

private :
  SM(){};

  static const int SCALE;               // num of pyramid
  static const int ORIENTATIONS;        // num of orientation
  static const int GAUSSIAN_K_SIZE;     // gaussian filter kernel size
  static const int GABOR_K_SIZE;        // gabor filter kernel size
  static Gaussian gaussian_kernel;
  static GaborBank gabor_bank;

  //static vector<Map2D> calcFlicker(Map2D, Map2D);
  //static vector<Map2D> calcMotion(vector<Map2D>, vector<Map2D>);

  static vector<Map2D> makePyramidIntense(const vector<Map2D>&, const vector<Map2D>&, const vector<Map2D>&);
  
  static vector<Map2D> makePyramidRed(vector<Map2D>, vector<Map2D>, vector<Map2D>);
  static vector<Map2D> makePyramidGreen(vector<Map2D>, vector<Map2D>, vector<Map2D>);
  static vector<Map2D> makePyramidBlue(vector<Map2D>, vector<Map2D>, vector<Map2D>);
  static vector<Map2D> makePyramidYellow(vector<Map2D>, vector<Map2D>, vector<Map2D>);

  static vector<Map2D> makePyramidOrientation(const vector<Map2D>&, int);
  static Map2D calcOrientation(Map2D, int);

  static vector<Map2D> makePyramid(Map2D);
  static Map2D downSampling(Map2D);
  static Map2D upSampling(Map2D);
};

const int SM::SCALE = 9;
const int SM::ORIENTATIONS = 4;
const int SM::GAUSSIAN_K_SIZE = 7;
const int SM::GABOR_K_SIZE = 7;
Gaussian SM::gaussian_kernel = Gaussian(GAUSSIAN_K_SIZE);
GaborBank SM::gabor_bank = GaborBank(GABOR_K_SIZE, ORIENTATIONS);

Map2D SM::calcSaliency(Map2D r, Map2D g, Map2D b)
{
  vector<Map2D> pyramid_r = makePyramid(r);
  vector<Map2D> pyramid_g = makePyramid(g);
  vector<Map2D> pyramid_b = makePyramid(b);

  vector<Map2D> I = makePyramidIntense(pyramid_r, pyramid_g, pyramid_b);

  //double threshold = 0.1 * I[0].max();
  //pyramid_r[0].clampZero(threshold);
  //pyramid_g[0].clampZero(threshold);
  //pyramid_b[0].clampZero(threshold);
  //pyramid_r[0] /= I[0];
  //pyramid_g[0] /= I[0];
  //pyramid_b[0] /= I[0];

  vector<Map2D> R = makePyramidRed(pyramid_r, pyramid_g, pyramid_b);
  vector<Map2D> G = makePyramidGreen(pyramid_r, pyramid_g, pyramid_b);
  vector<Map2D> B = makePyramidBlue(pyramid_r, pyramid_g, pyramid_b);
  vector<Map2D> Y = makePyramidYellow(pyramid_r, pyramid_g, pyramid_b);

  vector<vector<Map2D> > O;
  for (int i=0;i<ORIENTATIONS;++i) {
    vector<Map2D> tmp = makePyramidOrientation(I, i);
    O.push_back(tmp);
  }
  
  return R[0];
}

vector<Map2D> SM::makePyramidIntense(const vector<Map2D>& r, const vector<Map2D>& g, const vector<Map2D>& b)
{
  vector<Map2D> I;

  for (int n=0;n<SCALE;++n) {
    Map2D intense = r[n] + g[n] + b[n] / 3.0;
    I.push_back(intense);
  }

  return I;
}

vector<Map2D> SM::makePyramidRed(vector<Map2D> r, vector<Map2D> g, vector<Map2D> b)
{
  vector<Map2D> R;

  Map2D red;
  for (int n=0;n<SCALE;++n) {
    red = r[n] - (g[n] + b[n]) / 2.0;
    red.clampZero(0);

    R.push_back(red);
  }

  return R;
}

vector<Map2D> SM::makePyramidGreen(vector<Map2D> r, vector<Map2D> g, vector<Map2D> b)
{
  vector<Map2D> G;

  for (int n=0;n<SCALE;++n) {
    Map2D green = g[n] - (b[n] + r[n]) / 2.0;
    green.clampZero(0);

    G.push_back(green);
  }

  return G;
}

vector<Map2D> SM::makePyramidBlue(vector<Map2D> r, vector<Map2D> g, vector<Map2D> b)
{
  vector<Map2D> B;

  for (int n=0;n<SCALE;++n) {
    Map2D blue = b[n] - (r[n] + g[n]) / 2.0;
    blue.clampZero(0);

    B.push_back(blue);
  }

  return B;
}

vector<Map2D> SM::makePyramidYellow(vector<Map2D> r, vector<Map2D> g, vector<Map2D> b)
{
  vector<Map2D> Y;

  for (int n=0;n<SCALE;++n) {
    Map2D yellow = (r[n] + g[n]) / 2.0 - Map2D::abs(r[n] - g[n]) / 2.0 - b[n];
    yellow.clampZero(0);

    Y.push_back(yellow);
  }

  return Y;
}

Map2D SM::calcOrientation(Map2D intense, int o_number)
{
  if (o_number < 0 || ORIENTATIONS <= o_number) {
    cerr << "ERROR: SM::calcOrientation(Map2D, int)  o_nubmer out of range" << endl;
    exit(1);
  }

  Map2D orientation = gabor_bank.kernel[o_number].filter2D(intense);
  
  return orientation;
}

vector<Map2D> SM::makePyramidOrientation(const vector<Map2D>& I, int o_number)
{
  if (o_number < 0 || ORIENTATIONS <= o_number) {
    cerr << "ERROR: SM::calcOrientation(const Map2D&, int)  o_nubmer out of range" << endl;
    exit(1);
  }

  vector<Map2D> O;

  for (int n=0;n<SCALE;++n) {
    Map2D orientation = gabor_bank.kernel[o_number].filter2D(I[n]);
    O.push_back(orientation);
  }

  return O;
}

vector<Map2D> SM::makePyramid(Map2D image)
{
  vector<Map2D> pyramid;

  pyramid.push_back(image); // original image

  for (int i=0;i<SCALE-1;++i) {
    Map2D tmp = downSampling(pyramid[i]);
    pyramid.push_back(tmp);
  }

  return pyramid;
}


Map2D SM::downSampling(Map2D image)
{
  int rows = (int)(image.rows / 2.0);
  int cols = (int)(image.cols / 2.0);
  
  if (rows < 1 || cols < 1) {
    cerr << "ERROR: downSampling image too small" << endl;
    exit(1);
  }

  Map2D _image(rows, cols);
 
  Map2D tmp = gaussian_kernel.filter2D(image);

  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
      _image.data[y][x] = tmp.data[(int)(2*y)][(int)(2*x)];
    }
  }

  return _image;
}


Map2D SM::upSampling(Map2D image) //  bilinear
{
  int rows = image.rows * 2;
  int cols = image.cols * 2;

  Map2D _image(rows, cols);

  for (int y=0;y<image.rows;++y) {
    for (int x=0;x<image.cols;++x) {
      double value = image.data[y][x];
      _image.data[2*y][2*x] = value;
      _image.data[2*y][2*x+1] = value;
      _image.data[2*y+1][2*x] = value;
      _image.data[2*y+1][2*x+1] = value;
      
      /*_image.data[2*y][2*x] = value;
      _image.data[2*y][2*x+1] = (value * image.data[y][x+1]) / 2.0;
      _image.data[2*y+1][2*x] = (value * image.data[y+1][x]) / 2.0;
      _image.data[2*y+1][2*x+1] = (value * image.data[y+1][x+1]) / 2.0;*/
    }
  }

  return _image;
}


int main(int argc, char* argv[])
{
  if (argc != CMD_LINE_NUM) {
    cerr << "ERROR: usage " << argv[0] << " <input directory>" << endl;
    return -1;
  }

  Map2D frame_r;
  Map2D frame_g;
  Map2D frame_b;
  
  int frame_cnt = 1;
  char pic_name[64];
  
  bool flag = true;

  for(;flag;) {
    sprintf(pic_name, "%s/R/frame%05d.pgm", argv[1], frame_cnt);
    flag = getImage(pic_name, frame_r);
    sprintf(pic_name, "%s/G/frame%05d.pgm", argv[1], frame_cnt);
    flag = getImage(pic_name, frame_g);
    sprintf(pic_name, "%s/B/frame%05d.pgm", argv[1], frame_cnt);
    flag = getImage(pic_name, frame_b);   

    if (!flag) {
      break;
    }
    
    sprintf(pic_name, "output/S/frame%05d.pgm", frame_cnt);
    flag = saveImage(pic_name, SM::calcSaliency(frame_r, frame_g, frame_b));

    frame_cnt++;
  }

  return 0;
}
