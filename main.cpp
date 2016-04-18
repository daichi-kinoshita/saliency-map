#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <vector>
#include <algorithm>

using namespace std;

static const int CMD_LINE_NUM = 2;

// get image (.pgn)
Vec2D getImage(char*);

// save image (.pgn) 
bool saveImage(char*, Vec2D);

/* vector<vector<double>> class */
class Vec2D
{
public :
  vector<vector<double> > data;
  int rows;
  int cols;

  Vec2D();
  Vec2D(int, int);
  Vec2D(const Vec2D&);

  void set(int, int);
  bool empty();
  void clampZero(double);
  double min();
  double max();
  static Vec2D abs(Vec2D);
  void normalizeRange(void);

  template<class T> friend Vec2D operator + (Vec2D, T);
  template<class T> friend Vec2D operator + (T, Vec2D);
  friend Vec2D operator + (Vec2D, Vec2D);
  template<class T> friend Vec2D operator - (Vec2D, T);
  template<class T> friend Vec2D operator - (T, Vec2D);
  friend Vec2D operator - (Vec2D, Vec2D);
  template<class T> friend Vec2D operator * (Vec2D, T);
  template<class T> friend Vec2D operator * (T, Vec2D);
  friend Vec2D operator * (Vec2D, Vec2D);
  template<class T> friend Vec2D operator / (Vec2D, T);
  template<class T> friend Vec2D operator / (T, Vec2D);
  friend Vec2D operator / (Vec2D, Vec2D);

  Vec2D& operator = (const Vec2D&);
  template<class T> Vec2D& operator += (const T&);
  Vec2D& operator += (const Vec2D&);
  template<class T> Vec2D& operator -= (const T&);
  Vec2D& operator -= (const Vec2D&);
  template<class T> Vec2D& operator *= (const T&);
  Vec2D& operator *= (const Vec2D&);
  template<class T> Vec2D& operator /= (const T&);
  Vec2D& operator /= (const Vec2D&);

};

Vec2D::Vec2D()
{
  rows = 0;
  cols = 0;
}

Vec2D::Vec2D(int rows, int cols)
{
  set(rows, cols);
}

Vec2D::Vec2D(const Vec2D& obj)
{
  this->set(obj.rows, obj.cols);

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      this->data[i][j] = obj.data[i][j];
    }
  }
}

void Vec2D::set(int rows, int cols)
{
  data = vector<vector<double> >(rows, vector<double>(cols, 0));
  this->rows = rows;
  this->cols = cols;
}

void Vec2D::clampZero(double x)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {

      double tmp = this->data[i][j];
      if (tmp < x) {
	this->data[i][j] = 0;
      }

    }
  }
}

double Vec2D::min(void)
{
  double min = DBL_MAX;
  for (int y=0;y<this->rows;++y) {
    double tmp = *min_element(this->data[y].begin(), this->data[y].end());
    if (tmp < min) {
      min = tmp;
    }
  }

  return min;
}

double Vec2D::max(void)
{
  double max = 0.0;
  for (int y=0;y<rows;++y) {
    double tmp = *max_element(this->data[y].begin(), this->data[y].end());
    if (tmp > max) {
      max = tmp;
    }
  }

  return max;
}

Vec2D Vec2D::abs(Vec2D obj)
{
  Vec2D result(obj.rows, obj.cols);

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      
      double tmp = obj.data[i][j];
      if (tmp < 0) {
	tmp = std::abs(tmp);
      }
      result.data[i][j] = tmp;
    }
  }

  return result;
}

void Vec2D::normalizeRange(void)
{
  double minimum = this->min();
  double maximum = this->max();
  *this -= minimum;
  if (minimum < maximum) {
    *this /= maximum - minimum;
  }

}

template<class T> Vec2D operator + (Vec2D obj, T n)
{
  Vec2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] += n;
    }
  }

  return result;
}

template<class T> Vec2D operator + (T n, Vec2D obj)
{
  Vec2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] += n;
    }
  }

  return result;
}

Vec2D operator + (Vec2D left, Vec2D right)
{
  if (left.rows != left.rows || right.cols != right.cols) {
    cerr << "ERROR : Vec2D class operator \"+\" failed" << endl;
    exit(1);
  }

  Vec2D result = left;

  for (int i=0;i<left.rows;++i) {
    for (int j=0;j<left.cols;++j) {
      result.data[i][j] += right.data[i][j];
    }
  }

  return result;
}

template<class T> Vec2D operator - (Vec2D obj, T n)
{
 Vec2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] -= n;
    }
  }

  return result;
}

template<class T> Vec2D operator - (T n, Vec2D obj)
{
 Vec2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] -= n;
    }
  }

  return result;
}

Vec2D operator - (Vec2D left, Vec2D right)
{
  if (left.rows != right.rows || left.cols != right.cols) {
    cerr << "ERROR : Vec2D class operator \"-\" failed" << endl;
    exit(1);
  }

  Vec2D result = left;

  for (int i=0;i<left.rows;++i) {
    for (int j=0;j<left.cols;++j) {
      result.data[i][j] -= right.data[i][j];
    }
  }

  return result;
}

template<class T> Vec2D operator * (Vec2D obj, T n)
{
 Vec2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] *= n;
    }
  }

  return result;
}

template<class T> Vec2D operator * (T n, Vec2D obj)
{
 Vec2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] *= n;
    }
  }

  return result;
}

Vec2D operator * (Vec2D left, Vec2D right)
{
  if (left.rows != right.rows || left.cols != right.cols) {
    cerr << "ERROR : Vec2D class operator \"*\" failed" << endl;
    exit(1);
  }

  Vec2D result = left;

  for (int i=0;i<left.rows;++i) {
    for (int j=0;j<left.cols;++j) {
      result.data[i][j] *= right.data[i][j];
    }
  }

  return result;
}

template<class T> Vec2D operator / (Vec2D obj, T n)
{
 Vec2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] /= n;
    }
  }

  return result;
}

template<class T> Vec2D operator / (T n, Vec2D obj)
{
 Vec2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] /= n;
    }
  }

  return result;
}

Vec2D operator / (Vec2D left, Vec2D right)
{
  if (left.rows != right.rows || left.cols != right.cols) {
    cerr << "ERROR : Vec2D class operator \"/\" failed" << endl;
    exit(1);
  }

  Vec2D result = left;

  for (int i=0;i<left.rows;++i) {
    for (int j=0;j<left.cols;++j) {
      result.data[i][j] /= right.data[i][j];
    }
  }

  return result;
}

Vec2D& Vec2D::operator = (const Vec2D& obj)
{
  this->set(obj.rows, obj.cols);

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      this->data[i][j] = obj.data[i][j];
    }
  }

  return *this;
}

template<class T> Vec2D& Vec2D::operator += (const T& n)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] += n;
    }
  }

  return *this;
}

Vec2D& Vec2D::operator += (const Vec2D& obj)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] += obj.data[i][j];
    }
  }

  return *this;
}

template<class T> Vec2D& Vec2D::operator -= (const T& n)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] -= n;
    }
  }

  return *this;
}

Vec2D& Vec2D::operator -= (const Vec2D& obj)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] -= obj.data[i][j];
    }
  }

  return *this;
}

template<class T> Vec2D& Vec2D::operator *= (const T& n)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] *= n;
    }
  }

  return *this;
}

Vec2D& Vec2D::operator *= (const Vec2D& obj)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] *= obj.data[i][j];
    }
  }

  return *this;
}

template<class T> Vec2D& Vec2D::operator /= (const T& n)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] /= n;
    }
  }

  return *this;
}

Vec2D& Vec2D::operator /= (const Vec2D& obj)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] /= obj.data[i][j];
    }
  }

  return *this;
}

bool Vec2D::empty()
{
  if (this->rows == 0 || this->cols == 0) {
    return true;
  }
  else {
    return false;
  }
}

/* Filter class */
class Filter
{
public :
  Filter(int);

  void filter2D(Vec2D&);

protected :
  Vec2D kernel;
  double sigma;   // filter sigma
  int k_size;     // filter kernel size
};

Filter::Filter(int k_size)
{
  if (k_size % 2 == 0) {
    k_size++;
  }

  this->k_size = k_size;
  kernel = Vec2D(k_size, k_size);
}

/* image * kernel */
void Filter::filter2D(Vec2D& image)
{
  Vec2D tmp = image;

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

	  sum += tmp.data[yy][xx] * kernel.data[dy+h_k_size][dx+h_k_size];
	}
      }

      image.data[y][x] = sum;
    }
  }

  image.clampZero(0);

}

/* Gaussian class */
class GaussianFilter : public Filter
{
public :
  GaussianFilter(int);
  GaussianFilter(int, double);

private :
  void makeGaussianKernel(void);
};

GaussianFilter::GaussianFilter(int k_size) : Filter(k_size)
{
  sigma = 0.3 * (k_size / 2.0 - 1) + 0.8;

  makeGaussianKernel();
}

GaussianFilter::GaussianFilter(int k_size, double sigma) : Filter(k_size)
{
  this->sigma = sigma;

  makeGaussianKernel();
}

void GaussianFilter::makeGaussianKernel(void)
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
class GaborFilter : public Filter
{
public :
  GaborFilter(int, double, double, double, double);
  //size, theta, sigma, gamma, psi, band_width
  GaborFilter(int, double, double, double, double, double);

private :
  static const double _GAMMA;       // gamma default
  static const double _PSI;         // psi default
  static const double _BAND_W;      // band width default

  double lambda;      // gabor filter lambda
  double gamma;       // gabor filter gamma
  double psi;         // gabor filter psi
  double band_w;      // band width
  double theta;       // orientation rad
  
  void makeGaborKernel();
};

GaborFilter::GaborFilter(int k_size, double theta, double gamma=1.0, double psi=0.0, double band_w=1.0) : Filter(k_size)
{
  lambda = k_size / 1.8;

  this->gamma = gamma;
  this->psi = psi;
  this->band_w = band_w;
  this->theta = theta;

  sigma = sqrt(log(2.0) / 2.0) / M_PI * lambda;
  sigma *= (pow(2, band_w) + 1.0 ) / (pow(2, band_w) - 1.0);

  makeGaborKernel();
}

GaborFilter::GaborFilter(int k_size, double theta, double sigma, double gamma=1.0, double psi=0.0, double band_w=1.0) : Filter(k_size)
{
  lambda = (int)(k_size / 2.0);

  this->gamma = gamma;
  this->psi = psi;
  this->band_w = band_w;
  this->theta = theta;
  this->sigma = sigma;

  makeGaborKernel();
}

void GaborFilter::makeGaborKernel(void)
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

/* calculate Saliency Map class */
class SM
{
public :
  static const int SCALE;               // num of gaussian pyramid
  static const int ORIENTATION;         // num of orientation
  static const int GAUSSIAN_K_SIZE;     // gaussian filter kernel size
  static const int GABOR_K_SIZE;        // gabor filter kernel size

  Vec2D calcSaliency(Vec2D, Vec2D, Vec2D);
  //static Vec2D calcSaliency(Vec2D, Vec2D, Vec2D, Vec2D, vector<Vec2D>);
  static Vec2D calcIntense(Vec2D, Vec2D, Vec2D);
  //static vector<Vec2D> calcOrientation(Vec2D);


private :
  vector<Vec2D> Color;    //R,G,B,Y
  vector<Vec2D> O;        //Orientation
  Vec2D I;                //Intense

  vector<Vec2D> calcColor(Vec2D, Vec2D, Vec2D, Vec2D);
  //static vector<Vec2D> calcFlicker(Vec2D, Vec2D);
  //static vector<Vec2D> calcMotion(vector<Vec2D>, vector<Vec2D>);
  SM(){};
};

const int SM::SCALE = 9;
const int SM::ORIENTATION = 4;
const int SM::GAUSSIAN_K_SIZE = 9;
const int SM::GABOR_K_SIZE = 17;

Vec2D SM::calcSaliency(Vec2D r, Vec2D g, Vec2D b)
{
  Vec2D saliency(r.rows, r.cols);

  I = calcIntense(r, g, b);

  Color = calcColor(r, g, b, I);

  //GaborFilter gabor000(GABOR_K_SIZE);
  //filter2D(I, gabor_kernel);

  return I;
}

Vec2D SM::calcIntense(Vec2D r, Vec2D g, Vec2D b)
{
  Vec2D I = (r + g + b) / 3.0;

  return I;
}

vector<Vec2D> SM::calcColor(Vec2D r, Vec2D g, Vec2D b, Vec2D I)
{
  vector<Vec2D> color;

  double threshold = 0.1 * I.max();
  r.clampZero(threshold);
  g.clampZero(threshold);
  b.clampZero(threshold);

  r = r / I;
  g = g / I;
  b = b / I;

  Vec2D R = r - (g + b) / 2.0;
  Vec2D G = g - (b + r) / 2.0;
  Vec2D B = b - (r + g) / 2.0;
  Vec2D Y = (r + g) / 2.0 - Vec2D::abs(r - g) / 2.0 - b;

  R.clampZero(0);
  G.clampZero(0);
  B.clampZero(0);
  Y.clampZero(0);

  color.push_back(R);
  color.push_back(G);
  color.push_back(B);
  color.push_back(Y);

  return color;
}

/*vector<Vec2D> SM::calcOrientation(Vec2D I)
{
  vector<Vec2D> O;

  return O;
  }*/

bool getImage(char* pic_name, Vec2D& image)
{
  ifstream pic_stream(pic_name);
  if (pic_stream == NULL) {
    return false;
  }

  //header
  string tmp;
  int rows;
  int cols;
  pic_stream >> tmp;
  pic_stream >> cols >> rows;

  image.set(rows, cols);

  //data (0.0-1.0)
  double data;
  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
      pic_stream >> data;
      image.data[y][x] = data / 255.0;
    }
  }
  
  return true;
}

bool saveImage(char* file_name, Vec2D image)
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

  //header
  pic_stream << "P2" << endl;
  pic_stream << cols << " " << rows << endl;
  pic_stream << 255 << endl;
  
  //data (0-255)
  image.normalizeRange();
  image *= 255; 
  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
      pic_stream << (int)image.data[y][x] << " ";
    }
    pic_stream << endl;
  }

  return true;
}


int main(int argc, char* argv[])
{
  if (argc != CMD_LINE_NUM) {
    cerr << "ERROR: usage " << argv[0] << " <input directory>" << endl;
    return -1;
  }

  Vec2D frame_r;
  Vec2D frame_g;
  Vec2D frame_b;
  
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

    //sprintf(pic_name, "output/I/frame%05d.pgm", frame_cnt);
    //flag = saveImage(pic_name, SM::calcIntense(frame_r, frame_g, frame_b));

    sprintf(pic_name, "output/S/frame%05d.pgm", frame_cnt);
    Vec2D frame_s = SM::calcSaliency(frame_r, frame_g, frame_b);
    flag = saveImage(pic_name, frame_s);
    
    frame_cnt++;
  }

  return 0;
}
