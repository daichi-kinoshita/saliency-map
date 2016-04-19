#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>  // DBL_MAX
#include <vector>
#include <algorithm>

using namespace std;

static const int CMD_LINE_NUM = 2;

class Map2D;

// get image (.pgn)
Map2D getImage(char*);

// save image (.pgn) 
bool saveImage(char*, Map2D);

/* vector<vector<double>> class */
class Map2D
{
public :
  vector<vector<double> > data;
  int rows;
  int cols;

  Map2D();
  Map2D(int, int);
  Map2D(const Map2D&);

  void set(int, int);
  bool empty();
  void clampZero(double);
  double min();
  double max();
  static Map2D abs(Map2D);
  void resize(int, int);
  void normalizeRange(void);

  template<class T> friend Map2D operator + (Map2D, T);
  template<class T> friend Map2D operator + (T, Map2D);
  friend Map2D operator + (Map2D, Map2D);
  template<class T> friend Map2D operator - (Map2D, T);
  template<class T> friend Map2D operator - (T, Map2D);
  friend Map2D operator - (Map2D, Map2D);
  template<class T> friend Map2D operator * (Map2D, T);
  template<class T> friend Map2D operator * (T, Map2D);
  friend Map2D operator * (Map2D, Map2D);
  template<class T> friend Map2D operator / (Map2D, T);
  template<class T> friend Map2D operator / (T, Map2D);
  friend Map2D operator / (Map2D, Map2D);

  Map2D& operator = (const Map2D&);
  template<class T> Map2D& operator += (const T&);
  Map2D& operator += (const Map2D&);
  template<class T> Map2D& operator -= (const T&);
  Map2D& operator -= (const Map2D&);
  template<class T> Map2D& operator *= (const T&);
  Map2D& operator *= (const Map2D&);
  template<class T> Map2D& operator /= (const T&);
  Map2D& operator /= (const Map2D&);

};

Map2D::Map2D()
{
  rows = 0;
  cols = 0;
}

Map2D::Map2D(int rows, int cols)
{
  set(rows, cols);
}

Map2D::Map2D(const Map2D& obj)
{
  this->set(obj.rows, obj.cols);

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      this->data[i][j] = obj.data[i][j];
    }
  }
}

void Map2D::set(int rows, int cols)
{
  data = vector<vector<double> >(rows, vector<double>(cols, 0));
  this->rows = rows;
  this->cols = cols;
}

void Map2D::clampZero(double x)
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

double Map2D::min(void)
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

double Map2D::max(void)
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

Map2D Map2D::abs(Map2D obj)
{
  Map2D result(obj.rows, obj.cols);

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

void Map2D::resize(int rows, int cols)
{
  Map2D tmp = *this;

  this->data = vector<vector<double> >(rows, vector<double>(cols, 0));
  this->rows = rows;
  this->cols = cols;

  int _rows;
  int _cols;
  if (tmp.rows < this->rows) {
    _rows = tmp.rows;
  }
  else {
    _rows = this->rows;
  }

  if (tmp.cols < this->cols) {
    _cols = tmp.cols;
  }
  else {
    _cols = this->cols;
  }

  for (int y=0;y<_rows;++y) {
    for (int x=0;x<_cols;++y) {
      this->data[y][x] = tmp.data[y][x];
    }
  }
}

void Map2D::normalizeRange(void)
{
  double minimum = this->min();
  double maximum = this->max();
  *this -= minimum;
  if (minimum < maximum) {
    *this /= maximum - minimum;
  }

}

template<class T> Map2D operator + (Map2D obj, T n)
{
  Map2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] += n;
    }
  }

  return result;
}

template<class T> Map2D operator + (T n, Map2D obj)
{
  Map2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] += n;
    }
  }

  return result;
}

Map2D operator + (Map2D left, Map2D right)
{
  if (left.rows != left.rows || right.cols != right.cols) {
    cerr << "ERROR : Map2D class operator \"+\" failed" << endl;
    exit(1);
  }

  Map2D result = left;

  for (int i=0;i<left.rows;++i) {
    for (int j=0;j<left.cols;++j) {
      result.data[i][j] += right.data[i][j];
    }
  }

  return result;
}

template<class T> Map2D operator - (Map2D obj, T n)
{
 Map2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] -= n;
    }
  }

  return result;
}

template<class T> Map2D operator - (T n, Map2D obj)
{
 Map2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] -= n;
    }
  }

  return result;
}

Map2D operator - (Map2D left, Map2D right)
{
  if (left.rows != right.rows || left.cols != right.cols) {
    cerr << "ERROR : Map2D class operator \"-\" failed" << endl;
    exit(1);
  }

  Map2D result = left;

  for (int i=0;i<left.rows;++i) {
    for (int j=0;j<left.cols;++j) {
      result.data[i][j] -= right.data[i][j];
    }
  }

  return result;
}

template<class T> Map2D operator * (Map2D obj, T n)
{
 Map2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] *= n;
    }
  }

  return result;
}

template<class T> Map2D operator * (T n, Map2D obj)
{
 Map2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] *= n;
    }
  }

  return result;
}

Map2D operator * (Map2D left, Map2D right)
{
  if (left.rows != right.rows || left.cols != right.cols) {
    cerr << "ERROR : Map2D class operator \"*\" failed" << endl;
    exit(1);
  }

  Map2D result = left;

  for (int i=0;i<left.rows;++i) {
    for (int j=0;j<left.cols;++j) {
      result.data[i][j] *= right.data[i][j];
    }
  }

  return result;
}

template<class T> Map2D operator / (Map2D obj, T n)
{
 Map2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] /= n;
    }
  }

  return result;
}

template<class T> Map2D operator / (T n, Map2D obj)
{
 Map2D result = obj;

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      result.data[i][j] /= n;
    }
  }

  return result;
}

Map2D operator / (Map2D left, Map2D right)
{
  if (left.rows != right.rows || left.cols != right.cols) {
    cerr << "ERROR : Map2D class operator \"/\" failed" << endl;
    exit(1);
  }

  Map2D result = left;

  for (int i=0;i<left.rows;++i) {
    for (int j=0;j<left.cols;++j) {
      result.data[i][j] /= right.data[i][j];
    }
  }

  return result;
}

Map2D& Map2D::operator = (const Map2D& obj)
{
  this->set(obj.rows, obj.cols);

  for (int i=0;i<obj.rows;++i) {
    for (int j=0;j<obj.cols;++j) {
      this->data[i][j] = obj.data[i][j];
    }
  }

  return *this;
}

template<class T> Map2D& Map2D::operator += (const T& n)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] += n;
    }
  }

  return *this;
}

Map2D& Map2D::operator += (const Map2D& obj)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] += obj.data[i][j];
    }
  }

  return *this;
}

template<class T> Map2D& Map2D::operator -= (const T& n)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] -= n;
    }
  }

  return *this;
}

Map2D& Map2D::operator -= (const Map2D& obj)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] -= obj.data[i][j];
    }
  }

  return *this;
}

template<class T> Map2D& Map2D::operator *= (const T& n)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] *= n;
    }
  }

  return *this;
}

Map2D& Map2D::operator *= (const Map2D& obj)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] *= obj.data[i][j];
    }
  }

  return *this;
}

template<class T> Map2D& Map2D::operator /= (const T& n)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] /= n;
    }
  }

  return *this;
}

Map2D& Map2D::operator /= (const Map2D& obj)
{
  for (int i=0;i<this->rows;++i) {
    for (int j=0;j<this->cols;++j) {
      this->data[i][j] /= obj.data[i][j];
    }
  }

  return *this;
}

bool Map2D::empty()
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
  Filter(const Filter&);

  void filter2D(Map2D&);

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

Filter::Filter(const Filter& obj)
{
  this->kernel = obj.kernel;
  this->sigma = obj.sigma;
  this->k_size = obj.k_size;
}

/* image * kernel */
void Filter::filter2D(Map2D& image)
{
  Map2D tmp = image;

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
  sigma = 0.3 * (k_size / 2.0 - 1) + 0.8;

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
  double lambda = k_size / 1.8;

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

/* calculate Saliency Map class */
class SM
{
public :
  static const int SCALE;               // num of gaussian pyramid
  static const int ORIENTATION;         // num of orientation
  static const int GAUSSIAN_K_SIZE;     // gaussian filter kernel size
  static const int GABOR_K_SIZE;        // gabor filter kernel size

  SM(Map2D, Map2D, Map2D);
  
  Map2D Saliency;         //Saliency map
  Map2D I;
  //vector<Map2D> I;        //Intense

  vector<Map2D> Ori000;   //Orientation
  vector<Map2D> Ori045;
  vector<Map2D> Ori090;
  vector<Map2D> Ori135;
  
private :
  static Gaussian gaussian_kernel;
  static Gabor gabor_kernel000;
  static Gabor gabor_kernel045;
  static Gabor gabor_kernel090;
  static Gabor gabor_kernel135;

  vector<Map2D> calcColor(Map2D, Map2D, Map2D, Map2D);
  //static vector<Map2D> calcFlicker(Map2D, Map2D);
  //static vector<Map2D> calcMotion(vector<Map2D>, vector<Map2D>);

  void calcSaliency(Map2D, Map2D, Map2D);
  //void calcSaliency(Map2D, Map2D, Map2D, Map2D, vector<Map2D>);
  void calcIntense(vector<Map2D>, vector<Map2D>, vector<Map2D>);
  //void calcOrientation(Map2D);

  static Map2D downSampling(Map2D);
  static Map2D upSampling(Map2D);
  static vector<Map2D> makePyramid(Map2D);
};

const int SM::SCALE = 9;
const int SM::ORIENTATION = 4;
const int SM::GAUSSIAN_K_SIZE = 9;
const int SM::GABOR_K_SIZE = 17;

Gaussian SM::gaussian_kernel = Gaussian(GAUSSIAN_K_SIZE);
Gabor SM::gabor_kernel000 = Gabor(GABOR_K_SIZE, 0.0 * M_PI);
Gabor SM::gabor_kernel045 = Gabor(GABOR_K_SIZE, 0.25 * M_PI);
Gabor SM::gabor_kernel090 = Gabor(GABOR_K_SIZE, 0.5 * M_PI);
Gabor SM::gabor_kernel135 = Gabor(GABOR_K_SIZE, 0.75 * M_PI);

SM::SM(Map2D r, Map2D g, Map2D b)
{
  vector<Map2D> pyramid_r = makePyramid(r);
  vector<Map2D> pyramid_g = makePyramid(g);
  vector<Map2D> pyramid_b = makePyramid(b);

  I = (r + g + b) / 3.0;
  //vector<Map2D> intense = makePyramid(I);

  //calcIntense(pyramid_r, pyramid_g, pyramid_b);
  //calcColor(pyramid_r, pyramid_g, pyramid_b, I);
  //calcOrientation();
}

void SM::calcSaliency(Map2D r, Map2D g, Map2D b)
{

}

void SM::calcIntense(vector<Map2D> r, vector<Map2D> g, vector<Map2D> b)
{
  Map2D intense;
  
  for (int n=0;n<SCALE;++n) {
    intense = r[n] + g[n] + b[n] / 3.0;
    // I.push_back(intense);
  }

}

vector<Map2D> SM::calcColor(Map2D r, Map2D g, Map2D b, Map2D I)
{
  vector<Map2D> color;

  double threshold = 0.1 * I.max();
  r.clampZero(threshold);
  g.clampZero(threshold);
  b.clampZero(threshold);

  r /= I;
  g /= I;
  b /= I;

  Map2D R = r - (g + b) / 2.0;
  Map2D G = g - (b + r) / 2.0;
  Map2D B = b - (r + g) / 2.0;
  Map2D Y = (r + g) / 2.0 - Map2D::abs(r - g) / 2.0 - b;

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

/*vector<Map2D> SM::calcOrientation(Map2D I)
{
  vector<Map2D> O;

  return O;
}*/

Map2D SM::downSampling(Map2D image)
{
  int rows = (int)(image.rows / 2.0);
  int cols = (int)(image.cols / 2.0);

  if (rows < 1 || cols < 1) {
    cerr << "ERROR: downSampling image too small" << endl;
    exit(1);
  }

  Map2D _image(rows, cols);

  gaussian_kernel.filter2D(image);

  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
      _image.data[y][x] = image.data[2*y][2*x];
    }
  }

  //image.resize(rows, cols);

  return _image;
}

Map2D SM::upSampling(Map2D image)
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
    }
  }

  return _image;
}

vector<Map2D> SM::makePyramid(Map2D image)
{
  vector<Map2D> pyramid;

  pyramid.push_back(image); // original image

  for (int i=0;i<SCALE-1;++i) {
    Map2D tmp = downSampling(pyramid[i]);
    pyramid.push_back(tmp);
  }

}

bool getImage(char* pic_name, Map2D& image)
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

bool saveImage(char* file_name, Map2D image)
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

  Map2D frame_r;
  Map2D frame_g;
  Map2D frame_b;
  
  int frame_cnt = 1;
  char pic_name[64];
  
  bool flag = true;

  SM *frame_s;

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
    frame_s = new SM(frame_r, frame_g, frame_b);
    flag = saveImage(pic_name, frame_s->I);
    delete frame_s;

    frame_cnt++;
  }

  return 0;
}
