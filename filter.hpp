/*---------------------------------------------------------------------
It is a library for performing a filter processing on the image data.

Filter class
Gaussian:filter class
Gabor:filter class
GaborBank class

----------------------------------------------------------------------*/

#ifndef __FILTER_HPP
#define __FILTER_HPP

#include <cmath>
#include <vector>
#include "map.hpp"

/* Filter class */
class Filter
{
public :
  Filter(int);

  Map2D Filter2D(Map2D);

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

Map2D Filter::Filter2D(Map2D image)
{
  Map2D _image = image;

  int rows = image.Rows();
  int cols = image.Cols();
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

  //_image = Map2D::abs(_image);

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
  void MakeGaussianKernel(void);
};

Gaussian::Gaussian(int k_size) : Filter(k_size)
{
  sigma = 0.3 * (k_size / 2.0 - 1) + 0.5;

  MakeGaussianKernel();
}

Gaussian::Gaussian(int k_size, double sigma) : Filter(k_size)
{
  this->sigma = sigma;

  MakeGaussianKernel();
}

Gaussian::Gaussian(const Gaussian& obj) : Filter(obj){}

void Gaussian::MakeGaussianKernel(void)
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
  void MakeGaborKernel(double, double, double, double);
};

Gabor::Gabor(int k_size, double theta, double band_w=1.0, double gamma=1.0, double psi=0.0) : Filter(k_size)
{
  double lambda = k_size / 2.0;

  sigma = sqrt(log(2.0) / 2.0) / M_PI * lambda;
  sigma *= (pow(2, band_w) + 1.0 ) / (pow(2, band_w) - 1.0);

  MakeGaborKernel(theta, lambda, gamma, psi);
}

Gabor::Gabor(const Gabor& obj) : Filter(obj){}

void Gabor::MakeGaborKernel(double theta, double lambda, double gamma, double psi)
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

    if (i == 0 || i == 2) {
      Gabor tmp(k_size, (double)i / orientation * M_PI, 1.205);
      kernel.push_back(tmp);
    }
    else {
      Gabor tmp(k_size, (double)i / orientation * M_PI);
      kernel.push_back(tmp);
    }

  }
}

GaborBank::GaborBank(const GaborBank& obj)
{
  this->size = obj.size;
  this->kernel = obj.kernel;
}


#endif
