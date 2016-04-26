#ifndef INCLUDE_MAP_HPP
#define INCLUDE_MAP_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <algorithm>
#include "map.hpp"

using namespace std;

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
  Map2D result = obj;
  
  int rows = obj.rows;
  int cols = obj.cols;

  for (int i=0;i<rows;++i) {
    for (int j=0;j<cols;++j) {
      
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

#endif
