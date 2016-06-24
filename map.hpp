/*---------------------------------------------------------------------
It is a library that contains a class dealing with image

Map2D class
RGB class

----------------------------------------------------------------------*/

#ifndef __MAP_HPP
#define __MAP_HPP

#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <algorithm>

using namespace std;
  
class Map2D
{
private :
  int rows;
  int cols;
    
public :
  vector<vector<double> > data;
    
  Map2D(void);
  Map2D(int, int);
  Map2D(const Map2D&);
    
  int Rows(void) const;
  int Cols(void) const;
  void Resize(int, int);
  bool Empty(void);
  void ClampZero(double);
  double Min(void);
  double Max(void);
  static Map2D Abs(Map2D);
  void NormalizeRange(int);
    
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


Map2D::Map2D(void)
{
  rows = 0;
  cols = 0;
}

Map2D::Map2D(int rows, int cols)
{
  Resize(rows, cols);
}

Map2D::Map2D(const Map2D& obj)
{
  this->Resize(obj.Rows(), obj.Cols());

  for (int i=0;i<obj.Rows();++i) {
    for (int j=0;j<obj.Cols();++j) {
      this->data[i][j] = obj.data[i][j];
    }
  }
}

int Map2D::Rows(void) const
{
  return rows;
}

int Map2D::Cols(void) const
{
  return cols;
}

void Map2D::Resize(int rows, int cols)
{
  data = vector<vector<double> >(rows, vector<double>(cols, 0));
  this->rows = rows;
  this->cols = cols;
}

void Map2D::ClampZero(double x)
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

double Map2D::Min(void)
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

double Map2D::Max(void)
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

Map2D Map2D::Abs(Map2D obj)
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

void Map2D::NormalizeRange(int range)
{
  double min = this->Min();
  double max = this->Max();
  *this -= min;
  if (min < max) {
    *this /= (max - min);
    *this *= range;
  }

}

bool Map2D::Empty()
{
  if (this->rows == 0 || this->cols == 0) {
    return true;
  }
  else {
    return false;
  }
}

template<class T> Map2D operator + (Map2D obj, T n)
{
  Map2D result = obj;

  for (int i=0;i<obj.Rows();++i) {
    for (int j=0;j<obj.Cols();++j) {
      result.data[i][j] += n;
    }
  }

  return result;
}

template<class T> Map2D operator + (T n, Map2D obj)
{
  Map2D result = obj;

  for (int i=0;i<obj.Rows();++i) {
    for (int j=0;j<obj.Cols();++j) {
      result.data[i][j] += n;
    }
  }

  return result;
}

Map2D operator + (Map2D left, Map2D right)
{
  if (left.Rows() != left.Rows() || right.Cols() != right.Cols()) {
    cerr << "ERROR : Map2D class operator \"+\" failed" << endl;
    exit(1);
  }

  Map2D result = left;

  for (int i=0;i<left.Rows();++i) {
    for (int j=0;j<left.Cols();++j) {
      result.data[i][j] += right.data[i][j];
    }
  }

  return result;
}

template<class T> Map2D operator - (Map2D obj, T n)
{
  Map2D result = obj;

  for (int i=0;i<obj.Rows();++i) {
    for (int j=0;j<obj.Cols();++j) {
      result.data[i][j] -= n;
    }
  }

  return result;
}

template<class T> Map2D operator - (T n, Map2D obj)
{
  Map2D result = obj;

  for (int i=0;i<obj.Rows();++i) {
    for (int j=0;j<obj.Cols();++j) {
      result.data[i][j] -= n;
    }
  }

  return result;
}

Map2D operator - (Map2D left, Map2D right)
{
  if (left.Rows() != right.Rows() || left.Cols() != right.Cols()) {
    cerr << "ERROR : Map2D class operator \"-\" failed" << endl;
    exit(1);
  }

  Map2D result = left;

  for (int i=0;i<left.Rows();++i) {
    for (int j=0;j<left.Cols();++j) {
      result.data[i][j] -= right.data[i][j];
    }
  }

  return result;
}

template<class T> Map2D operator * (Map2D obj, T n)
{
  Map2D result = obj;

  for (int i=0;i<obj.Rows();++i) {
    for (int j=0;j<obj.Cols();++j) {
      result.data[i][j] *= n;
    }
  }

  return result;
}

template<class T> Map2D operator * (T n, Map2D obj)
{
  Map2D result = obj;

  for (int i=0;i<obj.Rows();++i) {
    for (int j=0;j<obj.Cols();++j) {
      result.data[i][j] *= n;
    }
  }

  return result;
}

Map2D operator * (Map2D left, Map2D right)
{
  if (left.Rows() != right.Rows() || left.Cols() != right.Cols()) {
    cerr << "ERROR : Map2D class operator \"*\" failed" << endl;
    exit(1);
  }

  Map2D result = left;

  for (int i=0;i<left.Rows();++i) {
    for (int j=0;j<left.Cols();++j) {
      result.data[i][j] *= right.data[i][j];
    }
  }

  return result;
}

template<class T> Map2D operator / (Map2D obj, T n)
{
  Map2D result = obj;

  if (n == 0) {
    return result;
  }

  for (int i=0;i<obj.Rows();++i) {
    for (int j=0;j<obj.Cols();++j) {
      result.data[i][j] /= n;
    }
  }

  return result;
}

Map2D operator / (Map2D left, Map2D right)
{
  if (left.Rows() != right.Rows() || left.Cols() != right.Cols()) {
    cerr << "ERROR : Map2D class operator \"/\" failed" << endl;
    exit(1);
  }

  Map2D result = left;

  for (int i=0;i<left.Rows();++i) {
    for (int j=0;j<left.Cols();++j) {
      double value = right.data[i][j];
      if (value != 0.0) {
	result.data[i][j] /= value;
      }
    }
  }

  return result;
}

Map2D& Map2D::operator = (const Map2D& obj)
{
  this->Resize(obj.Rows(), obj.Cols());

  for (int i=0;i<obj.Rows();++i) {
    for (int j=0;j<obj.Cols();++j) {
      this->data[i][j] = obj.data[i][j];
    }
  }

  return *this;
}

template<class T> Map2D& Map2D::operator += (const T& n)
{
  for (int i=0;i<this->Rows();++i) {
    for (int j=0;j<this->Cols();++j) {
      this->data[i][j] += n;
    }
  }

  return *this;
}

Map2D& Map2D::operator += (const Map2D& obj)
{
  for (int i=0;i<this->Rows();++i) {
    for (int j=0;j<this->Cols();++j) {
      this->data[i][j] += obj.data[i][j];
    }
  }

  return *this;
}

template<class T> Map2D& Map2D::operator -= (const T& n)
{
  for (int i=0;i<this->Rows();++i) {
    for (int j=0;j<this->Cols();++j) {
      this->data[i][j] -= n;
    }
  }

  return *this;
}

Map2D& Map2D::operator -= (const Map2D& obj)
{
  for (int i=0;i<this->Rows();++i) {
    for (int j=0;j<this->Cols();++j) {
      this->data[i][j] -= obj.data[i][j];
    }
  }

  return *this;
}

template<class T> Map2D& Map2D::operator *= (const T& n)
{
  for (int i=0;i<this->Rows();++i) {
    for (int j=0;j<this->Cols();++j) {
      this->data[i][j] *= n;
    }
  }

  return *this;
}

Map2D& Map2D::operator *= (const Map2D& obj)
{
  for (int i=0;i<this->Rows();++i) {
    for (int j=0;j<this->Cols();++j) {
      this->data[i][j] *= obj.data[i][j];
    }
  }

  return *this;
}

template<class T> Map2D& Map2D::operator /= (const T& n)
{
  for (int i=0;i<this->Rows();++i) {
    for (int j=0;j<this->Cols();++j) {
      this->data[i][j] /= n;
    }
  }

  return *this;
}

Map2D& Map2D::operator /= (const Map2D& obj)
{
  for (int i=0;i<this->Rows();++i) {
    for (int j=0;j<this->Cols();++j) {
      double value = obj.data[i][j];
      if (value != 0.0) {
	this->data[i][j] /= obj.data[i][j];
      }
    }
  }  
  return *this;
}

class Map2D_C3
{
public :
  vector<Map2D> color;

  Map2D_C3(Map2D, Map2D, Map2D);
    
  Map2D& operator[] (int i);
    
};

Map2D_C3::Map2D_C3(Map2D red, Map2D green, Map2D blue)
{
  color.push_back(red);
  color.push_back(green);
  color.push_back(blue);
}
  
Map2D& Map2D_C3::operator[](int i)
{
  if (i < 0 || color.size() <= i) {
    cerr << "ERROR: out of array range" << endl;
    cerr << "in Map2D, operator[](int)" << endl;
      exit(1);
  }
    
  return color[i];
}

#endif
