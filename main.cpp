#include "all.hpp"
#include "io_image.hpp"
#include "map.hpp"
#include "filter.hpp"
#include "progress.hpp"

using namespace std;
using namespace io;

static const int CMD_LINE_NUM = 4;

/* calculate SaliencyMap class */
class SM
{
public :
  static Map2D CalcSaliency(Map2D, Map2D, Map2D);
  //static Map2D calcSaliency(Map2D, Map2D, Map2D, Map2D, Map2D, Map2D);

private :
  SM(){};

  static const int SCALE;               // num of pyramid
  static const int ORIENTATIONS;        // num of orientation
  static const int GAUSSIAN_K_SIZE;     // gaussian filter kernel size
  static const int GABOR_K_SIZE;        // gabor filter kernel size
  static Gaussian gaussian_kernel;
  static GaborBank gabor_bank;

  //static vector<Map2D> CalcFlicker(Map2D, Map2D);
  //static vector<Map2D> CalcMotion(vector<Map2D>, vector<Map2D>);
  
  static vector<Map2D> MakePyramidIntense(Map2D, Map2D, Map2D);
  
  static vector<Map2D> MakePyramidRed(Map2D, Map2D, Map2D);
  static vector<Map2D> MakePyramidGreen(Map2D, Map2D, Map2D);
  static vector<Map2D> MakePyramidBlue(Map2D, Map2D, Map2D);
  static vector<Map2D> MakePyramidYellow(Map2D, Map2D, Map2D);

  static vector<Map2D> MakePyramidOrientation(const vector<Map2D>&, int);
  static Map2D CalcOrientation(Map2D, int);

  static vector<Map2D> MakePyramid(Map2D);
  static Map2D DownSampling(Map2D);
  static Map2D UpSampling(Map2D);
};

const int SM::SCALE = 9;
const int SM::ORIENTATIONS = 4;
const int SM::GAUSSIAN_K_SIZE = 7;
const int SM::GABOR_K_SIZE = 9;
Gaussian SM::gaussian_kernel = Gaussian(GAUSSIAN_K_SIZE);
GaborBank SM::gabor_bank = GaborBank(GABOR_K_SIZE, ORIENTATIONS);

Map2D SM::CalcSaliency(Map2D r, Map2D g, Map2D b)
{
  vector<Map2D> I = MakePyramidIntense(r, g, b);

  double threshold = 0.1 * I[0].Max();

  //r.clampZero(threshold);
  //g.clampZero(threshold);
  //b.clampZero(threshold);
  //r /= I[0];
  //g /= I[0];
  //b /= I[0];
 
  vector<Map2D> R = MakePyramidRed(r, g, b);
  vector<Map2D> G = MakePyramidGreen(r, g, b);
  vector<Map2D> B = MakePyramidBlue(r, g, b);
  vector<Map2D> Y = MakePyramidYellow(r, g, b);

  vector<vector<Map2D> > O;
  for (int i=0;i<ORIENTATIONS;++i) {
    vector<Map2D> tmp = MakePyramidOrientation(I, i);
    O.push_back(tmp);
  }
  
  return O[0][2];
}

vector<Map2D> SM::MakePyramidIntense(Map2D r, Map2D g, Map2D b)
{
  vector<Map2D> I;

  Map2D intense = (r + g + b) / 3.0;

  I = MakePyramid(intense);

  return I;
}

vector<Map2D> SM::MakePyramidRed(Map2D r, Map2D g, Map2D b)
{
  vector<Map2D> R;

  Map2D red = r - ((g + b) / 2.0);
  red.ClampZero(0);

  R = MakePyramid(red);

  return R;
}

vector<Map2D> SM::MakePyramidGreen(Map2D r, Map2D g, Map2D b)
{
  vector<Map2D> G;

  Map2D green = g - ((b + r) / 2.0);
  green.ClampZero(0);
  
  G = MakePyramid(green);

  return G;
}

vector<Map2D> SM::MakePyramidBlue(Map2D r, Map2D g, Map2D b)
{
  vector<Map2D> B;

  Map2D blue = b - ((r + g) / 2.0);
  blue.ClampZero(0);

  B = MakePyramid(blue);

  return B;
}

vector<Map2D> SM::MakePyramidYellow(Map2D r, Map2D g, Map2D b)
{
  vector<Map2D> Y;

  Map2D yellow = ((r + g) / 2.0) - (Map2D::Abs(r - g) / 2.0) - b;
  yellow.ClampZero(0);
  
  Y = MakePyramid(yellow);

  return Y;
}

Map2D SM::CalcOrientation(Map2D intense, int o_number)
{
  if (o_number < 0 || ORIENTATIONS <= o_number) {
    cerr << "ERROR: SM::calcOrientation(Map2D, int)  o_nubmer out of range" << endl;
    exit(1);
  }

  Map2D orientation = gabor_bank.kernel[o_number].Filter2D(intense);
  
  return orientation;
}

vector<Map2D> SM::MakePyramidOrientation(const vector<Map2D>& I, int o_number)
{
  if (o_number < 0 || ORIENTATIONS <= o_number) {
    cerr << "ERROR: SM::calcOrientation(const Map2D&, int)  o_nubmer out of range" << endl;
    exit(1);
  }

  vector<Map2D> O;

  for (int n=0;n<SCALE;++n) {
    Map2D orientation = gabor_bank.kernel[o_number].Filter2D(I[n]);
    O.push_back(orientation);
  }

  return O;
}

vector<Map2D> SM::MakePyramid(Map2D image)
{
  vector<Map2D> pyramid;

  pyramid.push_back(image); // original image

  for (int i=0;i<SCALE-1;++i) {
    Map2D tmp = DownSampling(pyramid[i]);
    pyramid.push_back(tmp);
  }

  return pyramid;
}


Map2D SM::DownSampling(Map2D image)
{
  int rows = (int)(image.Rows() / 2.0);
  int cols = (int)(image.Cols() / 2.0);
  
  if (rows < 1 || cols < 1) {
    cerr << "ERROR: downSampling image too small" << endl;
    cerr << "in SM class, downSampling(Map2D)" << endl;
    exit(1);
  }

  Map2D _image(rows, cols);
 
  Map2D tmp = gaussian_kernel.Filter2D(image);

  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
      _image.data[y][x] = tmp.data[(int)(2*y)][(int)(2*x)];
    }
  }

  return _image;
}


Map2D SM::UpSampling(Map2D image) //  bilinear
{
  int rows = image.Rows() * 2;
  int cols = image.Cols() * 2;

  Map2D _image(rows, cols);

  for (int y=0;y<rows;++y) {
    for (int x=0;x<cols;++x) {
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
    cerr << "ERROR: usage " << argv[0] << " <input> <output> <frames>" << endl;
    cerr << "in main(int, char**)" << endl;
    return -1;
  }

  vector<Map2D> frame; // color
  
  int frame_num = atoi(argv[3]);
  int frame_cnt = 1;
  char pic_name[64];
  
  bool flag = true;


  cerr << "processing saliency map..." << endl;
  cerr << "input directory : " << argv[1] << " directory" << endl;
  cerr << "output directory: " << argv[2] << " directory" << end;
  cerr << "frames: " << argv[3] << " frames" << endl;

  for(;flag;) {
    sprintf(pic_name, "%s/frame%05d.pnm", argv[1], frame_cnt);
    
    if (!Image::GetPNMImage(pic_name, frame)) {
      break;
    }  

    if (frame.empty()) {
      break;
    }
    
    sprintf(pic_name, "%s/frame%05d.pgm", argv[2], frame_cnt);
    if (!Image::SavePGMImage(pic_name, SM::CalcSaliency(frame[2], frame[1], frame[0]))) {
      break;
    }

    frame_cnt++;

    DispProgressBar(((double)frame_cnt / frame_num) * 100.0);
  }

  return 0;
}
