#ifndef __PROGRESS_HPP
#define __PROGRESS_HPP

#include <iostream>
#include <sstream>
#include <cstdlib>

using namespace std;

namespace io 
{

  static const int MAX_BAR_LEN = 50;
  
  void DispProgressBar(double percent)
  {
    string progress_bar;
    
    int bar_len = percent / (100 / MAX_BAR_LEN);

    for (int i=0;i<MAX_BAR_LEN;++i) {
      if (i < bar_len) {
	progress_bar += '=';
      }
      else if (i == bar_len) {
	progress_bar += '>';
      }
      else {
	progress_bar += ' ';
      }
    }

    printf("\r");
    cerr << "[" << progress_bar << "] " << (int)(percent + 0.5) << "%";
  }

}

#endif
