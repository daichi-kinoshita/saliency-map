#ifndef __PROGRESS_HPP
#define __PROGRESS_HPP

#include <iostream>
#include <sstream>
#include <cstdlib>

using namespace std;

void DispProgressBar(double percent)
{
  int max_bar_len = 25;
  string progress_bar;

  int bar_len = percent / (100 / max_bar_len);

  for (int i=0;i<max_bar_len;++i) {
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

#endif
