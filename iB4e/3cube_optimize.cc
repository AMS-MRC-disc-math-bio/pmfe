#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>

#ifdef GMP
  #include <gmpxx.h>
#endif

using namespace std;

int main(const int argc, const char * argv[]) 
{
  int a;

  for(int i = 0; i < 3; i++){ 
    cin >> a;
    if(a > 0) cout << 1; else cout << 0;
    cout << " ";
  }

  return 0;
}

