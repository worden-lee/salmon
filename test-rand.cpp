#include "rand.h"
#include <time.h>
#include <math.h>
#include <iostream>

int main()
{
  const double variance = 0.1; //1e-7;
  sgenrand2(time(0));
  for (int i = 0; i < 10000; ++i)
    cout << norml_distn(1,sqrt(variance)) - 1 << endl;
}
