//******************************************************************************
// random.h
// 2005.01.03
// Lango, Trevor M.
//
//******************************************************************************

#ifndef RANDOM_H
#define RANDOM_H

#include <iostream>
using std::ostream;

class Random;
class Random {

  friend ostream &operator<<( ostream &,
                              const Random * );

 private:
  double table[ 100 ][ 100 ];
  int m;
  int n;

 protected:

 public:
  Random( );
  ~Random( );

  double getRandom( );

};

#endif
