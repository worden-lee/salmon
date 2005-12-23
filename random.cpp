//******************************************************************************
// random.cpp
// 2005.01.03
// Lango, Trevor M.
//
//******************************************************************************

#include "random.h"

#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

ostream &operator<<( ostream &out,
                     const Random *someRandom ) {

  return out;

}

Random::Random( ) {

  srand( time( NULL ) + getpid( ) );

  for( int i = 0; i < 100; ++i ) {

    for( int j = 0; j < 100; ++j ) {

      table[ i ][ j ] = ( ( double )( rand( ) % RAND_MAX ) / RAND_MAX );

    }

  }

}

Random::~Random( ) {
}

double Random::getRandom( ) {

  m = ( rand( ) % 99 );
  n = ( rand( ) % 99 );

  double r = table[ m ][ n ];

  table[ m ][ n ] = ( ( double )( rand( ) % RAND_MAX ) / RAND_MAX );

  return ( r );

}

