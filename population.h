//******************************************************************************
// population.h
// 2004.12.27
// Lango, Trevor M.
//
//******************************************************************************

#ifndef POPULATION_H
#define POPULATION_H

#include "random.h"

#include "vnl_matrix.h"
#include "vnl_vector.h"

#include <iostream>
using std::ostream;

#include <fstream>
using std::fstream;

#include <vector>
using std::vector;

/*
#include <complex>
using std::complex;
*/

#include <new>

class Population;
class Population {

  friend ostream &operator<<( ostream &,
                              const Population * );

 private:

  // for bookkeeping convenience
  int year;

  // 'salmon' is a matrix where:
  // - rows represent years
  // - columns represent age classes
  vector< vector< double > > salmon;

  // 'n' is the max age of a salmon:
  const int maxAgePerSalmon;

  // 'beta' is the nest site competition intensity factor
  const double nestSiteCompetitionFactor;

  // 'alpha' is the maximum number of offspring for a given salmon:
  const double maxOffspringPerSalmon;

  // special case if set, no delta_l
  bool do_3d;

  // special case if set, just do impulsive noise
  bool do_impulse;

  // special case if set, only early ocean mortality varies
  bool early_ocean;

  // statistics:
  double s_mean;
  double s_variance;
  vector< double > s;

  double ac_mean;
  double ac_variance;
  vector< double > ac;

  double sigma;

  double delta_e_mean;
  double delta_l_mean;
  
  //double a; //u(t)=A*u(t-1)+B*v(t)
  //double b; //variance=(b^2)/(1-a^2)

  //now we make a and b matrices.
  bool do_color;
  vnl_matrix< double > a;
  vnl_vector< double > b;
  vnl_vector< double > u;

  vector< vnl_vector< double > > u_history;
  vector< double > v_history;

  vnl_vector< double > xstar;
  vnl_matrix< double > jacobian;
/*   vnl_vector< double > H_delta_e; */
/*   vnl_vector< double > H_delta_l; */
/*   vnl_vector< double > H_ac; */
/*   vnl_vector< double > H_s; */
  vnl_vector< double > H;
  vnl_vector< double > H1;

  Random *random;

 protected:

 public:

  Population( int n = 0,
              double beta = 0.0,
              double alpha = 0,
              double s_mean = 0.0,
              double s_variance = 0.0,
	      double ac_mean = 0.0,
	      double ac_variance = 0.0,
	      double sigma = 0.0,
              double r = 0.0,
              double theta = 0.0,
	      int T = 0,
	      bool threeD = false,
	      bool imp = false,
	      bool earlyOcean = false );
  ~Population( );

  // overloaded operator(s):
  void operator++( );

  // function functions:
  void doCohortExtinction( double thresh );

  // set functions:

  // get functions:
  const vector< vector< double > > &getSalmon( ) const;
  const int getMaxAgePerSalmon( ) const;
  double getNestSiteCompetitionFactor( ) const;
  const double getMaxOffspringPerSalmon( ) const;

  double getAc_Mean( ) const;
  double getS_Mean( ) const;

  double getAcVariance( ) const;
  double getSVariance( ) const;

  double getSigma( ) const;

  vnl_matrix< double > get_a( ) const;
  vnl_vector< double > get_b( ) const;

  vcl_complex< double > get_a_as_complex( ) const;
  
  vector< double > &getAc( );
  vector< double > &getS( );
  double getXStar( const int ) const;
  const vnl_vector< double> &getXStar() const;

  const vnl_matrix< double > &getJacobian( ) const;
  const vnl_vector< double > &getH( ) const;
  const vnl_vector< double > &getH1( ) const; // this is often 0

  double getMeanDelta_e( ) const;
  double getMeanDelta_l( ) const;

  const vector< vnl_vector< double > > &get_u_history( ) const;

  // helper functions:
  vector< double > getAgeClass( const int ) const;
  //  vector< vcl_complex< double > > getYear( const int = year ) const;
  vnl_vector< double > getYear( const int ) const;

  int getCurrentYear( ) const;

  // queries
  bool doColor( ) const;
  bool doThreeD( ) const;
  bool doImpulse( ) const;
};

#endif
