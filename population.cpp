//******************************************************************************
// population.cpp
// 2004.12.27
// Lango, Trevor M.
//
//******************************************************************************

#include "population.h"
#define BINOMIAL 0
#if BINOMIAL
#include "random.h"
#else
#include "rand.h"
#endif

#include <iostream>
using std::cout;
using std::cerr;

#include <new>

#include <cmath>
#include <unistd.h>
#include <stdlib.h>

ostream &operator<<( ostream &out,
                     const Population *somePopulation ) {

  out << "salmon matrix: \nyear\t";

  for( int i = 0; i < somePopulation->maxAgePerSalmon; ++i ) {

    out << "age" << i << '\t';

  }

  out << '\n';

  for( int i = 0; i < somePopulation->salmon[ 0 ].size( ); ++i ) {

    out << i << '\t';

    for( int j = 0; j < somePopulation->maxAgePerSalmon; ++j ) {

      out << somePopulation->salmon[ j ][ i ] << '\t';

    }

    out << '\n';

  }

  out << "maxAgePerSalmon: "
      << somePopulation->maxAgePerSalmon
      << '\n';

  out << "nestSiteCompetitionFactor: "
      << somePopulation->nestSiteCompetitionFactor
      << '\n';

  out << "maxOffspringPerSalmon: "
      << somePopulation->maxOffspringPerSalmon
      << '\n';

  return out;

}

typedef enum { NO_IMP, IMP_UP, IMP_DOWN } imp_dir_t;

static double getRandom( double mean, double variance, double min, double max,
			 imp_dir_t imp, int t )
{
  if (variance == 0)
    return mean;
//#define PER_3
#ifndef PER_3
  if (imp == IMP_UP)
    return ((t == 1)? mean + sqrt(variance) : mean);
  else if (imp == IMP_DOWN)
    return ((t == 1)? mean - sqrt(variance) : mean);
#else
  if (imp == IMP_UP)
    return ((t%3 == 1)? mean + sqrt(variance) : mean);
  else if (imp == IMP_DOWN)
    return ((t%3 == 1)? mean - sqrt(variance) : mean);
#endif
#if BINOMIAL
  if( random->getRandom( ) < 0.5 ) {
    return mean - sqrt( variance );
  }
  else {
    return mean + sqrt( variance );
  }
#else // !BINOMIAL
  double ss; 
  //    do
  //    { ss = norml_distn( s_mean, sqrt( s_variance ) );
  //    }while (ss < 0 || ss > 1);
  ss = norml_distn( mean, sqrt( variance ) );
  if (ss > max)      ss = max;
  else if (ss < min) ss = min;
  return ss;
#endif
}

double clip(double ss, double min, double max)
{
  if (ss > max)      ss = max;
  else if (ss < min) ss = min;
  return ss;
}

Population::Population( int someMaxAgePerSalmon,
                        double someNestSiteCompetitionFactor,
                        double someMaxOffspringPerSalmon,
                        double someS_mean,
                        double someS_variance,
                        double someAc_mean,
                        double someAc_variance,
                        double someSigma,
                        double r,
                        double theta,
			int T,
			bool threeD,
			bool imp,
			bool early )
: year( 0 ),
  maxAgePerSalmon( someMaxAgePerSalmon ),
  nestSiteCompetitionFactor( someNestSiteCompetitionFactor ),
  maxOffspringPerSalmon( someMaxOffspringPerSalmon ),
  do_3d( threeD ),
  do_impulse( imp ),
  early_ocean( early ),
  s_mean( someS_mean ),
  s_variance( someS_variance ),
  ac_mean( someAc_mean ),
  ac_variance( someAc_variance ),
  sigma( someSigma ),
  delta_e_mean( 0.5 * ( erfc( ( ac_mean - maxAgePerSalmon + 1.5 ) /
			      ( sigma * sqrt( 2 ) ) ) ) ),
  delta_l_mean( threeD ? 0 :
		(0.5 * ( erfc( ( maxAgePerSalmon - 0.5 - ac_mean ) /
			       ( sigma * sqrt( 2 ) ) ) )) ),
  xstar( maxAgePerSalmon, 0 ),
  jacobian( maxAgePerSalmon, maxAgePerSalmon, 0 ),
  H( maxAgePerSalmon, 0 ),
  H1( maxAgePerSalmon, 0 ),
  random( new Random( ) ) {

  double out_var;
  if ( s_variance > 0 )
    out_var = s_variance;
  else
    out_var = ac_variance;
  
  if (r != 0)
  {
    if (sin(theta) == 0)// 1 dimension of color
    {
      a = vnl_matrix< double >(1,1,0);
      b = vnl_vector< double >(1,0);
      u = vnl_vector< double >(1,0);
      a[0][0] = r;
    }
    else // 2 dimensions of color
    {
      a = vnl_matrix< double >(2,2,0);
      b = vnl_vector< double >(2,0);
      u = vnl_vector< double >(2,0);
      a[0][0] = a[1][1] = r*cos(theta);
      a[1][0] = r*sin(theta);
      a[0][1] = - a[1][0];
    }
    // either way
#if 0 // scale stationary variance
    b[0] = sqrt(out_var*(1 - r*r));
#else // scale expected variance for T years
    b[0] = sqrt( out_var*(1 - r*r)*(T - 1) /
		 (T - ( (2 + 2*r + r*r - pow(r,2*T)) /
			(1 - r*r) )
		  + ( (1 - pow(r,T)) * (1 + 2*r - pow(r,T)) /
		      (T * (1 - r) * (1 - r)) ) ) );
#endif // Wichmann et al. 2005
  }
  else // no color
  {
    a = vnl_matrix< double >(1,1,0);
    b = vnl_vector< double >(1,0);
    u = vnl_vector< double >(1,0);
    // leave a = 0
    b[0] = sqrt(out_var);
  }
  if (do_impulse)
  {
    imp_dir_t imp_dir = (s_variance>0 ? IMP_DOWN : IMP_UP);
    u[0] =  getRandom( 0, sqrt(out_var), 0, 1, imp_dir, 0 );
  }
  else
    u[0] = norml_distn(0.0, sqrt(out_var));
  
#if !BINOMIAL
  // time(0) not sufficient for we can run multiple times in a second
  sgenrand2( time(0)+getpid() );
#endif

  u_history.push_back(u);
  if (s_variance > 0)
    s.push_back(clip(s_mean + u[0],0,1));
  else
    s.push_back(s_mean);
  if (ac_variance > 0)
    ac.push_back(ac_mean + u[0]);
  else
    ac.push_back(ac_mean);

// set xstar(s):

  double c = 1.0 /
       ( ( ( delta_e_mean ) +
     ( ( 1 - delta_e_mean ) *
       ( 1 - delta_l_mean ) *
       ( s_mean ) ) +
     ( ( delta_l_mean ) *
       ( 1 - delta_e_mean ) *
       ( s_mean ) * ( s_mean ) ) ) *
         ( pow( s_mean, maxAgePerSalmon - 3 ) ) );

  xstar[0] =
    ( maxOffspringPerSalmon - c ) / nestSiteCompetitionFactor;

  for( int i = 1; i <= maxAgePerSalmon - 3; ++i )
    xstar[i] = pow( s_mean, i ) * xstar[ 0 ];

  xstar[maxAgePerSalmon - 2] =
    pow( s_mean, maxAgePerSalmon - 2 ) * ( 1 - delta_e_mean ) * xstar[ 0 ];

  xstar[maxAgePerSalmon - 1] =
    pow( s_mean, maxAgePerSalmon - 1 ) * delta_l_mean *
      ( 1 - delta_e_mean ) * xstar[ 0 ];

  // set jacobian:

  jacobian[ 0 ][ maxAgePerSalmon - 3 ] =
    delta_e_mean * c * c / maxOffspringPerSalmon;

  jacobian[ 0 ][ maxAgePerSalmon - 2 ] =
    ( 1 - delta_l_mean ) * c * c / maxOffspringPerSalmon;

  jacobian[ 0 ][ maxAgePerSalmon - 1 ] = c * c / maxOffspringPerSalmon;

  for( int i = 0; i < maxAgePerSalmon - 3; ++i )
    jacobian[ i + 1 ][ i ] = s_mean;

  jacobian[ maxAgePerSalmon - 2][ maxAgePerSalmon - 3 ] =
    ( 1 - delta_e_mean ) * s_mean;

  jacobian[ maxAgePerSalmon - 1 ][ maxAgePerSalmon - 2 ] =
    delta_l_mean * s_mean;


  double rstar = ( delta_e_mean * xstar[ maxAgePerSalmon - 3 ] +
		   ( 1.0 - delta_l_mean ) * xstar[ maxAgePerSalmon - 2 ] +
		   xstar[ maxAgePerSalmon - 1 ] );

  double dprime = ( maxOffspringPerSalmon /
		    ( ( 1.0 + nestSiteCompetitionFactor * rstar ) *
		      ( 1.0 + nestSiteCompetitionFactor * rstar ) ) );

  if ( ac_variance > 0 )
  {
    vnl_vector< double >
      H_delta_e( maxAgePerSalmon, 0 ),
      H_delta_l( maxAgePerSalmon, 0 );

    // set H_delta_e:
    H_delta_e[ 0 ] =
      xstar[ maxAgePerSalmon - 3 ] * dprime;
    H_delta_e[ maxAgePerSalmon - 2 ] =
      - s_mean * xstar[ maxAgePerSalmon - 3 ];
    
    // set H_delta_l;
    H_delta_l[ 0 ] =
      - xstar[ maxAgePerSalmon - 2 ] * dprime;
    H_delta_l[ maxAgePerSalmon - 1 ] =
      s_mean * xstar[ maxAgePerSalmon - 2 ];

    // set H_ac:
    double Pe = - ( ( 1 /  ( sqrt( 2 ) * sqrt( M_PI ) * sigma ) ) *
		    exp( - pow( ( maxAgePerSalmon - 1.5 - ac_mean ), 2 ) /
			 ( 2 * pow( sigma, 2 ) ) ) );
    double Pl = threeD ? 0 :
      ( 1 / ( sqrt( 2 ) * sqrt( M_PI ) * sigma ) ) *
      exp( - pow( ( maxAgePerSalmon - 0.5 - ac_mean ), 2 ) /
	   ( 2 * pow( sigma, 2 ) ) );
    // goes to nan
    if (sigma==0)
      Pe = Pl = 0;

#if 0 // old ac    
    H = ( Pe * H_delta_e ) + ( Pl * H_delta_l );
#else // new ac
    H = ( Pe * H_delta_e );
    H1 = ( Pl * H_delta_l );
#endif
  }
  else
  {
    // set H_s:
    H[ 0 ] = 0.0;
    H[ 1 ] = xstar[ 0 ];

    if (!early_ocean)
    {
      for( int i = 2; i <= maxAgePerSalmon - 3; ++i )
	H[ i ] = xstar[ i - 1 ];
    }
  
    if (!early_ocean || maxAgePerSalmon == 3)
      H[ maxAgePerSalmon - 2 ] =
	( 1 - delta_e_mean ) * xstar[ maxAgePerSalmon - 3 ];
    
    if (!early_ocean || maxAgePerSalmon == 2)
      H[ maxAgePerSalmon - 1 ] =
	delta_l_mean * xstar[ maxAgePerSalmon - 2 ];
  }
  
    // salmon has one entry per age class, each a time series
  salmon.resize( maxAgePerSalmon );
  // first entry, initial conditions, is xstar
  for( int i = 0; i < maxAgePerSalmon; ++i )
    salmon[ i ].push_back( xstar[ i ] );
}

Population::~Population( ) {

  delete random;

}

void Population::operator++( )
{
  ++year;

  if ( year % 1000 == 0 )
  { cout << ".";
    cout.flush();
  }

  double v;
  if (do_impulse)
  {
//    imp_dir_t imp_dir = (s_variance>0 ? IMP_DOWN : IMP_UP);
    imp_dir_t imp_dir = IMP_DOWN;
    v =  getRandom( 0, 1, 0, 1, imp_dir, year );
  }
  else
  {
    v = norml_distn( 0.0, 1.0 );
  }
  v_history.push_back( v );
  u = a*u + b*v; //gives u[t] = a*u[t-1]+b*v[t]
  u_history.push_back( u );

  if (s_variance > 0)
    s.push_back(clip(s_mean + u[0],0,1));
  else
    s.push_back(s_mean);
  if (ac_variance > 0)
    ac.push_back(ac_mean + u[0]);
  else
    ac.push_back(ac_mean);

  // set delta_e:
  double delta_e = ( 0.5 * ( erfc( ( ac[ year ] - maxAgePerSalmon + 1.5 ) /
				     ( sigma * sqrt( 2 ) ) ) ) );

  // set delta_l:
#if 0
  double delta_l =
    do_3d ? 0 : ( 0.5 * ( erfc( ( maxAgePerSalmon - 0.5 - ac[ year ] ) /
				( sigma * sqrt( 2 ) ) ) ) );
#else
  double delta_l =
    do_3d ? 0 : ( 0.5 * ( erfc( ( maxAgePerSalmon - 0.5 - ac[ year - 1 ] ) /
				( sigma * sqrt( 2 ) ) ) ) );
#endif

  double s_early = s.back();
  double s_late = early_ocean ? s_mean : s.back();
  
  //this won't work right if N <= 3

  salmon[ maxAgePerSalmon - 1 ].push_back( ( s_late *
					     delta_l *
					     salmon[ maxAgePerSalmon - 2 ]
					     [ year - 1 ] ) );

  salmon[ maxAgePerSalmon - 2 ].push_back( s_late *
					   ( 1 - delta_e ) *
					   salmon[ maxAgePerSalmon - 3 ]
					   [ year - 1 ] );

  for( int i = maxAgePerSalmon - 3; i > 1; --i ) {

    salmon[ i ].push_back( ( s_late *
                           salmon[ i - 1 ][ year - 1 ] ) );

  }

  salmon[ 1 ].push_back( ( s_early *
                           salmon[ 0 ][ year - 1 ] ) );

  double returningSalmon = ( ( delta_e *
			       salmon[ maxAgePerSalmon - 3 ][ year - 1 ] ) +
			     ( ( 1 - delta_l ) *
			       salmon[ maxAgePerSalmon - 2 ][ year - 1 ] ) +
			     salmon[ maxAgePerSalmon - 1 ][ year - 1 ] );

  salmon[ 0 ].push_back( ( maxOffspringPerSalmon * returningSalmon ) /
			 ( 1 + nestSiteCompetitionFactor * returningSalmon ) );

}

void Population::doCohortExtinction( double thresh )
{
  vector< vector<double> >::iterator class_i;
  for (class_i = salmon.begin(); class_i != salmon.end(); ++class_i)
  {
    // x_i(t) is at the end of the age class i time series
    double &xit = class_i->back();
    // wipe the cohort out if it's too small
    if (xit < thresh)
      xit = 0;
  }
}

const vector< vector< double > > &Population::getSalmon( ) const {

  return salmon;

}

const int Population::getMaxAgePerSalmon( ) const {

  return maxAgePerSalmon;

}

double Population::getNestSiteCompetitionFactor( ) const {

  return nestSiteCompetitionFactor;

}

const double Population::getMaxOffspringPerSalmon( ) const {

  return maxOffspringPerSalmon;

}

double Population::getAc_Mean( ) const {

  return ac_mean;

}

double Population::getS_Mean( ) const {

  return s_mean;

}

vnl_matrix< double > Population::get_a( ) const {

  return a;

}

vcl_complex< double > Population::get_a_as_complex( ) const
{
  if (a.rows() == 1)
    return a[0][0];
  else if (a.rows() == 2)
    return vcl_complex<double>(a[0][0], a[1][0]);
  //else
  cerr << "A has no complex counterpart\n";
  exit(-1);
}

vnl_vector< double > Population::get_b( ) const {

  return b;

}

double Population::getAcVariance( ) const {

  return ac_variance;

}

double Population::getSVariance( ) const {

  return s_variance;

}

double Population::getSigma( ) const
{
  return sigma;
}

double Population::getMeanDelta_e( ) const
{
  return delta_e_mean;
}

double Population::getMeanDelta_l( ) const
{
  return delta_l_mean;
}

vector< double > &Population::getAc( ) {

  return ac;

}

vector< double > &Population::getS( ) {

  return s;

}

const vnl_matrix< double > &Population::getJacobian( ) const
{
  return jacobian;
}

const vnl_vector< double > &Population::getH( ) const
{
  return H;
}

const vnl_vector< double > &Population::getH1( ) const
{
  return H1;
}

double Population::getXStar( const int someAgeClass ) const
{
  return xstar[ someAgeClass ];
}

const vnl_vector< double> &Population::getXStar() const
{
  return xstar;
}

const vector< vnl_vector< double > > &Population::get_u_history( ) const
{
  return u_history;
}

vector< double > Population::getAgeClass( const int someAgeClass ) const
{
  return ( salmon[ someAgeClass ] );
}

// vector< vcl_complex< double > > Population::getYear( const int someYear ) const
// {
//   vector< vcl_complex< double > > year;
//   for( int i = 0; i < maxAgePerSalmon; ++i ) {
//     year.push_back( salmon[ i ][ someYear ] );
//   }
//   return year;
// }

vnl_vector< double > Population::getYear( const int someYear ) const
{
  vnl_vector< double > y(maxAgePerSalmon);
  for( int i = 0; i < maxAgePerSalmon; ++i ) {
    y[ i ] = salmon[ i ][ someYear ];
  }
  return y;
}

int Population::getCurrentYear( ) const
{
  return year;
}

bool Population::doColor( ) const
{
  return ! a.is_zero( );
}

bool Population::doThreeD( ) const
{
  return do_3d;
}

