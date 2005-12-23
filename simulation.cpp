//******************************************************************************
// simulation.cpp
// 2005-01-03
// Lango, Trevor M.
// Worden, Lee
//
//******************************************************************************

#include "simulation.h"

#include "vnl_real_eigensystem.h"
#include "vnl_svd.h"
#include "vnl_matrix_inverse.h"
#include "vnl_complexify.h"
#include "vnl_fft_1d.h"
#include "vnl_determinant.h"
#include "vnl_math.h"
#include "vcl_complex.h"
#include "vnl_vector.h"

#include <complex>
#include <cmath>

#include <new>

#include <iostream>
using std::ofstream;
using std::fstream;
using std::ios;
using std::cout;
using std::cerr;
using std::endl;
using std::form;

#ifdef NEED_STL_FUNCTION
#include <stl_function.h>
#endif
#include <algorithm>
using std::sort;
using std::binary_function;
using std::copy;

// function for whether a complex number is actually positive real
template<typename T>
bool is_positive(vcl_complex<T> c)
{ return imag(c)==0 && real(c)>0; }

// little function object, does the compare for sorting eigenvalues by size
struct compareEigenvalues : public binary_function< int, int, bool > {
  vnl_diag_matrix< vcl_complex< double > > L;
  // constructor
  compareEigenvalues( vnl_diag_matrix< vcl_complex< double > > &LL )
    : L( LL ) {}
  // function call, used like (a < b)
  bool operator()( int a, int b ) const {
    return ( abs(L[a]) > abs(L[b]) ||
	     ((abs(L[a])==abs(L[b]))&&(arg(L[a])<arg(L[b]))) ||
// brute force: positive eigenvalue first!
	     (is_positive(L[a]) && !is_positive(L[b])) );
  }
};

// convert a vector<double> to vnl_vector<double>
template<class T>
vnl_vector<T> vnlize_vector(const vector<T> &before)
{
  vnl_vector<T> after(before.size());
  copy(before.begin(), before.end(), after.begin());
  return after;
}

// complex conjugate a vector. why isn't this in the library?
template<class T>
vnl_vector<T> conjugate(const vnl_vector<T> &source)
{
  vnl_vector<T> dest(source);
  vnl_c_vector<T>::conjugate(source.begin(),
                             dest.begin(),
                             source.size());  // size of block
  return dest;
}

// interspecies multiply (complex matrix * real vector)
template<class T>
vnl_vector< vcl_complex<T> > operator*( const vnl_matrix< vcl_complex<T> > &M,
					const vnl_vector< T > &V )
{
  return M * vnl_complexify(V);
}

ostream &operator<<( ostream &out,
                     const Simulation *someSimulation ) 
{
  out << "year\t";

  for( int i = 0; i < someSimulation->population->getMaxAgePerSalmon( ); ++i )
  {
    out << "age" << i << '\t';
  }
  out << '\n';

  for( int i = 0; i < someSimulation->maxTime; ++i )
  {
    out << i << '\t';
    for( int j = 0; j < someSimulation->population->getYear( i ).size( ); ++j )
    {
      out << ( someSimulation->population->getYear( i ) )[ j ] << '\t';
    }
    out << '\n';
  }

  return out;
}

Simulation::Simulation( const int someMaxTime,
                        Population *somePopulation,
                        char *some_output_quantity,
			const string someDirectory )
: maxTime( someMaxTime ),
  population( somePopulation ),
  output_quantity( some_output_quantity ),
  output_weights( somePopulation->getMaxAgePerSalmon( ) ),
  directory( someDirectory ),
  do_extinction( false )
{
  // set output_weights:
  switch( output_quantity[ 0 ] ) {

  case 'r': // recruits
    output_weights.fill( 0 );
    output_weights[ 0 ] = 1;
    break;

  case 't': // total population
    output_weights.fill( 1 );
    //output_weights.normalize();
    break;

  case 'd': // difference between this year's recruits and last year's
    output_weights.fill( 0 );
    output_weights[ 0 ] = population->getS_Mean();
    output_weights[ 1 ] = 1;
    output_weights.normalize();
    break;

  case 'c': // catch
    output_weights.fill( 0 );
    output_weights[ population->getMaxAgePerSalmon() - 3 ] =
      population->getMeanDelta_e();
    output_weights[ population->getMaxAgePerSalmon() - 2 ] =
      1 - population->getMeanDelta_l();
    output_weights[ population->getMaxAgePerSalmon() - 1 ] = 1;
    break;

  case 'x': // direction of x*, akin to total
    output_weights = population->getXStar();
    output_weights.normalize();
    break;
    
  case 'l': // last age class
    output_weights.fill( 0 );
    output_weights[ population->getMaxAgePerSalmon() - 1 ] = 1;
    break;
    
  default:
    break;
  }
  //cout << "weights " << output_weights << endl;
}

Simulation::~Simulation( )
{
}

void Simulation::doExtinction(double someThreshold, bool cohort_whether)
{
  do_extinction = true;
  cohort_extinction = cohort_whether;
  threshold = someThreshold;
}

const int Simulation::getMaxTime( ) const
{
  return maxTime;
}

Population *Simulation::getPopulation( ) const
{
  return population;
}

const vnl_vector< double > &Simulation::getOutput_Weights( )
{
  return output_weights;
}

void Simulation::run( )
{
  cout << "generating population versus time.";

  //extinction_time = -1;
  extinction_time = maxTime;

  //for( int i = 0; i < maxTime; ++i ) {
  // year 0 is initial conditions
  for( int i = 1; i < maxTime; ++i )
  {
    // generate salmon population versus time
    ++( *population );

    // check against extinction threshold
    if (do_extinction)
    {
      bool extinct;
      if (cohort_extinction)
      {
	population->doCohortExtinction( threshold );
	extinct = (population->getYear( population->getCurrentYear() ).sum()
		   == 0);
      }
      else
      {
	double quantity =
	  dot_product( output_weights,
		       population->getYear( population->getCurrentYear( ) ) );
	extinct = (quantity <= threshold);
      }
      if (extinct)
      {
	extinction_time = population->getCurrentYear();
	break;
      }
    }
  }

  if (do_extinction)
  {
    if (extinction_time >= 0 && extinction_time < maxTime)
      cout << "extinct at " << extinction_time << endl;
    else
      cout << "no extinction by " << population->getCurrentYear() + 1 << "\n";
  }
  else
    cout << "done!\n";
}

void Simulation::report( )
{
  /****** first set up linear algebra data objects ******/

  int nYears = population->getCurrentYear( ) + 1;

  const vnl_matrix< double > &J = population->getJacobian( );

  //construct M and K matrices for colored noise.

  //if no color set M=J and K=H sans complexify.

  vnl_matrix< double > M;
  vnl_vector< double > K;
  vnl_vector< double > q_hat;
  vnl_vector< double > q = output_weights;

  if ( population->doColor( ) )
  {
    //if we're doing colored noise, construct M and K bigger than J and H.
    const vnl_matrix< double > &A = population->get_a();
    const vnl_vector< double > &B = population->get_b();
    int M_dim =  J.rows() + A.rows();
     
    vnl_matrix< double > H_matrix( J.rows(), A.rows(), 0 );
    H_matrix.set_column( 0, population->getH() );

    // LW messing with this here, see if it works
    vnl_matrix< double > H1( J.rows(), A.rows(), 0 );
    H1.set_column( 0, population->getH1() );

    vnl_matrix< double > HA = H_matrix * A;
    vnl_vector< double > HB = H_matrix * B;

    //cout << "Here is H and H1:\n" << H_matrix << H1 << endl;
    
    vnl_matrix< double  > M_temp( M_dim, M_dim, 0 );
    M_temp.update( J );
    // LW messing with this here, see if it works
    M_temp.update( HA + H1, 0, J.cols() );
    M_temp.update( A, J.rows(), J.cols() );

    vnl_vector< double > K_temp( M_dim, 0 );
    K_temp.update( HB );
    K_temp.update( B, J.rows() );

    M = M_temp;
    K = K_temp;
    vnl_vector< double > q_hat_tmp( M_dim, 0 );
    q_hat_tmp.update( q );
    q_hat = q_hat_tmp;

  }
  else 
  {
    //if no color, i.e. no -c or -cp option, then M=J and K=H.
    M = J;
    K = population->getH();
    q_hat = q;
    
  }

  vnl_real_eigensystem J_eigensystem( J );
  vnl_diag_matrix< vcl_complex< double > > &J_Lambda = J_eigensystem.D;
  vnl_matrix< vcl_complex< double > > J_U = J_eigensystem.V;
  vnl_diag_matrix< vcl_complex<double> > J_sgn( J_U.cols() );
  for ( int i=0; i<J_U.cols(); ++i ) 
    J_sgn[i] = 1.0/(J_U[0][i]/abs(J_U[0][i]));
  J_U = J_U * J_sgn;
  J_U.normalize_columns();
  vnl_matrix< vcl_complex< double > > J_V( J_U.rows(), J_U.cols(), 0 );

  if (!J_U.is_finite())
  { cerr << "\nJ_U not finite! fudging J_V.\n";
  // leave J_V = 0, this screws things up, but they're screwed anyway.
  }
  else
  {
    vnl_matrix_inverse< vcl_complex<double> > J_Uinv( J_U );
    J_V = J_Uinv;
  }
  vnl_vector< vcl_complex< double > > J_H( vnl_complexify( population->getH() ) );
  vnl_vector< vcl_complex<double> > J_VH = J_V * J_H;

  vnl_vector< vcl_complex< double > > J_qU = vnl_complexify(q) * J_U;

  vector< int > J_indices( J_Lambda.size() );
  for( int i = 0; i < J_indices.size( ); ++i )
    J_indices[ i ] = i;
  sort( J_indices.begin(), J_indices.end(), compareEigenvalues( J_Lambda ) );

  //cout << "Here is M:\n" << M << endl;
  
  //now do the eigensystem stuff for M and K.
  vnl_real_eigensystem eigensystem( M );

  vnl_matrix< vcl_complex< double > > U = eigensystem.V;

  vnl_diag_matrix< vcl_complex<double> > sgn( U.cols() );
  for (int i = 0; i < U.cols(); ++i)
    sgn[i] = 1.0 / (U[0][i] / abs(U[0][i]));
  U = U * sgn;
  U.normalize_columns();
  
  vnl_matrix< vcl_complex< double > > V( U.rows(), U.cols(), 0 );
  if (!U.is_finite())
  { cerr << "\nU not finite! fudging V.\n";
  // leave V = 0, this screws things up, but they're screwed anyway.
  }
  else
  {
    vnl_matrix_inverse< vcl_complex< double > > Uinv( U );
    V = Uinv;
  }

  vnl_diag_matrix< vcl_complex< double > > &Lambda = eigensystem.D;

  vnl_vector< vcl_complex< double > > H = vnl_complexify( K );

  vnl_vector< vcl_complex< double > > VH = V * H;
  
  vnl_vector< vcl_complex< double > > qU = vnl_complexify(q_hat) * U;

  vector< int > indices( Lambda.size() );
  for( int i = 0; i < indices.size( ); ++i )
    indices[ i ] = i;
  sort( indices.begin(), indices.end(), compareEigenvalues( Lambda ) );

  // this changes a few things
  bool do_3d = population->doThreeD();

  /************** now begin reporting **************************/
  
  if (!do_extinction)
  {
    cout << "\ndoing observed transfer";
    observed_transfer( );
    cout << "done!\n\ndoing analytic transfer";
    vnl_vector< vcl_complex< double > > VH1 =
      J_V * vnl_complexify(population->getH1());
    analytic_transfer( nYears, J_qU, J_Lambda,
		       J_VH, VH1, indices );
    cout << "done!\n";
    cout << "\ndoing other output files";
  }

  if (do_extinction)
  {
    // possibly change this variance in case of color
    double input_variance
      = ( population->getSVariance( ) > 0 ?
	  population->getSVariance( ) : population->getAcVariance( ) );
    analytic_extinction_time( q_hat, U, Lambda, VH, indices,
			      input_variance, population->get_b()[0],
			      population->get_a_as_complex(), threshold );
    //analytic_extinction_time( q, J_U, J_Lambda, J_VH, J_indices,
    //                          input_variance, threshold );

    ofstream ext( (directory + "/extinction-time.out").c_str( ) );
    if (extinction_time >= 0)
    {
      ext << extinction_time << endl;
    }
    cout << "\ndoing output files";
  }
  
  { // output the noise signal
    ofstream u_file( (directory + "/u.out").c_str() );
    // output age class population values versus time to file:
    ofstream history_file( ( directory + "/population.out" ).c_str( ) );
    // output summary quantity of population versus time as well:
    ofstream qfile( ( directory + "/output-quantity.out" ).c_str( ) );

/*
    history_file << "year\t";

    for( int i = 0; i < population->getMaxAgePerSalmon( ); ++i ) {
      history_file << "age" << i << '\t';
    }
    history_file << '\n';
*/
    const vector< vnl_vector< double > > &u = population->get_u_history( );
    for( int i = 0; i <= population->getCurrentYear( ); ++i ) 
    {
      u_file << i;
      for( int j = 0; j < u[i].size( ); ++j )
	u_file << ' ' << u[i][j];
      u_file << '\n' ;

      const vnl_vector< double > &pop = population->getYear( i );
      //cout << '.';
      history_file << i << '\t';
      for( int j = 0; j < pop.size( ); ++j ) {
	//cout << '.';
#if GNU_OSTREAM
	history_file.form("%.12g\t", pop[ j ]);
#else
	history_file << pop[ j ] << '\t';
#endif
      }
      //cout << '.';
      history_file << '\n';

      qfile << dot_product( output_weights, pop ) << '\n';
      //qfile.form("%.12g\n", q);
      if ( ( i % 1000 ) == 0 )
      {
	cout << ".";
	cout.flush();
      }
    }
  }
  
  {
    ofstream deltas( ( directory + "/deltas.out" ).c_str( ) );
    deltas << population->getMeanDelta_e( ) << '\t'
	   << population->getMeanDelta_l( ) << '\n';
  }

  // write out the eigensystem, to check 
  
  {
    ofstream eigensystem( ( directory + "/eigensystem.txt" ).c_str( ) );
    eigensystem << "xstar:";
    for (int i = 0; i < J_indices.size(); ++i)
      eigensystem << population->getXStar(i) << ' ';
    eigensystem << "\n\nM:\n" << M;
    eigensystem << "\n\nK:\n" << K;
    eigensystem << "\n\nJ:\n" << J;
    eigensystem << "\n\nindices:";
    for (int i = 0; i < indices.size(); ++i)
      eigensystem << ' ' << indices[i];
    //eigensystem << "\n\nU before:\n" << U_before;
    eigensystem << "\n\nU:\n" << U;
    //eigensystem << "\n\nabs(U before):\n" << U_before.apply(cabs);
    //eigensystem << "\n\nabs(U):\n" << U.apply(cabs);
    eigensystem << "\nLambda:\n" << Lambda;
    
    eigensystem << "\n\nv_0 H u_0:\n" << (VH[0]*U.get_column(0)) << "\n\n";

    eigensystem << "\n\nv_1 H u_1 + v_2 H u_2:\n"
		<< (VH[1]*U.get_column(1) + VH[2]*U.get_column(2)) << "\n\n";

    /* Lee says pointless code.
    {
      vnl_vector< vcl_complex<double> > u0 = U.get_column(0);
      const vnl_vector< vcl_complex<double> > &xstar =
	vnl_complexify(population->getXStar());
      vcl_complex<double> angle = 
	dot_product(u0,xstar) / (u0.two_norm()*xstar.two_norm());
      eigensystem << "cos(angle(u0,xstar)) = "
		  << real(angle) << " = " << real(angle/M_PI) << " pi"
		  << "\n\n";
    }
    */
    
#define VAR(i,j) \
	 ( qU[ indices[ i ] ] * VH[ indices[ i ] ] * \
	   conj( qU[ indices[ j ] ] * VH[ indices[ j ] ] ) / \
	   ( 1.0 - Lambda[ indices[ i ] ] * conj( Lambda[ indices[ j ] ] ) ) )
    for ( int i = 0; i < H.size(); ++i )
    { if (imag(Lambda[ indices[ i ] ]) == 0)
      {
	eigensystem << "variance in q in subsystem " << indices[i]
		    << ": " << VAR(i,i) << endl;
      }
      if (i < Lambda.size() - 1 &&
	       (Lambda[ indices[ i ] ] == conj(Lambda[ indices[ i + 1 ] ]) ||
		i == 3)) // 3,4 special case
      {
	eigensystem << "variance in q in subsystem ("
		    << indices[ i ] << "," << indices[ i+1 ] << "): "
		    << VAR(i,i) + VAR(i,i+1) + VAR(i+1,i) + VAR(i+1,i+1)
		    << endl;
      }
    }
  } 

  
  {
    ofstream eigenvector;

    for( int i = 0; i < indices.size( ); ++i ) {

      eigenvector.open( ( directory + "/u" + ( char )( i + '0' ) + ".mag.out" ).c_str( ) );

      for( int j = 0; j < indices.size( ); ++j ) {
	eigenvector << abs( U.get( j, indices[ i ] ) ) << '\t';
//    *eigenvector << abs(real( U.get( j, indices[ i ] ) ) ) << '\t';
      }
      eigenvector << '\n';
      eigenvector.close( );

      eigenvector.open( ( directory + "/v" + ( char )( i + '0' ) + ".mag.out" ).c_str( ) );

      for( int j = 0; j < indices.size( ); ++j ) {
	eigenvector << abs( V.get( indices[ i ], j ) ) << '\t';
//    *eigenvector << abs(real( V.get( indices[ i ], j ) ) ) << '\t';
      }
      eigenvector << '\n';
      eigenvector.close( );
    }
  }
  
  
  { // this is for table{tbl:eigenvectors-coho-chinook} of paper
    ofstream eigenvector( (directory+"/eigenvectors.tex").c_str() );

    // DO_PM=1 means count eigenvalues with im >= 0 only
#define DO_PM 1
    vector<int> toprint;
    bool skip_next = false;
    for( int i = 0; i < indices.size( ); ++i )
    {
      if (skip_next)
      { skip_next = false; continue; }
      if (do_3d && i == indices.size() - 1)
	break;
      toprint.push_back(i);
      for( int j = 0; j < indices.size( ); ++j )
      {
	vcl_complex<double> uij = U.get( j, indices[ i ] );
	if ( fabs(imag(uij)) > 0.0001*fabs(real(uij)) && DO_PM )
	  skip_next = true;
      }
    }

    eigenvector << "\\begin{tabular}{";
    for( int i0 = 0; i0 < toprint.size( ); ++i0 )
      eigenvector << 'l';
    eigenvector << "}\n";

    for( int i0 = 0; i0 < toprint.size( ); ++i0 )
    {
      int i = toprint[i0];
      if (i0 < toprint.size() - 1 && toprint[i0+1] != i+1)
	eigenvector << (i>0 ? " & ":"") << "$u_{" << i << ',' << i+1 << "}$";
      else
	eigenvector << (i>0 ? " & ":"") << "$u_{" << i << "}$";
    }
    eigenvector << " \\\\\n";

    for( int j = 0; j < indices.size( ); ++j )
    {
      if (do_3d && j == indices.size( ) - 1)
	break;
      for( int i0 = 0; i0 < toprint.size( ); ++i0 )
      {
	if (i0>0)
	  eigenvector << " & ";
	vcl_complex<double> uij = U.get( j, indices[ toprint[ i0 ] ] );
#if GNU_OSTREAM
	const char * const form_big = "%6.4g", * const form_small = "%4.9f";
	double small_value = 1e-4;
	if ( fabs(imag(uij)) <= 0.0001*fabs(real(uij)) )
	{
	  eigenvector << '$';
	  if (fabs(real(uij)) >= small_value)
	     eigenvector.form(form_big, real(uij));
	  else
	    eigenvector.form(form_small, real(uij));
	  eigenvector << '$';
	}
	else if (DO_PM)
	{
	  eigenvector << '$';
	  if (fabs(real(uij)) >= small_value)
	    eigenvector.form(form_big, real(uij));
	  else
	    eigenvector.form(form_small, real(uij));
	  if (imag(uij)>0)
	    eigenvector << "\\pm";
	  else
	    eigenvector << "\\mp";
	  if (fabs(imag(uij)) >= small_value)
	    eigenvector.form(form_big, fabs(imag(uij)));
	  else
	    eigenvector.form(form_small, fabs(imag(uij)));
	  eigenvector << "i$";
	}
	else
	{
	  eigenvector << '$';
	  if (fabs(real(uij)) >= small_value)
	    eigenvector.form(form_big, real(uij));
	  else
	    eigenvector.form(form_small, real(uij));
	  if (imag(uij)>0)
	    eigenvector << '+';
	  if (fabs(imag(uij)) >= small_value)
	    eigenvector.form(form_big, imag(uij));
	  else
	    eigenvector.form(form_small, imag(uij));
	  eigenvector << "i$";
	}
#else
	if ( fabs(imag(uij)) <= 0.0001*fabs(real(uij)) )
	{
	  eigenvector << '$' << real(uij) << '$';
	}
	else if (DO_PM)
	{
	  eigenvector << '$' << real(uij)
		      << (imag(uij)>0 ? "\\pm":"\\mp")
		      << fabs(imag(uij)) << "i$";
	}
	else
	{
	  eigenvector << '$' << real(uij)
		      << (imag(uij)>0 ? "+":"")
		      << imag(uij) << "i$";
	}
#endif
      }
      if (do_3d ? (j < indices.size() - 2) : (j < indices.size() - 1))
	eigenvector << " \\\\";
      eigenvector << "\n";
    }
    eigenvector << "\\end{tabular}\n";
  }
  
if (0)
  { // output the 'modes' of population dynamics
    // these are defined using the eigenstructure of J not M,
    //  maybe should change
#undef DO_PM
#define DO_PM 1
    vector<int> toprint;
    bool skip_next = false;
    for( int i = 0; i < J_indices.size( ); ++i )
    {
      if (skip_next)
      { skip_next = false; continue; }
      toprint.push_back(i);
      for( int j = 0; j < J_indices.size( ); ++j )
      {
	vcl_complex<double> uij = J_U.get( j, J_indices[ i ] );
	if ( fabs(imag(uij)) > 0.0001*fabs(real(uij)) && DO_PM )
	  skip_next = true;
      }
    }

    const vnl_vector< double> &xstar = population->getXStar();

    vector<ofstream*> modes(J_U.cols());
    vector<ofstream*> rawmodes(J_U.cols());
    vector< vnl_vector< vcl_complex< double > > > u_columns(J_U.cols());
    for (vector<int>::iterator ii=toprint.begin(); ii!=toprint.end(); ++ii)
    {
      int i = *ii;
      modes[i] = new ofstream( (directory+"/mode"+char(i+'0')+".out").c_str() );
      rawmodes[i] =
	new ofstream( (directory+"/raw-mode"+char(i+'0')+".out").c_str() );
      u_columns[i] = J_U.get_column(i);
    }
    for (int t=0; t < population->getCurrentYear(); ++t)
    {
      if ( ( t % 1000 ) == 0 )
      {
	cout << ".";
	cout.flush();
      }
      const vnl_vector< double > y = population->getYear(t) - xstar;
      const vnl_vector< vcl_complex< double > > w = J_V*y;
      for (vector<int>::iterator ii=toprint.begin(); ii!=toprint.end(); ++ii)
      {
	int i = *ii;
        vnl_vector< vcl_complex< double > > modei = w[i]*u_columns[i];
	if ( fabs(imag(Lambda[i])) > 0.0001*fabs(real(Lambda[i])) )
	  //modei += conjugate(modei);
	  modei *= 2; // we'll just use real part
	// at this point modei should always be real
	*rawmodes[i] << t;
	for (int j = 0; j < modei.size(); ++j)
	  *rawmodes[i] << ' ' << modei[j].real();
	*rawmodes[i] << endl;

	modei += vnl_complexify(xstar);
	*modes[i] << t;
	for (int j = 0; j < modei.size(); ++j)
	  *modes[i] << ' ' << modei[j].real();
	*modes[i] << endl;
      }
    } // should destruct the ofstream's here
  }

  if(population->getSigma() > 0)
  { // this is for table{tbl:variance-components-coho-chinook} of paper
    ofstream components( (directory+"/variance-components.tex").c_str() );

    // for formatting we must include 4 columns:
    //  S_0, S_12, S_3, S_4
    // even if S_4 doesn't exist (as in the coho case)
    vcl_complex<double> entries[4] = 
      { VAR(0,0),
	VAR(1,1) + VAR(1,2) + VAR(2,1) + VAR(2,2),
	VAR(3,3),
      };
    if (indices.size() > 4)
    {
      entries[3] = VAR(4,4);
      //entries[4] = VAR(3,3) + VAR(3,4) + VAR(4,3) + VAR(4,4);
    }
    int i;
    for (i = 0; i < (indices.size()==5 ? 4 : 3); ++i)
    {
      if (i > 0) components << " \\\\ ";
      vcl_complex<double> c = entries[i];
      if (imag(entries[i]) != 0)
	cout << "alert: variance component " << c << " is not pure real!\n";
      components << real(c);
    }
    for ( ; i < 4; ++i )
      components << " \\\\ \\ ";
    components << '\n';
  }

  
  {
    ofstream eigenvalues( ( directory + "/eigenvalues-mag.out" ).c_str( ) );
    
    /*
    for( int i = 0; i < population->getMaxAgePerSalmon( ); ++i ) {
      eigenvalues << eigenSystem.D( i, i ).real( ) << '\t'
      << eigenSystem.D( i, i ).imag( ) << '\n'; 
    }
    */
    
    for( int i = 0; i < indices.size( ); ++i )
      eigenvalues << abs( Lambda( indices[ i ], indices[ i ] ) ) << '\t';
    eigenvalues << '\n';
    eigenvalues.close( );
    eigenvalues.open( ( directory + "/eigenvalues.out" ).c_str( ) );
    for( int i = 0; i < indices.size( ); ++i ) {
      if (do_3d && i == indices.size( ) - 1)
	break;
//      eigenvalues << real( Lambda[ indices[ i ] ] ) << '\t'
//		  << imag( Lambda[ indices[ i ] ] ) << '\n';
      eigenvalues << Lambda[ indices[ i ] ] <<  ' ';
    }
    eigenvalues << '\n';
    eigenvalues.close( );
    eigenvalues.open( ( directory + "/frequencies.out" ).c_str( ) );
    for( int i = 0; i < indices.size( ); ++i ) {
      if (do_3d && i == indices.size( ) - 1)
	break;
//      eigenvalues << real( Lambda[ indices[ i ] ] ) << '\t'
//		  << imag( Lambda[ indices[ i ] ] ) << '\n';
      eigenvalues << 2*M_PI / vcl_arg( Lambda[ indices[ i ] ] ) <<  ' ';
    }
    eigenvalues << '\n';
  } 

  
  {
    ofstream res( (directory + "/resonances.out").c_str( ) );
    ofstream resmag( (directory + "/resonances.mag.out").c_str( ) );
    for( int i = 0; i < population->getMaxAgePerSalmon( ); ++i ) {
      if (i > 0) {
	resmag << ' ';
	res << ' ';
      }
      vcl_complex<double> r =
	1.0 / ( 1.0 - Lambda[ indices[ i ] ] * conj( Lambda[ indices[ i ] ] ) );
      res << r;
      resmag << abs(r);
    }
    resmag << '\n';
    res << '\n';
  }
 
  
  {
    ofstream h( ( directory + "/forcing.out" ).c_str( ) );
    ofstream vh( ( directory + "/transformed-forcing.out" ).c_str( ) );
    ofstream vhmag( ( directory + "/transformed-forcing.mag.out" ).c_str( ) );

    for( int i = 0; i < H.size( ); ++i ) {
      if (i>0) { h << ' '; vh << ' '; vhmag << ' '; } 
      h << H[ i ].real( ); // assume H is real
      vh << VH[ indices[ i ] ];
      vhmag << abs( VH[ indices[ i ] ] );
    }
    h << '\n'; vh << '\n'; vhmag << '\n';
  }
  
  
  {
    ofstream qu( (directory + "/transformed-weights.out").c_str( ) );
    ofstream qumag( (directory + "/transformed-weights.mag.out").c_str( ) );
    for( int i = 0; i < qU.size( ); ++i ) {
      if (i>0) { qu << ' '; qumag << ' '; } 
      qu << qU[ indices[ i ] ];
      qumag << abs( qU[ indices[ i ] ] );
    }
    qu << '\n'; qumag << '\n';
  }
  
  
  {
    ofstream mf( (directory + "/mifi.out").c_str( ) );
    ofstream mfmag( (directory + "/mifi.mag.out").c_str( ) );
    // mi, fi, resonance = size of peak of xfer fn due to only the
    // most relevant mode
    ofstream mfr( (directory + "/mfr.out").c_str( ) );
    ofstream mfrmag( (directory + "/mfr.mag.out").c_str( ) );
    for( int i = 0; i < qU.size( ); ++i ) {
      if (i>0) { mf << ' '; mfmag << ' ';  mfr << ' '; mfrmag << ' ';} 
      vcl_complex< double > val = qU[ indices[ i ] ] * VH[ indices[ i ] ];
      mf << val;
      mfmag << abs( val );
      val /= (1 - abs( Lambda[ indices[ i ] ] ));
      mfr << val;
      mfrmag << abs( val );
    }
    mf << '\n'; mfmag << '\n'; mfr << '\n'; mfrmag << '\n';
  }
  

  // predicted variance of the output quantity
  
  {
    ofstream anal_var( ( directory + "/analytic-variance.out" ).c_str( ) );
    double input_variance
      = ( population->getSVariance( ) > 0 ?
	    population->getSVariance( ) : population->getAcVariance( ) );
    vcl_complex< double > var = 0;
    for ( int i = 0; i < H.size(); ++i )
      for ( int j = 0; j < H.size(); ++j )
      {
	var += qU[ i ] * VH[ i ] *
	  conj( qU[ j ] * VH[ j ] ) /
	  ( 1.0 - Lambda[ i ] * conj( Lambda[ j ] ) );
      }
    var *= input_variance;
    // trust this to come out approximately real
    anal_var << real( var ) << "\n";
    ofstream varobjs( ( directory + "/variance-objects.txt" ).c_str( ) );
    vnl_matrix< vcl_complex< double > > inner_matrix( H.size( ), H.size( ) );
    for ( int i = 0; i < H.size(); ++i )
      for ( int j = 0; j < H.size(); ++j )
	inner_matrix[ i ][ j ] =
	  1.0 / ( 1.0 - Lambda[ i ] * conj( Lambda[ j ] ) );

    varobjs << "q = \n" << q << "\n"
	    << "m = \n" << qU << "\n"
	    << "H = \n" << H << "\n"
	    << "f = \n" << VH << "\n"
	    << "Lambda = \n" << Lambda
	    << "( 1/(1 - li \\bar lj) ) =\n" << inner_matrix
	    << "E(u*u) = " << input_variance << "\n";
  } 

  // components of variance gain of the output quantity

  {
    ofstream var_comps( ( directory +
			  "/analytic-variance-components.out" ).c_str( ) );

    // first variances
if (0)
    for ( int i = 0; i < H.size(); ++i )
    {
      if (i > 0) var_comps << '\t';
      var_comps << qU[ i ] * VH[ i ] * conj( qU[ i ] * VH[ i ] ) /
	( 1.0 - Lambda[ i ] * conj( Lambda[ i ] ) );
    }
    
    // then covariances
    for ( int i = 0; i < H.size(); ++i )
    {
      for ( int j = 0; j < H.size(); ++j )
	//if (j != i)
      {
	vcl_complex<double> v_ij =
	  qU[ indices[ i ] ] * VH[ indices[ i ] ] *
	  conj( qU[ indices[ j ] ] * VH[ indices[ j ] ] ) /
	  ( 1.0 - Lambda[ indices[ i ] ] * conj( Lambda[ indices[ j ] ] ) );
	if (j > 0)
	  var_comps << ' ';
	if (fabs(imag(v_ij))<0.0001)
	  var_comps << real(v_ij);
	else
	  var_comps << v_ij;
      }
      var_comps << '\n';
    }
    } 
  
  {
    ofstream var_comps( ( directory +
			  "/analytic-variance-components-pairwise.out"
			  ).c_str( ) );

    for ( int i = 0; i < H.size(); ++i )
    {
      for ( int j = 0; j < H.size(); ++j )
      {
	if (j > 0)
	  var_comps << ' ';
	if (j < i)
	  var_comps << 0;
	else
	{
	  vcl_complex<double> v_ij =
	    qU[ indices[ i ] ] * VH[ indices[ i ] ] *
	    conj( qU[ indices[ j ] ] * VH[ indices[ j ] ] ) /
	    ( 1.0 - Lambda[ indices[ i ] ] * conj( Lambda[ indices[ j ] ] ) );
	  if (j == i)
	    var_comps << (fabs(imag(v_ij))<0.0001 ? real(v_ij) : v_ij);
	  else
	    var_comps << v_ij + conj(v_ij);
	}
      }
      var_comps << '\n';
    }
    } 

  // "contributions" to variance gain of the output quantity
  
  {
    ofstream var_comps( ( directory +
			  "/analytic-variance-contributions.out" ).c_str( ) );

    // the variance is sum over i and j of these quantities
    //  here I output how much each subspace i contributes to variance:
    //  all the components of variance involving i -- i.e. how much
    //  variance would disappear if subspace i were suppressed
    for ( int i = 0; i < H.size(); ++i )
    {
      vcl_complex<double> con = 0;
      for ( int j = 0; j < H.size(); ++j )
      {
	vcl_complex<double> v_ij =
	  qU[ indices[ i ] ] * VH[ indices[ i ] ] *
	    conj( qU[ indices[ j ] ] * VH[ indices[ j ] ] ) /
	    ( 1.0 - Lambda[ indices[ i ] ] * conj( Lambda[ indices[ j ] ] ) );
	if (j == i)
	  con += v_ij;
	else
	  con += v_ij + conj(v_ij);
      }
      var_comps << con;// indeed these are approx. real 
      //var_comps << abs(real(con)) << ' ';
    }
    var_comps << '\n';
  } 

  cout << "done!\n";
} 

//void Simulation::fft( ) {
void Simulation::observed_transfer( ) 
{
  
  int nYears = population->getCurrentYear( ) + 1;
  vnl_vector< double > average( population->getMaxAgePerSalmon( ), 0.0 );

  for( int i = 0; i < population->getMaxAgePerSalmon( ); ++i ) {
    const vector< double > &ageClass = population->getAgeClass( i );
    for( int j = 0; j < ageClass.size( ); ++j ) {
      average[ i ] += ageClass[ j ];
    }
    average[ i ] /= ageClass.size( );
  }

  vnl_vector< vcl_complex< double > > q_hat( nYears, 0.0 );
  vnl_fft_1d< double > fft( nYears );
  //vnl_fft_1d< double > u_hat( u.size( ) );

  for( int i = 0; i < population->getMaxAgePerSalmon( ); ++i ) {
    const vector< double > &yi = population->getAgeClass( i );
    for( int k = 0; k < nYears; ++k ) {
      q_hat[ k ] += output_weights[ i ] * ( yi[ k ] - average[ i ] );
//      output_fft[ k ] += output_weights[ i ] * ( y[ k ] - population->getXStar( i ) );
    }
  }

  fft.fwd_transform( q_hat );

  vector< vnl_vector< vcl_complex< double > > > u_hat;
#if 0
  vector<double> u_mean;
  if( population->getSVariance( ) > 0  ) {
    vector< double > &u = population->getS( );
    u_hat.push_back( vnl_complexify( vnlize_vector(u) ) );
    u_mean.push_back( population->getS_Mean( ) );
  }
  if( population->getAcVariance( ) > 0 ) {
    vector< double > &u = population->getAc( );
    u_hat.push_back( vnl_complexify( vnlize_vector(u) ) );
    u_mean.push_back( population->getAc_Mean( ) );
  }
  
  for( int i = 0; i < u_hat.size( ); ++i )
    for( int j = 0; j < u_hat[i].size( ); ++j )
      u_hat[i][j] -= u_mean[i];
#else
  const vector< vnl_vector< double > > &u = population->get_u_history();
  // has to be transposed, unfortunately
  for (int i = 0; i < u[0].size(); ++i )
  {
    u_hat.push_back( vnl_vector< vcl_complex<double> >(u.size()) );
    for (int j = 0; j < u.size(); ++j )
      u_hat[i][j] = u[j][i];
  }
#endif

  for( int i = 0; i < u_hat.size( ); ++i )
    fft.fwd_transform( u_hat[i] );
  
  ofstream population_fft_file( (directory +
				 "/population.fft.mag.out").c_str( ) );
  ofstream u_fft_file( (directory + "/u.fft.mag.out").c_str( ) );
  
  for( int i = 0; i < q_hat.size( ); ++i )
  {
//     history_fft_file << q_hat[ i ].real( )
// 		     << " "
// 		     << q_hat[ i ].imag( )
// 		     << '\n';
    population_fft_file << abs(q_hat[i]) << endl;
    if ( ( i % 1000 ) == 0 )
    {
      cout << ".";
      cout.flush();
    }
  }
  
  for( int i = 0; i < u_hat[0].size( ); ++i )
  {
//     u_fft_file << u_hat[ i ].real( )
// 	       << " "
// 	       << u_hat[ i ].imag( )
// 	       << '\n';
    for( int j = 0; j < u_hat.size( ); ++j )
      u_fft_file << abs(u_hat[j][i]) << ' ';
    u_fft_file << endl;
    if ( ( i % 1000 ) == 0 )
    {
      cout << ".";
      cout.flush();
    }
  }

  ofstream pointwise_file( (directory + "/observed-transfer.out").c_str( ) );

  for( int i = 0; i < q_hat.size( )/* && i < output_fft.size( )*/; ++i ) {
    for (int j = 0; j < u_hat.size( ); ++j )
      //if ( abs( u_hat[j][i] ) >= 0.1 )// censor big spikes
	pointwise_file << ( abs( q_hat[ i ] / u_hat[ j ][ i ] ) ) << ' ';
    //else
    //  pointwise_file << abs(q_hat[i]) << ' ';
    pointwise_file << '\n';
    if ( ( i % 1000 ) == 0 )
    {
      cout << ".";
      cout.flush();
    }
  } 
  
  }
  

vcl_complex<double>
  Simulation::Transfer(vcl_complex<double> z,
		       vnl_vector< vcl_complex<double> > &m,
		       vnl_diag_matrix< vcl_complex<double> > &Lambda,
		       vnl_vector< vcl_complex<double> > &f,
		       vnl_vector< vcl_complex<double> > &f1)
{
  
  vnl_diag_matrix< vcl_complex< double > > dz( Lambda.size( ) );

  for ( int i = 0; i < dz.size(); ++i )
    dz[i] = z / (z - Lambda[i]);

  vcl_complex< double > zinv(1.0/z);
  return dot_product( m, dz * (f + zinv*f1) );
 
}

vnl_vector< vcl_complex<double> >
  Simulation::Transfer(vcl_complex<double> z,
		       vnl_matrix< vcl_complex<double> > &U,
		       vnl_diag_matrix< vcl_complex<double> > &Lambda,
		       vnl_matrix< vcl_complex<double> > &V,
		       vnl_vector< vcl_complex<double> > &H)
#if 1
{ // I used 'time', I think this version is faster
  
  vnl_diag_matrix< vcl_complex< double > > dz( Lambda.size( ) );

  for ( int i = 0; i < dz.size(); ++i )
    dz[i] = z / (z - Lambda[i]);

  return U * (dz * (V * H));
  
}
#else
{
  
  vnl_vector< vcl_complex< double > > T( U.rows( ), 0 );
  const vnl_vector< vcl_complex< double > > & F( V * H );
  for( int k = 0; k < T.size( ); k++ ) {
    T += U.get_column( k ) * ( F[ k ] * z / (z - Lambda[ k ]) );
  }
  return T;
  
}
{
  
  vnl_vector< vcl_complex< double > > T( U.rows( ), 0 );

  for( int k = 0; k < T.size( ); k++ ) {

    T += U.get_n_columns( k, 1 ) * ( z / ( z - Lambda( k, k ) ) ) * V.get_n_rows( k, 1 ) * H;

  }
  
}
#endif

// once needed an abs that has complex type rather than double
//vcl_complex<double> cabs(vcl_complex<double>c)
//{
//  return abs(c);
//}

void Simulation::analytic_transfer
( int nYears, vnl_vector< vcl_complex< double > > &qU,
  vnl_diag_matrix< vcl_complex< double > > &Lambda,
  vnl_vector< vcl_complex< double > > &VH,
  vnl_vector< vcl_complex< double > > &VH1, vector< int > &indices)
{
  
  ofstream antr( (directory + "/analytic-transfer.out" ).c_str() ),
    tcomp( (directory + "/analytic-transfer-components.out").c_str() ),
    tpk( ( directory + "/transfer-peaks.out" ).c_str() ),
    tpar( ( directory +
	    "/analytic-transfer-components-parametric.out" ).c_str() );

  for( int w = 0; w < nYears; w++ ) {
    const vcl_complex<double> i( 0, 1 );
    vcl_complex<double> z = exp( ( 2 * M_PI * w * i ) / (double)nYears );
    
    //vnl_vector< vcl_complex<double> > T = Transfer(z, qU, Lambda, VH);
    //vcl_complex< double > output_quantity = dot_product( n, T );
    vcl_complex<double> output_quantity = Transfer(z, qU, Lambda, VH, VH1);
    
    antr << abs( output_quantity ) << '\n';

    vcl_complex< double > zinv(1.0/z);
    for ( int i = 0; i < Lambda.size(); i++ )
    {
      vcl_complex<double> ti = (qU[ indices[ i ] ] *
				(z / (z - Lambda[ indices[ i ] ])) *
				(VH[ indices[ i ] ] + zinv*VH1[ indices[ i ] ]) );
      tcomp << abs( ti ) << ' ';
      tpar << ti << ' ';
    }
    tcomp << '\n';
    tpar << '\n';
    
    if ( ( w % 1000 ) == 0 )
    {
      cout << ".";
      cout.flush();
    }
  }
  for (int i = 0; i < indices.size(); ++i)
  {
    double angle = arg( Lambda[ i ] );
    if ( angle < 0 )
      angle += 2 * M_PI;
    vcl_complex<double> z =
      Lambda[ indices[ i ] ] / abs( Lambda[ indices[ i ] ] );
    vcl_complex< double > zinv(1.0/z);
    tpk << angle << ' '
	<< abs( qU[ indices[ i ] ] * (z / (z - Lambda[ indices[ i ] ])) *
		( VH[ indices[ i ] ] + zinv * VH1[ indices[ i ] ] ) )
	<< '\n';
  }
  //    transfer_peaks << '\n';
  
}

void Simulation::analytic_extinction_time
( vnl_vector<double> &q,
  vnl_matrix< vcl_complex<double> > &U,
  vnl_diag_matrix< vcl_complex<double> > &Lambda,
  vnl_vector< vcl_complex<double> > &f, vector<int> &indices,
  double u_variance, double scaled_b, vcl_complex<double> a,
  double extinction_threshold)
{
  double T;
  const int n = Lambda.cols( );

  // the relevant quantity q here is always total population
  //  const vnl_vector<double> q(n, 1);

  vnl_vector<double> xstar(q.size(),0.0); 
  xstar.update(population->getXStar());

  const double qx = dot_product(q, xstar);
  const double alpha
    = extinction_threshold - qx;

  vcl_complex<double> sum = 0;
  vnl_vector< vcl_complex<double> > m = vnl_complexify(q) * U;
  for ( int i = 0; i < n; i++ )
    for ( int j = 0; j < n; j++ )
      sum += (f[i] * m[i] * conj(f[j] * m[j])
	      / (1.0 - Lambda[i] * conj(Lambda[j])));
  //sum *= u_variance;
  // sum == qVq*

  // scaled_b might not be sqrt(input_variance*(1-a*a))
  // due to scaling over in the Population class
  // here we correct for that.
  double unscaled_b = real( sqrt(u_variance*(1.0-a*conj(a))) );
  double variance_rescaling = unscaled_b / scaled_b;  

  //double argument = - alpha / (sqrt(2*u_variance*real(sum)));
  double argument = - alpha / (sqrt(2*variance_rescaling*real(sum)));
  double I = 0.5 * (1 - erf(argument));
  T = 1 / I;

  ofstream eobjs( (directory + "/extinction-objs.txt").c_str() );
  eobjs << "u_variance (sigma^2) = " << u_variance << endl;
  eobjs << "variance_rescaling = " << variance_rescaling << endl;
  eobjs << "q = " << q << endl;
  eobjs << "f = " << f << endl;
  eobjs << "m = " << m << endl;
  eobjs << "qVq* = " << sum << endl;
  eobjs << "alpha = " << alpha << endl;
  eobjs << "-alpha/(sqrt(2 qVq*) sigma) = " << argument << endl;
  eobjs << "1/2 (1 - erf( -alpha/(sqrt(2 qVq*) sigma) )) = " << I << endl;
  eobjs << "T = " << T << endl;
  eobjs << endl;

  if (vnl_math_isfinite(T))
  {
    ofstream etime( (directory + "/analytic-extinction-time.out").c_str() );
    etime << T << "\n";
  }
}
