//******************************************************************************
// simulation.h
// 2005.01.03
// Lango, Trevor M.
//
//******************************************************************************

#ifndef SIMULATION_H
#define SIMULATION_H

#include "population.h"

#include "vnl_matrix.h"
#include "vnl_diag_matrix.h"

#include <iostream>
using std::ostream;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <complex>
using std::complex;

#include <new>

class Simulation;
class Simulation {

  friend ostream &operator<<( ostream &,
                              const Simulation * );

 private:
  const int maxTime;
  Population *population;
  char *output_quantity;
  vnl_vector< double > output_weights;
  string directory;

  bool do_color;
  bool do_extinction;
  bool cohort_extinction;
  double threshold;
  int extinction_time;

/* 	// fft objects: */
/*   vector< vcl_complex< double > > output_fft; */
/*   vector< vcl_complex< double > > u_fft; */

/*   // file stream objects: */
/*   fstream *history_file; */

/*   // fft file stream objects: */
/*   fstream *history_transfer_file; */
/*   fstream *history_fft_file; */
/*   fstream *u_fft_file; */
/*   fstream *pointwise_file; */

/*   // */
/*   fstream *transformed_forcing_file; */
/*   fstream *forcing_file; */

 protected:

 public:
  Simulation( const int = 0,
              Population * = new Population,
              char *q = "r",
	      const string directory = "out" );
  ~Simulation( );

  // setup functions:
  void doExtinction(double thresh, bool cohort);
  void doColor();

  // get functions:
  const int getMaxTime( ) const;
  Population *getPopulation( ) const;
  const vnl_vector< double > &getOutput_Weights( );

  // simulation functions:
  void run( );
  void report( );

  // misc functions:
  //void fft( );
  //void pointwise_divide( );
  void observed_transfer( );
  vcl_complex<double>
    Transfer( vcl_complex<double> z,
	      vnl_vector< vcl_complex<double> > &m,
	      vnl_diag_matrix< vcl_complex<double> > &Lambda,
	      vnl_vector< vcl_complex<double> > &f,
	      vnl_vector< vcl_complex<double> > &f1);
  vnl_vector< vcl_complex<double> >
    Transfer( vcl_complex<double> z,
	      vnl_matrix< vcl_complex<double> > &U,
	      vnl_diag_matrix< vcl_complex<double> > &Lambda,
	      vnl_matrix< vcl_complex<double> > &V,
	      vnl_vector< vcl_complex<double> > &H );
  void analytic_transfer(int nYears,
			 vnl_vector< vcl_complex< double > > &nU,
			 vnl_diag_matrix< vcl_complex< double > > &Lambda,
			 vnl_vector< vcl_complex< double > > &VH,
			 vnl_vector< vcl_complex< double > > &VH1,
			 vector< int > &indices);
  void analytic_extinction_time
    (vnl_vector< double > &q,
     vnl_matrix< vcl_complex< double > > &U,
     vnl_diag_matrix< vcl_complex< double > > &Lambda,
     vnl_vector< vcl_complex< double > > &VH,
     vector< int > &indices,
     double u_variance, double scaled_b, vcl_complex<double> a,
     double extinction_threshold);
};

#endif

