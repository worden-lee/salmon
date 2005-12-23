//******************************************************************************
// salmon.cpp
// 2005.01.03
// Lango, Trevor M.
//
//******************************************************************************

#include <new>

#include <iostream>
using std::cout;
using namespace std;

#include <string>
using std::string;

#include <fstream>
using std::fstream;
using std::ios;

#include "population.h"
#include "simulation.h"

int usage_error()
{
  cout <<
    "usage: salmon 'n' 'beta' 'alpha' 't' 's_mean' 's_variance' 'ac_mean'\n"
    " 'ac_variance' 'sigma' 'output_quantity' 'directory'\n"
    " [-e 'threshold'|-ec 'cohort threshold] [-c 'color'] [-3|-imp]\n"
       << "\tn:\t\t\tmax age of a given salmon\n"
       << "\tbeta:\t\t\tnest site competition factor\n"
       << "\talpha:\t\t\tmaximum number of offspring for a given salmon\n"
       << "\tt:\t\t\ttime (in years) for the simulation to run\n"
       << "\ts_mean:\t\t\t\n"
       << "\ts_variance:\t\t\t\n"
       << "\tac_mean:\t\t\t\n"
       << "\tac_variance:\t\t\t\n"
       << "\tsigma:\t\t\t\n"
       << "\toutput_quantity:\t\t\t\n"
       << "\tdirectory:\t\t\t\n"
       << "\t-e threshold:\t\t\t(optional) extinction threshold\n"
       << "\t-ec threshold:\t\t\t(optional) per-cohort extinction threshold\n"
       << "\t-c color:\t\t\t(optional) real exponent for colored noise\n"
       << "\t-cp r theta:\t\t\t(optional) complex exponent for colored noise\n"
       << "\t-3:\t\t\t(special case) suppress delta_l for 3-d model\n"
       << "\t-imp:\t\t\t(special case) generate impulse response\n";
  return -1;
}


int main( int argc, char *argv[ ] ) {

  int n;
  double beta;
  int alpha;

  int simulationTime;

  double s_mean;
  double s_variance;

  double ac_mean;
  double ac_variance;
  double sigma;

  double r=0;
  double theta=0;  

  char *output_quantity;
  string directory;

  double ext_thresh = - HUGE;
  bool cohort_ext;

  bool do_color = false;
  bool threeD = false;
  bool imp = false;
  bool early = false;

  // parse the arguments:
  if( argc < 12 ) 
    return usage_error();
  else
  {
    n = atoi( *++argv );              // 1
    beta = atof( *++argv );           // 2
    alpha = atoi( *++argv );          // 3
    simulationTime = atoi( *++argv ); // 4
    s_mean = atof( *++argv );         // 5
    s_variance = atof( *++argv );     // 6
    ac_mean = atof( *++argv );        // 7
    ac_variance = atof( *++argv );    // 8
    sigma = atof( *++argv );          // 9
    output_quantity = *++argv;        // 10
    directory = *++argv;              // 11
    
    while (*++argv)
    {
      string arg(*argv);
      if (arg == "-e")
      {
	const char *carg = *++argv;
	if (!carg)
	{ cout << "-e must be followed by a number\n";
          return usage_error();
	}
	ext_thresh = atof( carg );
	cohort_ext = false;
      }
      else if (arg == "-ec")
      {
	const char *carg = *++argv;
	if (!carg)
	{ cout << "-ec must be followed by a number\n";
          return usage_error();
	}
	ext_thresh = atof( carg );
	cohort_ext = true;
      }
      else if (arg == "-c")
      {
	const char *carg = *++argv;
	if (!carg)
	{ cout << "-c must be followed by a number\n";
          return usage_error();
	}
	r = atof( carg );
	//do_color = true;
      }
      else if (arg == "-cp")
      {
        const char *carg = *++argv;
        if (!carg) {
           cout << "-cp must be followed by two numbers." << endl;
           return usage_error();
        }
        r = atof( carg );
        const char *carg2 = *++argv;
        if (!carg2) {
           cout << "-cp must be followed by two numbers." << endl;
           return usage_error();
        }
        theta = atof( carg2 );
        //do_color = true;
      }
      else if (arg == "-se")
	early = true;
      else if (arg == "-3")
	threeD = true;
      else if (arg == "-imp")
	imp = true;
      else
      {
	cout << "unrecognized argument " << arg << endl;
	return usage_error();
      }
    }
  }

  // set everything up:
  Population *myPopulation = new Population( n,
                                             beta,
                                             alpha,
                                             s_mean,
                                             s_variance,
					     ac_mean,
					     ac_variance,
					     sigma,
                                             r,
                                             theta,
					     simulationTime,
					     threeD,
					     imp,
					     early );

  Simulation *mySimulation = new Simulation( simulationTime,
                                             myPopulation,
                                             output_quantity,
					     directory );

  if ( ext_thresh  > - HUGE )
    mySimulation->doExtinction(ext_thresh, cohort_ext);
  
  if ( do_color ) 
  {
    //myPopulation->doColor();
    //mySimulation->doColor();
  }

  // run the simulation:
  mySimulation->run( );

  // output necessary info:
  mySimulation->report( );

  // clean up:
  delete myPopulation;
  delete mySimulation;

  return 0;

}

