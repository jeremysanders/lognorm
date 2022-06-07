// Copyright Jeremy Sanders 2022
// Released under the MIT licence

#include <xsTypes.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSstream.h>
#include <XSstreams.h>

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <valarray>
#include <string>

using namespace std;

#define T_MIN 0.0808
#define T_MAX 79.9

namespace
{
  double sqr(double x)
  {
    return x*x;
  }

  // get number of steps in cooling flow
  size_t get_num_steps()
  {
    size_t numsteps = 21;
    const std::string ntempsval =
      FunctionUtility::getModelString("LOGNORM_NTEMPS");
    if ( ntempsval != FunctionUtility::NOT_A_KEY() )
      {
	// get number of steps if set
	size_t tempstep;
	istringstream buffer(ntempsval);
	
	buffer >> tempstep;
	if( buffer.fail() )
	  {
	    tcerr << "Invalid value in LOGNORM_NTEMPS\n";
	  }
	else
	  {
	    numsteps = tempstep;
	  }
      }
    return numsteps;
  }

  // get max sigma
  double get_max_sigma()
  {
    double maxsigma = 3;
    const std::string maxsigmaval =
      FunctionUtility::getModelString("LOGNORM_MAXSIGMA");
    if ( maxsigmaval != FunctionUtility::NOT_A_KEY() )
      {
	// get number of steps if set
	size_t tempsigma;
	istringstream buffer(maxsigmaval);
	
	buffer >> tempsigma;
	if( buffer.fail() )
	  {
	    tcerr << "Invalid value in LOGNORM_MAXSIGMA\n";
	  }
	else
	  {
	    maxsigma = tempsigma;
	  }
      }
    return maxsigma;
  }

  extern "C"
  {
    void sumdem_(int* itype, int* swtch, float* ear, int* ne, float* abun,
		 float* dens, float* z, int* ninputt, float* intputt,
		 float* dem, int* ifl, int* qtherm, float* velocity,
		 float* photar, float* photerr, int* status);
  }

  // wrapper to get spectrum of thermal gas
  void calc_sumdem(const RealArray& energy,
		   const std::valarray<float>& abun,
		   float z,
		   const std::valarray<float>& temps,
		   const std::valarray<float>& dem,
		   RealArray& flux)
  {
    // copy energies to floats
    size_t esize = energy.size();
    valarray<float> ear(esize);
    for(size_t i = 0; i != esize; ++i)
      ear[i] = energy[i];

    // output array
    valarray<float> photar(esize-1);
    photar = 0.;

    int itype = 4; // apec
    int swtch = 2; // apec interpolate
    float dens = 1.;
    int numt = temps.size();
    int earsize = ear.size()-1;
    int ifl = 1;
    int qintherm = 1;
    float velocity = 0.;
    int status = 0;
    sumdem_( &itype, &swtch, &(ear[0]), &earsize,
	     const_cast<float*>(&(abun[0])), &dens, &z,
	     &numt, const_cast<float*>(&(temps[0])),
	     const_cast<float*>(&(dem[0])), &ifl, &qintherm, &velocity,
	     &(photar[0]), 0, &status);

    // write output into Real array
    flux.resize(esize-1);
    for(size_t i = 0; i != (esize-1); ++i)
      flux[i] = photar[i]; 
  }

  // compute the model given abundances
  void inner_calc(const RealArray& energy, RealArray& flux,
                  const valarray<float>& abun,
                  double kT, double logsigma, double redshift)
  {
    const size_t numsteps = get_num_steps();
    const double maxsigma = get_max_sigma();
    valarray<float> arry_dem(numsteps);
    valarray<float> arry_T(numsteps);

    const double deltalogt = logsigma * maxsigma * 2. / (numsteps-1);
    const double logkT = log(kT);
    const double minlogt = logkT - deltalogt*(numsteps-1)*0.5;

    double total_em = 0;
    for(size_t i = 0; i != numsteps; ++i)
      {
        const double logtstep = minlogt + i*deltalogt;
        const double tstep = exp(logtstep);
        if( tstep >= T_MIN && tstep <= T_MAX )
          {
            arry_T[i] = tstep;
            arry_dem[i] = exp( -0.5 * sqr((logtstep - logkT)/logsigma) );
            total_em += arry_dem[i];
          }
        else
          {
            arry_T[i] = 1;
            arry_dem[i] = 0;
          }

        if ( tpout.maxChatter() >= 15 )
          {
            tcout << "T=" << arry_T[i]
                  << " em=" << arry_dem[i]
                  << '\n';
          }

      }

    if(total_em > 1e-30)
      {
        // now normalize to give 1. as total
        arry_dem *= (1./total_em);
      }

    // have to copy arrays to floats to send to sumdem
    calc_sumdem(energy, abun, redshift, arry_T, arry_dem, flux);
  }

}

extern "C" void clclognorm(const RealArray& energy,
                           const RealArray& parameter,
                           int spectrum,
                           RealArray& flux,
                           RealArray& fluxError,
                           const string& init)

{
  const double kT = parameter[0];
  const double logsigma = parameter[1];
  const double Z = parameter[2];
  const double redshift = parameter[3];

  // need abundances as floats :-(
  valarray<float> abun(13);
  abun[0] = 1.0;  // He
  for(size_t i = 1; i != 13; ++i)
    abun[i] = Z;

  inner_calc(energy, flux, abun, kT, logsigma, redshift);
}

extern "C" void clcvlognorm(const RealArray& energy,
                            const RealArray& parameter,
                            int spectrum,
                            RealArray& flux,
                            RealArray& fluxError,
                            const string& init)

{
  const double kT = parameter[0];
  const double logsigma = parameter[1];

  valarray<float> abun(13);
  for(size_t i = 0; i != 13; ++i)
    abun[i] = parameter[i+2];

  const double redshift = parameter[15];

  inner_calc(energy, flux, abun, kT, logsigma, redshift);
}

