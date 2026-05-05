// Copyright Jeremy Sanders 2022
// Released under the MIT licence

#include <xsTypes.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSstream.h>
#include <XSstreams.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

using std::string;
using std::istringstream;
using std::log;
using std::exp;

#define T_MIN 0.0808
#define T_MAX 79.9

// Xspec/src/XSFunctions/calcMultiTempPlasma.cxx in XSpec
int calcMultiTempPlasma(const RealArray& energyArray, const int plasmaType,
                        const IntegerVector& Zarray, const RealArray& abun,
                        const Real dens, const Real z, const RealArray& Tarr,
                        const RealArray& Tbarr, const RealArray& DEMarr,
                        const int ifl, const bool qtherm, const Real velocity,
                        RealArray& fluxArray, RealArray& fluxErrArray);

namespace
{
  Real sqr(Real x)
  {
    return x*x;
  }

  // get number of steps in cooling flow
  size_t get_num_steps()
  {
    size_t numsteps = 21;
    const string ntempsval =
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
  Real get_max_sigma()
  {
    Real maxsigma = 3;
    const string maxsigmaval =
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

  // compute the model given abundances
  void inner_calc(const RealArray& energyArray, int spectrumNumber,
                  RealArray& fluxArray, RealArray& fluxErrArray,
                  const RealArray& abunin,
                  Real kT, Real logsigma, Real redshift, Real velocity)
  {
    const size_t numsteps = get_num_steps();
    const Real maxsigma = get_max_sigma();
    RealArray dem(numsteps);
    RealArray Tarr(numsteps);

    const Real deltalogt = logsigma * maxsigma * 2. / (numsteps-1);
    const Real logkT = log(kT);
    const Real minlogt = logkT - deltalogt*(numsteps-1)*0.5;

    Real total_em = 0;
    for(size_t i = 0; i != numsteps; ++i)
      {
        const Real logtstep = minlogt + i*deltalogt;
        const Real tstep = exp(logtstep);
        if( tstep >= T_MIN && tstep <= T_MAX )
          {
            Tarr[i] = tstep;
            dem[i] = exp( -0.5 * sqr((logtstep - logkT)/logsigma) );
            total_em += dem[i];
          }
        else
          {
            Tarr[i] = 1;
            dem[i] = 0;
          }

        if ( tpout.maxChatter() >= 15 )
          {
            tcout << "T=" << Tarr[i]
                  << " em=" << dem[i]
                  << '\n';
          }

      }

    if(total_em > 1e-30)
      {
        // now normalize to give 1. as total
        dem *= (1./total_em);
      }

    const int plasmaType = 6; // APEC

    // copy abundances (see bvapec.cxx)
    const size_t nAp(14);
    const int apelt[] = {1, 2, 6, 7, 8, 10, 12, 13, 14, 16, 18, 20, 26, 28};

    RealArray abun(nAp);
    IntegerVector Zarray(nAp);
    for (size_t i=0; i<nAp; i++) Zarray[i] = apelt[i];
    abun[0] = 1.0;
    for (size_t i=1; i<abun.size(); i++) abun[i] = abunin[i-1];

    calcMultiTempPlasma(energyArray, plasmaType, Zarray, abun,
                        (Real)1.0, redshift, Tarr, Tarr, dem,
                        spectrumNumber, true, velocity,
                        fluxArray, fluxErrArray);
  }

}

extern "C" void clclognorm(const RealArray& energy,
                           const RealArray& parameter,
                           int spectrum,
                           RealArray& flux,
                           RealArray& fluxError,
                           const string& init)
{
  const Real kT = parameter[0];
  const Real logsigma = parameter[1];
  const Real Z = parameter[2];
  const Real redshift = parameter[3];
  const Real velocity = 0.;

  // need abundances as floats :-(
  RealArray abun(13);
  abun[0] = 1.0;  // He
  for(size_t i = 1; i != 13; ++i)
    abun[i] = Z;

  inner_calc(energy, spectrum, flux, fluxError, abun, kT, logsigma, redshift, velocity);
}

extern "C" void clcvlognorm(const RealArray& energy,
                            const RealArray& parameter,
                            int spectrum,
                            RealArray& flux,
                            RealArray& fluxError,
                            const string& init)
{
  const Real kT = parameter[0];
  const Real logsigma = parameter[1];
  const Real velocity = 0.;

  RealArray abun(13);
  for(size_t i = 0; i != 13; ++i)
    abun[i] = parameter[i+2];

  const Real redshift = parameter[15];

  inner_calc(energy, spectrum, flux, fluxError, abun, kT, logsigma, redshift, velocity);
}

extern "C" void clcblognorm(const RealArray& energy,
                           const RealArray& parameter,
                           int spectrum,
                           RealArray& flux,
                           RealArray& fluxError,
                           const string& init)
{
  const Real kT = parameter[0];
  const Real logsigma = parameter[1];
  const Real Z = parameter[2];
  const Real redshift = parameter[3];
  const Real velocity = parameter[4];

  RealArray abun(13);
  abun[0] = 1.0;  // He
  for(size_t i = 1; i != 13; ++i)
    abun[i] = Z;

  inner_calc(energy, spectrum, flux, fluxError, abun, kT, logsigma, redshift, velocity);
}

extern "C" void clcbvlognorm(const RealArray& energy,
                             const RealArray& parameter,
                             int spectrum,
                             RealArray& flux,
                             RealArray& fluxError,
                             const string& init)
{
  const Real kT = parameter[0];
  const Real logsigma = parameter[1];

  RealArray abun(13);
  for(size_t i = 0; i != 13; ++i)
    abun[i] = parameter[i+2];

  const Real redshift = parameter[15];
  const Real velocity = parameter[16];

  inner_calc(energy, spectrum, flux, fluxError, abun, kT, logsigma, redshift, velocity);
}
