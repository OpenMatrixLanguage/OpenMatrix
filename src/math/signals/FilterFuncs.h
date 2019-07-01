/**
* @file FilterFuncs.h
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
* This file is part of the OpenMatrix Language ("OpenMatrix") software.
* Open Source License Information:
* OpenMatrix is free software. You can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* OpenMatrix is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
* You should have received a copy of the GNU Affero General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
* 
* Commercial License Information: 
* For a copy of the commercial license terms and conditions, contact the Altair Legal Department at Legal@altair.com and in the subject line, use the following wording: Request for Commercial License Terms for OpenMatrix.
* Altair's dual-license business model allows companies, individuals, and organizations to create proprietary derivative works of OpenMatrix and distribute them - whether embedded or bundled with other software - under a commercial license agreement.
* Use of Altair's trademarks and logos is subject to Altair's trademark licensing policies.  To request a copy, email Legal@altair.com and in the subject line, enter: Request copy of trademark and logo usage policy.
*/
#ifndef _Signals_WrapperFuncs_h
#define _Signals_WrapperFuncs_h

#include "SignalsExports.h"

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;
typedef hwTMatrix<int, hwTComplex<int> > hwMatrixI;

//------------------------------------------------------------------------------
//!
//! \brief Digital signal processing functions
//!
//------------------------------------------------------------------------------

// Filter generation

//!
//! Computes filter transfer function coefficients and returns status
//! \param order          Filter order
//! \param lowCutoffFreq
//! \param highCutoffFreq 
//! \param numerCoef
//! \param denomCoef
//! \param type           Optional argument
//!
SIGNALS_DECLS hwMathStatus Besself(int         order, 
                                   double      lowCutoffFreq, 
                                   double      highCutoffFreq,
                                   hwMatrix&   numerCoef, 
                                   hwMatrix&   denomCoef, 
                                   const char* type = "z");
//!
//! Computes filter transfer function coefficients and returns status
//! \param order          Filter order
//! \param lowCutoffFreq
//! \param highCutoffFreq 
//! \param numerCoef
//! \param denomCoef
//! \param type           Optional argument
//!
SIGNALS_DECLS hwMathStatus Besself3(int         order, 
                                    double      lowCutoffFreq, 
                                    double      highCutoffFreq,
                                    hwMatrix&   numerCoef, 
                                    hwMatrix&   denomCoef, 
                                    const char* type = "z");
//!
//! Computes Butterworth filter transfer function coefficients and returns status
//! \param order          Filter order
//! \param lowCutoffFreq
//! \param highCutoffFreq 
//! \param numerCoef
//! \param denomCoef
//! \param type           Optional argument
//!
SIGNALS_DECLS hwMathStatus Butter(int         order, 
                                  double      lowCutoffFreq, 
                                  double      highCutoffFreq,
                                  hwMatrix&   numerCoef, 
                                  hwMatrix&   denomCoef, 
                                  const char* type = "z");
//!
//! Computes Chebyshev type I filter transfer function coefficients and returns status
//! \param order          Filter order
//! \param lowCutoffFreq
//! \param highCutoffFreq 
//! \param passEdgeDb
//! \param numerCoef
//! \param denomCoef
//! \param type           Optional argument
//!
SIGNALS_DECLS hwMathStatus Cheby1(int         order, 
                                  double      lowCutoffFreq, 
                                  double      highCutoffFreq,
                                  double      passEdgeDb, 
                                  hwMatrix&   numerCoef, 
                                  hwMatrix&   denomCoef, 
                                  const char* type = "z");
//!
//! Computes Chebyshev type II filter transfer function coefficients and returns status
//! \param order          Filter order
//! \param lowCutoffFreq
//! \param highCutoffFreq 
//! \param stopEdgeDb
//! \param numerCoef
//! \param denomCoef
//! \param type           Optional argument
//!
SIGNALS_DECLS hwMathStatus Cheby2(int         order, 
                                  double      lowCutoffFreq, 
                                  double      highCutoffFreq,
                                  double      stopEdgeDb, 
                                  hwMatrix&   numerCoef, 
                                  hwMatrix&   denomCoef, 
                                  const char* type = "z");
//!
//! Computes Elliptic filter transfer function coefficients and returns status
//! \param order          Filter order
//! \param lowCutoffFreq
//! \param highCutoffFreq 
//! \param passEdgeDb
//! \param stopEdgeDb
//! \param numerCoef
//! \param denomCoef
//! \param type           Optional argument
//!
SIGNALS_DECLS hwMathStatus Ellip(int         order, 
                                 double      lowCutoffFreq, 
                                 double      highCutoffFreq,
                                 double      passEdgeDb, 
                                 double      stopEdgeDb,
                                 hwMatrix&   numerCoef, 
                                 hwMatrix&   denomCoef, 
                                 const char* type = "z");
//!
//! Computes FIR filter transfer function coefficients and returns status
//! \param order          Filter order
//! \param lowCutoffFreq
//! \param highCutoffFreq 
//! \param window
//! \param numerCoef
//! \param normalize      Optional argument
//!
SIGNALS_DECLS hwMathStatus Fir(int             order, 
                               double          lowCutoffFreq, 
                               double          highCutoffFreq,
                               const hwMatrix* window, 
                               hwMatrix&       numerCoef, 
                               bool            normalize = true);
//!
//! Computes FIR filter transfer function coefficients for multiband filters
//! \param order      Filter order
//! \param cutoffFreq 
//! \param numerCoef
//! \param type       Optional argument 
//! \param normalize  Optional argument
//!
SIGNALS_DECLS hwMathStatus Fir(int             order,
                               const hwMatrix& cutoffFreq,
                               const hwMatrix* window,
                               hwMatrix&       numerCoef,
                               const char*     type = "DC-0",
                               bool            normalize = true);
//!
//! Computes FIR filter transfer function coefficients for multiband filters
//! \param order      Filter order
//! \param cutoffFreq 
//! \param numerCoef
//! \param type       Optional argument 
//! \param normalize  Optional argument
//!
SIGNALS_DECLS hwMathStatus FirLS(int             order,
                                 const hwMatrix& freq,
                                 const hwMatrix& mag,
                                 const hwMatrix* weight,
                                 hwMatrix&       filterCoef);

// Filter order calculation

//!
//! Designs a Bessel filter and returns the status
//! \param passBandFreq
//! \param stopBandFreq
//! \param passEdgeDb
//! \param stopEdgeDb
//! \param order        Filter order
//! \param freqC
//! \param type         Optional argument
//!
SIGNALS_DECLS hwMathStatus BesselOrd(double      passBandFreq, 
                                     double      stopBandFreq,
                                     double      passEdgeDb, 
                                     double      stopEdgeDb, 
                                     int&        order,
                                     hwMatrix&   freqC,
                                     const char* type = "z");
//!
//! Designs a Bessel filter and returns the status
//! \param passBandFreq
//! \param stopBandFreq
//! \param passEdgeDb
//! \param stopEdgeDb
//! \param order        Filter order
//! \param freqC1
//! \param freqC2
//! \param type         Optional argument
//!
SIGNALS_DECLS hwMathStatus BesselOrd(const hwMatrix& passBandFreq, 
                                     const hwMatrix& stopBandFreq,
                                     double          passEdgeDb, 
                                     double          stopEdgeDb, 
                                     int&            order,
                                     hwMatrix&       freqC1,
                                     hwMatrix&       freqC2, 
                                     const char*     type = "z");
//!
//! Designs a Butterworth filter and returns the status
//! \param passBandFreq
//! \param stopBandFreq
//! \param passEdgeDb
//! \param stopEdgeDb
//! \param order
//! \param freqC
//! \param type         Optional argument
//!
SIGNALS_DECLS hwMathStatus ButterOrd(double      passBandFreq, 
                                     double      stopBandFreq,
                                     double      passEdgeDb,
                                     double      stopEdgeDb, 
                                     int&        order, 
                                     hwMatrix&   freqC, 
                                     const char* type = "z");
//!
//! Designs a Butterworth filter and returns the status
//! \param passBandFreq
//! \param stopBandFreq
//! \param passEdgeDb
//! \param stopEdgeDb
//! \param order        Filter order
//! \param freqC1
//! \param freqC2
//! \param type         Optional argument
//!
SIGNALS_DECLS hwMathStatus ButterOrd(const hwMatrix& passBandFreq, 
                                     const hwMatrix& stopBandFreq,
                                     double          passEdgeDb, 
                                     double          stopEdgeDb, 
                                     int&            order,
                                     hwMatrix&       freqC1, 
                                     hwMatrix&       freqC2, 
                                     const char*     type = "z");
//!
//! Designs a Chebyshev type I filter and returns the status
//! \param passBandFreq
//! \param stopBandFreq
//! \param passEdgeDb
//! \param stopEdgeDb
//! \param order        Filter order
//! \param freqC
//! \param type         Optional argument
//!
SIGNALS_DECLS hwMathStatus Cheby1Ord(double      passBandFreq, 
                                     double      stopBandFreq, 
                                     double      passEdgeDb,
                                     double      stopEdgeDb, 
                                     int&        order,
                                     hwMatrix&   freqC,
                                     const char* type = "z");
//!
//! Designs a Chebyshev type I filter and returns the status
//! \param passBandFreq
//! \param stopBandFreq
//! \param passEdgeDb
//! \param stopEdgeDb
//! \param order        Filter order
//! \param freqC1
//! \param freqC2
//! \param type         Optional argument
//!
SIGNALS_DECLS hwMathStatus Cheby1Ord(const hwMatrix& passBandFreq, 
                                     const hwMatrix& stopBandFreq,
                                     double          passEdgeDb, 
                                     double          stopEdgeDb, 
                                     int&            order,
                                     hwMatrix&       freqC1,
                                     hwMatrix&       freqC2,
                                     const char*     type = "z");
//!
//! Designs a Chebyshev type II filter and returns the status
//! \param passBandFreq
//! \param stopBandFreq
//! \param passEdgeDb
//! \param stopEdgeDb
//! \param order        Filter order
//! \param freqC
//! \param type         Optional argument
//!
SIGNALS_DECLS hwMathStatus Cheby2Ord(double      passBandFreq, 
                                     double      stopBandFreq, 
                                     double      passEdgeDb,
                                     double      stopEdgeDb, 
                                     int&        order, 
                                     hwMatrix&   freqC, 
                                     const char* type = "z");
//!
//! Designs a Chebyshev type II filter and returns the status
//! \param passBandFreq
//! \param stopBandFreq
//! \param passEdgeDb
//! \param stopEdgeDb
//! \param order        Filter order
//! \param freqC1
//! \param freqC2
//! \param type         Optional argument
//!
SIGNALS_DECLS hwMathStatus Cheby2Ord(const hwMatrix& passBandFreq, 
                                     const hwMatrix& stopBandFreq,
                                     double          passEdgeDb, 
                                     double          stopEdgeDb,
                                     int&            order,
                                     hwMatrix&       freqC1, 
                                     hwMatrix&       freqC2, 
                                     const char*     type = "z");
//!
//! Designs an Elliptic filter and returns the status
//! \param passBandFreq
//! \param stopBandFreq
//! \param passEdgeDb
//! \param stopEdgeDb
//! \param order
//! \param freqC
//! \param type         Optional argument
//!
SIGNALS_DECLS hwMathStatus EllipOrd(double      passBandFreq, 
                                    double      stopBandFreq, 
                                    double      passEdgeDb,
                                    double      stopEdgeDb, 
                                    int&        order, 
                                    hwMatrix&   freqC, 
                                    const char* type = "z");
//!
//! Designs an Elliptic filter and returns the status
//! \param passBandFreq
//! \param stopBandFreq
//! \param passEdgeDb
//! \param stopEdgeDb
//! \param order        Filter order
//! \param freqC1
//! \param freqC2
//! \param type         Optional argument
//!
SIGNALS_DECLS
hwMathStatus EllipOrd(const hwMatrix& passBandFreq, 
                      const hwMatrix& stopBandFreq,
                      double          passEdgeDb,
                      double          stopEdgeDb, 
                      int&            order,
                      hwMatrix&       freqC1, 
                      hwMatrix&       freqC2, 
                      const char*     type = "z");

// data filtering

//!
//! Filter a signal using the transfer function and return the status
//! \param numerCoef
//! \param denomCoef
//! \param inSignal
//! \param outSignal
//!
SIGNALS_DECLS hwMathStatus Filter(const hwMatrix& numerCoef, 
                                  const hwMatrix* denomCoef,
                                  const hwMatrix& inSignal, 
                                  hwMatrix&       outSignal);
//!
//! Filter a signal using the transfer function, making two passes to produce a
//! zero phase shift due to two-pass transer function
//! \param numerCoef
//! \param denomCoef
//! \param inSignal
//! \param outSignal
//!
SIGNALS_DECLS hwMathStatus FiltFilt(const hwMatrix& numerCoef, 
                                    const hwMatrix* denomCoef,
                                    const hwMatrix& inSignal, 
                                    hwMatrix&       outSignal);

// resampling

//!
//! Down sample a signal by an integer factor
//! \param inSignal  Signal to be downsampled
//! \param outSignal Re-sampled signal
//! \param k         Downsample period
//! \param phase     Offset to first sample to be retained within period length
//!
SIGNALS_DECLS hwMathStatus DownSample(const hwMatrix& inSignal, 
                                      hwMatrix&       outSignal,
                                      int             k, 
                                      int             phase = 0);
//!
//! Up sample a signal by an integer factor
//! \param inSignal  Signal to be upsampled
//! \param outSignal Re-sampled signal
//! \param k         Upsample factor
//! \param phase     Offset to the first sample, within period length
//!
SIGNALS_DECLS hwMathStatus UpSample(const hwMatrix& inSignal, 
                                    hwMatrix&       outSignal,
                                    int             k, 
                                    int             phase = 0);

// filter response

//!
//! Compute the response of a filter at specified frequencies using the transfer
//! function
//! \param numerCoef
//! \param denomCoef
//! \param freq
//! \param response
//! \param sampFreq  Optional argument
//!
SIGNALS_DECLS hwMathStatus Response(const hwMatrix& numerCoef,
                                    const hwMatrix* denomCoef,
                                    const hwMatrix& freq, 
                                    hwMatrix&       response,
                                    const double*   sampFreq = nullptr);
//!
//! Compute the magnitude response of a filter at specified frequencies using 
//! the transfer function
//! \param numerCoef
//! \param denomCoef
//! \param freq
//! \param response
//! \param sampFreq  Optional argument
//!
SIGNALS_DECLS hwMathStatus MagRes(const hwMatrix& numerCoef, 
                                  const hwMatrix* denomCoef,
                                  const hwMatrix& freq, 
                                  hwMatrix&       mag,
                                  const double*   sampFreq = nullptr);
//!
//! Compute the phase response of a filter at specified frequencies using 
//! the transfer function
//! \param numerCoef
//! \param denomCoef
//! \param freq
//! \param phase
//! \param sampFreq  Optional argument
//!
SIGNALS_DECLS hwMathStatus PhaseRes(const hwMatrix& numerCoef, 
                                    const hwMatrix* denomCoef,
                                    const hwMatrix& freq, 
                                    hwMatrix&       phase,
                                    const double*   sampFreq = nullptr);
//!
//! Compute the impulse response of a filter at specified frequencies using 
//! the transfer function
//! \param numerCoef
//! \param denomCoef
//! \param impulseRes
//! \param time       Optional argument
//! \param numPnts    Optional argument
//! \param sampFreq   Optional argument
//!
SIGNALS_DECLS hwMathStatus ImpulseRes(const hwMatrix& numerCoef, 
                                      const hwMatrix* denomCoef,
                                      hwMatrix&       impulseRes, 
                                      hwMatrix*       time = nullptr, 
                                      int             numPnts = 0, 
                                      double          sampFreq = 1.0);
//!
//! Estimate IIR transfer function coefficients from frequency response information 
//! \param response
//! \param freq
//! \param numerCoef
//! \param denomCoef
//! \param weight    Optional argument
//! \param type      Optional argument
//!
SIGNALS_DECLS hwMathStatus FilterID(const hwMatrix& response, 
                                    const hwMatrix& freq,
                                    hwMatrix&       numerCoef, 
                                    hwMatrix&       denomCoef,
                                    const hwMatrix* weight = nullptr, 
                                    const char*     type = "z");

// window functions

//!
//! Compute Chebyshev window weights and return status
//! \param weight
//! \param sideLobe
//! \param type     Optional argument
//!
SIGNALS_DECLS hwMathStatus ChebyWin(hwMatrix&   weight, 
                                    double      sideLobe,
                                    const char* type = "symmetric");
//!
//! Compute Hamming window weights and return status
//! \param weight
//! \param type   Optional argument
//!
SIGNALS_DECLS hwMathStatus HammWin(hwMatrix&   weight, 
                                   const char* type = "symmetric");
//!
//! Compute Hann window weights and return status
//! \param weight
//! \param type   Optional argument
//!
SIGNALS_DECLS hwMathStatus HannWin(hwMatrix&   weight, 
                                   const char* type = "symmetric");
//!
//! Compute Bartlett-Hann window weights and return status
//! \param weight
//! \param type   Optional argument
//!
SIGNALS_DECLS hwMathStatus BartHannWin(hwMatrix&   weight,
                                       const char* type = "symmetric");
//!
//! Compute Blackman window weights and return status
//! \param weight
//! \param type   Optional argument
//!
SIGNALS_DECLS hwMathStatus BlackmanWin(hwMatrix&   weight, 
                                       const char* type = "symmetric");
//!
//! Compute Welch window weights and return status
//! \param weight
//! \param type   Optional argument
//!
SIGNALS_DECLS hwMathStatus WelchWin(hwMatrix&   weight, 
                                    const char* type = "symmetric");
//!
//! Compute Parzen window weights and return status
//! \param weight
//! \param type   Optional argument
//!
SIGNALS_DECLS hwMathStatus ParzenWin(hwMatrix&   weight, 
                                     const char* type = "symmetric");
//!
//! Compute Kaiser-Bessel window weights and return status
//! \param weight
//! \param type   Optional argument
//!
SIGNALS_DECLS hwMathStatus KaiserBesselWin(hwMatrix&   weight, 
                                           double      beta,
                                           const char* type = "symmetric");

// acoustic weighted filters

//!
//! Acoustic A weight magnitude function
//! \param freq      Frequency
//! \param mag_in    Input magnitude
//! \param mag_out   Output magnitude
//! \param reference Optional argument
//!
SIGNALS_DECLS hwMathStatus dBa(const hwMatrix& freq, 
                               const hwMatrix& mag_in, 
                               hwMatrix&       mag_out,
                               double          reference = 1.0);
//!
//! Acoustic B weight magnitude function
//! \param freq      Frequency
//! \param mag_in    Input magnitude
//! \param mag_out   Output magnitude
//! \param reference Optional argument
//!
SIGNALS_DECLS hwMathStatus dBb(const hwMatrix& freq, 
                               const hwMatrix& mag_in, 
                               hwMatrix&       mag_out,
                               double          reference = 1.0);
//!
//! Acoustic C weight magnitude function
//! \param freq      Frequency
//! \param mag_in    Input magnitude
//! \param mag_out   Output magnitude
//! \param reference Optional argument
//!
SIGNALS_DECLS hwMathStatus dBc(const hwMatrix& freq, 
                               const hwMatrix& mag_in, 
                               hwMatrix&       mag_out,
                               double          reference = 1.0);
//!
//! Acoustic U weight magnitude function
//! \param freq      Frequency
//! \param mag_in    Input magnitude
//! \param mag_out   Output magnitude
//! \param reference Optional argument
//!
SIGNALS_DECLS hwMathStatus dBu(const hwMatrix& freq, 
                               const hwMatrix& mag_in, 
                               hwMatrix&       mag_out,
                               double          reference = 1.0);

// struct for findpeaks function
struct PeakInfo
{
    hwMatrix* parabol_pp;
    hwMatrix* parabol_x;
    hwMatrix* height;
    hwMatrix* baseline;
    hwMatrix* roots;
};

//!
//! Find peaks function
//! \param signal
//! \param peakprop
//! \param propval
//! \param twosided
//! \param peaks
//! \param loc
//! \param extra     Optional argument
//!
SIGNALS_DECLS hwMathStatus FindPeaks(const hwMatrix& signal,
                                     bool            twoSided,
                                     double          minPeakHeight,
                                     int             minPeakDistance,
                                     int             minPeakWidth,
                                     hwMatrix&       peaks,
                                     hwMatrixI&      index,
                                     int             indexOrigin,
                                     PeakInfo*       extra = nullptr);

#endif // _Signals_WrapperFuncs_h
