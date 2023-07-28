/**
* @file SignalsTboxFuncs.h
* @date January 2015
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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

#ifndef __SIGNALSTBOXFUNCS_OML_H__
#define __SIGNALSTBOXFUNCS_OML_H__

#pragma warning(disable: 4251)

#include "EvaluatorInt.h"
#include "SignalsTboxDefs.h"

template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

//------------------------------------------------------------------------------
//!
//! \brief oml Signals functions
//!
//------------------------------------------------------------------------------

extern "C" 
{
    //!
    //! Entry point which registers oml Signals functions with oml
    //! \param eval Evaluator interface
    //!
    SIGNALSOMLTBOX_DECLS int InitDll(EvaluatorInterface eval);
    //!
    //! Returns toolbox version
    //! \param eval Evaluator interface
    //!
    SIGNALSOMLTBOX_DECLS double GetToolboxVersion(EvaluatorInterface eval);
}

int determineFftSize(const hwMatrix *mtx);
int determineFftDim(const hwMatrix *mtx);

//!
//! Inverse Fast Fourier transform
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlIfft(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs);
//!
//! Fast Fourier transform
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFft(EvaluatorInterface           eval, 
            const std::vector<Currency>& inputs, 
            std::vector<Currency>&       outputs);
//!
//! Two dimensional Fast Fourier transform
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFft2(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs);
//!
//! Two dimensional inverse Fast Fourier transform
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlIfft2(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs,
              std::vector<Currency>&       outputs);
//!
//! N dimensional Fast Fourier transform
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFftN(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs,
             std::vector<Currency>&       outputs);
//!
//! N dimensional inverse Fast Fourier transform
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlIfftN(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs,
              std::vector<Currency>&       outputs);
//!
//! Short-time Fourier transform
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlStft(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs,
             std::vector<Currency>&       outputs);
//!
//! Inverse short-time Fourier transform
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlIstft(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs,
              std::vector<Currency>&       outputs);
//!
//! Spectrogram
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlSpectrogram(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
//!
//! Shift FFT related data to center DC component
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFftShift(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs);
//!
//! Reverse shift FFT related data to uncenter DC component
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlIFftShift(EvaluatorInterface           eval,
                  const std::vector<Currency>& inputs,
                  std::vector<Currency>&       outputs);
//!
//! Generates a vector of frequency locations
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFreq(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs);
//!
//! Creates one-sided spectrum amplitude vectors from the two-sided equivalents
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFold(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs);
//!
//! Performs cross correlation [xcorr command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlXcorr(EvaluatorInterface           eval, 
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs);
//!
//! Unwraps a vector of phase angles [unwrap command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlUnwrap(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs);
//!
//! Computes power spectral density [pwelch command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlPwelch(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs);
//!
//! Computes cross power spectral density [cpsd command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlCpsd(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs);
//!
//! Computes estimated transfer function [tfestimate command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlTFestimate(EvaluatorInterface           eval, 
				   const std::vector<Currency>& inputs, 
				   std::vector<Currency>&       outputs);
//!
//! Computes estimated mean square coherence [mscohere command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlMScohere(EvaluatorInterface           eval, 
				 const std::vector<Currency>& inputs, 
				 std::vector<Currency>&       outputs);
//!
//! Computes digital filter frequency response values [freqz command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFreqz(EvaluatorInterface           eval, 
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs);
//!
//! Computes analog filter frequency response values [freqs command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFreqs(EvaluatorInterface           eval, 
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs);
//!
//! Computes digital filter coefficients from frequency response values [invfreqz]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlInvFreqz(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs);
//!
//! Computes analog filter coefficients from frequency response values [infreqs]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlInvFreqs(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs);
//!
//! Computes the impulse response of a digital filter [impz command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlImpz(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs);
//!
//! Creates a digital FIR filter [fir1 command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFir1(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs);
//!
//! Creates a digital FIR filter [firls command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFirls(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs,
              std::vector<Currency>&       outputs);
//!
//! Creates a Bessel filter [besself command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlBesself(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs);
//!
//! Creates a Bessel filter [besself3 command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlBesself3(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs);
//!
//! Creates a Butterworth filter [butter command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlButter(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs);
//!
//! Creates a Chebyshev I filter [cheby1 command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlCheby1(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs);
//!
//! Creates a Chebyshev II filter [cheby2 command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlCheby2(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs);
//!
//! Creates a Elliptic filter [ellip command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlEllip(EvaluatorInterface           eval, 
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs);
//!
//! Designs a Butterworth filter [buttord command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlButtord(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs);
//!
//! Designs a Chebyshev I filter [cheb1ord command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlCheb1ord(EvaluatorInterface           eval, 
                  const std::vector<Currency>& inputs, 
                  std::vector<Currency>&       outputs);
//!
//! Designs a Chebyshev II filter [cheb2ord command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlCheb2ord(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs);
//!
//! Designs a Elliptic filter [ellipord command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlEllipord(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs);
//!
//! Returns Blackman window coefficients [blackman command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlBlackman(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs);
//!
//! Returns Bartlett-Hann window coefficients [barthann command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlBartHann(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs);
//!
//! Returns Hann window coefficients [hann command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlHann(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs,
             std::vector<Currency>&       outputs);
//!
//! Returns Hamming window coefficients [hamming command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlHamming(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs);
//!
//! Returns Parzen window coefficients [parzenwin command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlParzenwin(EvaluatorInterface           eval, 
                  const std::vector<Currency>& inputs, 
                  std::vector<Currency>&       outputs);
//!
//! Returns Welch window coefficients [welchwin command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlWelchwin(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs);
//!
//! Returns Dolph-Chebyshev window coefficients [chebwin command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlChebwin(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs);
//!
//! Returns Kaiser-Bessel window coefficients [kaiser command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlKaiser(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs);
//!
//! Filters a signal with an IIR or FIR filter [filter command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFilter(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs);
//!
//! Filter a signal forward and then backward, compensating for end effects [filtfilt command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFiltfilt(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs);
//!
//! Filters a signal with a 2D FIR filter [filter2 command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFilter2(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs);
//!
//! Computes the sinc function [sinc command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlSinc(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs,
             std::vector<Currency>&       outputs);
//!
//! Upsamples a signal [upsample command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlUpsample(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs);
//!
//! Downsamples a signal [downsample command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlDownsample(EvaluatorInterface           eval, 
                   const std::vector<Currency>& inputs, 
                   std::vector<Currency>&       outputs);
//!
//! Resamples a signal [resample command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlResample(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs);
//!
//! Resamples a signal [upfirdn command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlUpFirDn(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs);
//!
//! Find peaks of a signal [findpeaks command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFindPeaks(EvaluatorInterface           eval, 
                  const std::vector<Currency>& inputs, 
                  std::vector<Currency>&       outputs);
//!
//! Computes A-weighted acoustic function [dba command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmldbA(EvaluatorInterface           eval,
            const std::vector<Currency>& inputs,
            std::vector<Currency>&       outputs);
//!
//! Computes B-weighted acoustic function [dbb command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmldbB(EvaluatorInterface           eval,
            const std::vector<Currency>& inputs,
            std::vector<Currency>&       outputs);
//!
//! Computes C-weighted acoustic function [dbc command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmldbC(EvaluatorInterface           eval,
            const std::vector<Currency>& inputs,
            std::vector<Currency>&       outputs);
//!
//! Computes U-weighted acoustic function [dbu command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmldbU(EvaluatorInterface           eval,
            const std::vector<Currency>& inputs,
            std::vector<Currency>&       outputs);
//!
//! Computes rectangular pulse signal [rectpuls command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlRectPuls(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs);
//!
//! Computes triangular pulse signal [tripuls command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlTriPuls(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs);
//!
//! Computes Gaussian pulse signal [gauspuls command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlGausPuls(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs);
//!
//! Computes pulse train signal [gauspuls command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlPulsTran(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs);
//!
//! Computes square wave signal [square command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlSquare(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs,
               std::vector<Currency>&       outputs);
//!
//! Computes sawtooth wave signal [sawtooth command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlSawtooth(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs);
//!
//! Computes Dirichlet function signal [diric command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlDiric(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs,
              std::vector<Currency>&       outputs);
//!
//! Computes chirp signal [chirp command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlChirp(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs,
              std::vector<Currency>&       outputs);

#endif // __SIGNALSTBOXFUNCS_OML_H__         