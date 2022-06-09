/**
* @file FourierFuncs.h
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
#ifndef _Signals_FourierFuncs_h
#define _Signals_FourierFuncs_h

#include <vector>
#include "SignalsExports.h"

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
template <typename T1, typename T2> class hwTMatrixN;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;
typedef hwTMatrixN<double, hwTComplex<double> > hwMatrixN;

//------------------------------------------------------------------------------
//!
//! \brief Fourier transform functions
//!
//------------------------------------------------------------------------------

//!
//! Get the mean sample rate from a time vector and compute the standard 
//! deviation
//! \param time           Time vector
//! \param meanSampleRate Mean sample rate
//! \param stdDevRate     Standard deviation rate
//! \param scale          Optional parameter
//!
SIGNALS_DECLS hwMathStatus SampleRate(const hwMatrix& time, 
                                      double&         meanSampRate,
                                      double&         stdDevRate, 
                                      double          scale = 1.0);
//!
//! Get the sampling frequency from a time vector and returns the status
//! \param time     Time vector
//! \param sampFreq Sampling frequency
//!
SIGNALS_DECLS hwMathStatus SampleFreq(const hwMatrix& time,
                                      double&         sampFreq);
//!
//! Create a frequency vector from a sampling frequency  
//! \param numPnts  Number of points
//! \param sampFreq Sampling frequency
//! \param freq     Frequency vector
//!
SIGNALS_DECLS hwMathStatus Freq(int         numPnts, 
                                double      sampFreq, 
                                hwMatrix&   freq,
                                const char* option);
//!
//! Create a frequency vector from a time vector  
//! \param time Time vector
//! \param freq Frequency vector
//!
SIGNALS_DECLS hwMathStatus Freq(const hwMatrix& time, 
                                hwMatrix&       freq,
                                const char*     option);
//!
//! Create a 1-sided response vector by folding a 2-sided magnitude response vector  
//! \param mag_twosided 2-sided magnitude response vector
//! \param mag_onesided 1-sided response vector
//!
SIGNALS_DECLS hwMathStatus Fold(const hwMatrix& mag_twosided, 
                                hwMatrix&       mag_onesided);
//!
//! 2D FFT of a real or complex signal  
//! \param signal 
//! \param freqRes 
//!
SIGNALS_DECLS hwMathStatus Fft2(const hwMatrix& signal,
                                hwMatrix&       freqRes);
//!
//! 2D FFT of a real or complex signal  
//! \param signal 
//! \param m 
//! \param n 
//! \param freqRes 
//!
SIGNALS_DECLS hwMathStatus Fft2(const hwMatrix& signal,
                                int             m, 
                                int             n, 
                                hwMatrix&       freqRes);
//!
//! ND FFT of a real or complex signal  
//! \param signal 
//! \param freqRes 
//!
SIGNALS_DECLS hwMathStatus FftN(const hwMatrixN& signal,
                                hwMatrixN&       freqRes);
//!
//! ND FFT of a real or complex signal  
//! \param signal 
//! \param newdims 
//! \param freqRes 
//!
SIGNALS_DECLS hwMathStatus FftN(const hwMatrixN&        signal,
                                const std::vector<int>& newdims,
                                hwMatrixN&              freqRes);
//!
//! 2D inverse FFT of a real or complex signal  
//! \param freqRes 
//! \param signal 
//! \param assumeConjSym 
//!
SIGNALS_DECLS hwMathStatus Ifft2(const hwMatrix& freqRes,
                                 hwMatrix&       signal,
                                 bool            assumeConjSym);
//!
//! 2D inverse FFT of a real or complex signal
//! \param freqRes 
//! \param m 
//! \param n 
//! \param signal 
//! \param assumeConjSym 
//!
SIGNALS_DECLS hwMathStatus Ifft2(const hwMatrix& freqRes,
                                 int             m,
                                 int             n,
                                 hwMatrix&       signal,
                                 bool            assumeConjSym);
//!
//! ND inverse FFT of a real or complex signal  
//! \param freqRes 
//! \param signal 
//! \param assumeConjSym 
//!
SIGNALS_DECLS hwMathStatus IfftN(const hwMatrixN& freqRes,
                                 hwMatrixN&       signal,
                                 bool             assumeConjSym);
//!
//! FFT of a real or complex signal  
//! \param signal 
//! \param freqRes
//! \param fftSize Optional argument
//!
SIGNALS_DECLS hwMathStatus Fft(const hwMatrix& signal, 
                               hwMatrix&       freqRes,
                               int             fftSize = 0);
//!
//! FFT of a real or complex signal  
//! \param signal 
//! \param freqRes
//! \param fftSize Optional argument
//!
SIGNALS_DECLS hwMathStatus Fft(const hwMatrix& signal,
                               hwMatrix&       freqRes,
                               int             dim,
                               int             fftSize);
//!
//! FFT of a real or complex signal  
//! \param signal 
//! \param freqRes
//! \param fftSize Optional argument
//!
SIGNALS_DECLS hwMathStatus Fft(const hwMatrixN& signal,
                               hwMatrixN&       freqRes,
                               int              dim,
                               int              fftSize);
//!
//! Inverse FFT of a real or complex spectrum  
//! \param freqRes Real or complex spectrum
//! \param signal  Inverse FFT
//! \param fftSize Optional argument
//!
SIGNALS_DECLS hwMathStatus Ifft(const hwMatrix& freqRes,
                                hwMatrix&       signal, 
                                int             fftSize = 0);
//!
//! Inverse FFT of a real or complex signal  
//! \param signal 
//! \param freqRes
//! \param fftSize Optional argument
//!
SIGNALS_DECLS hwMathStatus Ifft(const hwMatrix& signal,
                                hwMatrix&       freqRes,
                                int             dim,
                                int             fftSize);
//!
//! Inverse FFT of a real or complex signal  
//! \param signal 
//! \param freqRes
//! \param fftSize Optional argument
//!
SIGNALS_DECLS hwMathStatus Ifft(const hwMatrixN& signal,
                                hwMatrixN&       freqRes,
                                int              dim,
                                int              fftSize);
//!
//! Short time FFT of a real or complex signal  
//! \param signal 
//! \param freqRes
//! \param fftSize Optional argument
//!
SIGNALS_DECLS hwMathStatus Stft(const hwMatrix& signal,
                                hwMatrix& freqRes,
                                int             fftSize = 0);
//!
//! Circular convolution of a real signal pair for periodic signals 
//! \param signal1 
//! \param signal2 
//! \param conv    Circular convolution
//! \param fftSize Optional argument
//!
SIGNALS_DECLS hwMathStatus ConvCirc(const hwMatrix& signal1, 
                                    const hwMatrix& signal2,
                                    hwMatrix&       conv, 
                                    int             fftSize = 0);
//!
//! Circular correlation of a real signal pair for periodic signals 
//! \param signal1 
//! \param signal2 
//! \param corr    Circular correlation
//! \param fftSize Optional argument
//!
SIGNALS_DECLS hwMathStatus CorrCirc(const hwMatrix& signal1, 
                                    const hwMatrix& signal2,
                                    hwMatrix&       corr,
                                    int             fftSize = 0);
//!
//! Coherence of a real signal pair with a window function for periodic signals
//! \param sysInput 
//! \param sysOutput
//! \param window             Window function
//! \param num_overlap_points Number of overlapping points 
//! \param cohere             Coherence
//! \param fftSize            Optional argument
//!
SIGNALS_DECLS hwMathStatus Coherence(const hwMatrix& sysInput, 
                                     const hwMatrix& sysOutput,
                                     const hwMatrix& window, 
                                     int             num_overlap_points,
                                     hwMatrix&       cohere, 
                                     int             fftSize = 0);
//!
//! Coherence of a real signal pair with a window function for periodic signals
//! \param sysInput 
//! \param sysOutput
//! \param blockSize          
//! \param num_overlap_points Number of overlapping points 
//! \param cohere             Coherence
//! \param fftSize            Optional argument
//!
SIGNALS_DECLS hwMathStatus Coherence(const hwMatrix& sysInput, 
                                     const hwMatrix& sysOutput,
                                     int             blockSize, 
                                     int             num_overlap_points,
                                     hwMatrix&       cohere, 
                                     int             fftSize = 0);
//!
//! Power spectral density of a real signal
//! \param signal   Input signal
//! \param sampFreq
//! \param density  Power spectral density        
//! \param fftSize  Optional argument
//!
SIGNALS_DECLS hwMathStatus PSD(const hwMatrix& signal, 
                               double          sampFreq, 
                               hwMatrix&       density,
                               int             fftSize = 0);
//!
//! Cross power spectral density of a real signal pair
//! \param signal1  Input signal 
//! \param signal2  Input signal 
//! \param sampFreq
//! \param density  Cross power spectral density        
//! \param fftSize  Optional argument
//!
SIGNALS_DECLS hwMathStatus CPSD(const hwMatrix& signal1, 
                                const hwMatrix& signal2, 
                                double          sampFreq,
                                hwMatrix&       density, 
                                int             fftSize = 0);
//!
//! Block power spectral density of a real signal with a window function
//! \param signal             Input signal 
//! \param window             Window function 
//! \param num_overlap_points Number of overlapping points 
//! \param sampFreq
//! \param density            Block power spectral density        
//! \param fftSize            Optional argument
//!
SIGNALS_DECLS hwMathStatus BlockPSD(const hwMatrix& signal, 
                                    const hwMatrix& window, 
                                    int             num_overlap_points,
                                    double          sampFreq, 
                                    hwMatrix&       density, 
                                    int             fftSize = 0);
//!
//! Block power spectral density of a real signal
//! \param signal             Input signal 
//! \param blockSize          
//! \param num_overlap_points Number of overlapping points 
//! \param sampFreq
//! \param density            Block power spectral density        
//! \param fftSize            Optional argument
//!
SIGNALS_DECLS hwMathStatus BlockPSD(const hwMatrix& signal, 
                                    int             blockSize, 
                                    int             num_overlap_points,
                                    double          sampFreq, 
                                    hwMatrix&       density, 
                                    int             fftSize = 0);
//!
//! Block cross power spectral density of a real signal pair with a window function
//! \param signal1            Input signal 
//! \param signal2            Input signal 
//! \param window             Window function
//! \param num_overlap_points Number of overlapping points 
//! \param sampFreq
//! \param density            Block power spectral density        
//! \param fftSize            Optional argument
//!
SIGNALS_DECLS hwMathStatus BlockCPSD(const hwMatrix& signal1, 
                                     const hwMatrix& signal2,
                                     const hwMatrix& window, 
                                     int             num_overlap_points,
                                     double          sampFreq, 
                                     hwMatrix&       density, 
                                     int             fftSize = 0);
//!
//! Block cross power spectral density of a real signal pair
//! \param signal1            Input signal 
//! \param signal2            Input signal 
//! \param blockSize
//! \param num_overlap_points Number of overlapping points 
//! \param sampFreq
//! \param density            Block power spectral density        
//! \param fftSize            Optional argument
//!
SIGNALS_DECLS hwMathStatus BlockCPSD(const hwMatrix& signal1, 
                                     const hwMatrix& signal2,
                                     int             blockSize, 
                                     int             num_overlap_points,
                                     double          sampFreq, 
                                     hwMatrix&       density, 
                                     int             fftSize = 0);

#endif // _Signals_FourierFuncs_h
