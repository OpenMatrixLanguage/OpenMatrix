/**
* @file WaveformFuncs.h
* @date February 2022
* Copyright (C) 2007-2022 Altair Engineering, Inc.  
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
#ifndef _Signals_WaveformFuncs_h
#define _Signals_WaveformFuncs_h

#include "SignalsExports.h"
#include <string>

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

//------------------------------------------------------------------------------
//!
//! \brief Digital signal waveform functions
//!
//------------------------------------------------------------------------------

//!
//! Computes a rectangular pulse and returns status
//! \param time          time vector
//! \param width         pulse width
//!
SIGNALS_DECLS hwMathStatus RectPulse(const hwMatrix& time, 
                                     double          width, 
                                     hwMatrix&       waveform);
//!
//! Computes a triangular pulse and returns status
//! \param time          time vector
//! \param width         pulse width
//! \param skew          pulse skew
//!
SIGNALS_DECLS hwMathStatus TriPulse(const hwMatrix& time,
                                    double          width,
                                    double          skew,
                                    hwMatrix&       waveform);
//!
//! Computes a Gaussian pulse and returns status
//! \param time          time vector
//! \param fc            center frequency
//! \param bw            fractional bandwidth in frequency domain
//! \param bwr           fractional bandwidth reference
//!
SIGNALS_DECLS hwMathStatus GausPulse(const hwMatrix& time,
                                     double          fc,
                                     double          bw,
                                     double          bwr,
                                     hwMatrix&       waveform);
//!
//! Computes a Gaussian pulse cutoff time and returns status
//! \param fc            center frequency
//! \param bw            fractional bandwidth in frequency domain
//! \param bwr           fractional bandwidth reference
//! \param tpr           trailing pulse rejection level
//!
SIGNALS_DECLS hwMathStatus GausPulse(double  fc,
                                     double  bw,
                                     double  bwr,
                                     double  tpr,
                                     double& tc);
//!
//! Computes a Pulse train and returns status
//! \param time          time vector
//! \param delay         vector of time delays
//! \param func          function name for the pulse
//! \param args          arguments for the pulse function
//!
SIGNALS_DECLS hwMathStatus PulsTran(const hwMatrix&    time,
                                    const hwMatrix&    delay,
                                    const std::string& func,
                                    const hwMatrix&    args,
                                    hwMatrix&          waveform);
//!
//! Computes a Pulse train and returns status
//! \param time          time vector
//! \param delay         vector of time delays
//! \param pulse         vector of pulse values
//! \param fs            sampling frequency
//! \param method        interpolation method
//!
SIGNALS_DECLS hwMathStatus PulsTran(const hwMatrix&    time,
                                    const hwMatrix&    delay,
                                    const hwMatrix&    pulse,
                                    double             fs,
                                    const std::string& method,
                                    hwMatrix&          waveform);
//!
//! Computes a square pulse and returns status
//! \param time          time vector
//! \param duty          pulse "on" time
//!
SIGNALS_DECLS hwMathStatus SquarePulse(const hwMatrix& time,
                                       double          duty,
                                       hwMatrix&       waveform);
//!
//! Computes a sawtooth pulse and returns status
//! \param time          time vector
//! \param duty          pulse "on" time
//!
SIGNALS_DECLS hwMathStatus SawToothPulse(const hwMatrix& time,
                                         double          width,
                                         hwMatrix&       waveform);
//!
//! Computes the Dirichlet function and returns status
//! \param time          time vector
//! \param duty          pulse "on" time
//!
SIGNALS_DECLS hwMathStatus Diric(const hwMatrix& time,
                                 int             n,
                                 hwMatrix&       waveform);
//!
//! Computes a chirp pulse and returns status
//! \param time          time vector
//! \param f0            frequency at time=0
//! \param t1            time t1
//! \param f1            frequency at time=t1
//! \param shape         pulse shape
//! \param phase         pulse phase shift at time=0
//!
SIGNALS_DECLS hwMathStatus ChirpPulse(const hwMatrix& time,
                                      double          f0,
                                      double          t1,
                                      double          f1,
                                      const char*     shape,
                                      double          phase,
                                      hwMatrix&       waveform);

#endif // _Signals_WaveformFuncs_h
