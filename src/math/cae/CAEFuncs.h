/**
* @file CAEFuncs.h
* @date February 2017
* Copyright (C) 2017-2018 Altair Engineering, Inc.  
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
#ifndef _CAE_Funcs_h
#define _CAE_Funcs_h

#include "CAEExports.h"

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

//------------------------------------------------------------------------------
//!
//! \brief CAE functions
//!
//------------------------------------------------------------------------------
//!
//! Returns hwMathStatus after executing the rainflow algorithm
//! \param time    
//! \param numBins  
//! \param minRange 
//! \param maxRange
//! \param output
//! \param hysteresis
//! \param result     Output matrix containing the result
//! 
CAEFUNCS_DECLS hwMathStatus RainFlowFunc(const hwMatrix& time, 
                                         int             numBins, 
                                         double          minRange, 
                                         double          maxRange,
                                         int             output,
                                         int             hysteresis,
                                         hwMatrix&       result);
//!
//! ISO 6487 function
//! \param inSignal  Input signal
//! \param sampFreq  Sampling frequency
//! \param cfc 
//! \param outSignal Output signal
//!
CAEFUNCS_DECLS hwMathStatus ISO6487(const hwMatrix& inSignal, 
                                    double          sampFreq, 
                                    double          cfc,
                                    hwMatrix&       outSignal);
//!
//! SAE filter
//! \param inSignal   Input signal 
//! \param sampFreq  
//! \param SAE_class
//! \param outSignal  
//! \param fftSize    Optional argument
//!
CAEFUNCS_DECLS hwMathStatus SAEFilter(const hwMatrix& inSignal, 
                                      double          sampFreq, 
                                      double          SAE_class,
                                      hwMatrix&       outSignal, 
                                      int             fftSize = 0);

//!
//! SAE filter 1995 function
//! \param inSignal  Input signal
//! \param sampFreq  Sampling frequency
//! \param cfc 
//! \param stdpad
//! \param direction
//! \param outSignal Output signal
//!
CAEFUNCS_DECLS hwMathStatus SAEFilt95(const hwMatrix& inSignal, 
                                      double          sampFreq, 
                                      double          cfc,
                                      int             stdpad, 
                                      int             direction, 
                                      hwMatrix&       outSignal);

#endif // _CAE_Funcs_h
