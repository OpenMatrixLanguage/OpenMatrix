/**
* @file SignalStatsFuncs.h
* @date February 2023
* Copyright (C) 2023 Altair Engineering, Inc.  
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
#ifndef _Signals_StatsFuncs_h
#define _Signals_StatsFuncs_h

#include "SignalsExports.h"
#include <string>

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;
template <typename T1, typename T2> class hwTMatrixN;
typedef hwTMatrixN<double, hwTComplex<double> > hwMatrixN;

//------------------------------------------------------------------------------
//!
//! \brief Digital signal statistic functions
//!
//------------------------------------------------------------------------------

//!
//! Computes peak to RMS values
//! \param data          signal vector or matrix
//! \param dim           dimension
//!
SIGNALS_DECLS hwMathStatus RSSQ(const hwMatrix& data,
                                int             dim,
                                hwMatrix&       result);
//!
//! Computes peak to RMS values
//! \param data          signal vector or matrix
//! \param dim           dimension
//!
SIGNALS_DECLS hwMathStatus RSSQ(const hwMatrixN& data,
                                int              dim,
                                hwMatrixN&       result);
//!
//!
//! Computes peak to RMS values
//! \param data          signal vector or matrix
//! \param dim           dimension
//!
SIGNALS_DECLS hwMathStatus Peak2RMS(const hwMatrix& data,
                                    int             dim,
                                    hwMatrix&       result);
//!
//! Computes peak to RMS values
//! \param data          signal vector or matrix
//! \param dim           dimension
//!
SIGNALS_DECLS hwMathStatus Peak2RMS(const hwMatrixN& data,
                                    int              dim,
                                    hwMatrixN&       result);
//!
//! Computes peak to peak values
//! \param data          signal vector or matrix
//! \param dim           dimension
//!
SIGNALS_DECLS hwMathStatus Peak2Peak(const hwMatrix& data,
                                     int             dim,
                                     hwMatrix&       result);
//!
//! Computes peak to peak values
//! \param data          signal vector or matrix
//! \param dim           dimension
//!
SIGNALS_DECLS hwMathStatus Peak2Peak(const hwMatrixN& data,
                                     int              dim,
                                     hwMatrixN&       result);
//!
//! Computes linear correlation of two vectors in the frequency domain
//! \param X            signal vector 1
//! \param Y            signal vector 2
//!
SIGNALS_DECLS hwMathStatus CorrLin(const hwMatrix& X,
                                   const hwMatrix& Y,
                                   int             maxlag,
                                   hwMatrix&       Rxy);
//!
//! Computes linear correlation of matrix columns in the frequency domain
//! \param X            signal matrix
//!
SIGNALS_DECLS hwMathStatus CorrLin(const hwMatrix& X,
                                   int             maxlag,
                                   hwMatrix&       Rxy);

#endif // _Signals_StatsFuncs_h
