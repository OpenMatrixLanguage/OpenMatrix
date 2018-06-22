/**
* @file StatUtilFuncs.h
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
* This file is part of the OpenMatrix Language (“OpenMatrix”) software.
* Open Source License Information:
* OpenMatrix is free software. You can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* OpenMatrix is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
* You should have received a copy of the GNU Affero General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
* 
* Commercial License Information: 
* For a copy of the commercial license terms and conditions, contact the Altair Legal Department at Legal@altair.com and in the subject line, use the following wording: Request for Commercial License Terms for OpenMatrix.
* Altair’s dual-license business model allows companies, individuals, and organizations to create proprietary derivative works of OpenMatrix and distribute them - whether embedded or bundled with other software - under a commercial license agreement.
* Use of Altair’s trademarks and logos is subject to Altair's trademark licensing policies.  To request a copy, email Legal@altair.com and in the subject line, enter: Request copy of trademark and logo usage policy.
*/
#ifndef _StatUtilFuncs_h
#define _StatUtilFuncs_h

#include "StatisticsExports.h"

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

//------------------------------------------------------------------------------
//!
//! \brief Statistics utility functions
//!
//------------------------------------------------------------------------------

//!
//! Returns the covariance of two vectors
//! \param x         Input vector
//! \param y         Input vector
//! \param n         Vector size
//! \param sampleCov Optional argument
//!
double Covariance(const double* x, 
                  const double* y,
                  int           n, 
                  bool          sampleCov = true);
//!
//! Computes the covariance of two vectors and returns status
//! \param X         Input
//! \param Y         Input
//! \param covar     Covariance of the input vectors
//! \param sampleCov Optional argument
//!
STATISTICS_DECLS hwMathStatus Covariance(const hwMatrix& X, 
                                         const hwMatrix& Y, 
                                         double&         covar, 
                                         bool            sampleCov = true);
//!
//! Computes the total sum of squares of a data vector and returns status
//! \param data Data matrix
//! \param sst  Total sum of squares
//!
STATISTICS_DECLS hwMathStatus TotalSumOfSquares(const hwMatrix& data, 
                                                double&         sst);

#endif // _StatUtilFuncs_h
