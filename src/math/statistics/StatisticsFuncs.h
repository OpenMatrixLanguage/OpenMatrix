/**
* @file StatisticsFuncs.h
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
#ifndef _Stats_StatisticsFuncs_h
#define _Stats_StatisticsFuncs_h

#include "StatisticsExports.h"

// forward declarations
class hwMathStatus;
class hwMersenneTwisterState;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;
typedef hwTMatrix<int, hwTComplex<int> > hwMatrixI;

//------------------------------------------------------------------------------
//!
//! \brief Statistics function
//!
//------------------------------------------------------------------------------

//!
//! Multiple linear regression
//! \param y
//! \param x
//! \param b
//! \param alpha
//! \param bci
//! \param r
//! \param rci
//! \param stats
//!
STATISTICS_DECLS hwMathStatus MultiRegress(const hwMatrix& y,
                                           const hwMatrix& X,
                                           hwMatrix&       b,
                                           double          alpha = 0.05,
                                           hwMatrix*       bci   = nullptr,
                                           hwMatrix*       r     = nullptr,
                                           hwMatrix*       rci   = nullptr,
                                           hwMatrix*       stats = nullptr);
//!
//! Exponential curve fit
//! \param x
//! \param y
//! \param coef
//! \param stats
//! \param yEst
//!
STATISTICS_DECLS hwMathStatus ExpCurveFit(const hwMatrix& x,
                                          const hwMatrix& y,
                                          hwMatrix&       coef,
                                          hwMatrix*       stats = nullptr,
                                          hwMatrix*       yEst = nullptr);
//!
//! Polynomial curve fit
//! \param x
//! \param y
//! \param order
//! \param coef
//! \param stats
//! \param yEst
//! \param scaled
//! \param mu
//! \param sigma
//!
STATISTICS_DECLS hwMathStatus PolyCurveFit(const hwMatrix& x,
                                           const hwMatrix& y,
                                           int             order,
                                           hwMatrix&       coef,
                                           hwMatrix*       stats  = nullptr,
                                           hwMatrix*       yEst   = nullptr,
                                           bool            scaled = false,
                                           double*         mu     = nullptr,
                                           double*         sigma  = nullptr);
//!
//! Remove the mean or regression line from a set of data
//! \param Y
//! \param str
//! \param Ymod
//!
STATISTICS_DECLS hwMathStatus Detrend(const hwMatrix& Y,
                                      const char*     str,
                                      hwMatrix&       Ymod);
//!
//! Error function
//! \param x
//!
STATISTICS_DECLS double ErrorFunc(double x);
//!
//! Compute the mean of a real data vector
//! \param data
//! \param xBar
//!
STATISTICS_DECLS hwMathStatus Mean(const hwMatrix& data,
                                   double&         xBar);
//!
//! Compute the means of the columns of a real matrix
//! \param A    Input
//! \param xBar
//!
STATISTICS_DECLS hwMathStatus Mean(const hwMatrix& A,
                                   hwMatrix&       xBar);
//!
//! Compute the median of a data vector
//! \param data   Input
//! \param median Median
//!
STATISTICS_DECLS hwMathStatus Median(const hwMatrix& data,
                                     double&         median);
//!
//! Compute the medians of the columns of a real matrix
//! \param A      Input
//! \param median Medians of the columns of A
//!
STATISTICS_DECLS hwMathStatus Median(const hwMatrix& A,
                                     hwMatrix&       median);
//!
//! Compute the average absolute deviation of a data vector
//! \param data   Input
//! \param avgDev Average absolute deviation
//!
STATISTICS_DECLS hwMathStatus AvgDev(const hwMatrix& data,
                                     double&         avgDev);
//!
//! Compute the standard deviation of a data vector
//! \param data      Input
//! \param stdDev    Standard deviation
//! \param sampleStd
//!
STATISTICS_DECLS hwMathStatus StdDev(const hwMatrix& data,
                                     double&         stdDev,
                                     bool            sampleStd = true);
//!
//! Compute the variance of a data vector
//! \param data      Input
//! \param variance  Variance
//! \param sampleVar
//!
STATISTICS_DECLS hwMathStatus Variance(const hwMatrix& data,
                                       double&         variance,
                                       bool            sampleVar = true);
//!
//!Compute the variances of the columns of a real matrix
//! \param A         Input
//! \param variance  Variance
//! \param sampleVar
//!
STATISTICS_DECLS hwMathStatus Variance(const hwMatrix& A,
                                       hwMatrix&       variance,
                                       bool            sampleVar = true);
//!
//! Compute the skewness of a data vector
//! \param data     Input
//! \param skewness Skewness
//!
STATISTICS_DECLS hwMathStatus Skewness(const hwMatrix& data,
                                       double&         skewness,
                                       bool            correctBias = true);
//!
//! Compute the kurtosis of a data vector
//! \param data     Input
//! \param skewness Skewness
//!
STATISTICS_DECLS hwMathStatus Kurtosis(const hwMatrix&  data,
                                       double&          kurtosis,
                                       bool             correctBias = true);
//!
//! Compute the root mean square value of a data vector
//! \param data Input
//! \param rms  Root mean square value
//!
STATISTICS_DECLS hwMathStatus RMS(const hwMatrix& data,
                                  double&         rms);
//!
//! Compute covariance matrix for the columns of a matrix
//! \param X         Input
//! \param cov       Covariance matrix
//! \param sampleVar
//!
STATISTICS_DECLS hwMathStatus Cov(const hwMatrix& X,
                                  hwMatrix&       cov,
                                  bool            sampleVar = true);
//!
//! Compute covariance matrix for the columns of two matrices
//! \param X
//! \param Y
//! \param cov         Covariance matrix
//! \param sampleCovar
//!
STATISTICS_DECLS hwMathStatus Cov(const hwMatrix& X,
                                  const hwMatrix& Y,
                                  hwMatrix&       cov,
                                  bool            sampleCovar = true);
//!
//! Compute correlation matrix for the columns of a matrix
//! \param X    Input
//! \param corr Correlation matrix  
//!
STATISTICS_DECLS hwMathStatus Corr(const hwMatrix& X,
                                   hwMatrix&       corr);
//!
//! Compute correlation matrix for the columns of two matrices
//! \param X    Input
//! \param Y    Input
//! \param corr Correlation matrix  
//!
STATISTICS_DECLS hwMathStatus Corr(const hwMatrix& X,
                                   const hwMatrix& Y,
                                   hwMatrix&       corr);
//!
//! Compute the geomtric mean of a data vector
//! \param data    Input
//! \param geoMean Geomtric mean
//!
STATISTICS_DECLS hwMathStatus GeoMean(const hwMatrix& data,
                                      double&         geoMean);
//!
//! Generate the histogram bin count for a data sample
//! \param data Input
//! \param bin
//!
STATISTICS_DECLS hwMathStatus PDF(const hwMatrix& data,
                                  hwMatrixI&      bin);
//!
//! Construct a full factorial design matrix
//! \param runs
//! \param matrix Full factorial design matrix
//!
STATISTICS_DECLS hwMathStatus FullDoe(const hwMatrixI& runs,
                                      hwMatrixI&       matrix);
//!
//! Construct a Box-Behnken design matrix
//! \param n
//! \param matrix Box-Behnken design matrix
//!
STATISTICS_DECLS hwMathStatus BBDoe(int        n,
                                    hwMatrixI& matrix);

#endif // _Stats_StatisticsFuncs_h
