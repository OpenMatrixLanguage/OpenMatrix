/**
* @file StatUtilFuncs.h
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
#ifndef _StatUtilFuncs_h
#define _StatUtilFuncs_h

#include <vector>
#include <string>
#include "StatisticsExports.h"

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

/*
namespace std
{
    template<typename T> class allocator;
    template<typename T, typename A = allocator<T> > class vector;
}
*/

//------------------------------------------------------------------------------
//!
//! \brief Statistics utility functions
//!
//------------------------------------------------------------------------------

//!
//! Returns the variance of a strided vector
//! \param x         Input vector
//! \param y         Input vector
//! \param n         Vector size
//! \param sampleCov Optional argument
//!
double Variance(const double* x,
                int           stride,
                int           n,
                bool          sampleVar = true);
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
//! Returns the covariance of a two strided vectors
//! \param x         Input vector
//! \param stride1   Spacing between elements
//! \param y         Input vector
//! \param stride2   Spacing between elements
//! \param n         Vector size
//! \param sampleCov Optional argument
//!
double Covariance(const double* pt1,
                  int           stride1,
                  const double* pt2,
                  int           stride2,
                  int           n,
                  bool          sampleVar = true);

enum MovingWindowEndType { SHRINK, DISCARD, SAME, USERVAL, PERIODIC };

//!
//! Compute the moving median of a strided vector
//! \param x           Input vector
//! \param n           Vector size
//! \param stride      Spacing between elements
//! \param window      Moving window
//! \param nb          Number of window points below the data index
//! \param unsortedIdx Indexing work vector
//! \param sortedIdx   Indexing work vector
//! \param median      Moving medians
//!
void MovingMedianHelper(const double*       data,
                        int                 n,
                        int                 stride,
                        hwMatrix&           window,
                        int                 nb,
                        bool                includeNaN,
                        MovingWindowEndType endtype,
                        double              userVal,
                        std::vector<int>&   unsortedIdx,
                        std::vector<int>&   sortedIdx,
                        double*             median);
//!
//! Compute the squared Euclidean distance between two ND points
//! \param pt1     Input point 1
//! \param stride1 Spacing between elements
//! \param pt2     Input point 2
//! \param stride2 Spacing between elements
//! \param n       Number of dimensions
//!
STATISTICS_DECLS double SqEuclidDistance(const double* pt1,
                                         int           stride1,
                                         const double* pt2,
                                         int           stride2,
                                         int           n);
//!
//! Compute the Euclidean distance between two ND points
//! \param pt1     Input point 1
//! \param stride1 Spacing between elements
//! \param pt2     Input point 2
//! \param stride2 Spacing between elements
//! \param n       Number of dimensions
//!
STATISTICS_DECLS double EuclidDistance(const double* pt1,
                                       int           stride1,
                                       const double* pt2,
                                       int           stride2,
                                       int           n);
//!
//! Compute the Standarized Euclidean distance between two ND points
//! \param pt1     Input point 1
//! \param stride1 Spacing between elements
//! \param pt2     Input point 2
//! \param stride2 Spacing between elements
//! \param n       Number of dimensions
//!
STATISTICS_DECLS double StdEuclidDistance(const double* pt1,
                                          int           stride1,
                                          const double* pt2,
                                          int           stride2,
                                          int           n);
//!
//! Compute the City Block distance between two ND points
//! \param pt1     Input point 1
//! \param stride1 Spacing between elements
//! \param pt2     Input point 2
//! \param stride2 Spacing between elements
//! \param n       Number of dimensions
//!
STATISTICS_DECLS double CityBlockDistance(const double* pt1,
                                          int           stride1,
                                          const double* pt2,
                                          int           stride2,
                                          int           n);
//!
//! Compute the Cosine distance between two ND points
//! \param pt1     Input point 1
//! \param stride1 Spacing between elements
//! \param pt2     Input point 2
//! \param stride2 Spacing between elements
//! \param n       Number of dimensions
//!
STATISTICS_DECLS double CosineDistance(const double* pt1,
                                       int           stride1,
                                       const double* pt2,
                                       int           stride2,
                                       int           n);
//!
//! Compute the Correlation distance between two ND points
//! \param pt1     Input point 1
//! \param stride1 Spacing between elements
//! \param pt2     Input point 2
//! \param stride2 Spacing between elements
//! \param n       Number of dimensions
//!
STATISTICS_DECLS double CorrelationDistance(const double* pt1,
                                            int           stride1,
                                            const double* pt2,
                                            int           stride2,
                                            int           n);
//!
//! Compute the Hamming distance between two ND points
//! \param pt1     Input point 1
//! \param stride1 Spacing between elements
//! \param pt2     Input point 2
//! \param stride2 Spacing between elements
//! \param n       Number of dimensions
//!
STATISTICS_DECLS double HammingDistance(const double* pt1,
                                        int           stride1,
                                        const double* pt2,
                                        int           stride2,
                                        int           n);
//!
//! Compute the Jaccard distance between two ND points
//! \param pt1     Input point 1
//! \param stride1 Spacing between elements
//! \param pt2     Input point 2
//! \param stride2 Spacing between elements
//! \param n       Number of dimensions
//!
STATISTICS_DECLS double JaccardDistance(const double* pt1,
                                        int           stride1,
                                        const double* pt2,
                                        int           stride2,
                                        int           n);
//!
//! Compute the Chebyshev distance between two ND points
//! \param pt1     Input point 1
//! \param stride1 Spacing between elements
//! \param pt2     Input point 2
//! \param stride2 Spacing between elements
//! \param n       Number of dimensions
//!
STATISTICS_DECLS double ChebyshevDistance(const double* pt1,
                                          int           stride1,
                                          const double* pt2,
                                          int           stride2,
                                          int           n);
//!
//! Compute the Minkowski distance between two ND points
//! \param pt1     Input point 1
//! \param stride1 Spacing between elements
//! \param pt2     Input point 2
//! \param stride2 Spacing between elements
//! \param n       Number of dimensions
//! \param stride  Spacing between dimensions
//!
STATISTICS_DECLS double MinkowskiDistance(const double* pt1,
                                          int           stride1,
                                          const double* pt2,
                                          int           stride2,
                                          int           n,
                                          int           p);
//!
//! Compute the L1 distance between two ND points
//! \param pt1     Input point 1
//! \param stride1 Spacing between elements
//! \param pt2     Input point 2
//! \param stride2 Spacing between elements
//! \param n       Number of dimensions
//!
STATISTICS_DECLS double L1Distance(const double* pt1,
                                   int           stride1,
                                   const double* pt2,
                                   int           stride2,
                                   int           n);
//!
//! Compute the Chi-squared distance between two ND points
//! \param pt1     Input point 1
//! \param stride1 Spacing between elements
//! \param pt2     Input point 2
//! \param stride2 Spacing between elements
//! \param n       Number of dimensions
//!
STATISTICS_DECLS double ChiSqDistance(const double* pt1,
                                      int           stride1,
                                      const double* pt2,
                                      int           stride2,
                                      int           n);
//!
//! Compute the Earth Mover's distance between two ND points
//! \param pt1     Input point 1
//! \param stride1 Spacing between elements
//! \param pt2     Input point 2
//! \param stride2 Spacing between elements
//! \param n       Number of dimensions
//!
STATISTICS_DECLS double EarthMoversDistance(const double* pt1,
                                            int           stride1,
                                            const double* pt2,
                                            int           stride2,
                                            int           n);
//!
//! Computes the total sum of squares of a data vector and returns status
//! \param data Data matrix
//! \param sst  Total sum of squares
//!
STATISTICS_DECLS hwMathStatus TotalSumOfSquares(const hwMatrix& data, 
                                                double&         sst);

#endif // _StatUtilFuncs_h
