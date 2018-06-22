/**
* @file MathUtilsFuncs.h
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
#ifndef _MathUtils_Funcs_h
#define _MathUtils_Funcs_h

#include "hwMatrix.h"
#include "MathUtilsExports.h"

//------------------------------------------------------------------------------
//!
//! \brief Miscellaneous math functions 
//!
//------------------------------------------------------------------------------

//!
//! Returns the status and the minimum value in a real matrix
//! \param A     Input matrix
//! \param value Pointer to the minimum value
//! \param index Optional pointer to the index of the minimum value 
//!
MATHUTILS_DECLS hwMathStatus Min(const hwMatrix& A, 
                                 double*         value, 
                                 int*            index = NULL);
//!
//! Returns the status and the minimum value in a complex matrix
//! \param A     Input matrix
//! \param value Pointer to the minimum value
//! \param index Optional pointer to the index of the minimum value 
//!
MATHUTILS_DECLS hwMathStatus Min(const hwMatrix& A, 
                                 hwComplex*      value, 
                                 int*            index = NULL);
//!
//! Returns the status and the minimum value of each column of a matrix
//! \param A   Input matrix
//! \param mag Pointer to the minimum values
//! \param row Optional pointer to the index of the minimum values 
//!
MATHUTILS_DECLS hwMathStatus Min(const hwMatrix& A, 
                                 hwMatrix*       mag, 
                                 hwMatrixI*      row = NULL);
//!
//! Returns the status and the minimum value from each vector of a matrix along 
//! a specified dimension
//! \param A     Input matrix
//! \param dim   Specified dimension
//! \param mag   Pointer to the minimum values
//! \param index Optional pointer to the index of the minimum values 
//!
MATHUTILS_DECLS hwMathStatus Min(const hwMatrix& A, 
                                 int             dim, 
                                 hwMatrix*       mag, 
                                 hwMatrixI*      index = NULL);
//!
//! Return the status and a matrix whose values are the minima of the arguments
//! \param A   Input matrix
//! \param b   Input value
//! \param mag Output matrix
//!
MATHUTILS_DECLS hwMathStatus Min(const hwMatrix& A, 
                                 double          b, 
                                 hwMatrix&       mag);
//!
//! Return the status and a matrix whose values are the minima of the arguments
//! \param A   Input matrix
//! \param b   Input value
//! \param mag Output matrix
//!
MATHUTILS_DECLS hwMathStatus Min(const hwMatrix&  A, 
                                 const hwComplex& b, 
                                 hwMatrix&        mag);
//!
//! Return the status and a matrix whose values are the minima of the arguments
//! \param A   Input matrix
//! \param B   Input matrix
//! \param mag Output matrix
//!
MATHUTILS_DECLS hwMathStatus Min(const hwMatrix& A,
                                 const hwMatrix& B,
                                 hwMatrix&       mag);
//!
//! Returns the status and maximum value of a real matrix
//! \param A     Input matrix
//! \param value Pointer to the maximum value
//! \param index Optional pointer to the index of the maximum value 
//!
MATHUTILS_DECLS hwMathStatus Max(const hwMatrix& A, 
                                 double*         value, 
                                 int*            index = NULL);
//!
//! Returns the status and maximum value of a complex matrix
//! \param A     Input matrix
//! \param value Pointer to the maximum value
//! \param index Optional pointer to the index of the maximum value 
//!
MATHUTILS_DECLS hwMathStatus Max(const hwMatrix& A, 
                                 hwComplex*      value, 
                                 int*            index = NULL);
//!
//! Returns the status and the maximum value of each column of a matrix
//! \param A   Input matrix
//! \param mag Output matrix with max values
//! \param row Optional argument giving the indices of the max values
//!
MATHUTILS_DECLS hwMathStatus Max(const hwMatrix& A, 
                                 hwMatrix*       mag, 
                                 hwMatrixI*      row = NULL);
//!
//! Returns the status and the maximum value from each vector of a matrix along 
//! the specified dimension
//! \param A     Input matrix
//! \param dim   Dimension
//! \param mag   Output matrix with max values
//! \param index Optional argument giving the indices of the max values
//!
MATHUTILS_DECLS hwMathStatus Max(const hwMatrix& A, 
                                 int             dim, 
                                 hwMatrix*       mag, 
                                 hwMatrixI*      index = NULL);
//!
//! Returns the status and a matrix whose values are the maxima of the arguments
//! \param A   Input matrix
//! \param b   Input value
//! \param mag Output matrix
//!
MATHUTILS_DECLS hwMathStatus Max(const hwMatrix& A, 
                                 double          b, 
                                 hwMatrix&       mag);
//!
//! Returns the status and a matrix whose values are the maxima of the arguments
//! \param A   Input matrix
//! \param b   Input value
//! \param mag Output matrix
//!
MATHUTILS_DECLS hwMathStatus Max(const hwMatrix&  A, 
                                 const hwComplex& b, 
                                 hwMatrix&        mag);
//!
//! Returns the status and a matrix whose values are the maxima of the arguments
//! \param A   Input matrix
//! \param B   Input matrix
//! \param mag Output matrix
//!
MATHUTILS_DECLS hwMathStatus Max(const hwMatrix& A, 
                                 const hwMatrix& B,
                                 hwMatrix&       mag);
//!
//! Returns the status and the minimum value of a vector in absolute value
//! \param A     Input matrix
//! \param value Output value
//! \param index Optional argument giving the index of the value in the matrix
//!
MATHUTILS_DECLS hwMathStatus AbsMin(const hwMatrix& A, 
                                    double*         value, 
                                    int*            index = NULL);
//!
//! Returns the status and the minimum value of each column of a matrix in 
//! absolute value
//! \param A     Input matrix
//! \param mag Output matrix
//! \param row Optional argument giving the indicies of the values
//!
MATHUTILS_DECLS hwMathStatus AbsMin(const hwMatrix& A,
                                    hwMatrix*       mag, 
                                    hwMatrixI*      row = NULL);
//!
//! Returns the status and the maximum value of a vector in absolute value
//! \param A     Input matrix
//! \param value Output matrix
//! \param index Optional argument giving the index
//!
MATHUTILS_DECLS hwMathStatus AbsMax(const hwMatrix& A,
                                    double*         value,
                                    int*            index = NULL);
//!
//! Returns the status and the maximum value of each column of a matrix in
//! absolute value
//! \param A   Input matrix
//! \param mag Output matrix
//! \param row Optional argument giving the index
//!
MATHUTILS_DECLS hwMathStatus AbsMax(const hwMatrix& A, 
                                    hwMatrix*       mag, 
                                    hwMatrixI*      row = NULL);
//!
//! Returns the status and peaks from a vector
//! \param A      Input matrix
//! \param mag    Output matrix
//! \param index  Gives the index of the peak(s)
//! \param option Optional argument
//!
MATHUTILS_DECLS hwMathStatus Peaks(const hwMatrix& A, 
                                   hwMatrix*       mag, 
                                   hwMatrixI*      index = NULL, 
                                   int             option = 0);
//!
//! Returns the status and the first index where A(index) = value
//! \param A     Input matrix
//! \param value Given value
//! \param index Index where value is found
//! \param tol   Optional argument giving the tolerance
//!
MATHUTILS_DECLS hwMathStatus GetIndex(const hwMatrix& A, 
                                      double          value, 
                                      int&            index, 
                                      double          tol = 0.0);
//!
//! Returns the status and the first pair, (i, j),  where A(i,j) = value
//! \param A     Input matrix
//! \param value Given value
//! \param ii    Row where the value is found
//! \param jj    Column where the value is found
//! \param tol   Optional argument giving the tolerance
//! \param byCol Optional argument specifying if columns are searched first
//!
MATHUTILS_DECLS hwMathStatus GetIndexPair(const hwMatrix& A, 
                                          double          value, 
                                          int&            ii, 
                                          int&            jj,
                                          double          tol   = 0.0, 
                                          bool            byCol = false);
//!
//! Returns the result of a bisection search on a sorted data array
//! \param data Sorted data array
//! \param npts Number of points
//! \param find Value to find
//!
MATHUTILS_DECLS int BinarySearch(const double* data,
                                 int           npts, 
                                 double        find);
//!
//! Returns the status and the bin locations to create a histogram
//! \param data Input matrix
//! \param bin  Bin locations
//!
MATHUTILS_DECLS hwMathStatus Bins(const hwMatrix& data,
                                  hwMatrix&       bin);
//!
//! Returns the status and the result of sorting a matrix by column in a given direction
//! \param unsorted Unsorted real or complex matrix
//! \param sorted   Sorted values
//! \param index    Optional argument for indices of sorted values
//! \param ascend   Optional argument, true if sorting is in ascending order
//!
MATHUTILS_DECLS hwMathStatus Sort(const hwMatrix& unsorted, 
                                  hwMatrix&       sorted, 
                                  hwMatrixI*      index = NULL,
                                  bool            ascend = true);
#endif // _MathUtils_Funcs_h
