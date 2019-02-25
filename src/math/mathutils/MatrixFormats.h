/**
* @file MatrixFormats.h
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
#ifndef _MathUtils_MatrixFormats_h
#define _MathUtils_MatrixFormats_h

//------------------------------------------------------------------------------
//!
//! \brief Functions to store LAPACK / Sundials special matrix formats 
//!
//------------------------------------------------------------------------------

#include "hwMatrix.h"
#include "MathUtilsExports.h"

//!
//! Returns the status after sizing a dense matrix
//! \param m        Number of rows
//! \param n        Number of columns
//! \param A        Dense matrix
//! \param initZero Optional argument if true, initializes elements of A to 0
//!
MATHUTILS_DECLS hwMathStatus SizeDenseMatrix(int         m, 
                                             int         n, 
                                             hwMatrix&   A, 
                                             bool        initZero = false);
//!
//! Returns the status after sizing a band matrix
//! \param n        
//! \param kl       
//! \param ku       
//! \param A        Given matrix to resize
//! \param initZero Optional argument if true, initializes elements of A to 0
//!
MATHUTILS_DECLS hwMathStatus SizeBandMatrix(int         n, 
                                            int         kl, 
                                            int         ku, 
                                            hwMatrix&   A, 
                                            bool        initZero = false);
//!
//! Returns the status after sizing a symmetric band matrix
//! \param n        
//! \param kd       
//! \param A        Given matrix to resize
//! \param initZero Optional argument if true, initializes elements of A to 0
//!
MATHUTILS_DECLS hwMathStatus SizeSymBandMatrix(int         n, 
                                               int         kd, 
                                               hwMatrix&   A, 
                                               bool        initZero = false);
//!
//! Returns status after setting an element in a dense matrix
//! \param i     Row index
//! \param j     Column index
//! \param value Value to set
//! \param A     Matrix to set value in
//!
MATHUTILS_DECLS hwMathStatus SetDenseMatrixElem(int       i, 
                                                int       j, 
                                                double    value,
                                                hwMatrix& A);
//!
//! Returns status after setting an element in a band matrix
//! \param i     Row index
//! \param j     Column index
//! \param value Value to set
//! \param A     Matrix to set value in
//! \param kl
//! \param ku
//!
MATHUTILS_DECLS
hwMathStatus SetBandMatrixElem(int       i, 
                               int       j, 
                               double    value,
                               hwMatrix& A, 
                               int       kl, 
                               int       ku);
//!
//! Returns status after setting an element in a symmetric band matrix
//! \param i     Row index
//! \param j     Column index
//! \param value Value to set
//! \param A     Matrix to set value in
//! \param kd
//!
MATHUTILS_DECLS hwMathStatus SetSymBandMatrixElem(int       i, 
                                                  int       j, 
                                                  double    value,
                                                  hwMatrix& A, 
                                                  int       kd);

#endif // _MathUtils_MatrixFormats_h
