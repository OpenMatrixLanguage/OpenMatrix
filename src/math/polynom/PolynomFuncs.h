/**
* @file PolynomFuncs.h
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
#ifndef _Polynom_WrapperFuncs_h
#define _Polynom_WrapperFuncs_h

#include "PolynomExports.h"

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTComplex<double> hwComplex;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

template <typename T1, typename T2> class hwTMatrixN;
typedef hwTMatrixN<double, hwTComplex<double> > hwMatrixN;


//------------------------------------------------------------------------------
//!
//! \brief Polynomial functions
//!
//------------------------------------------------------------------------------

//!
//! Computes polynomial roots and returns status
//! \param coef  Coefficients in descending order
//! \param roots Polynomial roots
//!
POLYNOM_DECLS hwMathStatus PolyRoots(const hwMatrix& coef, 
                                     hwMatrix&       roots);

//!
//! Evaluates a polynomial and returns status
//! \param coef Coefficients in ascending order
//! \param x
//! \param y    Real value
//!
POLYNOM_DECLS hwMathStatus PolyVal(const hwMatrix& coef, 
                                   double          x, 
                                   double&         y);
//!
//! Evaluates a polynomial and returns status
//! \param coef Coefficients in ascending order
//! \param x
//! \param y    Complex value
//!
POLYNOM_DECLS hwMathStatus PolyVal(const hwMatrix&  coef, 
                                   const hwComplex& x,
                                   hwComplex&       y);
//!
//! Evaluates a polynomial and returns status
//! \param coef Coefficients in ascending order
//! \param x
//! \param y    Output matrix
//!
POLYNOM_DECLS hwMathStatus PolyVal(const hwMatrix& coef, 
                                   const hwMatrix& x, 
                                   hwMatrix&       y);
//!
//! Performs polynomial division or deconvolution and returns status
//! \param A Input matrix
//! \param B Input matrix
//! \param Q Output matrix
//! \param R Output matrix
//!
POLYNOM_DECLS hwMathStatus PolyDivide(const hwMatrix& A, 
                                      const hwMatrix& B,
                                      hwMatrix&       Q, 
                                      hwMatrix&       R);
//!
//! Computes the derivative of a polynomial and returns status
//! \param A Input
//! \param D Output
//!
POLYNOM_DECLS hwMathStatus PolyDer(const hwMatrix& A, 
                                   hwMatrix&       D);
//!
//! Computes the derivative of a polynomial product and returns status
//! \param A Input
//! \param B 
//! \param D Output
//!
POLYNOM_DECLS hwMathStatus PolyDer(const hwMatrix& A, 
                                   const hwMatrix& B, 
                                   hwMatrix&       D);
//!
//! Computes the derivative of a quotient and returns status
//! \param A Input
//! \param B 
//! \param P 
//! \param Q
//!
POLYNOM_DECLS hwMathStatus PolyDer(const hwMatrix& A, 
                                   const hwMatrix& B, 
                                   hwMatrix&       P, 
                                   hwMatrix&       Q);
//!
//! Computes the integral of a polynomial and returns status
//! \param A Input
//! \param I Output integral
//! \param k Optional argument
//!
POLYNOM_DECLS hwMathStatus PolyInt(const hwMatrix& A,
                                   hwMatrix&       I, 
                                   double          k = 0.0);
//!
//! Performs linear interpolation and returns status
//! \param x_old 
//! \param y_old 
//! \param x_new 
//! \param y_new
//! \param extrap Optional argument
//!
POLYNOM_DECLS hwMathStatus LinearInterp(const hwMatrix& x_old, 
                                        const hwMatrix& y_old,
                                        const hwMatrix& x_new, 
                                        hwMatrix&       y_new, 
                                        bool            extrap = false);
//!
//! Performs piecewise cubic hermite interpolation and returns status
//! \param x_old 
//! \param y_old 
//! \param x_new 
//! \param y_new
//! \param extrap Optional argument
//!
POLYNOM_DECLS hwMathStatus PchipInterp(const hwMatrix& x_old, 
                                       const hwMatrix& y_old,
                                       const hwMatrix& x_new, 
                                       hwMatrix&       y_new,
                                       bool            extrap = false);
//!
//! Performs knot-a-not cubic spline interpolation and returns status
//! \param x_old 
//! \param y_old 
//! \param x_new 
//! \param y_new
//! \param extrap Optional argument
//!
POLYNOM_DECLS hwMathStatus Spline(const hwMatrix& x_old, 
                                  const hwMatrix& y_old,
                                  const hwMatrix& x_new, 
                                  hwMatrix&       y_new,
                                  bool            extrap = false);

//!
//! Computes knot-a-not cubic spline coefficients and returns status
//! \param x_old 
//! \param y_old 
//! \param coefs Coefficients
//!
POLYNOM_DECLS hwMathStatus Spline(const hwMatrix& x_old, 
                                  const hwMatrix& y_old,
                                  hwMatrix&       coefs);
//!
//! Performs clamped cubic spline interpolation and returns status
//! \param x_old 
//! \param y_old 
//! \param fp1
//! \param fp2
//! \param x_new
//! \param y_new
//! \param extrap Optional argument
//! 
POLYNOM_DECLS hwMathStatus Spline(const hwMatrix& x_old, 
                                  const hwMatrix& y_old,
                                  double          fp1,
                                  double          fp2, 
                                  const hwMatrix& x_new,
                                  hwMatrix&       y_new, 
                                  bool            extrap = false);
//!
//! Computes clamped cubic spline coefficients and returns status
//! \param x_old 
//! \param y_old 
//! \param fp1
//! \param fp2
//! \param coefs Coefficients
//!
POLYNOM_DECLS hwMathStatus Spline(const hwMatrix& x_old, 
                                  const hwMatrix& y_old,
                                  double          fp1, 
                                  double          fp2, 
                                  hwMatrix&       coefs);
//!
//! Performs bilinear interpolation and returns status
//! \param x_old 
//! \param y_old 
//! \param z_old 
//! \param x_new 
//! \param y_new
//! \param z_new 
//! \param extrap Optional argument
//!
POLYNOM_DECLS hwMathStatus BilinearInterp(const hwMatrix& x_old, 
                                          const hwMatrix& y_old, 
                                          const hwMatrix& z_old,
                                          const hwMatrix& x_new, 
                                          const hwMatrix& y_new, 
                                          hwMatrix&       z_new,
                                          bool            extrap = false);
//!
//! Performs spline surface interpolation and returns status
//! \param x_old 
//! \param y_old 
//! \param z_old 
//! \param x_new
//! \param y_new
//! \param z_new
//! \param extrap Optional argument
//! 
POLYNOM_DECLS hwMathStatus Spline2D(const hwMatrix& x_old, 
                                    const hwMatrix& y_old, 
                                    const hwMatrix& z_old,
                                    const hwMatrix& x_new, 
                                    const hwMatrix& y_new, 
                                    hwMatrix&       z_new,
                                    bool            extrap = false);
//!
//! Computes gradient of trilinear interpolation and returns status
//! \param x_old 
//! \param y_old 
//! \param z_old 
//! \param v_new
//! \param grad
//! \param extrap Optional argument
//! 

POLYNOM_DECLS hwMathStatus TrilinearInterpGrad(const hwMatrix&  x_old,
                                               const hwMatrix&  y_old,
                                               const hwMatrix&  z_old,
                                               const hwMatrixN& val_old,
                                               const hwMatrix&  v_new,
                                               hwMatrix&        grad);

#endif // _Polynom_WrapperFuncs_h
