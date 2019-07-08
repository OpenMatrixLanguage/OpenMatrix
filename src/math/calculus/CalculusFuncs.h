/**
* @file CalculusFuncs.h
* @date June, 2007
* Copyright (C) 2007-2019 Altair Engineering, Inc.  
* This file is part of the OpenMatrix Language ("OpenMatrix") software.
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
#ifndef _Calculus_WrapperFuncs_h
#define _Calculus_WrapperFuncs_h

#include "CalculusExports.h"
#include "hwGaussQuadrature.h"

typedef hwMathStatus (*QuadFunc2)(double x, double& y);

//------------------------------------------------------------------------------
//!
//! \brief Calculus toolbox functions 
//!
//------------------------------------------------------------------------------

//!
//! Returns hwMathStatus and gets the derivative using differences
//! \param x          Input x
//! \param y          Input y
//! \param derivative Output derivative
//! 
CALCULUS_DECLS hwMathStatus Derivative(const hwMatrix& x, 
                                       const hwMatrix& y, 
                                       hwMatrix&       derivative);
//!
//! Returns hwMathStatus and gets integral using trapezoidal method
//! \param x        Input x
//! \param y        Input y
//! \param integral Output integral
//! 
CALCULUS_DECLS hwMathStatus TrapZ(const hwMatrix& x, 
                                  const hwMatrix& y, 
                                  double&         integral);
//!
//! Returns hwMathStatus and gets cumulative integral using trapezoidal method
//! \param x        Input x
//! \param y        Input y
//! \param integral Output integral
//! 
CALCULUS_DECLS hwMathStatus CumTrapZ(const hwMatrix& x, 
                                     const hwMatrix& y, 
                                     hwMatrix&       integral);
//!
//! Returns hwMathStatus and gets integral using adaptive quadrature
//! \param pFunc    
//! \param a  
//! \param b 
//! \param area   Output 
//! \param count
//! \param reltol Optional relative tolerance
//! \param abs    Optional absolute tolerance
//!
CALCULUS_DECLS hwMathStatus Quad(const QuadFunc1 pFunc, 
                                 double          a, 
                                 double          b, 
                                 double&         area,
                                 int&            count, 
                                 double          reltol = 1.0e-3, 
                                 double          abstol = 1.0e-6);
//!
//! Returns hwMathStatus and gets integral using adaptive Simpson's rule
//! \param pFunc    
//! \param a  
//! \param b 
//! \param area  Output area 
//! \param count
//! \param tol   Optional tolerance
//!
CALCULUS_DECLS hwMathStatus QuadV(const QuadFunc2 pFunc, 
                                  double          a, 
                                  double          b, 
                                  double&         area,
                                  int&            count, 
                                  double          tol = 1.0e-6);

#endif // _Calculus_WrapperFuncs_h
