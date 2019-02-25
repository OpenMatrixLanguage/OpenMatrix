/**
* @file OptimizationFuncs.h
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
#ifndef __HWOPTIMIZATIONFUNCS_H__
#define __HWOPTIMIZATIONFUNCS_H__

//------------------------------------------------------------------------------
//!
//! \brief Optimization functions 
//!
//------------------------------------------------------------------------------

#include "hwOptimizationExports.h"

#include "hwGenericFuncFitter.h"
#include "hwTOMS748.h"
#include "hwUnconstrMinimizer.h"

typedef hwMathStatus (*BrentFunc)(double  x, 
                                  double& min);
typedef hwMathStatus (*NelderMeadFunc)(const hwMatrix& X, 
                                       double&         result);

//!
//! Nonlinear curve fit which works on surfaces and curves
//! \param pRespFunc  Response function
//! \param pJacFunc   Jacobian function
//! \param P
//! \param X
//! \param y
//! \param stats
//! \param yEst
//! \param maxIter     Max iterations
//! \param maxFuncEval
//! \param tolf
//! \param tolx
//! \param objHist
//! \param designHist
//! \param userData    User data
//!
HWOPTIMIZATION_DECLS hwMathStatus NLCurveFit(const LSqFitFunc pRespFunc,
                                             const LSqFitFunc pJacFunc,
                                             hwMatrix&        P,
                                             const hwMatrix&  X,
                                             const hwMatrix&  y,
                                             int&             maxIter,
                                             int&             maxFuncEval,
                                             hwMatrix*        stats       = nullptr,
                                             hwMatrix*        yEst        = nullptr,
                                             double           tolf        = 1.0e-6,
                                             double           tolx        = 1.0e-6,
                                             hwMatrix*        objHist     = nullptr,
                                             hwMatrix*        designHist  = nullptr,
                                             const hwMatrix*  userData    = nullptr);

//!
//! Solve a nonlinear system of equations
//! \param pRespFunc   Response function
//! \param pJacFunc    Jacobian function
//! \param P
//! \param minVal
//! \param numEqns_
//! \param maxIter     Max iterations
//! \param maxFuncEval
//! \param tolf
//! \param tolx
//! \param objHist
//! \param designHist
//! \param userData
//! 
HWOPTIMIZATION_DECLS hwMathStatus NLSolve(const LSqFitFunc pRespFunc,
                                          const LSqFitFunc pJacFunc,
                                          hwMatrix&        P,
                                          double&          minVal,
                                          int&             maxIter,
                                          int&             maxFuncEval,
                                          int              numEqns_    = -1,
                                          double           tolf        = 1.0e-6,
                                          double           tolx        = 1.0e-6,
                                          hwMatrix*        objHist     = nullptr,
                                          hwMatrix*        designHist  = nullptr,
                                          const hwMatrix*  userData    = nullptr);
//!
//! Find the minimum of a function without derivatives using Brent's method
//! \param f           Function which evaluates f(x) for any x in interval (a,b)
//! \param a           Left endpoint of initial interval
//! \param b           Right endpoint of initial interval
//! \param fa          f(a)
//! \param fb          f(b)
//! \param min
//! \param fmin        Abcissa approximating the point where f  attains a minimum
//! \param tol         Desired length of the interval of uncertainty of the final result
//! \param numIters    Number of iterations
//! \param numFunEvals
//!
HWOPTIMIZATION_DECLS hwMathStatus Brent(const BrentFunc f,
                                        double&         a,
                                        double&         b,
                                        double&         fa,
                                        double&         fb,
                                        double&         min,
                                        double&         fmin,
                                        int&            maxIter,
                                        int&            maxFuncEval, 
                                        double          tol,
                                        hwMatrix*       objHist    = nullptr,
                                        hwMatrix*       designHist = nullptr);
//!
//! Find zero of a function without derivatives using the TOMS 748 algorithm
//! \param f
//! \param a
//! \param b
//! \param fa           f(a)
//! \param fb           f(b)
//! \param root
//! \param froot 
//! \param tol
//! \param numIters
//! \param numFunEvals
//!
HWOPTIMIZATION_DECLS hwMathStatus F_zero(const TOMS748Func f,
                                         double&           a,
                                         double&           b,
                                         double&           fa,
                                         double&           fb,
                                         double&           root,
                                         double&           froot,
                                         int&              maxIter,
                                         int&              maxFuncEval,
                                         double            tol,
                                         hwMatrix*         objHist    = nullptr,
                                         hwMatrix*         designHist = nullptr);
//!
//! Finds minimum of a function without derivatives using Nelder-Mead simplex method
//! \param f
//! \param optPoint
//! \param minVal
//! \param numFunEvals
//! \param tolx
//!
HWOPTIMIZATION_DECLS hwMathStatus NelderMead(const NelderMeadFunc f,
                                             hwMatrix&            optPoint,
                                             double&              minVal,
                                             int&                 maxIter,
                                             int&                 maxFuncEval,
                                             double               tolx = 1.0e-6,
                                             hwMatrix*            objHist = nullptr,
                                             hwMatrix*            designHist = nullptr);

//!
//! Find the unconstrained minimum of a multivariate function
//! \param pRespFunc
//! \param pGradFunc
//! \param P
//! \param minVal
//! \param maxIter
//! \param maxFuncEval
//! \param tolf
//! \param tolx
//! \param objHist
//! \param designHist
//! \param userData
//!
HWOPTIMIZATION_DECLS hwMathStatus FMinUncon(const UnConMinObjFunc  pRespFunc,
                                            const UnConMinGradFunc pGradFunc,
                                            hwMatrix&              P,
                                            double&                minVal,
                                            int&                   maxIter,
                                            int&                   maxFuncEval,
                                            double                 tolf        = 1.0e-6,
                                            double                 tolx        = 1.0e-6,
                                            hwMatrix*              objHist     = nullptr,
                                            hwMatrix*              designHist  = nullptr, 
                                            const hwMatrix*        userData    = nullptr);


#endif // __HWOPTIMIZATIONFUNCS_H__
