/**
* @file hwGenericFuncFitter.h
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
#ifndef _Optimization_hwGenericFuncFitter_h
#define _Optimization_hwGenericFuncFitter_h

#include "hwGaussNewtLSq.h"

//!
//! \typedef LSqFitFunc
//! Usage: 1. For least squares curve/surface fitting, X is the domain variable 
//! to fit y=f(X). The goal is to find parameters(P) that the provide the best fit
//! 2. For solving nonlinear systems of equations of the form f(P)=0, pass a 
//! dummy variable for X
//!
typedef hwMathStatus (*LSqFitFunc)(const hwMatrix& P,
                                   const hwMatrix& X,
                                   const hwMatrix* userData,
                                   hwMatrix&       result);

//------------------------------------------------------------------------------
//!
//! \class hwGenericFuncFitter
//! \brief Generic function fitting class 
//!
//------------------------------------------------------------------------------
class hwGenericFuncFitter : public hwGaussNewtLSqFit
{
public:
    // \todo: The two constructors could be separated into two classes, 
    // hwNLCurveFit and hwNLSolve.

    //!
    //! Constructor - Find parameters when fitting equations of the form y=f(x)
    //! \param pObjFunc    Objective function
    //! \param pJacFunc    Jacobian function
    //! \param X_
    //! \param y_
    //! \param lowerBound
    //! \param upperBound
    //! \param maxIter     Max iterations
    //! \param maxFuncEval Max function evaluations
    //! \param tolf
    //! \param tolx
    //! \param userData
    //!
    hwGenericFuncFitter(const LSqFitFunc pObjFunc,
                        const LSqFitFunc pJacFunc,
                        const hwMatrix&  P,
                        const hwMatrix&  X_,
                        const hwMatrix&  y_,
                        const hwMatrix*  lowerBound,
                        const hwMatrix*  upperBound,
                        int              maxIter     = 200,
                        int              maxFuncEval = 400, 
                        double           tolf        = 1.0e-6,
                        double           tolx        = 1.0e-6,
                        const hwMatrix*  userData    = NULL);
    //!
    //! Constructor - Find parameters when solving f(P)=0
    //! \param pObjFunc    Objective function
    //! \param pJacFunc    Jacobian function
    //! \param P
    //! \param numEqns     Number of equations
    //! \param maxIter     Max iterations
    //! \param maxFuncEval
    //! \param tolf
    //! \param tolx
    //! \param userData
    //!
    hwGenericFuncFitter(const LSqFitFunc pObjFunc,
                        const LSqFitFunc pJacFunc,
                        const hwMatrix&  P,
                        int              numEqns,
                        int              maxIter     = 200,
                        int              maxFuncEval = 400, 
                        double           tolf        = 1.0e-6,
                        double           tolx        = 1.0e-6,
                        const hwMatrix*  userData    = NULL);
    //!
    //! Destructor
    //!
    virtual ~hwGenericFuncFitter();
    //!
    //! Compute parameters that minimize the objective function
    //!
    hwMathStatus Compute();

protected:
    //!
    //! Evaluate fitted function
    //! \param y_est
    //!
    void EvalFittedFunc(hwMatrix& y_est);
    //!
    //! Evaluate residuals
    //! \param residual
    //!
    void EvalResiduals(hwMatrix& residual);
    //!
    //! Evaluate Jacobian matrix
    //!
    void EvalJacobian();
    //!
    //! Evaluate approximate Hessian matrix
    //!
    void EvalHessian();
    //!
    //! Evaluate steepest descent step
    //! \param Gstep
    //!
    void EvalSteepDescentStep(hwMatrix& Gstep);
    //!
    //! Evaluate Newton step
    //! \param Nstep
    //!
    void EvalNewtonStep(hwMatrix& Nstep);
    //!
    //! Scale Matrix
    //!
    void ScaleMatrix();
    //!
    //! Enforce Bounds
    //! \param Pcand
    //! \param Pstep
    //!
    void EnforceBounds(const hwMatrix& Pcand, hwMatrix& Pstep);

private:
    LSqFitFunc      m_pObjFunc;    //!< objective function
    LSqFitFunc      m_pJacFunc;    //!< Jacobian - if NULL use numerical derivatices
    const hwMatrix* X;             //!< domain variable(s) for curve/surface fitting
    const hwMatrix* y;             //!< range variable for curve/surface fitting
    const hwMatrix* m_lowerBound;  //!< lower parameter bounds
    const hwMatrix* m_upperBound;  //!< upper parameter bounds
    const hwMatrix* m_userData;    //!< User data
    hwMatrix        D;             //!< Scale Matrix
    hwMatrix        DinvSqrt;      //!< Scale Matrix
    hwMatrix        M;             //!< Modified Hessian
};

#endif // _Optimization_hwGenericFuncFitter_h
