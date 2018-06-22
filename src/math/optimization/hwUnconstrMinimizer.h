/**
* @file hwUnconstrMinimizer.h
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
#ifndef _Optimization_hwUnconstrMinimizer_h
#define _Optimization_hwUnconstrMinimizer_h

#include "hwTrustRegionMinimizer.h"

//!
//! \typedef UnConMinObjFunc
//!
typedef hwMathStatus (*UnConMinObjFunc)(const hwMatrix& P, 
                                        const hwMatrix* userData, 
                                        double&         min);
//!
//! \typedef UnConMinGradFunc
//!
typedef hwMathStatus (*UnConMinGradFunc)(const hwMatrix& P, 
                                         const hwMatrix* userData, 
                                         hwMatrix&       grad);

//------------------------------------------------------------------------------
//!
//! \class hwUnconstrMinimizer
//! \brief Unconstrained minimization class
//!
//------------------------------------------------------------------------------
class hwUnconstrMinimizer : public hwTrustRegionMinimizer
{
public:
    //!
    //! Constructor
    //! \param pObjFunc    Objective function
    //! \param pGradFunc   Gradient
    //! \param P
    //! \param maxIter     Max iterations
    //! \param maxFuncEval
    //! \param tolf
    //! \param tolx
    //! \param userData    User data
    //! 
    hwUnconstrMinimizer(const UnConMinObjFunc  pObjFunc, 
                        const UnConMinGradFunc pGradFunc,
                        const hwMatrix&        P, 
                        int                    maxIter     = 200, 
                        int                    maxFuncEval = 400,
                        double                 tolf        = 1.0e-6, 
                        double                 tolx        = 1.0e-6,
                        const hwMatrix*        userData    = NULL);
    //!
    //! Destructor
    //!
    virtual ~hwUnconstrMinimizer();

    //!
    //! Reinitializes and called if object is used more than once
    //!
    void ReInitialize() { initStep = true; } 

protected:
    bool             initStep;      //!< True if this is the initial step
    hwMatrix         Gk;            //!< gradient at previous step
    hwMatrix         Gb;	        //!< backup of Gk, used when step is rejected
    hwMatrix         Pk;	        //!< paramter estimates at previous step
    hwMatrix         Bk;	        //!< Hessian approx at previous step
    const hwMatrix*  m_userData;    //!< user Data
    UnConMinObjFunc  m_pObjFunc;    //!< objective function
    UnConMinGradFunc m_pGradFunc;   //!< Gradient - if NULL use numerical derivatives

    //!
    //! Evaluates objective function
    //! \param result Result
    //!
    void EvalObjectiveFunc(double& result);
    //! 
    //! Evaluates gradient vector
    //!
    void EvalGradient();
    //!
    //! Evaluate approximate Hessian matrix
    //!
    void EvalHessian();
    //!
    //! Evaluates steepest descent strp
    //! \param GStep
    //!
    void EvalSteepDescentStep(hwMatrix& Gstep);
    //!
    //! Evaluate Newton step
    //! \param Nstep 
    //!
    void EvalNewtonStep(hwMatrix& Nstep);
    //!
    //! Evaluate trust region quadratic model
    //! \param f
    //! \param Pstep
    //! \param m_new
    //!
    void EvalQuadModelFunc(double          f, 
                           const hwMatrix& Pstep, 
                           double&         m_new);
    //!
    //! Resets gradient if step was rejected
    //!
    void ResetGradient() { G = Gb; }
};

#endif // _Optimization_hwUnconstrMinimizer_h
