/**
* @file hwGaussNewtLSqFit.h
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
#ifndef _Optimization_GaussNewtLSq_h
#define _Optimization_GaussNewtLSq_h

#include "hwTrustRegionMinimizer.h"

//------------------------------------------------------------------------------
//!
//! \class hwGaussNewtLSqFit
//! \brief Least squares fitting base class 
//!
//------------------------------------------------------------------------------
class hwGaussNewtLSqFit : public hwTrustRegionMinimizer
{
public:
    //!
    //! Constructor
    //! \param P
    //! \param numEqns
    //! \param maxIter
    //! \param maxFuncEval
    //! \param tolf
    //! \param tolx
    //!
    hwGaussNewtLSqFit(const hwMatrix& P, 
                      int             numEqns, 
                      int             maxIter     = 200,
                      int             maxFuncEval = 400, 
                      double          tolf        = 1.0e-6, 
                      double          tolx        = 1.0e-6);
    //!
    //! Destructor
    //!
    virtual ~hwGaussNewtLSqFit();

protected:
    int      m_numEqns; //!< Number of equations
    hwMatrix J;         //!< Jacobian of fitting functions

    //!
    //! Evaluate fitting function
    //! \param estimated y values Output
    //!
    virtual void EvalFittedFunc(hwMatrix& y_est) = 0;
    //!
    //! Evaluate residuals function
    //! \param estimated residual values Output
    //!
    virtual void EvalResiduals(hwMatrix& residual) = 0;
    //!
    //! Evaluate Jacobian function
    //!
    virtual void EvalJacobian() = 0;
    //!
    //! Evaluate objective function
    //! \param result Output
    //!
    void EvalObjectiveFunc(double& result);
    //!
    //! Evaluate gradient vector
    //!
    void EvalGradient();
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
    //! Evaluate trust region quadratic model
    //! \param f
    //! \param Pstep
    //! \param m_new
    //!
    void EvalQuadModelFunc(double          f, 
                           const hwMatrix& Pstep, 
                           double&         m_new);
private:
    hwMatrix F;         //!< Vector of fitting function residuals
    hwMatrix JT;        //!< Transposed Jacobian
};

#endif // _Optimization_GaussNewtLSq_h
