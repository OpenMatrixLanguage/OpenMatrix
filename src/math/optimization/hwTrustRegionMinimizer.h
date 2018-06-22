/**
* @file hwTrustRegionMinimizer.h
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
#ifndef _Optimization_hwTrustRegionMinimizer_h
#define _Optimization_hwTrustRegionMinimizer_h

#include "hwMatrix.h"

//------------------------------------------------------------------------------
//!
//! \class hwTrustRegionMinimizer
//! \brief Trust region base class
//!
//------------------------------------------------------------------------------
class hwTrustRegionMinimizer
{
public:
    //!
    //! Constructor
    //! \param P_
    //! \param maxIter     Max iterations
    //! \param maxFuncEval
    //! \param tolf        Function tolerance
    //! \param tolx        Step tolerance
    //!
    hwTrustRegionMinimizer(const hwMatrix& P_, 
                           int             maxIter, 
                           int             maxFuncEval,
                           double          tolf, 
                           double          tolx);
    //!
    //! Destructor
    //!
    virtual ~hwTrustRegionMinimizer();

    //!
    //! Get parameters as a vector
    //! \param Params Parameters
    //!
    virtual void GetParams(hwMatrix& Params);

    //!
    //! Returns status after checking object validity. \todo: Should be virtual?
    //!
    hwMathStatus ValidObjectCheck() const { return m_status; }
    //!
    //! Compute parameters that minimize the objective function
    //!
    hwMathStatus Compute();
    //!
    //! Gets number of iterations
    //!
    int Iterations() const { return m_numIter; }
    //!
    //! Gets objective function value
    double ObjFuncVal() const { return m_objFuncVal; }
    //!
    //! Sets design variable history 
    //! \param dv_history 
    //!
    void RequestDesignHist(hwMatrix*& dv_history) { DV_history = dv_history; }
    //!
    //! Sets objective history
    //! \param obj_history
    //!
    void RequestObjectiveHist(hwMatrix*& obj_history) { Obj_history = obj_history; }

protected:
    int          m_numParams;    //!<
    int          m_numFuncEvals; //!< number of function evals remaining
    double       m_numDeriv_eps; //!< used with numerical derivatives
    hwMatrix     P;              //!< parameters
    hwMatrix     G;    	         //!< gradient
    hwMatrix     B;	             //!< approximate Hessian
    hwMathStatus m_status;       //!< Status

    //!
    //! Get parameter at the specified index
    //! \param index Given index
    //!
    virtual double GetParam(int index) const;
    //!
    //! Set parameter at the specified index
    //! \param index Given index
    //! \param value Value to set
    //!
    virtual void SetParam(int    index, 
                          double value);
    //!
    //! Set parameters as a vector
    //! \param Params Parameters
    //!
    virtual void SetParams(const hwMatrix& Params);

    //!
    //!  Get initial estimate of the parameter vector
    //! \param Params Parameters
    //!
    virtual void GetInitialEstimate(hwMatrix& Params);
    //!
    //! Evaluate objective function
    //! \param result Output
    //!
    virtual void EvalObjectiveFunc(double& result) = 0;
    //!
    //! Evaluate gradient vector
    //!
    virtual void EvalGradient() = 0;
    //!
    //! Evaluate approximate Hessian matrix
    //!
    virtual void EvalHessian() = 0;
    //!
    //! Evaluate steepest descent step
    //! \param Gstep
    //!
    virtual void EvalSteepDescentStep(hwMatrix& Gstep) = 0;
    //!
    //! Evaluate Newton step
    //! \param Nstep
    //!
    virtual void EvalNewtonStep(hwMatrix& Nstep) = 0;
    //!
    //! Resets gradient vector if step is rejected
    //!
    virtual void ResetGradient() {}

    //!
    //! Update iteration history
    //! \param Pcand
    //! \param f
    //!
    void UpdateHistory(const hwMatrix& Pcand, 
                       double          f);
private:
    int       m_maxIter;     //!< Max number of iterations
    int       m_numIter;     //!< iteration number
    int       m_numHistPnts; //!< number of history points
    double    m_tolf;        //!< function tolerance
    double    m_tolx;        //!< step tolerance
    double    m_objFuncVal;  //!<
    hwMatrix* Obj_history;   //!< objective variable history
    hwMatrix* DV_history;    //!< design variable history

    //!
    //! Evaluate trust region quadratic model
    //! \param f
    //! \param Pstep
    //! \param m_new
    virtual void EvalQuadModelFunc(double          f, 
                                   const hwMatrix& Pstep,
                                   double&         m_new) = 0;
};

#endif // _Optimization_hwTrustRegionMinimizer_h
