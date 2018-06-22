/**
* @file DiffEqFuncs.h
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
#ifndef _DiffEq_WrapperFuncs_h
#define _DiffEq_WrapperFuncs_h

#include "DiffEqExports.h"

#include "hwRungeKutta.h"
#include "hwCvodeWrap.h"
#include "hwIdaWrap.h"
#include "cvode/cvode.h"
#include "cvode/cvode_spils.h"
#include "ida/ida.h"
#include "ida/ida_spils.h"

//------------------------------------------------------------------------------
//!
//! \brief ODE/DAE Solvers
//! When needed the following variables should be declared as follows:
//! hwMatrix y(numEqn);      // holds initial conditions
//! hwMatrix yp(numEqn);     // holds initial conditions
//! hwMatrix time(numTimes);
//! hwMatrix ySolution(numTimes, numEqn);
//!
//------------------------------------------------------------------------------

//!
//! Differential equation solver
//! \param sysfunc
//! \param tStart
//! \param tStop
//! \param numTimes
//! \param y
//! \param ySolution
//! \param relerr    Optional argument: relative error
//! \param abserr    Optional argument: absolute error
//! \param userData  Optional argument
//!
DIFFEQ_DECLS hwMathStatus RK45(RKF45Fn_client  sysfunc,
                               double          tStart, 
                               double          tStop,
                               int             numTimes, 
                               const hwMatrix& y, 
                               hwMatrix&       ySolution,
                               double          relerr   = 1.0e-3, 
                               const hwMatrix* abserr   = nullptr,
                               const hwMatrix* userData = nullptr);
//!
//! Differential equation solver
//! \param sysfunc
//! \param time
//! \param y
//! \param timeSolution
//! \param ySolution
//! \param relerr       Optional argument: relative error
//! \param abserr       Optional argument: absolute error
//! \param userData     Optional argument
//!
DIFFEQ_DECLS hwMathStatus RK45(RKF45Fn_client  sysfunc, 
                               const hwMatrix& time, 
                               const hwMatrix& y,
                               hwMatrix*       timeSolution, 
                               hwMatrix&       ySolution, 
                               double          relerr   = 1.0e-3,
                               const hwMatrix* abserr   = nullptr, 
                               const hwMatrix* userData = nullptr);
//!
//! Differential equation solver which wraps CVODE functions
//! \param sysfunc
//! \param rootfunc
//! \param jacDfunc
//! \param time
//! \param y
//! \param timeSolution
//! \param ySolution
//! \param job
//! \param relerr       Optional argument: relative error
//! \param reltol       Optional argument: absolute error
//! \param userData     Optional argument
//!
hwMathStatus ODE(CVRhsFn_client      sysfunc,
                 CVRootFn_client     rootfunc,
                 CVDenseJacFn_client jacDfunc,
                 const hwMatrix&     time,
                 const hwMatrix&     y,
                 hwMatrix*           timeSolution,
                 hwMatrix&           ySolution,
                 const char*         job,
                 double              reltol,
                 const               hwMatrix* abstol,
                 const               hwMatrix* userData);

//!
//! Differential equation solver which wraps CVODE functions
//! \param sysfunc
//! \param rootfunc
//! \param time
//! \param y
//! \param timeSolution
//! \param ySolution
//! \param relerr       Optional argument: relative error
//! \param reltol       Optional argument: absolute error
//! \param userData     Optional argument
//!
DIFFEQ_DECLS hwMathStatus ODE11(CVRhsFn_client  sysfunc,
                                CVRootFn_client rootfunc,
                                const hwMatrix& time,
                                const hwMatrix& y,
                                hwMatrix*       timeSolution,
                                hwMatrix&       ySolution,
                                double          reltol   = 1.0e-3,
                                const hwMatrix* abstol   = nullptr,
                                const hwMatrix* userData = nullptr);

//!
//! Differential equation solver which wraps CVODE functions
//! \param sysfunc
//! \param rootfunc
//! \param jacDfunc
//! \param time
//! \param y
//! \param timeSolution
//! \param ySolution
//! \param relerr       Optional argument: relative error
//! \param reltol       Optional argument: absolute error
//! \param userData     Optional argument
//!
DIFFEQ_DECLS hwMathStatus ODE22a(CVRhsFn_client      sysfunc,
                                 CVRootFn_client     rootfunc,
                                 CVDenseJacFn_client jacDfunc,
                                 const hwMatrix&     time,
                                 const hwMatrix&     y,
                                 hwMatrix*           timeSolution,
                                 hwMatrix&           ySolution,
                                 double              reltol   = 1.0e-3,
                                 const hwMatrix*     abstol   = nullptr,
                                 const hwMatrix*     userData = nullptr);

//!
//! Differential algrbraic Equation solver which wraps IDA functions
//! \param sysfunc
//! \param rootfunc
//! \param jacDfunc
//! \param time
//! \param y
//! \param yp
//! \param timeSolution
//! \param ySolution
//! \param job
//! \param reltol       Optional argument: relative error
//! \param abstol       Optional argument: absolute error
//! \param userData     Optional argument
//!
hwMathStatus DAE(IDAResFn_client      sysfunc,
                 IDARootFn_client     rootfunc,
                 IDADenseJacFn_client jacDfunc,
                 const hwMatrix&      time,
                 const hwMatrix&      y,
                 const hwMatrix&      yp,
                 hwMatrix*            timeSolution,
                 hwMatrix&            ySolution,
                 const char*          job,
                 double               reltol,
                 const hwMatrix*      abstol,
                 const hwMatrix*      userData);

//!
//! Differential equation solver which wraps CVODE functions
//! \param sysfunc
//! \param rootfunc
//! \param jacDfunc
//! \param time
//! \param y
//! \param timeSolution
//! \param ySolution
//! \param relerr       Optional argument: relative error
//! \param reltol       Optional argument: absolute error
//! \param userData     Optional argument
//!
DIFFEQ_DECLS hwMathStatus DAE11a(IDAResFn_client      sysfunc,
                                 IDARootFn_client     rootfunc,
                                 IDADenseJacFn_client jacDfunc,
                                 const hwMatrix&      time,
                                 const hwMatrix&      y,
                                 const hwMatrix&      yp,
                                 hwMatrix*            timeSolution,
                                 hwMatrix&            ySolution,
                                 double               reltol   = 1.0e-3,
                                 const hwMatrix*      abstol   = nullptr,
                                 const hwMatrix*      userData = nullptr);

#endif // _DiffEq_WrapperFuncs_h
