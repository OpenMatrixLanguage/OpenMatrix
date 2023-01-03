/**
* @file DiffEqFuncs.h
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
#ifndef _DiffEq_WrapperFuncs_h
#define _DiffEq_WrapperFuncs_h

#include "DiffEqExports.h"

#include "hwArkodeWrap.h"
#include "hwCvodeWrap.h"
#include "hwIdaWrap.h"

//------------------------------------------------------------------------------
//!
//! \brief ODE/DAE Solvers
//! When needed the following variables should be declared as follows:
//! hwMatrix y(numEqn);      // holds initial conditions
//! hwMatrix yp(numEqn);     // holds initial conditions
//!
//------------------------------------------------------------------------------

DIFFEQ_DECLS hwMathStatus RK45(ARKRhsFn_client      sysfunc,
                               ARKRootFn_client     rootfunc,
                               int                  nrtfn,
                               ARKDenseJacFn_client jacDfunc,
                               const hwMatrix&      time,
                               const hwMatrix&      y,
                               hwMatrix*            timeSolution, 
                               hwMatrix&            ySolution, 
                               double               relerr      = 0.001,
                               const hwMatrix*      abserr      = nullptr, 
                               double               maxstep     = -999.0,
                               const hwMatrix*      userData    = nullptr,
                               hwMatrix*            pEventTime  = nullptr,
                               hwMatrix*            pEventFnVal = nullptr,
                               hwMatrix*            pEventIndx  = nullptr);
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
                 int                 nrtfn,
                 CVDenseJacFn_client jacDfunc,
                 const hwMatrix&     time,
                 const hwMatrix&     y,
                 hwMatrix*           timeSolution,
                 hwMatrix&           ySolution,
                 const char*         job,
                 double              reltol,
                 const hwMatrix*     abstol,
                 double              maxstep,
                 const hwMatrix*     userData,
                 hwMatrix*           pEventTime,
                 hwMatrix*           pEventFnVal,
                 hwMatrix*           pEventIndx);

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
DIFFEQ_DECLS hwMathStatus ODE113(CVRhsFn_client  sysfunc,
                                 CVRootFn_client rootfunc,
                                 int             nrtfn,
                                 const hwMatrix& time,
                                 const hwMatrix& y,
                                 hwMatrix*       timeSolution,
                                 hwMatrix&       ySolution,
                                 double          reltol   = 0.001,
                                 const hwMatrix* abstol   = nullptr,
                                 double          maxstep  = -999.0,
                                 const hwMatrix* userData = nullptr,
                                 hwMatrix*       pEventTime  = nullptr,
                                 hwMatrix*       pEventFnVal = nullptr,
                                 hwMatrix*       pEventIndx  = nullptr);

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
DIFFEQ_DECLS hwMathStatus ODE15s(CVRhsFn_client      sysfunc,
                                 CVRootFn_client     rootfunc,
                                 int                 nrtfn,
                                 CVDenseJacFn_client jacDfunc,
                                 const hwMatrix&     time,
                                 const hwMatrix&     y,
                                 hwMatrix*           timeSolution,
                                 hwMatrix&           ySolution,
                                 double              reltol   = 0.001,
                                 const hwMatrix*     abstol   = nullptr,
                                 double              maxstep  = -999.0,
                                 const hwMatrix*     userData = nullptr,
                                 hwMatrix*           pEventTime  = nullptr,
                                 hwMatrix*           pEventFnVal = nullptr,
                                 hwMatrix*           pEventIndx  = nullptr);

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
                 int                  nrtfn,
                 IDADenseJacFn_client jacDfunc,
                 const hwMatrix&      time,
                 const hwMatrix&      y,
                 const hwMatrix&      yp,
                 hwMatrix*            timeSolution,
                 hwMatrix&            ySolution,
                 const char*          job,
                 double               reltol,
                 const hwMatrix*      abstol,
                 double               maxstep,
                 const hwMatrix*      userData,
                 hwMatrix*            pEventTime,
                 hwMatrix*            pEventFnVal,
                 hwMatrix*            pEventIndx);

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
DIFFEQ_DECLS hwMathStatus DAE15i(IDAResFn_client      sysfunc,
                                 IDARootFn_client     rootfunc,
                                 int                  nrtfn,
                                 IDADenseJacFn_client jacDfunc,
                                 const hwMatrix&      time,
                                 const hwMatrix&      y,
                                 const hwMatrix&      yp,
                                 hwMatrix*            timeSolution,
                                 hwMatrix&            ySolution,
                                 double               reltol   = 0.001,
                                 const hwMatrix*      abstol   = nullptr,
                                 double               maxstep  = -999.0,
                                 const hwMatrix*      userData = nullptr,
                                 hwMatrix*            pEventTime  = nullptr,
                                 hwMatrix*            pEventFnVal = nullptr,
                                 hwMatrix*            pEventIndx  = nullptr);

#endif // _DiffEq_WrapperFuncs_h
