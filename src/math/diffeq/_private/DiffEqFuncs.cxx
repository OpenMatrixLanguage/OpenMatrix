/**
* @file DiffEqFuncs.cxx
* @date June, 2007
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

#include "DiffEqFuncs.h"

//#include "hwCvodeWrap.h"
//#include "hwIdaWrap.h"

//------------------------------------------------------------------------------
// Differential equation solver
//------------------------------------------------------------------------------
/*
hwMathStatus RK45(ARKRhsFn_client      sysfunc,
                  ARKRootFn_client     rootfunc,
                  ARKDenseJacFn_client jacDfunc,
                  double               tStart,
                  double               tStop,
                  int                  numTimes, 
                  const hwMatrix&      y, 
                  hwMatrix&            ySolution,
                  double               relerr, 
                  const hwMatrix*      abserr, 
                  double               maxstep,
                  const hwMatrix*      userData)
{
    hwArkWrap rkf45(sysfunc, rootfunc, jacDfunc, tStart, y, relerr, abserr,
        maxstep, userData);

    hwMathStatus status = rkf45.GetStatus();

    if (!status.IsOk())
    {
        if (status.GetArg2() == 4)
        {
            status.SetArg2(8);
        }
        if (status.GetArg1() == 1)
        {
        }
        else if (status.GetArg1() == 2)
        {
            status.SetArg1(5);
        }
        else if (status.GetArg1() == 3)
        {
            status.SetArg1(7);
        }
        else if (status.GetArg1() == 4)
        {
            status.SetArg1(8);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    status = rkf45.FillMatrix(tStart, tStop, numTimes, ySolution);

    if (!status.IsOk())
    {
        if (status.GetArg2() == 2)
        {
            status.SetArg2(3);
        }

        switch (status.GetArg1())
        {
            case 111: status.SetArg1(1);  break;
            case 1:   status.SetArg1(2);  break;
            case 2:   status.SetArg1(4);  break;
            case 3:   status.SetArg1(4);  break;
            case 4:   status.SetArg1(6);  break;
            default:  status.ResetArgs(); break;
        }

        if (status.GetArg2() == 2)
        {
            status.SetArg2(4);
        }

        return status;
    }

    return status;
}
*/
//------------------------------------------------------------------------------
// Differential equation solver
//------------------------------------------------------------------------------
hwMathStatus RK45(ARKRhsFn_client      sysfunc,
                  ARKRootFn_client     rootfunc,
                  int                  nrtfn,
                  ARKDenseJacFn_client jacDfunc,
                  const hwMatrix&      time,
                  const hwMatrix&      y,
                  hwMatrix*            timeSolution, 
                  hwMatrix&            ySolution, 
                  double               relerr,
                  const hwMatrix*      abserr, 
                  double               maxstep,
                  const hwMatrix*      userData,
                  hwMatrix*            pEventTime,
                  hwMatrix*            pEventFnVal,
                  hwMatrix*            pEventIndx)
{
    if (!time.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 5);
    }

    if (time.IsEmpty())
    {
        ySolution.Dimension(0, 0, hwMatrix::REAL);
        return hwMathStatus();
    }

    hwArkWrap rkf45(sysfunc, rootfunc, nrtfn, jacDfunc, time(0),
                    y, relerr, abserr, maxstep, userData,
                    pEventTime, pEventFnVal, pEventIndx);

    hwMathStatus status = rkf45.GetStatus();

    if (!status.IsOk())
    {
        if (status.GetArg2() == 8)
        {
            status.SetArg2(10);
        }

        switch (status.GetArg1())
        {
            case 7:  status.SetArg1(9);  break;
            case 8:  status.SetArg1(10);  break;
            case 9:  status.SetArg1(11);  break;
            default: status.ResetArgs(); break;
        }

        return status;
    }

    status = rkf45.FillMatrix(time, timeSolution, ySolution);

    if (!status.IsOk())
    {
        switch (status.GetArg1())
        {
            case 111: status.SetArg1(1);  break;
            case 222: status.SetArg1(4);  break;
            case 333: status.SetArg1(2);  break;
            case 1:   status.SetArg1(5);  break;
            case 2:   status.SetArg1(7);  break;
            case 3:   status.SetArg1(8);  break;
            default:  status.ResetArgs(); break;
        }

        if (status.GetArg2() == 2)
        {
            status.SetArg2(7);
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Differential equation solver
//------------------------------------------------------------------------------
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
                 hwMatrix*           pEventIndx)
{
    if (!time.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 5);
    }

    if (time.IsEmpty())
    {
        ySolution.Dimension(0, 0, hwMatrix::REAL);
        return hwMathStatus();
    }

    hwCvodeWrap cvode(sysfunc, rootfunc, nrtfn, jacDfunc, time(0),
                      y, job, reltol, abstol, maxstep, userData,
                      pEventTime, pEventFnVal, pEventIndx);

    hwMathStatus status = cvode.GetStatus();

    if (!status.IsOk())
    {
        int arg1 = status.GetArg1();

        if (arg1 > 6)
        {
            status.SetArg1(arg1+2);
        }
        return status;
    }

    status = cvode.FillMatrix(time, timeSolution, ySolution);

    if (!status.IsOk())
    {
        switch (status.GetArg1())
        {
            case 111: status.SetArg1(1);  break;
            case 222: status.SetArg1(4);  break;
            case 333: status.SetArg1(2);  break;
            case 1:   status.SetArg1(5);  break;
            case 2:   status.SetArg1(7);  break;
            case 3:   status.SetArg1(8);  break;
            default:  status.ResetArgs(); break;
        }

        if (status.GetArg2() == 2)
        {
            status.SetArg2(7);
        }
        else if (status.GetArg2() == 3)
        {
            status.SetArg2(8);
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Differential equation solver, uses Adams method, functional iteration
//------------------------------------------------------------------------------
hwMathStatus ODE11(CVRhsFn_client  sysfunc,
                   CVRootFn_client rootfunc,
                   int             nrtfn,
                   const hwMatrix& time,
                   const hwMatrix& y,
                   hwMatrix*       timeSolution,
                   hwMatrix&       ySolution,
                   double          reltol,
                   const hwMatrix* abstol,
                   double          maxstep,
                   const hwMatrix* userData,
                   hwMatrix*       pEventTime,
                   hwMatrix*       pEventFnVal,
                   hwMatrix*       pEventIndx)
{
    const char*  job    = "11";
    hwMathStatus status = ODE(sysfunc, rootfunc, nrtfn, (CVDenseJacFn_client) nullptr,
                              time, y, timeSolution, ySolution, job, reltol, abstol,
                              maxstep, userData, pEventTime, pEventFnVal, pEventIndx);

    if (!status.IsOk())
    {
        int arg1 = status.GetArg1();
        int arg2 = status.GetArg2();

        if (arg1 > 4)
        {
            status.SetArg1(arg1-1);
        }

        if (arg2 > 9)
        {
            status.SetArg2(arg2-2);
        }
        else if (arg2 > 4)
        {
            status.SetArg2(arg2-1);
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Differential equation solver, uses BDF method, Newton iteration, dense linear solver
//------------------------------------------------------------------------------
hwMathStatus ODE22a(CVRhsFn_client      sysfunc,
                    CVRootFn_client     rootfunc,
                    int                 nrtfn,
                    CVDenseJacFn_client jacDfunc,
                    const hwMatrix&     time,
                    const hwMatrix&     y,
                    hwMatrix*           timeSolution,
                    hwMatrix&           ySolution,
                    double              reltol,
                    const hwMatrix*     abstol,
                    double              maxstep,
                    const hwMatrix*     userData,
                    hwMatrix*           pEventTime,
                    hwMatrix*           pEventFnVal,
                    hwMatrix*           pEventIndx)
{
    const char*  job    = "22a";
    hwMathStatus status = ODE(sysfunc, rootfunc, nrtfn, jacDfunc, time, y,
                              timeSolution, ySolution, job, reltol, abstol,
                              maxstep, userData, pEventTime, pEventFnVal, pEventIndx);

    if (!status.IsOk())
    {
        int arg1 = status.GetArg1();
        int arg2 = status.GetArg2();

        if (arg1 > 9)
        {
            status.SetArg1(arg1-1);
        }
        if (arg2 > 9)
        {
            status.SetArg2(arg2-1);
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Differential algebraic equation solver, wraps IDA functions
//------------------------------------------------------------------------------
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
                 hwMatrix*            pEventIndx)
{
    if (!time.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 5);
    }
    if (time.IsEmpty())
    {
        ySolution.Dimension(0, 0, hwMatrix::REAL);
        return hwMathStatus();
    }

    hwIdaWrap ida(sysfunc, rootfunc, nrtfn, jacDfunc, time(0), y,
                  yp, job, reltol, abstol, maxstep, userData,
                  pEventTime, pEventFnVal, pEventIndx);

    hwMathStatus status = ida.GetStatus();

    if (!status.IsOk())
    {
        int arg1 = status.GetArg1();
        int arg2 = status.GetArg2();

        if (arg1 > 7)
        {
            status.SetArg1(arg1+2);
        }
        if (arg2 > 7)
        {
            status.SetArg2(arg2+2);
        }
        return status;
    }

    status = ida.FillMatrix(time, timeSolution, ySolution);

    if (!status.IsOk())
    {
        switch (status.GetArg1())
        {
            case 111: status.SetArg1(1);  break;
            case 222: status.SetArg1(4);  break;
            case 333: status.SetArg1(2);  break;
            case 1:   status.SetArg1(5);  break;
            case 2:   status.SetArg1(8);  break;
            case 3:   status.SetArg1(9);  break;
            default:  status.ResetArgs(); break;
        }

        if (status.GetArg2() == 2)
        {
            status.SetArg2(8);
        }
        else if (status.GetArg2() == 3)
        {
            status.SetArg2(9);
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Differential algebraic equation solver, uses BDF method, Newton iteration, dense linear solver
//------------------------------------------------------------------------------
hwMathStatus DAE11a(IDAResFn_client      sysfunc,
                    IDARootFn_client     rootfunc,
                    int                  nrtfn,
                    IDADenseJacFn_client jacDfunc,
                    const hwMatrix&      time,
                    const hwMatrix&      y,
                    const hwMatrix&      yp,
                    hwMatrix*            timeSolution,
                    hwMatrix&            ySolution,
                    double               reltol,
                    const hwMatrix*      abstol,
                    double               maxstep,
                    const hwMatrix*      userData,
                    hwMatrix*            pEventTime,
                    hwMatrix*            pEventFnVal,
                    hwMatrix*            pEventIndx)
{
    const char*  job    = "11a";
    hwMathStatus status = DAE(sysfunc, rootfunc, nrtfn, jacDfunc, time, y, yp,
                              timeSolution, ySolution, job, reltol, abstol,
                              maxstep, userData, pEventTime, pEventFnVal, pEventIndx);

    if (!status.IsOk())
    {
        int arg1 = status.GetArg1();
        int arg2 = status.GetArg2();

        if (arg1 > 10)
        {
            status.SetArg1(arg1-1);
        }
        if (arg2 > 10)
        {
            status.SetArg2(arg2-1);
        }
    }

    return status;
}
