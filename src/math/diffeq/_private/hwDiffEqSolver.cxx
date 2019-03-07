/**
* @file hwDiffEqSolver.cxx
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

#include "hwDiffEqSolver.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwDiffEqSolver::hwDiffEqSolver(const hwMatrix& y) 
{
    if (y.IsEmpty())
    {
        m_status(HW_MATH_ERR_EMPTYMATRIX, 1);
    }
    else if (!y.IsReal())
    {
        m_status(HW_MATH_ERR_COMPLEX, 1);
    }
    else if (!y.IsVector())
    {
        m_status(HW_MATH_ERR_VECTOR, 1);
    }
    else
    {
        m_status = m_y.Dimension(y.M(), y.N(), hwMatrix::REAL);

        if (!m_status.IsOk())
        {
            m_status.ResetArgs();
        }
        else
        {
            m_y = y;
        }
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwDiffEqSolver::~hwDiffEqSolver()
{
}
//------------------------------------------------------------------------------
// Fills output matrix with results from each step
//------------------------------------------------------------------------------
hwMathStatus hwDiffEqSolver::FillMatrix(double    tStart, 
                                        double    tStop,
                                        int       numTimes, 
                                        hwMatrix& ySolution)
{
    if (tStart > tStop)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINTERVAL, 1, 2);
    }

    int numEqns = m_y.Size();
    
    m_status = ySolution.Dimension(numTimes, numEqns, hwMatrix::REAL);

    if (!m_status.IsOk())
    {
        if (m_status.GetArg1() == 0)
        {
            m_status.SetArg1(4);
        }
        else if (m_status.GetArg1() == 1)
        {
            m_status(HW_MATH_ERR_NONPOSITIVE, 3);
        }
        else
        {
            m_status.ResetArgs();
        }
        return m_status;
    }

    // store initial conditions for output
    for (int j = 0; j < numEqns; j++)
    {
        ySolution(0, j) = m_y(j);
    }

    // compute and store output at remaining times
    int    i;
    int    flag;
    double t = tStart;
    double tout;
    for (i = 1; i < numTimes; i++)
    {
        tout = tStart + ((double) i / (double) (numTimes-1)) * (tStop - tStart);
        TakeStep(t, tout, flag);

        if (Success(flag) == true)          // was the full step completed successfully?
        {
            for (int j = 0; j < numEqns; j++)
            {
                ySolution(i, j) = m_y(j);
            }
        }
        else if (Continue(flag) == true)    // can execution continue?
        {
            --i;
        }
        else
        {
            break;
        }
    }

    // report zeros for failed cases
    for ( ; i < numTimes; i++)
    {
        for (int j = 0; j < numEqns; j++)
        {
            ySolution(i, j) = 0.0;
        }
    }

    return m_status;
}
//------------------------------------------------------------------------------
// Fills output matrix with results from each step
//------------------------------------------------------------------------------
hwMathStatus hwDiffEqSolver::FillMatrix(const hwMatrix& time,
                                        hwMatrix*       timeSolution, 
                                        hwMatrix&       ySolution)
{
    if (!time.IsReal())
    {
        return m_status(HW_MATH_ERR_COMPLEX, 1);
    }
    
    if (!time.IsVector())
    {
        return m_status(HW_MATH_ERR_VECTOR, 1);
    }

    int numTimes = time.Size();
    int numEqns  = m_y.Size();

    if (numTimes == 2)
    {
        OneStepMode = true;
        SetStopTime(time(1));   // Sundials specific

        if (!timeSolution)
        {
            return m_status(HW_MATH_ERR_NULLPOINTER, 2);
        }

        m_status = timeSolution->Dimension(numTimes, 1, hwMatrix::REAL);

        if (!m_status.IsOk())
        {
            if (m_status.GetArg1() == 0)
            {
                m_status.SetArg1(2);
            }
            else
            {
                m_status.ResetArgs();
            }
            return m_status;
        }
    }
    else
    {
        OneStepMode  = false;
        timeSolution = NULL;    // not used
    }

    m_status = ySolution.Dimension(numTimes, numEqns, hwMatrix::REAL);

    if (!m_status.IsOk())
    {
        if (m_status.GetArg1() == 0)
        {
            m_status.SetArg1(3);
        }
        else
        {
            m_status.ResetArgs();
        }
        return m_status;
    }

    if (ySolution.IsEmpty())
    {
        return m_status;
    }

    int    i;
    int    flag;
    double t = time(0);
    double tout;

    // store initial conditions for output
    if (OneStepMode)
    {
        (*timeSolution)(0, 0) = time(0);
        flag = -1;  // RK45 specific
    }
    else
    {
        flag = 1;   // RK45 specific
    }

    for (int j = 0; j < numEqns; j++)
    {
        ySolution(0, j) = m_y(j);
    }

    // compute and store output at remaining times
    for (i = 1; i < numTimes; i++)
    {
        if (OneStepMode)
        {
            tout = time(1);
        }
        else
        {
            tout = time(i);
        }

        TakeStep(t, tout, flag);

        if (Success(flag) == true)          // full interval completed
        {
            if (OneStepMode)
            {
                m_status = timeSolution->Resize(i+1, 1);

                if (!m_status.IsOk())
                {
                    m_status.SetArg1(2);
                    return m_status;
                }

                m_status = ySolution.Resize(i+1, numEqns);

                if (!m_status.IsOk())
                {
                    m_status.SetArg1(3);
                    return m_status;
                }

                (*timeSolution)(i, 0) = t;
            }

            for (int j = 0; j < numEqns; j++)
            {
                ySolution(i, j) = m_y(j);
            }
        }
        else if (Continue(flag) == true)    // more steps to take
        {
            if (OneStepMode)
            {
                m_status = timeSolution->Resize(i+1, 1);

                if (!m_status.IsOk())
                {
                    m_status.SetArg1(2);
                    return m_status;
                }

                m_status = ySolution.Resize(i+1, numEqns);

                if (!m_status.IsOk())
                {
                    m_status.SetArg1(3);
                    return m_status;
                }

                (*timeSolution)(i, 0) = t;

                for (int j = 0; j < numEqns; j++)
                {
                    ySolution(i, j) = m_y(j);
                }
                ++numTimes;
            }
        }
        else
        {
            break;
        }
    }

    if (!m_status.IsOk())
    {
        return m_status;
    }

    // report zeros for failed cases
    for ( ; i < numTimes; i++)
    {
        if (OneStepMode)
        {
            (*timeSolution)(i, 0) = 0.0;
        }

        for (int j = 0; j < numEqns; j++)
        {
            ySolution(i, j) = 0.0;
        }
    }

    return m_status;
}
