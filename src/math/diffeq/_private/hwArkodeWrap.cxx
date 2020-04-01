/**
* @file hwArkodeWrap.cxx
* @date February 2019
* Copyright (C) 2007-2019 Altair Engineering, Inc.
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

#include "stdlib.h"                 // to free N_Vectors
#include "hwArkodeWrap.h"
#include "nvector/nvector_serial.h"
#include "arkode/arkode_direct.h"     // access ARKDls interface

// Client function pointers
static ARKRhsFn_client sysfunc_client       = nullptr;
static ARKRootFn_client rootfunc_client     = nullptr;
static ARKDenseJacFn_client jacDfunc_client = nullptr;

//------------------------------------------------------------------------------
// ARKODE functions called by Sundials, converting SUNDIALS to built-in types
//------------------------------------------------------------------------------
static int sysfunc_ARKODE(realtype t, N_Vector y, N_Vector yp, void* f_data)
{
    // system function - yp is the output
    return sysfunc_client(t, NV_DATA_S(y), NV_DATA_S(yp), f_data);
}
//------------------------------------------------------------------------------
// Root finding function
//------------------------------------------------------------------------------
static int rootfunc_ARKODE(realtype t, N_Vector y, realtype *gout, void* g_data)
{
    // not currently in use
    // root finding function - gout is the output
    return rootfunc_client(t, NV_DATA_S(y), gout, g_data);
}
//------------------------------------------------------------------------------
// Dense Jacobian function
//------------------------------------------------------------------------------
static int jacDfunc_ARKODE(realtype  t,
                           N_Vector  y,
                           N_Vector  yp,
                           SUNMatrix J,
                           void*     jac_data,
                           N_Vector  tmp1,
                           N_Vector  tmp2,
                           N_Vector  tmp3)
{
    // dense Jacobian - J is the output
    return jacDfunc_client(SM_COLUMNS_D(J), t, NV_DATA_S(y), NV_DATA_S(yp),
           jac_data, SM_COLS_D(J));
}

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwArkWrap::hwArkWrap(ARKRhsFn_client      sysfunc,
                     ARKRootFn_client     rootfunc,
                     ARKDenseJacFn_client jacDfunc,
                     double               tin,
                     const hwMatrix&      y_,
                     double               reltol_,
                     const hwMatrix*      abstol_,
                     double               maxstep,
                     const hwMatrix*      userData)
    : hwDiffEqSolver(y_)
    , arkode_mem(nullptr)
    , y(nullptr)
    , A(nullptr)
    , LS(nullptr)
{
    if (!m_status.IsOk())
    {
        if (m_status.GetArg1() == 1)
        {
            m_status.SetArg1(8);
        }
        else
        {
            m_status.ResetArgs();
        }
        return;
    }

    if (!sysfunc)
    {
        m_status(HW_MATH_ERR_NULLPOINTER, 1);
        return;
    }

    if (reltol_ < 0.0)
    {
        m_status(HW_MATH_ERR_NONPOSITIVE, 6);
        return;
    }

    // ARKodeCreate
    arkode_mem = ARKodeCreate();

    if (Check_flag((void *)arkode_mem, "ARKodeCreate", 0))
    {
        m_status(HW_MATH_ERR_ALLOCFAILED);
        return;
    }

    // Create serial vector of length NEQ for I.C.
    y = N_VNew_Serial(m_y.Size());

    if (Check_flag((void *)y, "N_VNew_Serial", 0))
    {
        m_status(HW_MATH_ERR_ALLOCFAILED);
        return;
    }

    free(NV_DATA_S(y));                 // Modified to use input data pointer
    NV_OWN_DATA_S(y) = false;
    NV_DATA_S(y) = m_y.GetRealData();

    // Initialize and populate ARKODE object
    sysfunc_client = sysfunc;

    int flag = ARKodeInit(arkode_mem, sysfunc_ARKODE, NULL, tin, y);
    // if (Check_flag(&flag, "ARKodeInit", 1)) return(1);

    int numEqns = m_y.Size();

    // Create dense SUNMatrix for use in linear solver
    A = SUNDenseMatrix(numEqns, numEqns);
    if (Check_flag((void *)A, "SUNDenseMatrix", 0))
    {
        m_status(HW_MATH_ERR_ALLOCFAILED);
        return;
    }

    // Create dense SUNLinearSolver object for use by ARKode
    LS = SUNDenseLinearSolver(y, A);
    if (Check_flag((void *)LS, "SUNDenseLinearSolver", 0))
    {
        m_status(HW_MATH_ERR_ALLOCFAILED);
        return;
    }

    // Call ARKDlsSetLinearSolver to attach the matrix and linear solver to ARKode
    flag = ARKDlsSetLinearSolver(arkode_mem, LS, A);
    // if(Check_flag(&flag, "ARKDlsSetLinearSolver", 1)) return(1);

    // Set the Jacobian routine to jacDfunc_ARKODE (user-supplied)
    if (jacDfunc)
    {
        jacDfunc_client = jacDfunc;
        flag = ARKDlsSetJacFn(arkode_mem, (ARKDlsJacFn)jacDfunc_ARKODE);
        // if(Check_flag(&flag, "ARKDlsJacFn", 1)) return(1);
    }
    else
    {
        flag = ARKDlsSetJacFn(arkode_mem, nullptr);
        // if(Check_flag(&flag, "ARKDlsJacFn", 1)) return(1);
    }

    // Call ARKodeRootInit to specify the root function rootfunc_ARKODE with 2 components
    if (rootfunc)
    {
        rootfunc_client = rootfunc;
        flag = ARKodeRootInit(arkode_mem, 2, rootfunc_ARKODE);
        // if (Check_flag(&flag, "ARKodeRootInit", 1)) return(1);
    }

    // Set the pointer to user-defined data
    if (userData)
    {
        flag = ARKodeSetUserData(arkode_mem, (void*)userData);
        // if (Check_flag(&flag, "ARKodeSetFdata", 1)) return(1);
    }

    // Manage tolerances
    if (abstol_)
    {
        if (abstol_->IsEmpty())
        {
            double abstol = 1.0e-6;

            flag = ARKodeSStolerances(arkode_mem, reltol_, abstol);
            // if (Check_flag(&flag, "ARKodeSStolerances", 1)) return(1);
        }
        else
        {
            if (!abstol_->IsReal())
            {
                m_status(HW_MATH_ERR_COMPLEX, 7);
                return;
            }

            if (!abstol_->IsVector())
            {
                m_status(HW_MATH_ERR_VECTOR, 7);
                return;
            }

            if (abstol_->Size() == 1)
            {
                double abstol = (*abstol_)(0);

                if (abstol <= 0.0)
                {
                    m_status(HW_MATH_ERR_NONPOSITIVE, 7);
                    return;
                }

                flag = ARKodeSStolerances(arkode_mem, reltol_, abstol);
                // if (Check_flag(&flag, "ARKodeSStolerances", 1)) return(1);
            }
            else if (abstol_->Size() == numEqns)
            {
                for (int i = 0; i < numEqns; i++)
                {
                    if ((*abstol_)(i) <= 0.0)
                    {
                        m_status(HW_MATH_ERR_NONPOSITIVE, 7);
                        return;
                    }
                }

                // Create serial vector of length NEQ for abstol
                N_Vector abstol = N_VNew_Serial(numEqns);

                if (Check_flag((void *)abstol, "N_VNew_Serial", 0))
                {
                    m_status(HW_MATH_ERR_ALLOCFAILED);
                    return;
                }

                free(NV_DATA_S(abstol));             // Modify to use input data pointer
                NV_OWN_DATA_S(abstol) = false;
                NV_DATA_S(abstol) = const_cast<double*>(abstol_->GetRealData());

                flag = ARKodeSVtolerances(arkode_mem, reltol_, abstol);
                // if (Check_flag(&flag, "ARKodeSVtolerances", 1)) return(1);

                N_VDestroy_Serial(abstol);           // Free abstol vector
            }
            else
            {
                m_status(HW_MATH_ERR_ARRAYSIZE, 5, 8);
                return;
            }
        }
    }
    else
    {
        double abstol = 1.0e-6;

        flag = ARKodeSStolerances(arkode_mem, reltol_, abstol);
        // if (Check_flag(&flag, "ARKodeSStolerances", 1)) return(1);
    }

    if (maxstep > 0.0)
    {
        ARKodeSetMaxStep(arkode_mem, maxstep);
    }
    else if (maxstep != -999.0)
    {
        m_status(HW_MATH_ERR_NONPOSITIVE, 8);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwArkWrap::~hwArkWrap()
{
    if (arkode_mem)
    {
        ARKodeFree(&arkode_mem); // Free integrator memory
    }

    if (y)
    {
        N_VDestroy_Serial(y);    // Free y vector
    }

    if (A)
    {
        SUNMatDestroy(A);        // Free the matrix memory
    }

    if (LS)
    {
        SUNLinSolFree(LS);       // Free the linear solver memory
    }
}
//------------------------------------------------------------------------------
// Perform an integration step
//------------------------------------------------------------------------------
void hwArkWrap::TakeStep(double& t, double tout, int& flag)
{
    if (!arkode_mem)
    {
        flag = ARK_MEM_NULL;
        return;
    }

    // Advance y from t to tout if possible
    if (OneStepMode)
    {
        flag = ARKode(arkode_mem, tout, y, &t, ARK_ONE_STEP);
    }
    else
    {
        flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
    }
}
//------------------------------------------------------------------------------
// Returns true if successful
//------------------------------------------------------------------------------
bool hwArkWrap::Success(int flag)
{
    if (OneStepMode)
    {
        if (flag == ARK_TSTOP_RETURN)
        {
            return true;
        }
    }
    else if (flag == ARK_SUCCESS)
    {
        return true;
    }

    return false;
}
//------------------------------------------------------------------------------
// Returns true if execution can continue
//------------------------------------------------------------------------------
bool hwArkWrap::Continue(int flag)
{
    if (flag == ARK_SUCCESS)
    {
        if (OneStepMode)
        {
            return true;
        }
        return false;
    }

    if (flag < 0)
    {
        switch (flag)
        {
            case -1: m_status(HW_MATH_ERR_CVODE_WORKLOAD);  break;
            case -2: m_status(HW_MATH_ERR_CVODE_ACCURACY);  break;
            case -3: m_status(HW_MATH_ERR_CVODE_ERR);       break;
            case -4: m_status(HW_MATH_ERR_CVODE_CONV);      break;
            case -8: m_status(HW_MATH_ERR_USERFUNCFAIL, 1); break;
            case -6: m_status(HW_MATH_ERR_USERFUNCFAIL, 2); break;
            default: m_status(HW_MATH_ERR_USERFUNCFAIL);    break;
        }
        return false;
    }

    return true;
}
//------------------------------------------------------------------------------
// Set stop time for one step mode
//------------------------------------------------------------------------------
void hwArkWrap::SetStopTime(double tstop)
{
    int flag = ARKodeSetStopTime(arkode_mem, tstop);
}
//------------------------------------------------------------------------------
// Check ARKODE flag
//------------------------------------------------------------------------------
int hwArkWrap::Check_flag(void* flagvalue, const char* funcname, int opt)
{
    if (opt == 0 && !flagvalue) // opt == 0 => Sundials allocated memory
    {
        return 1;               // No memory was allocated by Sundials
    }
    else if (opt == 1)          // opt == 1 => Sundials returned a flag
    {
        int* errflag = static_cast<int *>(flagvalue);
        if (errflag && *errflag < 0)
        {
            return 1;
        }
    }
    else if (opt == 2 && !flagvalue) // opt == 2 => Sundials allocated memory
    {
        return 1;                    // No memory was allocated by Sundials
    }

    return 0;
}
