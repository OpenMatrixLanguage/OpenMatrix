/**
* @file hwIdaWrap.cxx
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

#include "hwIdaWrap.h"

#include "stdlib.h"                 // to free N_Vectors

#include "nvector/nvector_serial.h"
#include "ida/ida_direct.h"         // access CVDls interface

// Client function pointers
static IDAResFn_client sysfunc_client       = nullptr;
static IDARootFn_client rootfunc_client     = nullptr;
static IDADenseJacFn_client jacDfunc_client = nullptr;

//------------------------------------------------------------------------------
// IDA functions called by SUNDIALS, converting from SUNDIALS to built-in types
//------------------------------------------------------------------------------
static int sysfunc_IDA(realtype t, 
                       N_Vector y, 
                       N_Vector yp, 
                       N_Vector r,
                       void*    r_data)
{
    // system function - r is the output
    return sysfunc_client(t, NV_DATA_S(y), NV_DATA_S(yp), NV_DATA_S(r), r_data);
}
//------------------------------------------------------------------------------
// Root finding function - currently not in use
//------------------------------------------------------------------------------
static int rootfunc_IDA(realtype  t,
                        N_Vector  y, 
                        N_Vector  yp, 
                        realtype* gout,
                        void*     g_data)
{
    // root finding function - gout is the output
    return rootfunc_client(t, NV_DATA_S(y), NV_DATA_S(yp), gout, g_data);
}
//------------------------------------------------------------------------------
// Dense Jacobian function
//------------------------------------------------------------------------------
static int jacDfunc_IDA(realtype  t, 
                        realtype  cj, 
                        N_Vector  y,
                        N_Vector  yp, 
                        N_Vector  r, 
                        SUNMatrix J, 
                        void*     jac_data,
                        N_Vector  tempv1, 
                        N_Vector  tempv2, 
                        N_Vector  tempv3)
{
    // dense Jacobian - J is the output
    return jacDfunc_client(SM_COLUMNS_D(J), t, cj, NV_DATA_S(y), NV_DATA_S(yp),
                           NV_DATA_S(r), jac_data, SM_COLS_D(J));
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwIdaWrap::hwIdaWrap(IDAResFn_client      sysfunc, 
                     IDARootFn_client     rootfunc, 
                     IDADenseJacFn_client jacDfunc,
                     double               tin, 
                     const hwMatrix&      y_, 
                     const hwMatrix&      yp_, 
                     const char*          job,
                     double               reltol_, 
                     const hwMatrix*      abstol_, 
                     const hwMatrix*      userData)
    : hwDiffEqSolver(y_)
    , ida_mem       (nullptr)
    , y             (nullptr)
    , yp            (nullptr)
    , A             (nullptr)
    , LS            (nullptr)
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

    if (!job)
    {
        m_status(HW_MATH_ERR_NULLPOINTER, 7);
        return;
    }

    if (reltol_ < 0.0)
    {
        m_status(HW_MATH_ERR_NONPOSITIVE, 8);
        return;
    }

    if (yp_.Size() != m_y.Size())
    {
        m_status(HW_MATH_ERR_ARRAYSIZE, 5, 6);
        return;
    }

    m_status = m_yp.Dimension(yp_.M(), yp_.N(), hwMatrix::REAL);

    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
    }
    else
    {
        m_yp = yp_;
    }

    // IDACreate
    if (job[0] == '1' && job[1] == '1')
    {
        ida_mem = IDACreate();
    }
    else
    {
        m_status(HW_MATH_ERR_INVALIDINPUT);
        return;
    }

    if (Check_flag((void *)ida_mem, "IDACreate", 0))
    {
        m_status(HW_MATH_ERR_ALLOCFAILED);
        return;
    }

    // Create serial vectors of length NEQ for I.C.
    y = N_VNew_Serial(m_y.Size());

    if (Check_flag((void *)y, "N_VNew_Serial", 0))
    {
        m_status(HW_MATH_ERR_ALLOCFAILED);
        return;
    }

    free(NV_DATA_S(y));                 // Modified to use input data pointer
    NV_OWN_DATA_S(y) = false;
    NV_DATA_S(y) = m_y.GetRealData();

    yp = N_VNew_Serial(m_yp.Size());

    if (Check_flag((void *)yp, "N_VNew_Serial", 0))
    {
        m_status(HW_MATH_ERR_ALLOCFAILED);
        return;
    }

    free(NV_DATA_S(yp));                 // Modified to use input data pointer
    NV_OWN_DATA_S(yp) = false;
    NV_DATA_S(yp) = m_yp.GetRealData();

    // Initialize and populate IDA object
    sysfunc_client = sysfunc;

    int flag = IDAInit(ida_mem, sysfunc_IDA, tin, y, yp);
    // if (Check_flag(&flag, "IDAInit", 1)) return(1);

    int numEqns = m_y.Size();

    if (job[2] == 'a')
    {
        // Create dense SUNMatrix for use in linear solver
        A = SUNDenseMatrix(numEqns, numEqns);
        if (Check_flag((void *)A, "SUNDenseMatrix", 0))
        {
            m_status(HW_MATH_ERR_ALLOCFAILED);
            return;
        }

        // Create dense SUNLinearSolver object for use by CVode
        LS = SUNDenseLinearSolver(y, A);
        if (Check_flag((void *)LS, "SUNDenseLinearSolver", 0))
        {
            m_status(HW_MATH_ERR_ALLOCFAILED);
            return;
        }

        // Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode
        flag = IDADlsSetLinearSolver(ida_mem, LS, A);
        // if(Check_flag(&flag, "IDADlsSetLinearSolver", 1)) return(1);

        // Set the Jacobian routine to Jac (user-supplied)
        if (jacDfunc)
        {
            jacDfunc_client = jacDfunc;
            flag = IDADlsSetJacFn(ida_mem, (IDADlsJacFn) jacDfunc_IDA);
            // if(Check_flag(&flag, "IDADlsSetJacFn", 1)) return(1);
        }
        else
        {
            flag = IDADlsSetJacFn(ida_mem, nullptr);
            // if(Check_flag(&flag, "IDADlsSetJacFn", 1)) return(1);
        }
    }

    // Call IDARootInit to specify the root function rootfunc_IDA with 2 components
    if (rootfunc)
    {
        rootfunc_client = rootfunc;
        flag = IDARootInit(ida_mem, 2, rootfunc_IDA);
        // if (Check_flag(&flag, "IDARootInit", 1)) return(1);
    }

    // Set the pointer to user-defined data
    if (userData)
    {
        flag = IDASetUserData(ida_mem, (void*) userData);
        // if (Check_flag(&flag, "IDASetUserData", 1)) return(1);
    }

    // Manage tolerances
    if (abstol_)
    {
        if (abstol_->IsEmpty())
        {
            double abstol = 1.0e-6;
            flag = IDASStolerances(ida_mem, reltol_, abstol);
            // if (Check_flag(&flag, "IDASStolerances", 1)) return(1);
        }
        else
        {
            if (!abstol_->IsReal())
            {
                m_status(HW_MATH_ERR_COMPLEX, 9);
                return;
            }

            if (!abstol_->IsVector())
            {
                m_status(HW_MATH_ERR_VECTOR, 9);
                return;
            }

            if (abstol_->Size() == 1)
            {
                double abstol = (*abstol_)(0);

                if (abstol <= 0.0)
                {
                    m_status(HW_MATH_ERR_NONPOSITIVE, 9);
                    return;
                }

                flag = IDASStolerances(ida_mem, reltol_, abstol);
                // if (Check_flag(&flag, "IDASStolerances", 1)) return(1);
            }
            else if (abstol_->Size() == numEqns)
            {
                int numEqns = abstol_->Size();

                for (int i = 0; i < numEqns; i++)
                {
                    if ((*abstol_)(i) <= 0.0)
                    {
                        m_status(HW_MATH_ERR_NONPOSITIVE, 9);
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

                flag = IDASVtolerances(ida_mem, reltol_, abstol);
                // if (Check_flag(&flag, "IDASVtolerances", 1)) return(1);

                N_VDestroy_Serial(abstol);           // Free abstol vector
            }
            else
            {
                m_status(HW_MATH_ERR_ARRAYSIZE, 5, 9);
                return;
            }
        }
    }
    else
    {
        double abstol = 1.0e-6;

        flag = IDASStolerances(ida_mem, reltol_, abstol);
        // if (Check_flag(&flag, "IDASStolerances", 1)) return(1);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwIdaWrap::~hwIdaWrap()
{
    if (ida_mem)
    {
        IDAFree(&ida_mem);     // Free integrator memory
    }

    if (y)
    {
        N_VDestroy_Serial(y);  // Free y vector
    }

    if (yp)
    {
        N_VDestroy_Serial(yp); // Free yp vector
    }

    if (A)
    {
        SUNMatDestroy(A);     // Free the matrix memory
    }

    if (LS)
    {
        SUNLinSolFree(LS);   // Free the linear solver memory
    }
}
//------------------------------------------------------------------------------
// Perform an integration step
//------------------------------------------------------------------------------
void hwIdaWrap::TakeStep(double& t, double tout, int& flag)
{
    if (!ida_mem)
    {
        flag = IDA_MEM_NULL;
        return;
    }

    // Advance y from t to tout if possible
    if (OneStepMode)
    {
        flag = IDASolve(ida_mem, tout, &t, y, yp, IDA_ONE_STEP);
    }
    else
    {
        flag = IDASolve(ida_mem, tout, &t, y, yp, IDA_NORMAL);
    }
}
//------------------------------------------------------------------------------
// Returns true if successful
//------------------------------------------------------------------------------
bool hwIdaWrap::Success(int flag)
{
    if (OneStepMode)
    {
        if (flag == IDA_TSTOP_RETURN)
        {
            return true;
        }
    }
    else if (flag == IDA_SUCCESS)
    {
        return true;
    }

    return false;
}
//------------------------------------------------------------------------------
// Returns true if execution can continue
//------------------------------------------------------------------------------
bool hwIdaWrap::Continue(int flag)
{
    if (flag == IDA_SUCCESS)
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
            case -4:  m_status(HW_MATH_ERR_IDA_WORKLOAD);    break;
            case -5:  m_status(HW_MATH_ERR_IDA_ACCURACY);    break;
            case -6:  m_status(HW_MATH_ERR_IDA_ERR);         break;
            case -7:  m_status(HW_MATH_ERR_IDA_CONV);        break;
            case -11: m_status(HW_MATH_ERR_USERFUNCFAIL, 1); break;
            case -9:  m_status(HW_MATH_ERR_USERFUNCFAIL, 2); break;
            default:  m_status(HW_MATH_ERR_USERFUNCFAIL);    break;
        }
        return false;
    }

    return true;
}
//------------------------------------------------------------------------------
// Set stop time for one step mode
//------------------------------------------------------------------------------
void hwIdaWrap::SetStopTime(double tstop)
{
    int flag = IDASetStopTime(ida_mem, tstop);
}
//------------------------------------------------------------------------------
// Check IDA flag
//------------------------------------------------------------------------------
int hwIdaWrap::Check_flag(void* flagvalue, const char* funcname, int opt)
{
    /*
    * Check function return value...
    *   opt == 0 means SUNDIALS function allocates memory so check if
    *            returned NULL pointer
    *   opt == 1 means SUNDIALS function returns a flag so check if
    *            flag >= 0
    *   opt == 2 means function allocates memory so check if returned
    *            NULL pointer 
    */

    int *errflag;

    if (opt == 0 && !flagvalue) // Sundials allocated memory
    {
        return 1;
    }
    else if (opt == 1)          // Sundials function returned flag
    {
        int* errflag = static_cast<int *>(flagvalue);
        if (errflag && *errflag < 0)
        {
            return 1;
        }
    }
    else if (opt == 2 && !flagvalue) // Sundials allocated memory
    {
        return 1;
    }

    return 0;
}
