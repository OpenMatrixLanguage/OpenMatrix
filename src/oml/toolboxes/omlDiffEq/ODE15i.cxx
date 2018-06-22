/**
* @file ODE15i.cxx
* @date January 2015
* Copyright (C) 2017-2018 Altair Engineering, Inc.  
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

#include "OmlDiffEqTboxFuncs.h"

#include "BuiltInFuncsUtils.h"
#include "FunctionInfo.h"
#include "OML_Error.h"
#include "StructData.h"

#include "DiffEqFuncs.h"

// file scope variables and functions
static int                 ODE15i_sys_size;
static FunctionInfo*       ODE15i_sys_func;
static FunctionInfo*       ODE15i_JAC_func;
static std::string         ODE15i_sys_name;
static FUNCPTR             ODE15i_sys_pntr;
static FUNCPTR             ODE15i_JAC_pntr;
static EvaluatorInterface* ODE15i_eval_ptr;

//------------------------------------------------------------------------------
// Wrapper for dae system function called by ode algorithm from oml scripts
//------------------------------------------------------------------------------
static int ODE15i_file_func(double  t, 
                            double* y, 
                            double* yp, 
                            double* r, 
                            void*   userData)
{
    std::vector<Currency> inputs;
    inputs.push_back(t);

    hwMatrix* yMatrix = nullptr;

    if (ODE15i_sys_func && ODE15i_sys_func->Nargin() > 1)
    {
        yMatrix = new hwMatrix(ODE15i_sys_size, 1, hwMatrix::REAL);

        for (int i = 0; i < ODE15i_sys_size; ++i)
        {
            (*yMatrix)(i) = y[i];
        }
        inputs.push_back(yMatrix);
    }

    hwMatrix* ypMatrix = nullptr;

    if (ODE15i_sys_func && ODE15i_sys_func->Nargin() > 2)
    {
        ypMatrix = new hwMatrix(ODE15i_sys_size, 1, hwMatrix::REAL);

        for (int i = 0; i < ODE15i_sys_size; ++i)
        {
            (*ypMatrix)(i) = yp[i];
        }
        inputs.push_back(ypMatrix);
    }

    Currency  result;

    if (ODE15i_sys_func)
    {
        result = ODE15i_eval_ptr->CallInternalFunction(ODE15i_sys_func, inputs);
    }
    else if (ODE15i_sys_pntr)
    {
        result = ODE15i_eval_ptr->CallFunction(ODE15i_sys_name, inputs);
    }
    else
    {
        return -1;
    }

    if (result.IsMatrix())
    {
        const hwMatrix* rMatrix = result.Matrix();
        if (!rMatrix || !rMatrix->IsReal() || rMatrix->Size() != ODE15i_sys_size)
        {
            return -1;
        }
        for (int i = 0; i < ODE15i_sys_size; ++i)
        {
            r[i] = (*rMatrix)(i);
        }
    }
    else if (result.IsScalar())
    {
        if (ODE15i_sys_size != 1)
        {
            return -1;
        }
        r[0] = result.Scalar();
    }
    else
    {
        return -1;
    }

    return 0;
}
//------------------------------------------------------------------------------
// Wrapper for Jacobian of ode system function and is called by dae algorithm
//------------------------------------------------------------------------------
static int ODE15i_jac_file_func(long int N, 
                                double   t, 
                                double   cj, 
                                double*  y,
                                double*  yp, 
                                double*  r, 
                                void*    jac_data, 
                                double** Jac)
{
    std::vector<Currency> inputs;
    inputs.push_back(t);

    hwMatrix* yMatrix = nullptr;

    if (ODE15i_JAC_func && ODE15i_JAC_func->Nargin() > 1)
    {
        yMatrix = new hwMatrix(ODE15i_sys_size, 1, hwMatrix::REAL);

        for (int i = 0; i < ODE15i_sys_size; ++i)
        {
            (*yMatrix)(i) = y[i];
        }
        inputs.push_back(yMatrix);
    }

    hwMatrix* ypMatrix = nullptr;

    if (ODE15i_JAC_func && ODE15i_JAC_func->Nargin() > 2)
    {
        ypMatrix = new hwMatrix(ODE15i_sys_size, 1, hwMatrix::REAL);

        for (int i = 0; i < ODE15i_sys_size; ++i)
        {
            (*ypMatrix)(i) = yp[i];
        }
        inputs.push_back(ypMatrix);
    }

    std::vector<Currency> outputs;

    if (ODE15i_JAC_func)
    {
        outputs = ODE15i_eval_ptr->DoMultiReturnFunctionCall(ODE15i_JAC_func, 
            inputs, static_cast<int>(inputs.size()), 2, true);
    }
    else if (ODE15i_JAC_pntr)
    {
        outputs = ODE15i_eval_ptr->DoMultiReturnFunctionCall(ODE15i_JAC_pntr, 
            inputs, static_cast<int>(inputs.size()), 2, true);
    }
    if (outputs.size() != 2)
    {
        return -1;
    }

    Currency dFdy  = outputs[0];
    Currency dFdyp = outputs[1];

    if (dFdy.IsScalar() && dFdyp.IsScalar())
    {
        if (N != 1)
        {
            return -1;
        }
        Jac[0][0] = dFdy.Scalar() + cj * dFdyp.Scalar();
    }
    else if (dFdy.IsMatrix() && dFdyp.IsMatrix())
    {
        const hwMatrix* dFdyMtx  = dFdy.Matrix();
        const hwMatrix* dFdypMtx = dFdyp.Matrix();

        if (dFdyMtx->M() != N || dFdyMtx->N() != N)
        {
            return -1;
        }
        if (dFdypMtx->M() != N || dFdypMtx->N() != N)
        {
            return -1;
        }
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Jac[i][j] = (*dFdyMtx)(j, i) + cj * (*dFdypMtx)(j, i);
            }
        }
    }
    else
    {
        return -1;
    }

    return 0;
}
//------------------------------------------------------------------------------
// Solves a system of stiff differential algebraic equations [ode15i]
//------------------------------------------------------------------------------
bool OmlOde15i(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    int nargin = eval.GetNarginValue();

    if (nargin < 4)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    std::string JacFuncName;
    
    FunctionInfo* funcInfo2 = nullptr;
    FUNCPTR       funcPntr2 = nullptr;

    if (!inputs[0].IsFunctionHandle())
    {
        throw OML_Error(OML_ERR_HANDLE, 1, OML_VAR_TYPE);
    }

    std::string   funcName  = inputs[0].FunctionHandle()->FunctionName();
    FunctionInfo* funcInfo1 = nullptr;
    FUNCPTR       funcPntr1 = nullptr;

    if (funcName == "anonymous")
    {
        funcInfo1 = inputs[0].FunctionHandle();
    }
    else if (!eval.FindFunctionByName(funcName, &funcInfo1, &funcPntr1))
    {
        throw OML_Error(OML_ERR_FUNCNAME, 1);
    }

    if (funcInfo1 && (funcInfo1->Parameters().size() < 1 || 
        funcInfo1->Parameters().size() > 3))
    {
        throw OML_Error(OML_ERR_NUMARGIN, 1);
    }

    if (!inputs[1].IsScalar() && !inputs[1].IsMatrix())
    {
        throw OML_Error(OML_ERR_SCALARVECTOR, 2);
    }
    if (!inputs[2].IsScalar() && !inputs[2].IsMatrix())
    {
        throw OML_Error(OML_ERR_SCALARVECTOR, 3);
    }
    if (!inputs[3].IsScalar() && !inputs[3].IsMatrix())
    {
        throw OML_Error(OML_ERR_SCALARVECTOR, 4);
    }
    const hwMatrix*  time         = inputs[1].ConvertToMatrix();
    const hwMatrix*  y            = inputs[2].ConvertToMatrix();
    double           reltol       = 1.0e-3;
    const hwMatrix*  abstol       = nullptr;
    bool             deleteAbsTol = false;
    const hwMatrix*  yp           = inputs[3].ConvertToMatrix();

    if (nargin > 4)
    {
        if (inputs[4].IsStruct())
        {
            StructData* options = inputs[4].Struct();

            if (options->N() != 1)
            {
                throw OML_Error(OML_ERR_OPTION, 5);
            }
            const Currency& reltol_C   = options->GetValue(0, -1, "RelTol");
            const Currency& abstol_C   = options->GetValue(0, -1, "AbsTol");
            const Currency& jacobian_C = options->GetValue(0, -1, "Jacobian");

            if (!reltol_C.IsEmpty())
            {
                if (reltol_C.IsScalar())
                {
                    reltol = reltol_C.Scalar();
                }
                else
                {
                    throw OML_Error(OML_ERR_SCALARVECTOR, 5, OML_VAR_ABSTOL);
                }
            }

            if (abstol_C.IsEmpty())
            {
                abstol = EvaluatorInterface::allocateMatrix(1, 1, 1.0e-6);
                deleteAbsTol = true;
            }
            else if (abstol_C.IsScalar() || abstol_C.IsMatrix())
            {
                abstol = abstol_C.ConvertToMatrix();
            }
            else
            {
                throw OML_Error(OML_ERR_SCALARVECTOR, 4, OML_VAR_ABSTOL);
            }

            if (!jacobian_C.IsEmpty())
            {
                if (jacobian_C.IsFunctionHandle())
                {
                    JacFuncName = jacobian_C.FunctionHandle()->FunctionName();

                    if (funcName == "anonymous")
                    {
                        funcInfo2 = jacobian_C.FunctionHandle();
                    }
                    else
                    {
                        if (!eval.FindFunctionByName(JacFuncName, &funcInfo2, &funcPntr2))
                        {
                            throw OML_Error(HW_ERROR_INVFUNCNAME);
                        }
                    }

                    if (funcInfo2 && (funcInfo2->Parameters().size() < 1 || funcInfo2->Parameters().size() > 3))
                    {
                        throw OML_Error(OML_ERR_NUMARGIN, 5);
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_OPTIONVAL, 5, OML_VAR_JACOBIAN);
                }
            }
            else
            {
                funcInfo2 = nullptr;
            }
        }
        else if (inputs[4].IsMatrix())
        {
            if (inputs[4].Matrix()->M() != 0 || inputs[4].Matrix()->M() != 0)
            {
                throw OML_Error(OML_ERR_HANDLE_EMPTY, 5);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_HANDLE_EMPTY, 5);
        }
    }

    if (nargin > 5)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    // set file scope variables
    ODE15i_sys_func = funcInfo1;
    ODE15i_JAC_func = funcInfo2;
    ODE15i_sys_name = funcName;
    ODE15i_sys_pntr = funcPntr1;
    ODE15i_JAC_pntr = funcPntr2;
    ODE15i_sys_size = y->Size();
    ODE15i_eval_ptr = &eval;

    // call algorithm
    hwMatrix*        timeSolution = nullptr;
    hwMatrix*        ySolution    = new hwMatrix;
    hwMatrix*        userData     = nullptr;
    IDARootFn_client rootFunc     = nullptr;
    hwMathStatus     status;

    if (time->Size() == 2)
    {
        // one step mode
        timeSolution = new hwMatrix;

        if (!ODE15i_JAC_func)
        {
            status = DAE11a(ODE15i_file_func, rootFunc, 
                (IDADenseJacFn_client) NULL, *time, *y, *yp, timeSolution, 
                *ySolution, reltol, abstol, userData);
        }
        else
        {
            status = DAE11a(ODE15i_file_func, rootFunc, ODE15i_jac_file_func,
                *time, *y, *yp, timeSolution, *ySolution, reltol, abstol, userData);
        }
    }
    else
    {
        // interval mode
        timeSolution = (hwMatrix*) time;
        timeSolution->IncrRefCount();

        if (!ODE15i_JAC_func)
        {
            status = DAE11a(ODE15i_file_func, rootFunc, (IDADenseJacFn_client) nullptr, 
                *time, *y, *yp, nullptr, *ySolution, reltol, abstol, userData);
        }
        else
        {
            status = DAE11a(ODE15i_file_func, rootFunc, ODE15i_jac_file_func,
                *time, *y, *yp, nullptr, *ySolution, reltol, abstol, userData);
        }
    }

    if (deleteAbsTol)
    {
        delete abstol;
        abstol = nullptr;
    }

    if (!status.IsOk())
    {
        if (time->Size() == 2)
        {
            delete timeSolution;
            timeSolution = nullptr;
        }
        else
        {
            timeSolution->DecrRefCount();
        }

        if (status.GetArg1() == 111)
        {
            status.SetUserFuncName(funcName);
            status.SetArg1(1);
        }
        else if (status.GetArg1() == 222)
        {
            status.SetUserFuncName(JacFuncName);
            status.SetArg1(3);
        }
        else if (status.GetArg1() == 4)
        {
            status.SetArg1(2);
        }
        else if (status.GetArg1() == 5)
        {
            status.SetArg1(3);
        }
        else if (status.GetArg1() == 6)
        {
            status.SetArg1(4);
        }
        else if (status.GetArg1() > 8 && status.GetArg1() < 11)
        {
            status.SetArg1(5);
        }
        else
        {
            status.ResetArgs();
        }
        if (status.GetArg2() == 5)
        {
            status.SetArg2(3);
        }
        else if (status.GetArg2() == 6)
        {
            status.SetArg2(4);
        }
        else if (status.GetArg2() > 8 && status.GetArg2() < 11)
        {
            status.SetArg2(5);
        }
        if (status.IsWarning())
        {
            BuiltInFuncsUtils::SetWarning(eval, status.GetMessageString());
        }
        else
        {
            delete ySolution;
            ySolution = nullptr;
            throw OML_Error(status);
        }
    }

    // pack outputs
    outputs.push_back(timeSolution);
    outputs.push_back(ySolution);

    return true;
}
