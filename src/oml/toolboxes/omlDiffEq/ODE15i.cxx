/**
* @file ODE15i.cxx
* @date January 2015
* Copyright (C) 2017-2018 Altair Engineering, Inc.  
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

#include "OmlDiffEqTboxFuncs.h"

#include "BuiltInFuncsUtils.h"
#include "FunctionInfo.h"
#include "OML_Error.h"
#include "StructData.h"
#include "hwMathException.h"

#include "DiffEqFuncs.h"

// file scope variables and functions
static EvaluatorInterface* ODE15i_eval_ptr;
static int                 ODE15i_sys_size;
static FunctionInfo*       ODE15i_sys_func;
static FunctionInfo*       ODE15i_JAC_func;
static std::string         ODE15i_sys_name;
static FUNCPTR             ODE15i_sys_pntr;
static FUNCPTR             ODE15i_JAC_pntr;
static int                 ODE15i_event_size;
static FunctionInfo*       ODE15i_event_func;
static FUNCPTR             ODE15i_event_pntr;

//------------------------------------------------------------------------------
// Wrapper for DAE system function called by DAE algorithm from oml scripts
//------------------------------------------------------------------------------
static int ODE15i_file_func(double  t, 
                            double* y, 
                            double* yp, 
                            double* r, 
                            void*   userData)
{
    std::vector<Currency> inputs;
    inputs.push_back(t);

    if (ODE15i_sys_func)
    {
        if (ODE15i_sys_func->Nargin() > 1)
        {
            hwMatrix* yMatrix = new hwMatrix(ODE15i_sys_size, 1, hwMatrix::REAL);

            for (int i = 0; i < ODE15i_sys_size; ++i)
            {
                (*yMatrix)(i) = y[i];
            }
            inputs.push_back(yMatrix);
        }

        if (ODE15i_sys_func->Nargin() > 2)
        {
            hwMatrix* ypMatrix = new hwMatrix(ODE15i_sys_size, 1, hwMatrix::REAL);

            for (int i = 0; i < ODE15i_sys_size; ++i)
            {
                (*ypMatrix)(i) = yp[i];
            }
            inputs.push_back(ypMatrix);
        }
    }

    Currency result;

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
// Wrapper for Jacobian of DAE system function and is called by DAE algorithm
//------------------------------------------------------------------------------
static int ODE15i_jac_file_func(int64_t  N,
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

    if (ODE15i_JAC_func)
    {
        if (ODE15i_JAC_func->Nargin() > 1)
        {
            hwMatrix* yMatrix = new hwMatrix(ODE15i_sys_size, 1, hwMatrix::REAL);

            for (int i = 0; i < ODE15i_sys_size; ++i)
            {
                (*yMatrix)(i) = y[i];
            }

            inputs.push_back(yMatrix);
        }

        if (ODE15i_JAC_func->Nargin() > 2)
        {
            hwMatrix* ypMatrix = new hwMatrix(ODE15i_sys_size, 1, hwMatrix::REAL);

            for (int i = 0; i < ODE15i_sys_size; ++i)
            {
                (*ypMatrix)(i) = yp[i];
            }

            inputs.push_back(ypMatrix);
        }
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
// Wrapper for DAE event function in oml scripts and is called by DAE algorithm
//------------------------------------------------------------------------------
static int ODE15i_event_file_func(double  t,
                                  double* y,
                                  double* yp,
                                  double* gout,
                                  void*   userData)
{
    std::vector<Currency> inputs;
    inputs.push_back(t);
    std::vector<Currency> outputs;

    if (ODE15i_event_func)
    {
        if (ODE15i_event_func->Nargin() > 1)
        {
            hwMatrix* yMatrix = new hwMatrix(ODE15i_sys_size, 1, hwMatrix::REAL);

            for (int i = 0; i < ODE15i_sys_size; ++i)
            {
                (*yMatrix)(i) = y[i];
            }

            inputs.push_back(yMatrix);
        }

        if (ODE15i_event_func->Nargin() > 2)
        {
            hwMatrix* ypMatrix = new hwMatrix(ODE15i_sys_size, 1, hwMatrix::REAL);

            for (int i = 0; i < ODE15i_sys_size; ++i)
            {
                (*ypMatrix)(i) = yp[i];
            }

            inputs.push_back(ypMatrix);
        }
    }

    try
    {
        if (ODE15i_event_func)
        {
            outputs = ODE15i_eval_ptr->DoMultiReturnFunctionCall(ODE15i_event_func, inputs, (int)inputs.size(), 3, true);
        }
        else if (ODE15i_event_pntr)
        {
            outputs = ODE15i_eval_ptr->DoMultiReturnFunctionCall(ODE15i_event_pntr, inputs, (int)inputs.size(), 3, true);
        }
        else
        {
            return -1;
        }
    }
    catch (OML_Error&)
    {
        // ODE15iResetStaticVars();
        throw;
    }
    catch (hwMathException&)
    {
        // ODE15iResetStaticVars();
        throw;
    }

    // manage events
    if (outputs.size() != 3)
        return -1;

    if (outputs[0].IsScalar() && outputs[1].IsScalar() && outputs[2].IsScalar())
    {
        gout[0] = outputs[0].Scalar();
        double isTerminal = outputs[1].Scalar();
        double direction = outputs[2].Scalar();

        if (isTerminal != 0.0 && isTerminal != 1.0)
        {
            return -1;
        }

        if (direction != 0.0 && direction != 1.0 && direction != -1.0)
        {
            return -1;
        }

        WriteEventDataIDA(&isTerminal, &direction, 1);
    }
    else if (outputs[0].IsMatrix() && outputs[1].IsMatrix() && outputs[2].IsMatrix())
    {
        hwMatrix* value = outputs[0].GetWritableMatrix();
        hwMatrix* isTerminal = outputs[1].GetWritableMatrix();
        hwMatrix* direction = outputs[2].GetWritableMatrix();

        if (!value->IsVector() || !value->IsReal() ||
            !isTerminal->IsVector() || !isTerminal->IsReal() ||
            !direction->IsVector() || !direction->IsReal())
        {
            return -1;
        }

        if (value->Size() != ODE15i_event_size || isTerminal->Size() != ODE15i_event_size ||
            direction->Size() != ODE15i_event_size)
        {
            return -1;
        }

        for (int i = 0; i < ODE15i_event_size; ++i)
        {
            gout[i] = (*value)(i);

            if ((*isTerminal)(i) != 0.0 && (*isTerminal)(i) != 1.0)
            {
                return -1;
            }

            if ((*direction)(i) != 0.0 && (*direction)(i) != 1.0 && (*direction)(i) != -1.0)
            {
                return -1;
            }
        }

        WriteEventDataIDA(isTerminal->GetRealData(), direction->GetRealData(), ODE15i_event_size);
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
    int nargout = eval.GetNargoutValue();

    if (nargin < 4)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (nargout > 2 && nargout != 5)
    {
        throw OML_Error(OML_ERR_NUMARGOUT);
    }

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
    else if (!eval.FindFunctionByName(funcName, &funcInfo1, &funcPntr1, nullptr))
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
    const hwMatrix* time         = inputs[1].ConvertToMatrix();
    const hwMatrix* y            = inputs[2].ConvertToMatrix();
    double          reltol       = 0.001;
    const hwMatrix* abstol       = nullptr;
    double          maxstep      = -999.0;
    bool            deleteAbsTol = false;
    const hwMatrix* yp           = inputs[3].ConvertToMatrix();
    std::string     JacFuncName;
    FunctionInfo*   funcInfo2    = nullptr;
    FUNCPTR         funcPntr2    = nullptr;
    std::string     eventsFuncName;
    FunctionInfo*   funcInfo3    = nullptr;
    FUNCPTR         funcPntr3    = nullptr;

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
            const Currency& maxstep_C  = options->GetValue(0, -1, "MaxStep");
            const Currency& jacobian_C = options->GetValue(0, -1, "Jacobian");
            const Currency& events_C   = options->GetValue(0, -1, "Events");

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
                throw OML_Error(OML_ERR_SCALARVECTOR, 5, OML_VAR_ABSTOL);
            }

            if (!maxstep_C.IsEmpty())
            {
                if (maxstep_C.IsScalar())
                {
                    maxstep = maxstep_C.Scalar();
                }
                else
                {
                    throw OML_Error(OML_ERR_SCALAR, 5, OML_VAR_MAXSTEP);
                }
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
                        if (!eval.FindFunctionByName(JacFuncName, &funcInfo2, &funcPntr2, nullptr))
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

            if (!events_C.IsEmpty())
            {
                if (events_C.IsFunctionHandle())
                {
                    eventsFuncName = events_C.FunctionHandle()->FunctionName();

                    if (!eval.FindFunctionByName(eventsFuncName, &funcInfo3, &funcPntr3, nullptr))
                        throw OML_Error(OML_ERR_FUNCNAME, 5);

                    if (funcInfo3 && funcInfo3->Parameters().size() < 2)
                        throw OML_Error(OML_ERR_NUMARGIN, 5);
                }
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

    // call event user function to determine how many events exist
    ODE15i_event_func = funcInfo3;
    ODE15i_eval_ptr = &eval;

    if (ODE15i_event_func)
    {
        try
        {
            hwMatrix* y_temp = new hwMatrix(*y);
            hwMatrix* yp_temp = new hwMatrix(*yp);

            if (y_temp->N() != 1)
                y_temp->Transpose();

            if (yp_temp->N() != 1)
                yp_temp->Transpose();

            std::vector<Currency> temp_inputs;
            std::vector<Currency> temp_outputs;
            Currency result;
            temp_inputs.push_back((*time)(0));
            temp_inputs.push_back(y_temp);
            temp_inputs.push_back(yp_temp);

            // find number of constraint returns
            temp_outputs = ODE15i_eval_ptr->DoMultiReturnFunctionCall(ODE15i_event_func, temp_inputs, 1, 3, true);

            if (temp_outputs.size() != 3)
            {
                throw OML_Error(OML_ERR_FUNCNAME, 5);
            }

            Currency& V = temp_outputs[0];
            Currency& T = temp_outputs[1];
            Currency& D = temp_outputs[2];

            // TODO: verify output sizes
            if (V.IsScalar() && T.IsScalar() && D.IsScalar())
            {
                ODE15i_event_size = 1;
            }
            else if (V.IsMatrix() && T.IsMatrix() && D.IsMatrix())
            {
                hwMatrix* v = V.GetWritableMatrix();
                hwMatrix* t = T.GetWritableMatrix();
                hwMatrix* d = D.GetWritableMatrix();

                if (!v->IsVector() || !t->IsVector() || !d->IsVector())
                {
                    throw OML_Error(OML_ERR_FUNCNAME, 5);
                }

                ODE15i_event_size = v->Size();

                if (t->Size() != ODE15i_event_size || d->Size() != ODE15i_event_size)
                {
                    throw OML_Error(OML_ERR_FUNCNAME, 5);
                }
            }
        }
        catch (OML_Error&)
        {
            throw;
        }
    }

    std::unique_ptr<hwMatrix> pEventTime  = nullptr;
    std::unique_ptr<hwMatrix> pEventFnVal = nullptr;
    std::unique_ptr<hwMatrix> pEventIndx  = nullptr;

    if (ODE15i_event_size)
    {
        pEventTime.reset(EvaluatorInterface::allocateMatrix());
        pEventFnVal.reset(EvaluatorInterface::allocateMatrix());
        pEventIndx.reset(EvaluatorInterface::allocateMatrix());
    }

    // set file scope variables
    ODE15i_sys_func = funcInfo1;
    ODE15i_JAC_func = funcInfo2;
    ODE15i_sys_name = funcName;
    ODE15i_sys_pntr = funcPntr1;
    ODE15i_JAC_pntr = funcPntr2;
    ODE15i_sys_size = y->Size();

    // call algorithm
    hwMatrix*        timeSolution = nullptr;
    hwMatrix*        ySolution    = new hwMatrix;
    hwMatrix*        userData     = nullptr;
    IDARootFn_client rootFunc     = nullptr;
    hwMathStatus     status;

    if (ODE15i_event_size)
    {
        rootFunc = ODE15i_event_file_func;
    }

    if (time->Size() == 2)
    {
        // one step mode
        timeSolution = new hwMatrix;

        if (!ODE15i_JAC_func)
        {
            status = DAE15i(ODE15i_file_func, rootFunc, ODE15i_event_size,
                            (IDADenseJacFn_client) nullptr, *time, *y, *yp,
                            timeSolution, *ySolution, reltol, abstol, maxstep, userData,
                            pEventTime.get(), pEventFnVal.get(), pEventIndx.get());
        }
        else
        {
            status = DAE15i(ODE15i_file_func, rootFunc, ODE15i_event_size,
                            ODE15i_jac_file_func, *time, *y, *yp, timeSolution,
                            *ySolution, reltol, abstol, maxstep, userData,
                            pEventTime.get(), pEventFnVal.get(), pEventIndx.get());
        }
    }
    else
    {
        // interval mode
        timeSolution = (hwMatrix*) time;
        timeSolution->IncrRefCount();

        if (!ODE15i_JAC_func)
        {
            status = DAE15i(ODE15i_file_func, rootFunc, ODE15i_event_size,
                            (IDADenseJacFn_client) nullptr, *time, *y, *yp, nullptr,
                            *ySolution, reltol, abstol, maxstep, userData,
                            pEventTime.get(), pEventFnVal.get(), pEventIndx.get());
        }
        else
        {
            status = DAE15i(ODE15i_file_func, rootFunc, ODE15i_event_size,
                            ODE15i_jac_file_func, *time, *y, *yp, nullptr,
                            *ySolution, reltol, abstol, maxstep, userData,
                            pEventTime.get(), pEventFnVal.get(), pEventIndx.get());
        }
    }

    if (deleteAbsTol)
    {
        delete abstol;
        abstol = nullptr;
    }

    if (!status.IsOk())
    {
        if (status.GetArg1() == 1)
        {
            status.SetUserFuncName(funcName);
        }
        else if (status.GetArg1() == 2)
        {
            status.SetUserFuncName(eventsFuncName);
            status.SetArg1(5);
        }
        else if (status.GetArg1() == 4)
        {
            status.SetUserFuncName(JacFuncName);
            status.SetArg1(5);
        }
        else
        {
            switch (status.GetArg1())
            {
                case  5: status.SetArg1(2);  break;
                case  6: status.SetArg1(3);  break;
                case  7: status.SetArg1(4);  break;
                case 10: status.SetArg1(5);  break;
                case 11: status.SetArg1(5);  break;
                case 12: status.SetArg1(5);  break;
                default: status.ResetArgs(); break;
            }
        }

        if (status.GetArg2() == 7)
        {
            status.SetArg2(4);
        }
        else if (status.GetArg2() == 11)
        {
            status.SetArg2(4);
        }

        if (status.IsWarning())
        {
            BuiltInFuncsUtils::SetWarning(eval, status.GetMessageString());
        }
        else
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

            delete ySolution;
            ySolution = nullptr;
            throw OML_Error(status);
            return false;
        }
    }

    // pack outputs
    outputs.push_back(timeSolution);
    outputs.push_back(ySolution);

    if (nargout == 5)
    {
        outputs.push_back(pEventTime.release());
        outputs.push_back(pEventFnVal.release());
        outputs.push_back(pEventIndx.release());
    }

    return true;
}
