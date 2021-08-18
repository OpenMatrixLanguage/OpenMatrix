/**
* @file ODE113.cxx
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

#include "BuiltInFuncs.h"
#include "BuiltInFuncsUtils.h"
#include "FunctionInfo.h"
#include "OML_Error.h"
#include "StructData.h"

#include "DiffEqFuncs.h"

// file scope variables and functions
static EvaluatorInterface* ODE113_eval_ptr;
static int                 ODE113_sys_size;
static FunctionInfo*       ODE113_sys_func;
static std::string         ODE113_sys_name;
static FUNCPTR             ODE113_sys_pntr;
static int                 ODE113_event_size;
static FunctionInfo*       ODE113_event_func;
static FUNCPTR             ODE113_event_pntr;

//------------------------------------------------------------------------------
// Wrapper for the ODE system, used in ODE algorithms called from oml scripts
//------------------------------------------------------------------------------
static int ODE113_file_func(double t, double* y, double* yp, void* userData)
{
    std::vector<Currency> inputs;
    inputs.push_back(t);

    if (ODE113_sys_func && ODE113_sys_func->Nargin() > 1)
    {
        hwMatrix* yMatrix = new hwMatrix(ODE113_sys_size, 1, hwMatrix::REAL);

        for (int i = 0; i < ODE113_sys_size; ++i)
        {
            (*yMatrix)(i) = y[i];
        }

        inputs.push_back(yMatrix);
    }

    Currency result;

    if (ODE113_sys_func)
    {
        result = ODE113_eval_ptr->CallInternalFunction(ODE113_sys_func, inputs);
    }
    else if (ODE113_sys_pntr)
    {
        result = ODE113_eval_ptr->CallFunction(ODE113_sys_name, inputs);
    }
    else
    {
        return -1;
    }

    if (result.IsMatrix())
    {
        const hwMatrix* ypMatrix = result.Matrix();
        if (!ypMatrix || !ypMatrix->IsReal() || 
            ypMatrix->Size() != ODE113_sys_size)
        {
            return -1;
        }

        for (int i = 0; i < ODE113_sys_size; ++i)
        {
            yp[i] = (*ypMatrix)(i);
        }
    }
    else if (result.IsScalar())
    {
        if (ODE113_sys_size != 1)
        {
            return -1;
        }
        yp[0] = result.Scalar();
    }
    else
    {
        return -1;
    }

    return 0;
}
//------------------------------------------------------------------------------
// Wrapper for ODE event function in oml scripts and is called by ode algorithm
//------------------------------------------------------------------------------
static int ODE113_event_file_func(double  t,
                                  double* y,
                                  double* gout,
                                  void*   userData)
{
    std::vector<Currency> inputs;
    inputs.push_back(t);
    std::vector<Currency> outputs;

    if (ODE113_event_func && ODE113_event_func->Nargin() > 1)
    {
        hwMatrix* yMatrix = new hwMatrix(ODE113_sys_size, 1, hwMatrix::REAL);

        for (int i = 0; i < ODE113_sys_size; ++i)
        {
            (*yMatrix)(i) = y[i];
        }

        inputs.push_back(yMatrix);
    }

    try
    {
        if (ODE113_event_func)
        {
            outputs = ODE113_eval_ptr->DoMultiReturnFunctionCall(ODE113_event_func, inputs, (int)inputs.size(), 3, true);
        }
        else if (ODE113_event_pntr)
        {
            outputs = ODE113_eval_ptr->DoMultiReturnFunctionCall(ODE113_event_pntr, inputs, (int)inputs.size(), 3, true);
        }
        else
        {
            return -1;
        }
    }
    catch (OML_Error&)
    {
        // ODE113ResetStaticVars();
        throw;
    }
    catch (hwMathException&)
    {
        // ODE113ResetStaticVars();
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

        WriteEventDataCV(&isTerminal, &direction, 1);
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

        if (value->Size() != ODE113_event_size || isTerminal->Size() != ODE113_event_size ||
            direction->Size() != ODE113_event_size)
        {
            return -1;
        }

        for (int i = 0; i < ODE113_event_size; ++i)
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

        WriteEventDataCV(isTerminal->GetRealData(), direction->GetRealData(), ODE113_event_size);
    }
    else
    {
        return -1;
    }

    return 0;
}
//------------------------------------------------------------------------------
// Solves a system of non-stiff differential equations - Adams method [ode113]
//------------------------------------------------------------------------------
bool OmlOde113(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    int nargin = eval.GetNarginValue();
    int nargout = eval.GetNargoutValue();

    if (nargin < 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (nargout > 2 && nargout != 5)
    {
        throw OML_Error(OML_ERR_NUMARGOUT);
    }

    Currency cur1 = inputs[0];

    if (!cur1.IsFunctionHandle())
    {
        throw OML_Error(OML_ERR_HANDLE, 1, OML_VAR_TYPE);
    }    

    std::string   funcName = cur1.FunctionHandle()->FunctionName();
    FunctionInfo* funcInfo = nullptr;
    FUNCPTR       funcPntr = nullptr;

    if (funcName == "anonymous")
    {
        funcInfo = inputs[0].FunctionHandle();
    }
    else if (!eval.FindFunctionByName(funcName, &funcInfo, &funcPntr, nullptr))
    {
        throw OML_Error(OML_ERR_FUNCNAME, 1);
    }

    if (funcInfo && funcInfo->Parameters().size() != 1 && 
        funcInfo->Parameters().size() != 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN, 1);
    }

    Currency cur2 = inputs[1];
    if (!cur2.IsScalar() && !cur2.IsMatrix())
    {
        throw OML_Error(OML_ERR_SCALARVECTOR, 2);
    }

    Currency cur3 = inputs[2];
    if (!cur3.IsScalar() && !cur3.IsMatrix())
    {
        throw OML_Error(OML_ERR_SCALARVECTOR, 3);
    }

    const hwMatrix* time         = cur2.ConvertToMatrix();
    const hwMatrix* y            = cur3.ConvertToMatrix();
    double          reltol       = 0.001;
    const hwMatrix* abstol       = nullptr;
    double          maxstep      = -999.0;
    bool            deleteAbsTol = false;
    std::string     eventsFuncName;
    FunctionInfo*   funcInfo2    = nullptr;
    FUNCPTR         funcPntr2    = nullptr;

    if (nargin > 3)
    {
        if (inputs[3].IsStruct())
        {
            StructData* options = inputs[3].Struct();
            if (options->N() != 1)
            {
                throw OML_Error(OML_ERR_OPTION, 4);
            }
            const Currency& reltol_C  = options->GetValue(0, -1, "RelTol");
            const Currency& abstol_C  = options->GetValue(0, -1, "AbsTol");
            const Currency& maxstep_C = options->GetValue(0, -1, "MaxStep");
            const Currency& events_C  = options->GetValue(0, -1, "Events");

            if (!reltol_C.IsEmpty())
            {
                if (reltol_C.IsScalar())
                {
                    reltol = reltol_C.Scalar();
                }
                else
                {
                    throw OML_Error(OML_ERR_SCALARVECTOR, 4, OML_VAR_RELTOL);
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

            if (!maxstep_C.IsEmpty())
            {
                if (maxstep_C.IsScalar())
                {
                    maxstep = maxstep_C.Scalar();
                }
                else
                {
                    throw OML_Error(OML_ERR_SCALAR, 4, OML_VAR_MAXSTEP);
                }
            }

            if (!events_C.IsEmpty())
            {
                if (events_C.IsFunctionHandle())
                {
                    eventsFuncName = events_C.FunctionHandle()->FunctionName();

                    if (!eval.FindFunctionByName(eventsFuncName, &funcInfo2, &funcPntr2, nullptr))
                        throw OML_Error(OML_ERR_FUNCNAME, 4);

                    if (funcInfo2 && funcInfo2->Parameters().size() < 2)
                        throw OML_Error(OML_ERR_NUMARGIN, 4);
                }
            }
        }
        else if (inputs[3].IsMatrix())
        {
            if (inputs[3].Matrix()->M() != 0 || inputs[3].Matrix()->M() != 0)
            {
                throw OML_Error(OML_ERR_HANDLE_EMPTY, 3);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_HANDLE_EMPTY, 3);
        }
    }

    if (nargin > 4)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    // call event user function to determine how many events exist
    ODE113_event_func = funcInfo2;
    ODE113_eval_ptr = &eval;

    if (ODE113_event_func)
    {
        try
        {
            hwMatrix* y_temp = new hwMatrix(*y);

            if (y_temp->N() != 1)
                y_temp->Transpose();

            std::vector<Currency> temp_inputs;
            std::vector<Currency> temp_outputs;
            Currency result;
            temp_inputs.push_back((*time)(0));
            temp_inputs.push_back(y_temp);

            // find number of constraint returns
            temp_outputs = ODE113_eval_ptr->DoMultiReturnFunctionCall(ODE113_event_func, temp_inputs, 1, 3, true);

            if (temp_outputs.size() != 3)
            {
                throw OML_Error(OML_ERR_FUNCNAME, 4);
            }

            Currency& V = temp_outputs[0];
            Currency& T = temp_outputs[1];
            Currency& D = temp_outputs[2];

            // TODO: verify output sizes
            if (V.IsScalar() && T.IsScalar() && D.IsScalar())
            {
                ODE113_event_size = 1;
            }
            else if (V.IsMatrix() && T.IsMatrix() && D.IsMatrix())
            {
                hwMatrix* v = V.GetWritableMatrix();
                hwMatrix* t = T.GetWritableMatrix();
                hwMatrix* d = D.GetWritableMatrix();

                if (!v->IsVector() || !t->IsVector() || !d->IsVector())
                {
                    throw OML_Error(OML_ERR_FUNCNAME, 4);
                }

                ODE113_event_size = v->Size();

                if (t->Size() != ODE113_event_size || d->Size() != ODE113_event_size)
                {
                    throw OML_Error(OML_ERR_FUNCNAME, 4);
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

    if (ODE113_event_size)
    {
        pEventTime.reset(EvaluatorInterface::allocateMatrix());
        pEventFnVal.reset(EvaluatorInterface::allocateMatrix());
        pEventIndx.reset(EvaluatorInterface::allocateMatrix());
    }

    // set file scope variables
    ODE113_sys_func = funcInfo;
    ODE113_sys_name = funcName;
    ODE113_sys_pntr = funcPntr;
    ODE113_sys_size = y->Size();

    // call algorithm
    hwMatrix*       timeSolution = nullptr;
    hwMatrix*       ySolution    = new hwMatrix;
    hwMatrix*       userData     = nullptr;
    CVRootFn_client rootFunc     = nullptr;
    hwMathStatus    status;

    if (ODE113_event_size)
    {
        rootFunc = ODE113_event_file_func;
    }

    if (time->Size() == 2)
    {
        // one step mode
        timeSolution = new hwMatrix;
        status = ODE11(ODE113_file_func, rootFunc, ODE113_event_size, *time, *y,
                       timeSolution, *ySolution, reltol, abstol, maxstep, userData,
                       pEventTime.get(), pEventFnVal.get(), pEventIndx.get());
    }
    else
    {
        // interval mode
        timeSolution = (hwMatrix*) time;
        timeSolution->IncrRefCount();
        status = ODE11(ODE113_file_func, rootFunc, ODE113_event_size, *time, *y,
                       nullptr, *ySolution, reltol, abstol, maxstep, userData,
                       pEventTime.get(), pEventFnVal.get(), pEventIndx.get());
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
            status.SetArg1(4);
        }
        else if (status.GetArg1() == 4)
        {
            // status.SetUserFuncName(JacFuncName);  // not yet available
        }
        else
        {
            switch (status.GetArg1())
            {
                case  5:  status.SetArg1(2);  break;
                case  6:  status.SetArg1(3);  break;
                case  9:  status.SetArg1(4);  break;
                case 10:  status.SetArg1(4);  break;
                case 11:  status.SetArg1(4);  break;
                default: status.ResetArgs(); break;
            }
        }

        if (status.GetArg2() == 10)
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
