/**
* @file fsolve.cxx
* @date September 2017
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
#include "OmlOptimizationTboxFuncs.h"

#include "BuiltInFuncsUtils.h"
#include "FunctionInfo.h"
#include "OML_Error.h"
#include "StructData.h"
#include "hwMathException.h"

#include "hwOptimizationFuncs.h"

// File scope variables and functions
static std::vector<EvaluatorInterface*> FSOLVE_eval_ptr_stack;
static std::vector<FunctionInfo*>       FSOLVE_oml_func_stack;
static std::vector<FUNCPTR>             FSOLVE_oml_pntr_stack;
static std::vector<bool>                FSOLVE_oml_anon_stack;
static std::vector<bool>                FSOLVE_oml_Prow_stack; // row vector checks

//------------------------------------------------------------------------------
// System of equation to be solved. Wrapper function for the system function
//------------------------------------------------------------------------------
static hwMathStatus NLSolveFuncs(const hwMatrix& P, 
                                 const hwMatrix& X,
                                 const hwMatrix* userData, 
                                 hwMatrix&       residual)
{
    // X is not used - see lsqcurvefit.cxx
    if (residual.Size() == 0)
        return hwMathStatus();

    EvaluatorInterface* FSOLVE_eval_ptr        = FSOLVE_eval_ptr_stack.back();
    FunctionInfo*       FSOLVE_oml_func        = FSOLVE_oml_func_stack.back();
    FUNCPTR             FSOLVE_oml_pntr        = FSOLVE_oml_pntr_stack.back();
    bool                FSOLVE_oml_func_isanon = FSOLVE_oml_anon_stack.back();
    bool                FSOLVE_oml_func_isProw = FSOLVE_oml_Prow_stack.back();

    std::vector<Currency> inputs;
    std::vector<Currency> outputs;
    Currency result;
    hwMatrix* P_temp = EvaluatorInterface::allocateMatrix(&P);

    if (FSOLVE_oml_func_isProw && P_temp->M() > 1)
        P_temp->Transpose();

    inputs.push_back(P_temp);

    try
    {
        FSOLVE_eval_ptr->Mark();

        if (FSOLVE_oml_func_isanon)
        {
            result = FSOLVE_eval_ptr->CallInternalFunction(FSOLVE_oml_func, inputs);
            outputs.push_back(result);
        }
        else if (FSOLVE_oml_func)
        {
            outputs = FSOLVE_eval_ptr->DoMultiReturnFunctionCall(FSOLVE_oml_func,
                inputs, static_cast<int>(inputs.size()), 1, true);
        }
        else if (FSOLVE_oml_pntr)
        {
            outputs = FSOLVE_eval_ptr->DoMultiReturnFunctionCall(FSOLVE_oml_pntr, 
                inputs, static_cast<int>(inputs.size()), 1, true);
        }
        else
        {
            return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);
        }

        FSOLVE_eval_ptr->Unmark();
    }
    catch (OML_Error&)
    {
        FSOLVE_eval_ptr->Restore();
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);
    }
    catch (hwMathException&)
    {
        FSOLVE_eval_ptr->Restore();
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);
    }

    if (outputs.size() == 1 || outputs.size() == 2)
    {
        Currency& resid = outputs[0];

        if (resid.IsScalar() && residual.Size() == 1)
        {
            residual(0, 0) = resid.Scalar();
        }
        else if (resid.IsMatrix())
        {
            const hwMatrix* resMatrix = resid.Matrix();

            if (!resMatrix->IsReal())
                return hwMathStatus(HW_MATH_ERR_USERFUNCMATRIX, 111);
            if (resMatrix->Size() != residual.Size())
                return hwMathStatus(HW_MATH_ERR_USERFUNCSIZE, 111);
            for (int i = 0; i < resMatrix->Size(); ++i)
                residual(i) = (*resMatrix)(i);
        }
        else
        {
            return hwMathStatus(HW_MATH_ERR_USERFUNCREAL, 111);
        }
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Jacobian matrix function. Wrapper for the objective function Jacobian
//------------------------------------------------------------------------------
static hwMathStatus NLSolveJacobian(const hwMatrix& P, 
                                    const hwMatrix& X,
                                    const hwMatrix* userData, 
                                    hwMatrix&       J)
{
    // X is not used - see NLCurveFit.cxx
    EvaluatorInterface* FSOLVE_eval_ptr        = FSOLVE_eval_ptr_stack.back();
    FunctionInfo*       FSOLVE_oml_func        = FSOLVE_oml_func_stack.back();
    FUNCPTR             FSOLVE_oml_pntr        = FSOLVE_oml_pntr_stack.back();
    bool                FSOLVE_oml_func_isanon = FSOLVE_oml_anon_stack.back();
    bool                FSOLVE_oml_func_isProw = FSOLVE_oml_Prow_stack.back();

    std::vector<Currency> inputs;
    std::vector<Currency> outputs;
    hwMatrix* P_temp = EvaluatorInterface::allocateMatrix(&P);

    if (FSOLVE_oml_func_isProw && P_temp->M() > 1)
        P_temp->Transpose();

    inputs.push_back(P_temp);

    try
    {
        FSOLVE_eval_ptr->Mark();

        if (FSOLVE_oml_func)
        {
            outputs = FSOLVE_eval_ptr->DoMultiReturnFunctionCall(FSOLVE_oml_func, 
                inputs, static_cast<int>(inputs.size()), 2, true);
        }
        else if (FSOLVE_oml_pntr)
        {
            outputs = FSOLVE_eval_ptr->DoMultiReturnFunctionCall(FSOLVE_oml_pntr, 
                inputs, static_cast<int>(inputs.size()), 2, true);
        }
        else
        {
            return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 222);
        }

        FSOLVE_eval_ptr->Unmark();
    }
    catch (OML_Error&)
    {
        FSOLVE_eval_ptr->Restore();
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 222);
    }
    catch (hwMathException&)
    {
        FSOLVE_eval_ptr->Restore();
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 222);
    }

    if (outputs.size() == 2)
    {
        Currency& objJac = outputs[1];

        if (objJac.IsScalar())
        {
            if (P.Size() != 1 || J.Size() != 1)
                return hwMathStatus(HW_MATH_ERR_USERFUNCSIZE, 222);

            J(0) = objJac.Scalar();
        }
        else if (objJac.IsMatrix())
        {
            const hwMatrix* jacMatrix = objJac.Matrix();

            if (!jacMatrix->IsReal())
                return hwMathStatus(HW_MATH_ERR_USERFUNCREALMAT, 222);
            if (jacMatrix->M() != J.M() || jacMatrix->N() != J.N())
                return hwMathStatus(HW_MATH_ERR_USERFUNCSIZE, 222);

            J = (*jacMatrix);
        }
        else
        {
            return hwMathStatus(HW_MATH_ERR_USERFUNCREAL, 222);
        }
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 222);
    }
    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Finds a solution of a system of real, nonlinear equations [fsolve]
//------------------------------------------------------------------------------
bool OmlFsolve(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    int nargin = eval.GetNarginValue();
    int nargout = eval.GetNargoutValue();

    if (nargin < 2 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency cur1 = inputs[0];

    if (!cur1.IsFunctionHandle())
        throw OML_Error(OML_ERR_HANDLE, 1);

    std::string   funcName = cur1.FunctionHandle()->FunctionName();
    FunctionInfo* funcInfo = nullptr;
    FUNCPTR       funcPntr = nullptr;

    if (funcName == "anonymous")
    {
        funcInfo = cur1.FunctionHandle();
        (funcInfo->RedirectedFunction().empty()) ?
            FSOLVE_oml_anon_stack.push_back(true):  // true anonymous function
            FSOLVE_oml_anon_stack.push_back(false); // handle to an inner function
    }
    else
    {
        FSOLVE_oml_anon_stack.push_back(false);

        if (!eval.FindFunctionByName(funcName, &funcInfo, &funcPntr, nullptr))
        {
            throw OML_Error(OML_ERR_FUNCNAME, 1);
        }
    }

    if (funcInfo && funcInfo->Parameters().size() != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN, 1);
    }

    Currency cur2 = inputs[1];

    if (!cur2.IsScalar() && !cur2.IsMatrix())
        throw OML_Error(OML_ERR_SCALARVECTOR, 2);

    hwMatrix* optPoint = EvaluatorInterface::
                         allocateMatrix(cur2.ConvertToMatrix());

    if (!optPoint->IsFinite())
        throw OML_Error(OML_ERR_FINITE, 2);

    bool   displayHist        = false;
    bool   analyticalJacobian = false;
    int    maxIter            = 400;
    int    maxFuncEval        = 1000000;
    double tolf               = 1.0e-7;
    double tolx               = 1.0e-7;
    double minVal;

    if (nargin > 2)
    {
        if (inputs[2].IsStruct())
        {
            StructData* opt = inputs[2].Struct();

            if (opt->N() != 1)
                throw OML_Error(OML_ERR_OPTION, 3);
            
            Currency jacobianC = opt->GetValue(0, -1, "Jacobian");

            if (jacobianC.IsString())
            {
                std::string val (jacobianC.StringVal());

                if (val == "on")
                    analyticalJacobian = true;
                else if (val != "off")
                    throw OML_Error(OML_ERR_FUNCSWITCH, 3, OML_VAR_JACOBIAN);
            }
            else if (!jacobianC.IsEmpty())
            {
                throw OML_Error(OML_ERR_FUNCSWITCH, 3, OML_VAR_JACOBIAN);
            }

            Currency maxIterC = opt->GetValue(0, -1, "MaxIter");

            if (maxIterC.IsPositiveInteger())
                maxIter = static_cast<int>(maxIterC.Scalar());
            else if (!maxIterC.IsEmpty())
                throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_MAXITER);               

            Currency maxFunEvalC = opt->GetValue(0, -1, "MaxFunEvals");

            if (maxFunEvalC.IsPositiveInteger())
                maxFuncEval = static_cast<int>(maxFunEvalC.Scalar());
            else if (!maxFunEvalC.IsEmpty()) 
                throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_MAXFUNEVALS);

            Currency tolfC = opt->GetValue(0, -1, "TolFun");

            if (tolfC.IsScalar())
                tolf = tolfC.Scalar();
            else if (!tolfC.IsEmpty())
                throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_TOLFUN);

            Currency tolxC = opt->GetValue(0, -1, "TolX");

            if (tolxC.IsScalar())
                tolx = tolxC.Scalar();
            else if (!tolxC.IsEmpty())
                throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_TOLX);

            Currency displayC = opt->GetValue(0, -1, "Display");

            if (displayC.IsString())
            {
                std::string val (displayC.StringVal());

                if (val == "iter")
                    displayHist = true;
                else if (val != "off")
                    throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_DISPLAY);
            }
            else if (!displayC.IsEmpty())
            {
                throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_DISPLAY);
            }
        }
        else if (inputs[2].IsMatrix())
        {
            if (inputs[2].Matrix()->M() != 0 || inputs[2].Matrix()->M() != 0)
                throw OML_Error(OML_ERR_HANDLE_EMPTY, 3);
        }
        else
        {
            throw OML_Error(OML_ERR_HANDLE_EMPTY, 3);
        }
    }

    // set file scope variables
    FSOLVE_eval_ptr_stack.push_back(&eval);
    FSOLVE_oml_func_stack.push_back(funcInfo);
    FSOLVE_oml_pntr_stack.push_back(funcPntr);

    if (optPoint->N() > 1)
        FSOLVE_oml_Prow_stack.push_back(true);
    else
        FSOLVE_oml_Prow_stack.push_back(false);

    // call user function to figure out how many equations are in the system
    std::vector<Currency> temp_outputs;

    std::vector<Currency> temp_inputs;
    temp_inputs.push_back(EvaluatorInterface::allocateMatrix(optPoint));

    if (funcName == "anonymous")
    {
        std::string script_func = funcInfo->RedirectedFunction();

        if (!script_func.empty()) // handle to an inner function
        {
            FunctionInfo* script_func_fi;
            FUNCPTR       dummy;

            if (!eval.FindFunctionByName(script_func, &script_func_fi, &dummy, nullptr))
                throw OML_Error(OML_ERR_FUNCNAME, 1);

            FSOLVE_oml_anon_stack.push_back(false); // handle to an inner function
        }
        else
        {
            FSOLVE_oml_anon_stack.push_back(true);  // true anonymous function
        }
    }

    if (FSOLVE_oml_anon_stack.back())
    {
        Currency cur = eval.CallInternalFunction(funcInfo, temp_inputs);
        temp_outputs.push_back(cur);
    }
    else
    {
        temp_outputs = eval.DoMultiReturnFunctionCall(funcInfo,
             temp_inputs, static_cast<int>(temp_inputs.size()), 1, true);
    }

    // call algorithm
    int      numEqns = -1;
    Currency result  = temp_outputs[0];

    if (result.IsScalar())
        numEqns = 1;
    else if (result.IsMatrix())
        numEqns = result.Matrix()->Size();

    hwMatrix*    objHist    = (displayHist || nargout > 3) ?
                              EvaluatorInterface::allocateMatrix() : nullptr;
    hwMatrix*    designHist = (nargout > 3) ?
                              EvaluatorInterface::allocateMatrix() : nullptr;
    hwMatrix*    userData = nullptr;
    hwMathStatus status;

    if (analyticalJacobian)
    {
        status = NLSolve(NLSolveFuncs, NLSolveJacobian, *optPoint, minVal, 
            maxIter, maxFuncEval, numEqns, tolf, tolx, objHist, designHist, 
            userData);
    }
    else
    {
        status = NLSolve(NLSolveFuncs, (LSqFitFunc) NULL, *optPoint, minVal, 
            maxIter, maxFuncEval, numEqns, tolf, tolx, objHist, designHist, 
            userData);
    }

    // display history
    if (displayHist && (status.IsOk() || status.IsWarning() || status.IsInfoMsg()))
    {
        std::string line = "Iteration    f(x)\n";
        eval.PrintResult(line);

        for (int i = 0; i < objHist->Size(); ++i)
        {
            char intChar[8];
            char numChar[20];

#ifdef OS_WIN  // sprintf_s is a safer option
            sprintf_s(intChar, sizeof(intChar), "%5d",    i+1);
            sprintf_s(numChar, sizeof(numChar), "%12.5f", (*objHist)(i));
#else
            sprintf(intChar, "%5d", i+1);
            sprintf(numChar, "%12.5f", (*objHist)(i));
#endif
            line = std::string(intChar) + "  " + std::string(numChar) + "\n";
            eval.PrintResult(line);
        }

        line = "\n";
        eval.PrintResult(line);
    }

    if (!status.IsOk())
    {
        if (status == HW_MATH_WARN_MAXITERATE ||
            status == HW_MATH_WARN_MAXFUNCEVAL ||
            status == HW_MATH_WARN_LOCALMIN ||
            status == HW_MATH_WARN_NOTCONVERGE ||
            status.IsInfoMsg())
        {
            // warning message has been replaced by a return value
        }
        else
        {
            if (status.GetArg1() == 111)
            {
                status.SetUserFuncName(funcName);
                status.SetArg1(1);
            }
            else if (status.GetArg1() == 3)
            {
                status.SetArg1(2);
            }
            else if (status.GetArg1() > 4 && status.GetArg1() < 10)
            {
                status.SetArg1(3);
            }
            else
            {
                status.ResetArgs();
            }

            if (status.GetArg2() == 7)
            {
                status.SetArg2(status.GetArg1());
                status.SetArg1(1);
            }

            if (status.IsWarning())
            {
                BuiltInFuncsUtils::SetWarning(eval, status.GetMessageString());
            }
            else
            {
                FSOLVE_eval_ptr_stack.clear();
                FSOLVE_oml_func_stack.clear();
                FSOLVE_oml_pntr_stack.clear();
                FSOLVE_oml_anon_stack.clear();
                FSOLVE_oml_Prow_stack.clear();

                delete optPoint;
                optPoint = nullptr;

                if (objHist)
                {
                    delete objHist;
                    objHist = nullptr;
                }

                if (designHist)
                {
                    delete designHist;
                    designHist = nullptr;
                }

                throw OML_Error(status);
            }
        }
    }

    // pack outputs
    outputs.push_back(optPoint);

    if (nargout > 1)
    {
        hwMatrix dummy;
        hwMatrix* fval = EvaluatorInterface::allocateMatrix(numEqns, 1,
                         true);
        hwMathStatus status2 = NLSolveFuncs(*optPoint, dummy, userData, *fval);
        outputs.push_back(fval);
    }

    if (nargout > 2)
    {
        if (status == HW_MATH_INFO_TOLFCONV)
            outputs.push_back(1.0);
        else if (status == HW_MATH_INFO_TOLXCONV)
            outputs.push_back(2.0);
        else if (status == HW_MATH_INFO_TOLFCONV_R)
            outputs.push_back(3.0);
        else if (status == HW_MATH_INFO_TOLXCONV_R)
            outputs.push_back(4.0);
        else if (status == HW_MATH_WARN_MAXITERATE ||
                 status == HW_MATH_WARN_MAXFUNCEVAL ||
                 status == HW_MATH_WARN_NOTCONVERGE)
            outputs.push_back(0.0);
        else if (status == HW_MATH_WARN_NOSOLUTION)
            outputs.push_back(-2.0);
        else if (status == HW_MATH_INFO_SMALLTRUST)
            outputs.push_back(-3.0);
    }

    if (nargout > 3)
    {
        objHist->Transpose();
        designHist->Transpose();

        Currency out = EvaluatorInterface::allocateStruct();
        out.Struct()->SetValue(0, -1, "iterations", maxIter);
        out.Struct()->SetValue(0, -1, "nfev",       maxFuncEval);
        out.Struct()->SetValue(0, -1, "xiter",      designHist);
        out.Struct()->SetValue(0, -1, "fvaliter",   objHist);
        outputs.push_back(out);
    }
    else
    {
        if (objHist)
        {
            delete objHist;
            objHist = nullptr;
        }

        if (designHist)
        {
            delete designHist;
            designHist = nullptr;
        }
    }

    FSOLVE_eval_ptr_stack.pop_back();
    FSOLVE_oml_func_stack.pop_back();
    FSOLVE_oml_pntr_stack.pop_back();
    FSOLVE_oml_anon_stack.pop_back();
    FSOLVE_oml_Prow_stack.pop_back();

    return true;
}
