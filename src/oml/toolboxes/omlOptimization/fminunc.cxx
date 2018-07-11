/**
* @file fminunc.cxx
* @date September 2017
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

#include "OmlOptimizationTboxFuncs.h"

#include "BuiltInFuncsUtils.h"
#include "FunctionInfo.h"
#include "OML_Error.h"
#include "StructData.h"

#include "hwOptimizationFuncs.h"

// File scope variables and functions
static EvaluatorInterface* FMINUNC_eval_ptr        = nullptr;
static FunctionInfo*       FMINUNC_oml_func        = nullptr;
static FUNCPTR             FMINUNC_oml_pntr        = nullptr;
static bool                FMINUNC_oml_func_isanon(nullptr);

//------------------------------------------------------------------------------
// Wrapper for the objective function called in OmlFminunc by oml scripts
//------------------------------------------------------------------------------
static hwMathStatus FMinUnconFunc(const hwMatrix& P, 
                                  const hwMatrix* userData, 
                                  double&         min)
{
    std::vector<Currency> outputs;
    std::vector<Currency> inputs;

    inputs.push_back(EvaluatorInterface::allocateMatrix(&P));
    int numinputs = static_cast<int>(inputs.size()); 

    if (FMINUNC_oml_func_isanon)
    {
        Currency result = FMINUNC_eval_ptr->CallInternalFunction(
                          FMINUNC_oml_func, inputs);
        outputs.push_back(result);
    }
    else if (FMINUNC_oml_func)
    {
        outputs = FMINUNC_eval_ptr->DoMultiReturnFunctionCall(
            FMINUNC_oml_func, inputs, numinputs, 1, true);
    }
    else if (FMINUNC_oml_pntr)
    {
        outputs = FMINUNC_eval_ptr->DoMultiReturnFunctionCall(
            FMINUNC_oml_pntr, inputs, numinputs, 1, true);
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);
    }

    if (outputs.size() == 1 || outputs.size() == 2)
    {
        Currency obj = outputs[0];

        if (obj.IsScalar())
        {
            min = obj.Scalar();
        }
        else if (obj.IsMatrix())
        {
            const hwMatrix* minMatrix = obj.Matrix();

            if (minMatrix->Size() != 1)
            {
                return hwMathStatus(HW_MATH_ERR_USERFUNCSIZE, 111);
            }

            if (minMatrix->IsReal())
            {
                min = (*minMatrix)(0);
            }
            else if (minMatrix->z(0).IsReal())
            {
                min = minMatrix->z(0).Real();
            }
            else
            {
                return hwMathStatus(HW_MATH_ERR_USERFUNCREAL, 111);
            }
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
// Gradient vector function wrapper
//------------------------------------------------------------------------------
static hwMathStatus FMinUnconGradient(const hwMatrix& P, 
                                      const hwMatrix* userData, 
                                      hwMatrix&       grad)
{
    std::vector<Currency> outputs;
    std::vector<Currency> inputs;

    inputs.push_back(EvaluatorInterface::allocateMatrix(&P));
    int numinputs = static_cast<int>(inputs.size());    

    if (FMINUNC_oml_func)
    {
        outputs = FMINUNC_eval_ptr->DoMultiReturnFunctionCall(FMINUNC_oml_func, 
                  inputs, numinputs, 2, true);
    }
    else if (FMINUNC_oml_pntr)
    {
        outputs = FMINUNC_eval_ptr->DoMultiReturnFunctionCall(FMINUNC_oml_pntr, 
                  inputs, numinputs, 2, true);
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 222);
    }

    if (outputs.size() != 2)
    {
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 222);
    }

    Currency objGrad = outputs[1];

    if (objGrad.IsScalar())
    {
        if (P.Size() != 1 || grad.Size() != 1)
        {
            return hwMathStatus(HW_MATH_ERR_USERFUNCSIZE, 222);
        }
        grad(0) = objGrad.Scalar();
    }
    else if (objGrad.IsMatrix())
    {
        const hwMatrix* gradMatrix = objGrad.Matrix();

        if (!gradMatrix->IsReal())
        {
            return hwMathStatus(HW_MATH_ERR_USERFUNCREAL, 222);
        }
        if (gradMatrix->M() != grad.M() || gradMatrix->N() != grad.N())
        {
            return hwMathStatus(HW_MATH_ERR_USERFUNCSIZE, 222);
        }
        for (int i = 0; i < grad.M(); ++i)
        {
            grad(i) = (*gradMatrix)(i);
        }
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 222);
    }
      
    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Finds the unconstrained minimum of a real function [fminunc]
//------------------------------------------------------------------------------
bool OmlFminunc(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    int nargin = eval.GetNarginValue();

    if (nargin < 2 || nargin > 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
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
        funcInfo = cur1.FunctionHandle();
        std::string script_func = funcInfo->RedirectedFunction();
        FMINUNC_oml_func_isanon = 
            (script_func.size()) ? false: // Inner function
                                   true;  // true anonymous function
    }
    else if (!eval.FindFunctionByName(funcName, &funcInfo, &funcPntr))
    {
        throw OML_Error(OML_ERR_FUNCNAME, 1);
    }

    if (funcInfo && funcInfo->Parameters().size() != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN, 1);
    }

    Currency  cur2 = inputs[1];

    if (!cur2.IsScalar() && !cur2.IsMatrix())
        throw OML_Error(OML_ERR_SCALARVECTOR, 2);

    hwMatrix* optPoint = EvaluatorInterface::
                         allocateMatrix(cur2.ConvertToMatrix());

    if (!optPoint->IsFinite())
        throw OML_Error(OML_ERR_FINITE, 2);

    bool   displayHist        = false;
    bool   analyticalGradient = false;
    int    maxIter            = 100;
    int    maxFuncEval        = 400;
    double tolf               = 1.0e-7;
    double tolx               = 1.0e-7;

    if (nargin > 2)
    {
        if (inputs[2].IsStruct())
        {
            StructData* opt = inputs[2].Struct();

            if (opt->N() != 1)
            {
                throw OML_Error(OML_ERR_OPTION, 3);
            }           

            Currency gradientC = opt->GetValue(0, -1, "GradObj");

            if (!gradientC.IsEmpty())
            {
                if (!gradientC.IsString())
                {
                    throw OML_Error(OML_ERR_FUNCSWITCH, 3, OML_VAR_JACOBIAN);
                }

                std::string val (gradientC.StringVal());

                if (val != "on" && val != "off")
                {
                    throw OML_Error(OML_ERR_FUNCSWITCH, 3, OML_VAR_JACOBIAN);
                }
                analyticalGradient = (val == "on") ? true : false;
            }

            Currency maxIterC = opt->GetValue(0, -1, "MaxIter");

            if (!maxIterC.IsEmpty())
            {
                if (!maxIterC.IsPositiveInteger())
                {
                    throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_MAXITER);
                }
                maxIter = static_cast<int>(maxIterC.Scalar());                   
            }

            Currency maxFunEvalC = opt->GetValue(0, -1, "MaxFunEvals");

            if (!maxFunEvalC.IsEmpty())
            {
                if (!maxFunEvalC.IsPositiveInteger())
                {
                    throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_MAXFUNEVALS);
                }
                maxFuncEval = static_cast<int>(maxFunEvalC.Scalar());
            }

            Currency tolfC = opt->GetValue(0, -1, "TolFun");

            if (!tolfC.IsEmpty())
            {
                if (!tolfC.IsScalar())
                {
                    throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_TOLFUN);
                }
                tolf = tolfC.Scalar();
            }

            Currency tolxC = opt->GetValue(0, -1, "TolX");

            if (!tolxC.IsEmpty())
            {
                if (!tolxC.IsScalar())
                {
                    throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_TOLX);
                }
                tolx = tolxC.Scalar();
            }

            Currency displayC = opt->GetValue(0, -1, "Display");

            if (!displayC.IsEmpty())
            {
                if (!displayC.IsString())
                {
                    throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_DISPLAY);
                }

                std::string val (displayC.StringVal());

                if (val == "iter")
                {
                    displayHist = true;
                }
                else if (val == "off")
                {
                    displayHist = false;
                }
                else
                {
                    throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_DISPLAY);
                }
            }
        }
        else if (inputs[2].IsMatrix())
        {
            if (inputs[2].Matrix()->M() != 0 || inputs[2].Matrix()->M() != 0)
            {
                throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_DISPLAY);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_HANDLE_EMPTY, 3);
        }
    }

    // set file scope variables
    FMINUNC_eval_ptr = &eval;
    FMINUNC_oml_func = funcInfo;
    FMINUNC_oml_pntr = funcPntr;

    // call algorithm
    hwMatrix*    objHist    = (displayHist) ? EvaluatorInterface::allocateMatrix() :
                              nullptr;
    hwMatrix*    designHist = nullptr;
    hwMatrix*    userData   = nullptr;
    hwMathStatus status;
    double       minVal;

    if (analyticalGradient)
    {
        status = FMinUncon(FMinUnconFunc, FMinUnconGradient, *optPoint, minVal,
            maxIter, maxFuncEval, tolf, tolx, objHist, designHist, userData);
    }
    else
    {
        status = FMinUncon(FMinUnconFunc, (UnConMinGradFunc) NULL, *optPoint, minVal,
            maxIter, maxFuncEval, tolf, tolx, objHist, designHist, userData);
    }

    // display history
    if (displayHist && status.IsOk())
    {

        std::string line = "Iteration      f(x)\n";
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

    if (objHist)
    {
        delete objHist;
        objHist = nullptr;
    }

    if (!status.IsOk())
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
        else if (status.GetArg1() > 4 && status.GetArg1() < 7)
        {
            status.SetArg1(3);
        }
        else
        {
            status.ResetArgs();
        }

        if (status.IsWarning())
        {
            BuiltInFuncsUtils::SetWarning(eval, status.GetMessageString());
        }
        else
        {
            if (optPoint)
            {
                delete optPoint;
                optPoint = nullptr;
            }
            throw OML_Error(status);
        }
    }

    // pack outputs
    int nargout = eval.GetNargoutValue();

    if (nargout > 0)
        outputs.push_back(optPoint);

    if (nargout > 1)
        outputs.push_back(minVal);

    return true;
}
