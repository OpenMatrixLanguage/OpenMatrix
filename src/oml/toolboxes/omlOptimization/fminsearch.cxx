/**
* @file fminsearch.cxx
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
static std::vector<EvaluatorInterface*> FMINSEARCH_eval_ptr_stack;
static std::vector<FunctionInfo*>       FMINSEARCH_oml_func_stack;
static std::vector<FUNCPTR>             FMINSEARCH_oml_pntr_stack;
static std::vector<std::string>         FMINSEARCH_oml_name_stack;

//------------------------------------------------------------------------------
// Helper function for fminsearch algorithm
//------------------------------------------------------------------------------
static hwMathStatus FMINSEARCH_file_func(const hwMatrix& X, double& y)
{
    std::vector<Currency> inputs;
    inputs.push_back(EvaluatorInterface::allocateMatrix(&X));
    Currency result;

    EvaluatorInterface* FMINSEARCH_eval_ptr = FMINSEARCH_eval_ptr_stack.back();
    FunctionInfo*       FMINSEARCH_oml_func = FMINSEARCH_oml_func_stack.back();
    FUNCPTR             FMINSEARCH_oml_pntr = FMINSEARCH_oml_pntr_stack.back();
    std::string         FMINSEARCH_oml_name = FMINSEARCH_oml_name_stack.back();

    if (FMINSEARCH_oml_func)
    {
        result = FMINSEARCH_eval_ptr->CallInternalFunction(
                 FMINSEARCH_oml_func, inputs);
    }
    else if (FMINSEARCH_oml_pntr)
    {
        result = FMINSEARCH_eval_ptr->CallFunction(FMINSEARCH_oml_name, inputs);
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);
    }

    if (result.IsScalar())
        y = result.Scalar();
    else
        return hwMathStatus(HW_MATH_ERR_USERFUNCREALNUM, 111);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Finds the unconstrained minimum of a real function using the Nelder-Mead 
// simplex algorithm [fminsearch]
//------------------------------------------------------------------------------
bool OmlFminsearch(EvaluatorInterface           eval, 
                   const std::vector<Currency>& inputs, 
                   std::vector<Currency>&       outputs)
{
    int nargin = eval.GetNarginValue();
    int nargout = eval.GetNargoutValue();

    if (nargin < 2 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::string funcName;
    Currency    cur1 = inputs[0];

    if (cur1.IsFunctionHandle())
        funcName = cur1.FunctionHandle()->FunctionName();
    else if (cur1.IsString())
        funcName = cur1.StringVal();
    else
        throw OML_Error(OML_ERR_HANDLE, 1, OML_VAR_TYPE);

    FunctionInfo* funcInfo = nullptr;
    FUNCPTR       funcPntr = nullptr;

    if (funcName == "anonymous")
        funcInfo = cur1.FunctionHandle();
    else if (!eval.FindFunctionByName(funcName, &funcInfo, &funcPntr))
        throw OML_Error(OML_ERR_FUNCNAME, 1);

    if (funcInfo && funcInfo->Parameters().size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN, 1);

    Currency cur2 = inputs[1];

    if (!cur2.IsScalar() && !cur2.IsMatrix())
        throw OML_Error(OML_ERR_SCALARVECTOR, 2);

    hwMatrix* optPoint = EvaluatorInterface::
                         allocateMatrix(cur2.ConvertToMatrix());

    if (!optPoint->IsFinite())
        throw OML_Error(OML_ERR_FINITE, 2);

    if (optPoint->M() == 1)
        optPoint->Transpose();

    bool   displayHist  = false;
    int    maxIter      = 100;
    int    maxFunEval   = 400;
    double tolx         = 1.0e-7;

    if (nargin > 2 && inputs[2].IsStruct())
    {
        StructData* opt = inputs[2].Struct();
        if (opt->N() != 1)
            throw OML_Error(OML_ERR_OPTION, 3);

        Currency maxIterC = opt->GetValue(0, -1, "MaxIter");
        if (maxIterC.IsPositiveInteger())
            maxIter = static_cast<int>(maxIterC.Scalar());
        else if (!maxIterC.IsEmpty())

            throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_MAXITER);               

        Currency maxFunEvalC = opt->GetValue(0, -1, "MaxFunEvals");

        if (maxFunEvalC.IsPositiveInteger())
            maxFunEval = static_cast<int>(maxFunEvalC.Scalar());
        else if (!maxFunEvalC.IsEmpty()) 
            throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_MAXFUNEVALS);

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

    // set file scope variables
    FMINSEARCH_eval_ptr_stack.push_back(&eval);
    FMINSEARCH_oml_func_stack.push_back(funcInfo);
    FMINSEARCH_oml_pntr_stack.push_back(funcPntr);
    FMINSEARCH_oml_name_stack.push_back(funcName);

    // call algorithm
    hwMatrix* objHist    = (displayHist || nargout > 3) ?
                           EvaluatorInterface::allocateMatrix() : nullptr;
    hwMatrix* designHist = (nargout > 3) ?
                           EvaluatorInterface::allocateMatrix() : nullptr;
    double    minVal;

    hwMathStatus status = NelderMead(FMINSEARCH_file_func, *optPoint, minVal, 
                          maxIter, maxFunEval, tolx, objHist, designHist);

    // display history
    if (displayHist && (status.IsOk() || status.IsWarning() || status.IsInfoMsg()))
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

    if (!status.IsOk())
    {
        if (status == HW_MATH_WARN_MAXITERATE ||
            status == HW_MATH_WARN_MAXFUNCEVAL ||
            status.IsInfoMsg())
        {
            // warning message has been replaced by a return value
        }
        else
        {
            FMINSEARCH_eval_ptr_stack.clear();
            FMINSEARCH_oml_func_stack.clear();
            FMINSEARCH_oml_pntr_stack.clear();
            FMINSEARCH_oml_name_stack.clear();

            if (status.GetArg1() == 111)
            {
                status.SetUserFuncName(funcName);
                status.SetArg1(1);
            }

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

            // Set warning or throw error
            BuiltInFuncsUtils::CheckMathStatus(eval, status);
        }
    }

    // pack outputs
    if (cur2.Matrix()->M() == 1)
        optPoint->Transpose();

    outputs.push_back(optPoint);

    if (nargout > 1)
        outputs.push_back(minVal);

    if (nargout > 2)
    {
        if (status.IsOk())
            outputs.push_back(1.0);
        else if (status == HW_MATH_WARN_MAXITERATE ||
                 status == HW_MATH_WARN_MAXFUNCEVAL)
            outputs.push_back(0.0);
    }

    if (nargout > 3)
    {
        Currency out = EvaluatorInterface::allocateStruct();
        out.Struct()->SetValue(0, -1, "iterations", maxIter);
        out.Struct()->SetValue(0, -1, "nfev",       maxFunEval);
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

    FMINSEARCH_eval_ptr_stack.pop_back();
    FMINSEARCH_oml_func_stack.pop_back();
    FMINSEARCH_oml_pntr_stack.pop_back();
    FMINSEARCH_oml_name_stack.pop_back();

    return true;
}
