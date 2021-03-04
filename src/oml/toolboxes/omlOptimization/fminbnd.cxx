/**
* @file fminbnd.cxx
* @date September 2017
* Copyright (C) 2017-2019 Altair Engineering, Inc.  
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

#include "hwOptimizationFuncs.h"

// File scope variables and functions
static std::vector<EvaluatorInterface*> FMINBND_eval_ptr_stack;
static std::vector<FunctionInfo*>       FMINBND_oml_func_stack;
static std::vector<FUNCPTR>             FMINBND_oml_pntr_stack;
static std::vector<std::string>         FMINBND_oml_name_stack;

//------------------------------------------------------------------------------
// Helper function for OmlFminbnd
//------------------------------------------------------------------------------
static hwMathStatus FMINBND_file_func(double x, double& y)
{
    std::vector<Currency> inputs;
    inputs.push_back(x);
    Currency result;

    EvaluatorInterface* FMINBND_eval_ptr = FMINBND_eval_ptr_stack.back();
    FunctionInfo*       FMINBND_oml_func = FMINBND_oml_func_stack.back();
    FUNCPTR             FMINBND_oml_pntr = FMINBND_oml_pntr_stack.back();
    std::string         FMINBND_oml_name = FMINBND_oml_name_stack.back();

    if (FMINBND_oml_func)
        result = FMINBND_eval_ptr->CallInternalFunction(FMINBND_oml_func, inputs);
    else if (FMINBND_oml_pntr)
        result = FMINBND_eval_ptr->CallFunction(FMINBND_oml_name, inputs);
    else
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);

    if (result.IsScalar())
        y = result.Scalar();
    else
        return hwMathStatus(HW_MATH_ERR_USERFUNCREALNUM, 111);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Finds the minimum of a univariate real function within an interval [fminbnd]
//------------------------------------------------------------------------------
bool OmlFminbnd(EvaluatorInterface           eval, 
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
    else if (!eval.FindFunctionByName(funcName, &funcInfo, &funcPntr, nullptr))
        throw OML_Error(OML_ERR_FUNCNAME, 1, OML_VAR_PARAMETER);

    if (funcInfo && funcInfo->Parameters().size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN, 1);

    Currency cur2 = inputs[1];

    if (!cur2.IsMatrix())
        throw OML_Error(OML_ERR_VECTOR, 2, OML_VAR_VARIABLE);

    const hwMatrix* interval = cur2.Matrix();

    if (!interval->IsReal())
        throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);

    if (interval->Size() != 2)
        throw OML_Error(OML_ERR_VECTOR2, 2, OML_VAR_DATA);

    bool    displayHist = false;
    int     maxIter     = 400;
    int     maxFunEval  = 1000000;
    double  tolx        = 1.0e-7;

    if (nargin > 2 && inputs[2].IsStruct())
    {
        StructData* opt = inputs[2].Struct();

        if (opt->N() != 1)
            throw OML_Error(OML_ERR_OPTION, 3);

        Currency curMaxIter = opt->GetValue(0, -1, "MaxIter");

        if (!curMaxIter.IsEmpty())
        {
            if (!curMaxIter.IsPositiveInteger())
                throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_MAXITER);                   

            maxIter = static_cast<int>(curMaxIter.Scalar()); 
        }

        Currency curmaxFunEvals = opt->GetValue(0, -1, "MaxFunEvals");

        if (!curmaxFunEvals.IsEmpty())
        {
            if (!curmaxFunEvals.IsPositiveInteger())
                throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_MAXFUNEVALS);                    

            maxFunEval = static_cast<int>(curmaxFunEvals.Scalar());        
        }

        Currency curtolx = opt->GetValue(0, -1, "TolX");

        if (!curtolx.IsEmpty())
        {
            if (!curtolx.IsScalar())
                throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_TOLX);

            tolx = curtolx.Scalar();
        }

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
    FMINBND_eval_ptr_stack.push_back(&eval);
    FMINBND_oml_func_stack.push_back(funcInfo);
    FMINBND_oml_pntr_stack.push_back(funcPntr);
    FMINBND_oml_name_stack.push_back(funcName);

    // call algorithm
    hwMatrix* objHist    = (displayHist || nargout > 3) ?
                           EvaluatorInterface::allocateMatrix() : nullptr;
    hwMatrix* designHist = (nargout > 3) ?
                           EvaluatorInterface::allocateMatrix() : nullptr;
    double    a          = (*interval)(0);
    double    b          = (*interval)(1);
    double    fa;
    double    fb;
    double    min;
    double    fmin;

    hwMathStatus status = Brent(FMINBND_file_func, a, b, fa, fb, min, fmin, 
                          maxIter, maxFunEval, tolx, objHist, designHist);

    // display history
    if (displayHist && (status.IsOk() || status.IsWarning() || status.IsInfoMsg()))
    {
        std::string line = "Iteration    f(lb)         f(ub)\n";
        eval.PrintResult(line);

        for (int i = 0; i < objHist->M(); ++i)
        {
            char intChar[8];
            char numChar1[20];
            char numChar2[20];

#ifdef OS_WIN  // sprintf_s is a safer option
            sprintf_s(intChar, sizeof(intChar), "%5d",    i+1);
            sprintf_s(numChar1, sizeof(numChar1), "%12.5f", (*objHist)(i, 0));
            sprintf_s(numChar2, sizeof(numChar2), "%12.5f", (*objHist)(i, 1));
#else
            sprintf(intChar, "%5d", i+1);
            sprintf(numChar1, "%12.5f", (*objHist)(i, 0));
            sprintf(numChar2, "%12.5f", (*objHist)(i, 1));
#endif
            line = std::string(intChar) + "  " + std::string(numChar1)
                                        + "  " + std::string(numChar2) + "\n";

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
            FMINBND_eval_ptr_stack.clear();
            FMINBND_oml_func_stack.clear();
            FMINBND_oml_pntr_stack.clear();
            FMINBND_oml_name_stack.clear();

            if (status.GetArg1() == 111)
            {
                status.SetUserFuncName(funcName);
                status.SetArg1(1);
            }
            else if (status.GetArg1() == 3)
            {
                status.SetArg1(2);
            }
            else if (status.GetArg2() == 3)
            {
                status.SetArg2(2);
            }
            else
            {
                status.ResetArgs();
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

            // Sets warning or throws error
            BuiltInFuncsUtils::CheckMathStatus(eval, status);
        }
    }

    // pack outputs
    outputs.push_back(min);

    if (nargout > 1)
        outputs.push_back(fmin);

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

    FMINBND_eval_ptr_stack.pop_back();
    FMINBND_oml_func_stack.pop_back();
    FMINBND_oml_pntr_stack.pop_back();
    FMINBND_oml_name_stack.pop_back();

    return true;
}
