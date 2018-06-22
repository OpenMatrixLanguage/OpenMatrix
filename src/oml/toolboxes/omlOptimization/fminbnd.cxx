/**
* @file fminbnd.cxx
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
static FunctionInfo*       FMINBND_oml_sys_func = nullptr;
static std::string         FMINBND_oml_sys_name;
static FUNCPTR             FMINBND_oml_sys_pntr = nullptr;
static EvaluatorInterface* FMINBND_eval_ptr     = nullptr;

//------------------------------------------------------------------------------
// Helper function for OmlFminbnd
//------------------------------------------------------------------------------
static hwMathStatus FMINBND_file_func(double x, double& y)
{
    std::vector<Currency> inputs;
    inputs.push_back(x);

    Currency result;

    if (FMINBND_oml_sys_func)
        result = FMINBND_eval_ptr->CallInternalFunction(FMINBND_oml_sys_func, inputs);
    else if (FMINBND_oml_sys_pntr)
        result = FMINBND_eval_ptr->CallFunction(FMINBND_oml_sys_name, inputs);
    else
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);

    if (result.IsScalar())
        y = result.Scalar();
    else
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);

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
        throw OML_Error(OML_ERR_PLOT_LIMDATA, 2, OML_VAR_DATA);

    int     maxIter    = 100;
    int     maxFunEval = 400;
    double  tolx       = 1.0e-7;

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
    }

    // set file scope variables
    FMINBND_oml_sys_func = funcInfo;
    FMINBND_oml_sys_name = funcName;
    FMINBND_oml_sys_pntr = funcPntr;
    FMINBND_eval_ptr     = &eval;

    // call algorithm
    double a = (*interval)(0);
    double b = (*interval)(1);
    double fa;
    double fb;
    double min;
    double fmin;

    hwMathStatus status = Brent(FMINBND_file_func, a, b, fa, fb, min, fmin, 
                          tolx, maxIter, maxFunEval);

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
        else if (status.GetArg2() == 3)
        {
            status.SetArg2(2);
        }
        else
        {
            status.ResetArgs();
        }

        // Sets warning or throws error
        BuiltInFuncsUtils::CheckMathStatus(eval, status);
    }

    // pack outputs
    size_t nargout = eval.GetNargoutValue();

    if (nargout >= 0)
        outputs.push_back(min);

    if (nargout > 1)
        outputs.push_back(fmin);

    return true;
}
