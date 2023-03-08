/**
* @file lsqcurvefit.cxx
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
static std::vector<EvaluatorInterface*> LSQCURVEFIT_eval_ptr_stack;
static std::vector<FunctionInfo*>       LSQCURVEFIT_oml_func_stack;
static std::vector<FUNCPTR>             LSQCURVEFIT_oml_pntr_stack;
static std::vector<bool>                LSQCURVEFIT_oml_anon_stack;
static std::vector<bool>                LSQCURVEFIT_oml_Prow_stack; // row vector checks

//------------------------------------------------------------------------------
// Wrapper for the system function called by lsqcurvefit algorithm.
// Vector of functions to be fitted, typically they mostly have the same form
//---------------------------------------------, -------------------------------
static hwMathStatus NLCurveFitFunc(const hwMatrix& P, 
                                   const hwMatrix& X,
                                   const hwMatrix* userData, 
                                   hwMatrix&       residual)
{
    if (residual.Size() == 0)
        return hwMathStatus();

    EvaluatorInterface* LSQCURVEFIT_eval_ptr        = LSQCURVEFIT_eval_ptr_stack.back();
    FunctionInfo*       LSQCURVEFIT_oml_func        = LSQCURVEFIT_oml_func_stack.back();
    FUNCPTR             LSQCURVEFIT_oml_pntr        = LSQCURVEFIT_oml_pntr_stack.back();
    bool                LSQCURVEFIT_oml_func_isanon = LSQCURVEFIT_oml_anon_stack.back();
    bool                LSQCURVEFIT_oml_func_isProw = LSQCURVEFIT_oml_Prow_stack.back();

    std::vector<Currency> inputs;
    std::vector<Currency> outputs;
    Currency result;
    hwMatrix* P_temp = EvaluatorInterface::allocateMatrix(&P);

    if (LSQCURVEFIT_oml_func_isProw && P_temp->M() > 1)
        P_temp->Transpose();

    inputs.push_back(P_temp);
    inputs.push_back(EvaluatorInterface::allocateMatrix(&X));

    try
    {
        LSQCURVEFIT_eval_ptr->Mark();

        if (LSQCURVEFIT_oml_func_isanon)
        {
            result = LSQCURVEFIT_eval_ptr->CallInternalFunction(
                LSQCURVEFIT_oml_func, inputs);
        
            outputs.push_back(result);
        }
        else if (LSQCURVEFIT_oml_func)
        {
            outputs = LSQCURVEFIT_eval_ptr->DoMultiReturnFunctionCall(
                LSQCURVEFIT_oml_func, inputs, static_cast<int>(inputs.size()), 
                1, true);
        }
        else if (LSQCURVEFIT_oml_pntr)
        {
            outputs = LSQCURVEFIT_eval_ptr->DoMultiReturnFunctionCall(
                LSQCURVEFIT_oml_pntr, inputs, static_cast<int>(inputs.size()), 
                1, true);
        }
        else
        {
            return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);
        }

        LSQCURVEFIT_eval_ptr->Unmark();
    }
    catch (OML_Error&)
    {
        LSQCURVEFIT_eval_ptr->Restore();
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);
    }
    catch (hwMathException&)
    {
        LSQCURVEFIT_eval_ptr->Restore();
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
                return hwMathStatus(HW_MATH_ERR_USERFUNCREAL, 111);

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
// Jacobian matrix function, wrapper for the objective function Jacobian
//------------------------------------------------------------------------------
static hwMathStatus NLCurveFitJacobian(const hwMatrix& P, 
                                       const hwMatrix& X,
                                       const hwMatrix* userData, 
                                       hwMatrix&       J)
{
    EvaluatorInterface* LSQCURVEFIT_eval_ptr        = LSQCURVEFIT_eval_ptr_stack.back();
    FunctionInfo*       LSQCURVEFIT_oml_func        = LSQCURVEFIT_oml_func_stack.back();
    FUNCPTR             LSQCURVEFIT_oml_pntr        = LSQCURVEFIT_oml_pntr_stack.back();
    bool                LSQCURVEFIT_oml_func_isanon = LSQCURVEFIT_oml_anon_stack.back();
    bool                LSQCURVEFIT_oml_func_isProw = LSQCURVEFIT_oml_Prow_stack.back();

    std::vector<Currency> inputs;
    std::vector<Currency> outputs;
    hwMatrix* P_temp = EvaluatorInterface::allocateMatrix(&P);

    if (LSQCURVEFIT_oml_func_isProw && P_temp->M() > 1)
        P_temp->Transpose();

    inputs.push_back(P_temp);
    inputs.push_back(EvaluatorInterface::allocateMatrix(&X));

    try
    {
        LSQCURVEFIT_eval_ptr->Mark();

        if (LSQCURVEFIT_oml_func)
        {
            outputs = LSQCURVEFIT_eval_ptr->DoMultiReturnFunctionCall(
                LSQCURVEFIT_oml_func, inputs, static_cast<int>(inputs.size()), 2, true);
        }
        else if (LSQCURVEFIT_oml_pntr)
        {
            outputs = LSQCURVEFIT_eval_ptr->DoMultiReturnFunctionCall(
                LSQCURVEFIT_oml_pntr, inputs, static_cast<int>(inputs.size()), 2, true);
        }
        else
        {
            return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);
        }

        LSQCURVEFIT_eval_ptr->Unmark();
    }
    catch (OML_Error&)
    {
        LSQCURVEFIT_eval_ptr->Restore();
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);
    }
    catch (hwMathException&)
    {
        LSQCURVEFIT_eval_ptr->Restore();
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);
    }

    if (outputs.size() == 2)
    {
        Currency& objJac = outputs[1];

        if (objJac.IsScalar())
        {
            if (P.Size() != 1 || J.Size() != 1)
            {
                return hwMathStatus(HW_MATH_ERR_USERFUNCSIZE, 111);
            }

            J(0) = objJac.Scalar();
        }
        else if (objJac.IsMatrix())
        {
            const hwMatrix* jacMatrix = objJac.Matrix();

            if (!jacMatrix->IsReal())
            {
                return hwMathStatus(HW_MATH_ERR_USERFUNCREALMAT, 222);
            }
            if (jacMatrix->M() != J.M() || jacMatrix->N() != J.N())
            {
                return hwMathStatus(HW_MATH_ERR_USERFUNCSIZE, 222);
            }

            J = (*jacMatrix);
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
// Finds the equation parameters that produce the least squares best fit to a 
// data set [lsqcurvefit]
//------------------------------------------------------------------------------
bool OmlLsqcurvefit(EvaluatorInterface           eval, 
                    const std::vector<Currency>& inputs, 
                    std::vector<Currency>&       outputs)
{
    int nargin = eval.GetNarginValue();
    int nargout = eval.GetNargoutValue();

    if (nargin < 4 || nargin == 5 || nargin > 7)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency cur1 = inputs[0];

    if (!cur1.IsFunctionHandle())
        throw OML_Error(OML_ERR_HANDLE, 1);

    FunctionInfo* funcInfo = nullptr;
    FUNCPTR       funcPntr = nullptr;
    std::string   funcName = cur1.FunctionHandle()->FunctionName();

    if (funcName == "anonymous")
    {
        funcInfo = cur1.FunctionHandle();
        (funcInfo->RedirectedFunction().empty()) ?
            LSQCURVEFIT_oml_anon_stack.push_back(true):  // true anonymous function
            LSQCURVEFIT_oml_anon_stack.push_back(false); // handle to an inner function
    }
    else
    {
        LSQCURVEFIT_oml_anon_stack.push_back(false);

        if (!eval.FindFunctionByName(funcName, &funcInfo, &funcPntr, nullptr))
        {
            throw OML_Error(OML_ERR_FUNCNAME, 1);
        }
    }

    if (funcInfo && funcInfo->Parameters().size() != 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN, 1);
    }

    Currency cur2 = inputs[1];

    if (!cur2.IsScalar() && !cur2.IsMatrix())
        throw OML_Error(OML_ERR_SCALARVECTOR, 2);

    hwMatrix* optParam = EvaluatorInterface::
                         allocateMatrix(cur2.ConvertToMatrix());

    if (!optParam->IsFinite())
        throw OML_Error(OML_ERR_FINITE, 2);

    const hwMatrix* X    = nullptr;
    Currency        cur3 = inputs[2];

    if (cur3.IsScalar() || cur3.IsMatrix())
        X = cur3.ConvertToMatrix();
    else
        throw OML_Error(OML_ERR_SCALARVECTOR, 3);

    const hwMatrix* y    = nullptr;
    Currency        cur4 = inputs[3];

    if (cur4.IsScalar() || cur4.IsMatrix())
        y = cur4.ConvertToMatrix();
    else
        throw OML_Error(OML_ERR_SCALARVECTOR, 4);
       
    const hwMatrix* lowerBound = nullptr;
    const hwMatrix* upperBound = nullptr;

    if (nargin > 5)
    {
        if (inputs[4].IsScalar() || inputs[4].IsMatrix())
        {
            lowerBound = inputs[4].ConvertToMatrix();

            if (lowerBound->M() == 0 && lowerBound->N() == 0)   // default: unbounded
                lowerBound = nullptr;
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARVECTOR, 5);
        }

        if (inputs[5].IsScalar() || inputs[5].IsMatrix())
        {
            upperBound = inputs[5].ConvertToMatrix();

            if (upperBound->M() == 0 && upperBound->N() == 0)   // default: unbounded
                upperBound = nullptr;
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARVECTOR, 6);
        }

        if (lowerBound && !upperBound)
            throw OML_Error(OML_ERR_ARRAYSIZE, 5, 6);
        else if (!lowerBound && upperBound)
            throw OML_Error(OML_ERR_ARRAYSIZE, 5, 6);
    }

    bool   displayHist        = false;
    bool   analyticalJacobian = false;
    int    maxIter            = 400;
    int    maxFuncEval        = 1000000;
    double tolf               = 1.0e-7;
    double tolx               = 1.0e-7;

    if (nargin > 6)
    {
        if (!inputs[6].IsStruct())
            throw OML_Error(OML_ERR_STRUCT, 7);

        StructData* opt = inputs[6].Struct();

        if (opt->N() != 1)
            throw OML_Error(OML_ERR_OPTION, 7);

        Currency jacobianC = opt->GetValue(0, -1, "Jacobian");

        if (!jacobianC.IsEmpty())
        {
            if (!jacobianC.IsString())
                throw OML_Error(OML_ERR_FUNCSWITCH, 7, OML_VAR_JACOBIAN);

            std::string val (jacobianC.StringVal());

            if (val == "on")
                analyticalJacobian = true;
            else if (val != "off")
                throw OML_Error(OML_ERR_FUNCSWITCH, 7, OML_VAR_JACOBIAN);
       }

        Currency maxIterC = opt->GetValue(0, -1, "MaxIter");
        if (!maxIterC.IsEmpty())
        {
            if (!maxIterC.IsPositiveInteger())
                throw OML_Error(OML_ERR_OPTIONVAL, 7, OML_VAR_MAXITER);

            maxIter = static_cast<int>(maxIterC.Scalar());
        }

        Currency maxFunEvalC = opt->GetValue(0, -1, "MaxFunEvals");

        if (!maxFunEvalC.IsEmpty())
        {
            if (!maxFunEvalC.IsPositiveInteger())
                throw OML_Error(OML_ERR_OPTIONVAL, 7, OML_VAR_MAXFUNEVALS);

            maxFuncEval = static_cast<int>(maxFunEvalC.Scalar());
        }

        Currency tolfC = opt->GetValue(0, -1, "TolFun");

        if (!tolfC.IsEmpty())
        {
            if (!tolfC.IsScalar())
                throw OML_Error(OML_ERR_OPTIONVAL, 7, OML_VAR_TOLFUN);

            tolf = tolfC.Scalar();
        }

        Currency tolxC = opt->GetValue(0, -1, "TolX");

        if (!tolxC.IsEmpty())
        {
            if (!tolxC.IsScalar())
                throw OML_Error(OML_ERR_OPTIONVAL, 7, OML_VAR_TOLX);

            tolx = tolxC.Scalar();
        }

        Currency displayC = opt->GetValue(0, -1, "Display");

        if (!displayC.IsEmpty())
        {
            if (!displayC.IsString())
                throw OML_Error(OML_ERR_OPTIONVAL, 7, OML_VAR_DISPLAY);

            std::string val (displayC.StringVal());

            if (val == "iter")
                displayHist = true;
            else if (val != "off")
                throw OML_Error(OML_ERR_OPTIONVAL, 7, OML_VAR_DISPLAY);
        }
    }

    // set file scope variables
    LSQCURVEFIT_eval_ptr_stack.push_back(&eval);
    LSQCURVEFIT_oml_func_stack.push_back(funcInfo);
    LSQCURVEFIT_oml_pntr_stack.push_back(funcPntr);
    
    if (optParam->N() > 1)
        LSQCURVEFIT_oml_Prow_stack.push_back(true);
    else
        LSQCURVEFIT_oml_Prow_stack.push_back(false);

    // call algorithm
    hwMatrix*    stats      = nullptr;
    hwMatrix*    yEst       = nullptr;
    hwMatrix*    objHist    = (displayHist || nargout > 1) ? 
                              EvaluatorInterface::allocateMatrix() : nullptr;
    hwMatrix*    designHist = (nargout > 4) ?
                              EvaluatorInterface::allocateMatrix() : nullptr;
    hwMatrix*    userData   = nullptr;
    hwMathStatus status;

    if (analyticalJacobian)
    {
        status = NLCurveFit(NLCurveFitFunc, NLCurveFitJacobian, *optParam, *X, 
            *y, lowerBound, upperBound, maxIter, maxFuncEval, stats, yEst, tolf,
            tolx, objHist, designHist, userData);
    }
    else
    {
        status = NLCurveFit(NLCurveFitFunc, (LSqFitFunc) NULL, *optParam, *X, 
            *y, lowerBound, upperBound, maxIter, maxFuncEval, stats, yEst, tolf,
            tolx, objHist, designHist, userData);
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

    // clean up
    if (!status.IsOk())
    {
        if (status == HW_MATH_WARN_MAXITERATE ||
            status == HW_MATH_WARN_MAXFUNCEVAL ||
            status == HW_MATH_WARN_NOTCONVERGE ||
            status.IsInfoMsg())
        {
            // warning message has been replaced by a return value
        }
        else
        {
            if (status.GetArg1() == 111)
            {
                if (funcName == "anonymous")
                {
                    funcName = funcInfo->RedirectedFunction();

                    if (!funcName.size())
                        funcName = "anonymous";
                }

                status.SetUserFuncName(funcName);
                status.SetArg1(1);
            }
            else if (status.GetArg1() == 3)
            {
                status.SetArg1(2);
            }
            else if (status.GetArg1() == 4)
            {
                status.SetArg1(3);
            }
            else if (status.GetArg1() == 5)
            {
                status.SetArg1(4);
            }
            else if (status.GetArg1() == 6)
            {
                status.SetArg1(5);
            }
            else if (status.GetArg1() == 7)
            {
                status.SetArg1(6);
            }
            else if (status.GetArg1() > 7 && status.GetArg1() < 12)
            {
                status.SetArg1(7);
            }
            else
            {
                status.ResetArgs();
            }

            if (status.GetArg2() == 4)
            {
                status.SetArg2(3);
            }
            else if (status.GetArg2() == 5)
            {
                status.SetArg2(4);
            }
            else if (status.GetArg2() == 6)
            {
                status.SetArg2(5);
            }
            else if (status.GetArg2() == 7)
            {
                status.SetArg2(6);
            }

            if (status.IsWarning())
            {
                BuiltInFuncsUtils::SetWarning(eval, status.GetMessageString());
            }
            else
            {
                LSQCURVEFIT_eval_ptr_stack.clear();
                LSQCURVEFIT_oml_func_stack.clear();
                LSQCURVEFIT_oml_pntr_stack.clear();
                LSQCURVEFIT_oml_anon_stack.clear();
                LSQCURVEFIT_oml_Prow_stack.clear();

                delete optParam;
                optParam = nullptr;

                if (objHist)
                {
                    delete objHist;
                    objHist = nullptr;
                }

                delete optParam;
                optParam = nullptr;

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
    outputs.push_back(optParam);

    if (nargout > 1)
    {
        int size = objHist->Size();

        if (size)
            outputs.push_back((*objHist)(size-1));
    }

    if (nargout > 2)
    {
        hwMatrix dummy;
        hwMatrix* res = EvaluatorInterface::allocateMatrix(y->Size(), 1,
                         true);
        hwMathStatus status2 = NLCurveFitFunc(*optParam, *X, userData, *res);

        int size = res->Size();

        for (int i = 0; i < size; ++i)
            (*res)(i) = (*y)(i) - (*res)(i);

        outputs.push_back(res);
    }

    if (nargout > 3)
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
        else if (status == HW_MATH_INFO_SMALLTRUST)
            outputs.push_back(-3.0);
    }

    if (nargout > 4)
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

    LSQCURVEFIT_eval_ptr_stack.pop_back();
    LSQCURVEFIT_oml_func_stack.pop_back();
    LSQCURVEFIT_oml_pntr_stack.pop_back();
    LSQCURVEFIT_oml_anon_stack.pop_back();
    LSQCURVEFIT_oml_Prow_stack.pop_back();

    return true;
}
