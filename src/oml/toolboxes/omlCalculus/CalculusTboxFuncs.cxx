/**
* @file CalculusTboxFuncs.cxx
* @date January, 2015
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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

#include "CalculusTboxFuncs.h"

#include <memory>  // For std::unique_ptr

#include "CalculusFuncs.h"

#include "BuiltInFuncsUtils.h"
#include "FunctionInfo.h"
#include "OML_Error.h"

#define CALC "Calculus"
#define TBOXVERSION 2019.0

// file scope variables and functions
static std::vector<FunctionInfo*>       quad_oml_sys_func_stack;
static std::vector<std::string>         quad_oml_sys_name_stack;
static std::vector<FUNCPTR>             quad_oml_sys_pntr_stack;
static std::vector<EvaluatorInterface*> quad_eval_ptr_stack;

//------------------------------------------------------------------------------
// Entry point which registers oml Calculus functions with oml
//------------------------------------------------------------------------------
int InitDll(EvaluatorInterface eval)
{
    eval.RegisterBuiltInFunction("quad", &OmlQuad, 
                                 FunctionMetaData(-4, 2, CALC));
    eval.RegisterBuiltInFunction("quadv", &OmlQuadv, 
                                 FunctionMetaData(-4, 2, CALC));
    eval.RegisterBuiltInFunction("trapz", &OmlTrapz, 
                                 FunctionMetaData(2, 1, CALC));
    eval.RegisterBuiltInFunction("cumtrapz", &OmlCumtrapz, 
                                 FunctionMetaData(2, 1, CALC));

    return 1;
}
//------------------------------------------------------------------------------
// Helper method for performing integrations
//------------------------------------------------------------------------------
static hwMathStatus quad_file_func(const hwMatrix& x, hwMatrix& y)
{
    std::vector<Currency> inputs;
    inputs.push_back(EvaluatorInterface::allocateMatrix(&x));
    Currency result;

    FunctionInfo*       quad_oml_sys_func  = quad_oml_sys_func_stack.back();
    std::string         quad_oml_sys_name  = quad_oml_sys_name_stack.back();
    FUNCPTR             quad_oml_sys_pntr  = quad_oml_sys_pntr_stack.back();
    EvaluatorInterface* quad_eval_ptr      = quad_eval_ptr_stack.back();

    if (quad_oml_sys_func)
        result = quad_eval_ptr->CallInternalFunction(quad_oml_sys_func, inputs);
    else if (quad_oml_sys_pntr)
        result = quad_eval_ptr->CallFunction(quad_oml_sys_name, inputs);
    else
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);

    if (result.IsMatrix())
    {
        y = *(result.GetWritableMatrix());

        if (x.M() != y.M() || x.N() != y.N())
            return hwMathStatus(HW_MATH_ERR_USERFUNCSIZE, 111);
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Helper method for performing integrations
//------------------------------------------------------------------------------
static hwMathStatus quadv_file_func(double x, double& y)
{
    std::vector<Currency> inputs;
    inputs.push_back(x);
    Currency result;

    FunctionInfo*       quad_oml_sys_func  = quad_oml_sys_func_stack.back();
    std::string         quad_oml_sys_name  = quad_oml_sys_name_stack.back();
    FUNCPTR             quad_oml_sys_pntr  = quad_oml_sys_pntr_stack.back();
    EvaluatorInterface* quad_eval_ptr      = quad_eval_ptr_stack.back();

    if (quad_oml_sys_func)
        result = quad_eval_ptr->CallInternalFunction(quad_oml_sys_func, inputs);
    else if (quad_oml_sys_pntr)
        result = quad_eval_ptr->CallFunction(quad_oml_sys_name, inputs);
    else
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);

    if (result.IsScalar())
        y = result.Scalar();
    else
        return hwMathStatus(HW_MATH_ERR_USERFUNCFAIL, 111);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Performs adaptive numerical integration (quadrature) and return true
//------------------------------------------------------------------------------
bool OmlQuad(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    // unpack inputs
    int nargin = eval.GetNarginValue();

    if (nargin < 3 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::string funcName;

    if (inputs[0].IsFunctionHandle())
        funcName = inputs[0].FunctionHandle()->FunctionName();
    else if (inputs[0].IsString())
        funcName = inputs[0].StringVal();
    else
        throw OML_Error(OML_ERR_HANDLE, 1, OML_VAR_TYPE);

    FunctionInfo* funcInfo = NULL;
    FUNCPTR       funcPntr = NULL;

    if (funcName == "anonymous")
    {
        funcInfo = inputs[0].FunctionHandle();
    }
    else
    {
        if (!eval.FindFunctionByName(funcName, &funcInfo, &funcPntr))
            throw OML_Error(OML_ERR_FUNCNAME, 1);
    }

    if (funcInfo && funcInfo->Parameters().size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN, 1);

    if (!inputs[1].IsScalar() && !inputs[1].IsMatrix())
        throw OML_Error(OML_ERR_SCALARMATRIX, 2);

    if (!inputs[2].IsScalar() && !inputs[2].IsMatrix())
        throw OML_Error(OML_ERR_SCALARMATRIX, 3);

    double reltol = sqrt(MACHEP);
    double abstol = sqrt(MACHEP);

    if (nargin > 3)
    {
        if (inputs[3].IsScalar())
        {
            reltol = inputs[3].Scalar();
        }
        else if (inputs[3].IsMatrix())
        {
            const hwMatrix* tol = inputs[3].Matrix();

            if (tol->Size() != 2)
                throw OML_Error(OML_ERR_OPTIONVAL, 4, OML_VAR_RELTOL);

            reltol = (*tol)(0);
            abstol = (*tol)(1);
        }
        else
        {
            throw OML_Error(OML_ERR_NUMERIC, 4, OML_VAR_RELTOL);
        }
    }

    // set file scope variables
    quad_oml_sys_func_stack.push_back(funcInfo);
    quad_oml_sys_name_stack.push_back(funcName);
    quad_oml_sys_pntr_stack.push_back(funcPntr);
    quad_eval_ptr_stack.push_back(&eval);

    // call algorithm
    hwMathStatus status;

    if (inputs[1].IsScalar() && inputs[2].IsScalar())
    {
        double a = inputs[1].Scalar();
        double b = inputs[2].Scalar();
        double area;
        int count = 0;

        status = Quad(quad_file_func, a, b, area, count, reltol, abstol);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 111)
            {
                status.SetUserFuncName(funcName);
                status.SetArg1(1);
            }
            else if (status.GetArg1() == 6)
            {
                status.SetArg1(4);
            }
            else if (status.GetArg1() == 7)
            {
                status.SetArg1(4);
            }
            
            quad_oml_sys_func_stack.clear();
            quad_oml_sys_name_stack.clear();
            quad_oml_sys_pntr_stack.clear();
            quad_eval_ptr_stack.clear();

            BuiltInFuncsUtils::CheckMathStatus(eval, status);
        }

        // pack outputs
        size_t nargout = eval.GetNargoutValue();

        if (nargout >= 0)
            outputs.push_back(area);

        if (nargout > 1)
            outputs.push_back(count);
    }
    else if (inputs[1].IsScalar() && inputs[2].IsMatrix())
    {
        int i   = 0;
        int cnt = 0; 

        double a = inputs[1].Scalar();
        const hwMatrix* B = inputs[2].Matrix();

        if (!B->IsReal())
            throw OML_Error(OML_ERR_REAL, 3, OML_VAR_MATRIX);

        std::unique_ptr<hwMatrix> area(EvaluatorInterface::allocateMatrix(
                                       B->M(), B->N(), hwMatrix::REAL));
        std::unique_ptr<hwMatrix> count(EvaluatorInterface::allocateMatrix(
                                       B->M(), B->N(), hwMatrix::REAL));
        for (i = 0; i < B->Size(); ++i)
        {
            cnt = 0;
            status = Quad(quad_file_func, a, (*B)(i), (*area)(i), cnt, reltol, abstol);
            (*count)(i) = static_cast<double>(cnt);

            if (!status.IsOk())
            {
                if (status.GetArg1() == 111)
                {
                    status.SetUserFuncName(funcName);
                    status.SetArg1(1);
                }
                else if (status.GetArg1() == 6)
                {
                    status.SetArg1(4);
                }
                else if (status.GetArg1() == 7)
                {
                    status.SetArg1(5);
                }

                quad_oml_sys_func_stack.clear();
                quad_oml_sys_name_stack.clear();
                quad_oml_sys_pntr_stack.clear();
                quad_eval_ptr_stack.clear();

                BuiltInFuncsUtils::CheckMathStatus(eval, status);
            }
        }

        // pack outputs
        size_t nargout = eval.GetNargoutValue();

        if (nargout >= 0)
            outputs.push_back(area.release());

        if (nargout > 1)
            outputs.push_back(count.release());
    }
    else if (inputs[1].IsMatrix() && inputs[2].IsScalar())
    {
        int i   = 0;
        int cnt = 0;

        const hwMatrix* A = inputs[1].Matrix();

        if (!A->IsReal())
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_MATRIX);

        double b = inputs[2].Scalar();
        std::unique_ptr<hwMatrix> area(EvaluatorInterface::allocateMatrix(A->M(), A->N(), hwMatrix::REAL));
        std::unique_ptr<hwMatrix> count(EvaluatorInterface::allocateMatrix(A->M(), A->N(), hwMatrix::REAL));

        for (i = 0; i < A->Size(); ++i)
        {
            cnt = 0;
            status = Quad(quad_file_func, (*A)(i), b, (*area)(i), cnt, reltol, abstol);
            (*count)(i) = (double) cnt;

            if (!status.IsOk())
            {
                if (status.GetArg1() == 111)
                {
                    status.SetUserFuncName(funcName);
                    status.SetArg1(1);
                }
                else if (status.GetArg1() == 6)
                {
                    status.SetArg1(4);
                }
                else if (status.GetArg1() == 7)
                {
                    status.SetArg1(5);
                }

                quad_oml_sys_func_stack.clear();
                quad_oml_sys_name_stack.clear();
                quad_oml_sys_pntr_stack.clear();
                quad_eval_ptr_stack.clear();

                BuiltInFuncsUtils::CheckMathStatus(eval, status);
            }
        }

        // pack outputs
        size_t nargout = eval.GetNargoutValue();

        if (nargout >= 0)
            outputs.push_back(area.release());

        if (nargout > 1)
            outputs.push_back(count.release());
    }
    else if (inputs[1].IsMatrix() && inputs[2].IsMatrix())
    {
        const hwMatrix* A = inputs[1].Matrix();
        const hwMatrix* B = inputs[2].Matrix();

        if (!A->IsReal())
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_MATRIX);

        if (!B->IsReal())
            throw OML_Error(OML_ERR_REAL, 3, OML_VAR_MATRIX);

        if (A->M() != B->M() || A->M() != B->M())
            throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

        std::unique_ptr<hwMatrix> area(EvaluatorInterface::allocateMatrix(
                                       A->M(), A->N(), hwMatrix::REAL));
        std::unique_ptr<hwMatrix> count(EvaluatorInterface::allocateMatrix(
                                       A->M(), A->N(), hwMatrix::REAL));
        int i   = 0;
        int cnt = 0;

        for (i = 0; i < A->Size(); ++i)
        {
            cnt = 0;
            status = Quad(quad_file_func, (*A)(i), (*B)(i), (*area)(i), cnt, reltol, abstol);
            (*count)(i) = (double) cnt;

            if (!status.IsOk())
            {
                if (status.GetArg1() == 111)
                {
                    status.SetUserFuncName(funcName);
                    status.SetArg1(1);
                }
                else if (status.GetArg1() == 6)
                {
                    status.SetArg1(4);
                }
                else if (status.GetArg1() == 7)
                {
                    status.SetArg1(5);
                }

                quad_oml_sys_func_stack.clear();
                quad_oml_sys_name_stack.clear();
                quad_oml_sys_pntr_stack.clear();
                quad_eval_ptr_stack.clear();

                BuiltInFuncsUtils::CheckMathStatus(eval, status);
            }
        }

        // pack outputs
        size_t nargout = eval.GetNargoutValue();

        if (nargout >= 0)
            outputs.push_back(area.release());

        if (nargout > 1)
            outputs.push_back(count.release());
    }

    quad_oml_sys_func_stack.pop_back();
    quad_oml_sys_name_stack.pop_back();
    quad_oml_sys_pntr_stack.pop_back();
    quad_eval_ptr_stack.pop_back();

    return true;
}
//------------------------------------------------------------------------------
// Performs numerical integration using Simpson's rule and return true
//------------------------------------------------------------------------------
bool OmlQuadv(EvaluatorInterface           eval, 
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs)
{
    // unpack inputs
    int nargin = eval.GetNarginValue();

    if (nargin < 3 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::string funcName;

    if (inputs[0].IsFunctionHandle())
        funcName = inputs[0].FunctionHandle()->FunctionName();
    else if (inputs[0].IsString())
        funcName = inputs[0].StringVal();
    else
        throw OML_Error(OML_ERR_HANDLE, 1, OML_VAR_TYPE);

    FunctionInfo* funcInfo = NULL;
    FUNCPTR       funcPntr = NULL;

    if (funcName == "anonymous")
    {
        funcInfo = inputs[0].FunctionHandle();
    }
    else
    {
        if (!eval.FindFunctionByName(funcName, &funcInfo, &funcPntr))
            throw OML_Error(OML_ERR_FUNCNAME, 1);
    }

    if (funcInfo && funcInfo->Parameters().size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN, 1);

    if (!inputs[1].IsScalar() && !inputs[1].IsMatrix())
        throw OML_Error(OML_ERR_SCALARMATRIX, 2);

    if (!inputs[2].IsScalar() && !inputs[2].IsMatrix())
        throw OML_Error(OML_ERR_SCALARMATRIX, 3);

    double abstol = sqrt(MACHEP);

    if (nargin > 3)
    {
        if (inputs[3].IsScalar())
            abstol = inputs[3].Scalar();
        else
            throw OML_Error(OML_ERR_NUMERIC, 4, OML_VAR_ABSTOL);
    }

    // set file scope variables
    quad_oml_sys_func_stack.push_back(funcInfo);
    quad_oml_sys_name_stack.push_back(funcName);
    quad_oml_sys_pntr_stack.push_back(funcPntr);
    quad_eval_ptr_stack.push_back(&eval);

    // call algorithm
    hwMathStatus status;

    if (inputs[1].IsScalar() && inputs[2].IsScalar())
    {
        double a = inputs[1].Scalar();
        double b = inputs[2].Scalar();
        double area;
        int count = 0;

        status = QuadV(quadv_file_func, a, b, area, count, abstol);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 111)
            {
                status.SetUserFuncName(funcName);
                status.SetArg1(1);
            }
            else if (status.GetArg1() == 6)
            {
                status.SetArg1(4);
            }
            else if (status.GetArg1() == 7)
            {
                status.SetArg1(4);
            }

            quad_oml_sys_func_stack.clear();
            quad_oml_sys_name_stack.clear();
            quad_oml_sys_pntr_stack.clear();
            quad_eval_ptr_stack.clear();

            BuiltInFuncsUtils::CheckMathStatus(eval, status);
        }

        // pack outputs
        size_t nargout = eval.GetNargoutValue();

        if (nargout >= 0)
            outputs.push_back(area);

        if (nargout > 1)
            outputs.push_back(count);
    }
    else if (inputs[1].IsScalar() && inputs[2].IsMatrix())
    {
        int i;
        double a = inputs[1].Scalar();
        const hwMatrix* B = inputs[2].Matrix();
        std::unique_ptr<hwMatrix> area(EvaluatorInterface::allocateMatrix(B->M(), B->N(), hwMatrix::REAL));
        std::unique_ptr<hwMatrix> count(EvaluatorInterface::allocateMatrix(B->M(), B->N(), hwMatrix::REAL));
        int cnt;

        if (!B->IsReal())
            throw OML_Error(OML_ERR_REAL, 3, OML_VAR_MATRIX);

        for (i = 0; i < B->Size(); ++i)
        {
            cnt = 0;
            status = QuadV(quadv_file_func, a, (*B)(i), (*area)(i), cnt, abstol);
            (*count)(i) = (double) cnt;

            if (!status.IsOk())
            {
                if (status.GetArg1() == 111)
                {
                    status.SetUserFuncName(funcName);
                    status.SetArg1(1);
                }
                else if (status.GetArg1() == 6)
                {
                    status.SetArg1(4);
                }
                else if (status.GetArg1() == 7)
                {
                    status.SetArg1(5);
                }

                quad_oml_sys_func_stack.clear();
                quad_oml_sys_name_stack.clear();
                quad_oml_sys_pntr_stack.clear();
                quad_eval_ptr_stack.clear();

                BuiltInFuncsUtils::CheckMathStatus(eval, status);
            }
        }

        // pack outputs
        size_t nargout = eval.GetNargoutValue();

        if (nargout >= 0)
            outputs.push_back(area.release());

        if (nargout > 1)
            outputs.push_back(count.release());
    }
    else if (inputs[1].IsMatrix() && inputs[2].IsScalar())
    {
        int i;
        const hwMatrix* A = inputs[1].Matrix();
        double b = inputs[2].Scalar();
        std::unique_ptr<hwMatrix> area(EvaluatorInterface::allocateMatrix(A->M(), A->N(), hwMatrix::REAL));
        std::unique_ptr<hwMatrix> count(EvaluatorInterface::allocateMatrix(A->M(), A->N(), hwMatrix::REAL));
        int cnt;

        if (!A->IsReal())
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_MATRIX);

        for (i = 0; i < A->Size(); ++i)
        {
            cnt = 0;
            status = QuadV(quadv_file_func, (*A)(i), b, (*area)(i), cnt, abstol);
            (*count)(i) = (double) cnt;

            if (!status.IsOk())
            {
                if (status.GetArg1() == 111)
                {
                    status.SetUserFuncName(funcName);
                    status.SetArg1(1);
                }
                else if (status.GetArg1() == 6)
                {
                    status.SetArg1(4);
                }
                else if (status.GetArg1() == 7)
                {
                    status.SetArg1(5);
                }

                quad_oml_sys_func_stack.clear();
                quad_oml_sys_name_stack.clear();
                quad_oml_sys_pntr_stack.clear();
                quad_eval_ptr_stack.clear();

                BuiltInFuncsUtils::CheckMathStatus(eval, status);
            }
        }

        // pack outputs
        size_t nargout = eval.GetNargoutValue();

        if (nargout >= 0)
            outputs.push_back(area.release());

        if (nargout > 1)
            outputs.push_back(count.release());
    }
    else if (inputs[1].IsMatrix() && inputs[2].IsMatrix())
    {
        int i;
        const hwMatrix* A = inputs[1].Matrix();
        const hwMatrix* B = inputs[2].Matrix();

        if (!A->IsReal())
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_MATRIX);

        if (!B->IsReal())
            throw OML_Error(OML_ERR_REAL, 3, OML_VAR_MATRIX);

        if (A->M() != B->M() || A->M() != B->M())
            throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

        std::unique_ptr<hwMatrix> area(EvaluatorInterface::allocateMatrix(A->M(), A->N(), hwMatrix::REAL));
        std::unique_ptr<hwMatrix> count(EvaluatorInterface::allocateMatrix(A->M(), A->N(), hwMatrix::REAL));
        int cnt;

        for (i = 0; i < A->Size(); ++i)
        {
            cnt = 0;
            status = QuadV(quadv_file_func, (*A)(i), (*B)(i), (*area)(i), cnt, abstol);
            (*count)(i) = (double) cnt;

            if (!status.IsOk())
            {
                if (status.GetArg1() == 111)
                {
                    status.SetUserFuncName(funcName);
                    status.SetArg1(1);
                }
                else if (status.GetArg1() == 6)
                {
                    status.SetArg1(4);
                }
                else if (status.GetArg1() == 7)
                {
                    status.SetArg1(5);
                }

                quad_oml_sys_func_stack.clear();
                quad_oml_sys_name_stack.clear();
                quad_oml_sys_pntr_stack.clear();
                quad_eval_ptr_stack.clear();

                BuiltInFuncsUtils::CheckMathStatus(eval, status);
            }
        }

        // pack outputs
        size_t nargout = eval.GetNargoutValue();

        if (nargout >= 0)
            outputs.push_back(area.release());

        if (nargout > 1)
            outputs.push_back(count.release());
    }

    quad_oml_sys_func_stack.pop_back();
    quad_oml_sys_name_stack.pop_back();
    quad_oml_sys_pntr_stack.pop_back();
    quad_eval_ptr_stack.pop_back();

    return true;
}
//------------------------------------------------------------------------------
// Performs numerical integration of discrete data using the trapezoid rule
//------------------------------------------------------------------------------
bool OmlTrapz(EvaluatorInterface           eval, 
              const std::vector<Currency>& inputs,
              std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_VARIABLE);

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_VECTOR, 2, OML_VAR_VARIABLE);

    const hwMatrix* mtx1 = inputs[0].ConvertToMatrix();
    const hwMatrix* mtx2 = inputs[1].ConvertToMatrix();
    double result;

    BuiltInFuncsUtils::CheckMathStatus(eval, TrapZ(*mtx1, *mtx2, result));
    outputs.push_back(result);

    return true;
}
//------------------------------------------------------------------------------
// Cumulative numerical integration of discrete data using the trapezoid rule
//------------------------------------------------------------------------------
bool OmlCumtrapz(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_VARIABLE);

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_VECTOR, 2, OML_VAR_VARIABLE);

    const hwMatrix* mtx1 = inputs[0].ConvertToMatrix();
    const hwMatrix* mtx2 = inputs[1].ConvertToMatrix();
    hwMatrix* result = EvaluatorInterface::allocateMatrix();

    BuiltInFuncsUtils::CheckMathStatus(eval, CumTrapZ(*mtx1, *mtx2, *result));
    outputs.push_back(result);

    return true;
}
//------------------------------------------------------------------------------
// Returns toolbox version
//------------------------------------------------------------------------------
double GetToolboxVersion(EvaluatorInterface eval)
{
    return TBOXVERSION;
}
