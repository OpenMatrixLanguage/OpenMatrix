/**
* @file MathUtilsTboxFuncs.cxx
* @date January 2015
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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

#include "MathUtilsTboxFuncs.h"

#include <memory>  // For std::unique_ptr

#include "BuiltInFuncs.h"
#include "BuiltInFuncsUtils.h"
#include "OML_Error.h"
#include "MathUtilsFuncs.h"
#include "SpecialFuncs.h"
#include "hwMatrix.h"
#include "MatrixNUtils.h"

#if defined(_DARWIN) || defined(LINUX)
  #include <stdlib.h>   // for _fcvt in omlRat
#endif

#define ELEM   "ElementaryMath"
#define STATAN "StatisticalAnalysis"
#define TBOXVERSION 2019.0

//------------------------------------------------------------------------------
// Entry point which registers toolbox with oml
//------------------------------------------------------------------------------
int InitDll(EvaluatorInterface eval)
{
    eval.RegisterBuiltInFunction("beta", OmlBeta, 
                                 FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("gamma", OmlGamma, 
                                 FunctionMetaData(1, 1, STATAN));
    eval.RegisterBuiltInFunction("factorial", OmlFactorial,
                                 FunctionMetaData(1, 1, STATAN));
    eval.RegisterBuiltInFunction("bins", OmlBins,
                                 FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("rat", OmlRat, 
                                 FunctionMetaData(-2, -2, ELEM));
    return 1;
}
//------------------------------------------------------------------------------
// Beta function
//------------------------------------------------------------------------------
bool OmlBeta(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsScalar() && inputs[1].IsScalar())
        outputs = doBeta(eval, inputs);
    else
        outputs.push_back(mtxFun(eval, inputs, 1, &doBeta)[0]);

    return true;
}
//------------------------------------------------------------------------------
// Executes the beta function and returns outputs
//------------------------------------------------------------------------------
std::vector<Currency> doBeta(EvaluatorInterface&          eval, 
                             const std::vector<Currency>& inputs)
{
    if (!inputs[0].IsScalar())
        throw OML_Error(OML_ERR_SCALAR, 1, OML_VAR_TYPE);

    if (!inputs[1].IsScalar())
        throw OML_Error(OML_ERR_SCALAR, 2, OML_VAR_TYPE);

    double d1 = inputs[0].Scalar();
    double d2 = inputs[1].Scalar();
    std::vector<Currency> result;

    if (IsNaN_T<double>(d1) || isinfinity(d1) || 
        IsNaN_T<double>(d2) || isinfinity(d2))
    {
        result.push_back(std::numeric_limits<double>::quiet_NaN());
    }
    else
    {
        result.push_back(BetaFunc(d1, d2));
    }

    return result;
}
//------------------------------------------------------------------------------
//  Gamma function
//------------------------------------------------------------------------------
bool OmlGamma(EvaluatorInterface           eval, 
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsScalar())
        outputs = doGamma(eval, inputs);
    else
        outputs.push_back(mtxFun(eval, inputs, 1, &doGamma)[0]);

    return true;
}
//------------------------------------------------------------------------------
// Executes the gamma function and returns outputs
//------------------------------------------------------------------------------
std::vector<Currency> doGamma(EvaluatorInterface&          eval, 
                              const std::vector<Currency>& inputs)
{
    if (!inputs[0].IsScalar())
        throw OML_Error(OML_ERR_SCALAR, 1, OML_VAR_TYPE);

    std::vector<Currency> result;
    result.push_back(GammaFunc(inputs[0].Scalar()));
    return result;
}
//------------------------------------------------------------------------------
//  Factorial function
//------------------------------------------------------------------------------
bool OmlFactorial(EvaluatorInterface           eval,
                  const std::vector<Currency>& inputs,
                  std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsScalar())
        outputs = doFactorial(eval, inputs);
    else
        outputs.push_back(mtxFun(eval, inputs, 1, &doFactorial)[0]);

    return true;
}
//------------------------------------------------------------------------------
// Executes the factorial function and returns outputs
//------------------------------------------------------------------------------
std::vector<Currency> doFactorial(EvaluatorInterface&          eval,
                                  const std::vector<Currency>& inputs)
{
    if (!inputs[0].IsInteger())
        throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_TYPE);

    double k = inputs[0].Scalar();
    double factorial;

    if (k < 0.0)
        throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_TYPE);

    if (k < 10.0)
    {
        switch (static_cast<int> (k))
        {
        case 0:  factorial = 1.0;      break;
        case 1:  factorial = 1.0;      break;
        case 2:  factorial = 2.0;      break;
        case 3:  factorial = 6.0;      break;
        case 4:  factorial = 24.0;     break;
        case 5:  factorial = 120.0;    break;
        case 6:  factorial = 720.0;    break;
        case 7:  factorial = 5040.0;   break;
        case 8:  factorial = 40320.0;  break;
        case 9:  factorial = 362880.0; break;
        default: break;
        }
    }
    else
    {
        factorial = floor(GammaFunc(k + 1.0) + 0.5);
    }

    std::vector<Currency> result;
    result.push_back(factorial);
    return result;
}
//------------------------------------------------------------------------------
// Divides range of data into given number of equal bins aand returns true
//------------------------------------------------------------------------------
bool OmlBins(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_TYPE);

    if (!inputs[1].IsPositiveInteger())
        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_TYPE);

    const hwMatrix* data = inputs[0].ConvertToMatrix();
    int numBins          = static_cast<int>(inputs[1].Scalar());
    std::unique_ptr<hwMatrix> result(EvaluatorInterface::allocateMatrix(1, 
                                     numBins, hwMatrix::REAL));

    BuiltInFuncsUtils::CheckMathStatus(eval, Bins(*data, *result));

    outputs.push_back(result.release());
    return true;
}
//------------------------------------------------------------------------------
// Returns rational fraction approximation
//------------------------------------------------------------------------------
bool OmlRat(EvaluatorInterface           eval, 
            const std::vector<Currency>& inputs, 
            std::vector<Currency>&       outputs)
{
    int nargin  = static_cast<int> (inputs.size());
    int nargout = eval.GetNargoutValue();

    if (nargin != 1 && nargin != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargout > 2)
        throw OML_Error(OML_ERR_NUMARGOUT);

    double       tol = 1.0e-6;
    hwMathStatus status;

    if (nargin == 2)
    {
        if (!inputs[1].IsScalar())
            throw OML_Error(OML_ERR_POSITIVE_SCALAR, 2, OML_VAR_VALUE);

        tol = inputs[1].Scalar();
    }

    if (inputs[0].IsScalar())
    {
        double   value = inputs[0].Scalar();
        double   num;
        double   den;
        hwMatrix cfTerms;

        if (nargin == 1)
        {
            tol = fabs(value) * tol;
        }

        BuiltInFuncsUtils::CheckMathStatus(eval, ContFrac(value, tol, num, den, cfTerms));

        if (nargout < 2)
        {
            int  cfsize = cfTerms.Size();
            int  decimal;
            int  sign;
            std::string str;

            #if defined(_DARWIN) || defined(LINUX)
              char* argChar;
            #else
              char argChar[25];
            #endif

            if (cfTerms(0) != 0.0)
            {
                #if defined(_DARWIN) || defined(LINUX)
                  argChar = fcvt(cfTerms(0), 0, &decimal, &sign);
                #else
                  _fcvt_s(argChar, 25, cfTerms(0), 0, &decimal, &sign);
                #endif

                if (!sign)
                    str = std::string(argChar);
                else
                    str = "-" + std::string(argChar);
            }
            else
            {
                str = "0";
            }

            for (int i = 1; i < cfsize-1; ++i)
            {
                #if defined(_DARWIN) || defined(LINUX)
                  argChar = fcvt(cfTerms(i), 0, &decimal, &sign);
                #else
                  _fcvt_s(argChar, 25, cfTerms(i), 0, &decimal, &sign);
                #endif

                if (!sign)
                    str += " + 1/(" + std::string(argChar);
                else
                    str += " + 1/(-" + std::string(argChar);
            }

            if (cfsize > 1)
            {
                #if defined(_DARWIN) || defined(LINUX)
                  argChar = fcvt(cfTerms(cfsize - 1), 0, &decimal, &sign);
                #else
                  _fcvt_s(argChar, 25, cfTerms(cfsize - 1), 0, &decimal, &sign);
                #endif

                if (!sign)
                    str += " + 1/" + std::string(argChar);
                else
                    str += " + 1/(-" + std::string(argChar) + ")";

                for (int i = 1; i < cfsize - 1; ++i)
                    str += ")";
            }

            outputs.push_back(str);
        }
        else
        {
            outputs.push_back(num);
            outputs.push_back(den);
        }
    }
    else if (inputs[0].IsMatrix())
    {
        const hwMatrix* data = inputs[0].ConvertToMatrix();
        int             size = data->Size();
        hwMatrix        cfTerms;

        if (nargin == 1)
        {
            double norm;
            BuiltInFuncsUtils::CheckMathStatus(eval, data->Norm(norm, 1));
            tol = norm * tol;
        }

        if (nargout < 2)
        {
            double num;
            double den;
            std::unique_ptr<hwMatrix> result(EvaluatorInterface::allocateMatrix());

            for (int i = 0; i < size; ++i)
            {
                BuiltInFuncsUtils::CheckMathStatus(eval, ContFrac((*data)(i), tol, num, den, cfTerms));

                int  cfsize = cfTerms.Size();
                int  decimal;
                int  sign;
                std::string str;

                #if defined(_DARWIN) || defined(LINUX)
                  char* argChar;
                #else
                  char argChar[25];
                #endif

                if (cfTerms(0) != 0.0)
                {
                    #if defined(_DARWIN) || defined(LINUX)
                      argChar = fcvt(cfTerms(0), 0, &decimal, &sign);
                    #else
                      _fcvt_s(argChar, 25, cfTerms(0), 0, &decimal, &sign);
                    #endif

                    if (!sign)
                        str = std::string(argChar);
                    else
                        str = "-" + std::string(argChar);
                }
                else
                {
                    str = "0";
                }

                for (int j = 1; j < cfsize-1; ++j)
                {
                    #if defined(_DARWIN) || defined(LINUX)
                      argChar = fcvt(cfTerms(j), 0, &decimal, &sign);
                    #else
                      _fcvt_s(argChar, 25, cfTerms(j), 0, &decimal, &sign);
                    #endif

                    if (!sign)
                        str += " + 1/(" + std::string(argChar);
                    else
                        str += " + 1/(-" + std::string(argChar);
                }

                if (cfsize > 1)
                {
                    #if defined(_DARWIN) || defined(LINUX)
                      argChar = fcvt(cfTerms(cfsize-1), 0, &decimal, &sign);
                    #else
                      _fcvt_s(argChar, 25, cfTerms(cfsize-1), 0, &decimal, &sign);
                    #endif

                    if (!sign)
                        str += " + 1/" + std::string(argChar);
                    else
                        str += " + 1/(-" + std::string(argChar) + ")";

                    for (int j = 1; j < cfsize - 1; ++j)
                        str += ")";
                }

                int strlen = static_cast<int> (str.size());
                int n      = result->N();

                if (i == 0)
                {
                    Currency strC(str);
                    (*result) = *(strC.Matrix());
                }
                else
                {
                    if (n > strlen)
                    {
                        std::string spaces(n-strlen, ' ');
                        str = str + spaces;
                    }
                    else if (n < strlen)
                    {
                        hwMatrix filler(i, strlen-n, hwMatrix::REAL);
                        filler.SetElements(32.0);   // ASCII space
                        hwMatrix copy(*result);
                        status = result->InsertColumns(copy, n, filler);
                        n = strlen;
                    }

                    hwMatrix copy(*result);
                    Currency strC(str);
                    status = result->InsertRows(copy, i, *(strC.Matrix()));
                }
            }

            Currency out(result.release());
            out.SetMask(Currency::MASK_STRING);

            outputs.push_back(out);
        }
        else
        {
            double    numI;
            double    denI;
            int       m    = data->M();
            int       n    = data->N();
            hwMatrix* num  = EvaluatorInterface::allocateMatrix(m, n, 0.0);
            hwMatrix* den  = EvaluatorInterface::allocateMatrix(m, n, 0.0);

            for (int i = 0; i < size; ++i)
            {
                BuiltInFuncsUtils::CheckMathStatus(eval,
                    ContFrac((*data)(i), tol, numI, denI, cfTerms));

                (*num)(i) = static_cast<double>(numI);
                (*den)(i) = static_cast<double>(denI);
            }

            outputs.push_back(num);
            outputs.push_back(den);
        }
    }
    else if (inputs[0].IsNDMatrix())
    {
        if (nargout < 2)
        {
            // convert from ND to vector
            const hwMatrixN*      matrix  = inputs[0].MatrixN();
            const double*         data    = matrix->GetRealData();
            hwMatrix*             slice2D = new hwMatrix(matrix->Size(), (void*) data, hwMatrix::REAL);
            std::vector<Currency> inputs2;

            inputs2.push_back(slice2D);

            for (int i = 1; i < inputs.size(); ++i)
            {
                inputs2.push_back(inputs[i]);
            }

            return OmlRat(eval, inputs2, outputs);
        }
        else
        {
            return oml_MatrixNUtil1(eval, inputs, outputs, OmlRat);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_TYPE);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns toolbox version
//------------------------------------------------------------------------------
double GetToolboxVersion(EvaluatorInterface eval)
{
    return TBOXVERSION;
}
