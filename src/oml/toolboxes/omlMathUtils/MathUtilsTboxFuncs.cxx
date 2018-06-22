/**
* @file MathUtilsTboxFuncs.cxx
* @date January 2015
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

#include "MathUtilsTboxFuncs.h"

#include <memory>  // For std::unique_ptr

#include "BuiltInFuncs.h"
#include "BuiltInFuncsUtils.h"
#include "OML_Error.h"
#include "MathUtilsFuncs.h"
#include "SpecialFuncs.h"

#define STATAN "StatisticalAnalysis"

//------------------------------------------------------------------------------
// Entry point which registers toolbox with oml
//------------------------------------------------------------------------------
int InitDll(EvaluatorInterface eval)
{
    eval.RegisterBuiltInFunction("beta", OmlBeta, 
                                 FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("gamma", OmlGamma, 
                                 FunctionMetaData(1, 1, STATAN));
    eval.RegisterBuiltInFunction("bins", OmlBins, 
                                 FunctionMetaData(2, 1, STATAN));
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
