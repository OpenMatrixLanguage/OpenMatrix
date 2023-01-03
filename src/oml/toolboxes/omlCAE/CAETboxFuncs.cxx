/**
* @file CAETboxFuncs.cxx
* @date February 2017
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

#include "CAETboxFuncs.h"

#include <memory>  // For std::unique_ptr

#include "CAEFuncs.h"
#include "hwMatrix.h"

#include "BuiltInFuncsUtils.h"
#include "OML_Error.h"

#define CAE "CAE"
#define TBOXVERSION 2019.0

//------------------------------------------------------------------------------
// Entry point which registers omlCAE toolbox with oml
//------------------------------------------------------------------------------
int InitDll(EvaluatorInterface eval)
{
    eval.RegisterBuiltInFunction("rainflow",  oml_rainflow, 
                                 FunctionMetaData(5, 1, CAE));
    eval.RegisterBuiltInFunction("iso6487",    OmlISO6487,
                                 FunctionMetaData(3, 1, CAE));
    eval.RegisterBuiltInFunction("saefilter",  OmlSAEfilter,
                                 FunctionMetaData(-4, 2, CAE));
    eval.RegisterBuiltInFunction("saefilt95",  OmlSAEfilt95,
                                 FunctionMetaData(5, 1, CAE));
    return 1;
}
//------------------------------------------------------------------------------
// Rainflow counting fatigue analysis function
//------------------------------------------------------------------------------
bool oml_rainflow(EvaluatorInterface           eval, 
                  const std::vector<Currency>& inputs, 
                  std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);
    
    if (!inputs[0].IsMatrix())
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);

    const hwMatrix* signal = inputs[0].Matrix();

    // defaults
    int numBins = 64;
    int hysteresis = 0;     // place holder
    std::string method = "3pt-ASTM";

    if (nargin > 1)
    {
        if (inputs[1].IsPositiveInteger())
        {
            numBins = static_cast<int>(inputs[1].Scalar());
        }
        else if (!(inputs[1].IsMatrix() && inputs[1].Matrix()->Is0x0()))
        {
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);
        }
    }

#if 0 // Unused code - hysteris input is not used
    if (!inputs[4].IsInteger())
        throw OML_Error(OML_ERR_NATURALNUM, 5, OML_VAR_VALUE);

    int hysteresis = (int) inputs[4].Scalar();
#endif 

    if (nargin > 2)
    {
        if (inputs[2].IsString())
        {
            method = inputs[2].StringVal();
        }
        else if (!(inputs[2].IsMatrix() && inputs[2].Matrix()->Is0x0()))
        {
            throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);
        }
    }

    std::unique_ptr<hwMatrix> H(EvaluatorInterface::allocateMatrix());
    std::unique_ptr<hwMatrix> R(EvaluatorInterface::allocateMatrix());
    std::unique_ptr<hwMatrix> M(EvaluatorInterface::allocateMatrix());

    if (method == "3pt-serial")     // HG method
    {
        std::unique_ptr<hwMatrix> C(EvaluatorInterface::allocateMatrix());  // output

        hwMathStatus status = RainFlowFunc(*signal, numBins, hysteresis,
                                           *C, *H, *R, *M);

        BuiltInFuncsUtils::CheckMathStatus(eval, status);
        outputs.push_back(C.release());
    }
    else
    {
        // get the cycle count matrix from a script
        std::vector<Currency> inputs2;
        inputs2.push_back(inputs[0]);
        Currency C;

        if (method == "3pt-ASTM")
        {
            // call Beck Chen 4 point algorithm script 
            C = eval.CallFunction("rainflow3_ASTM", inputs2);
        }
        else if (method == "3pt-recursive")
        {
            // call Beck Chen 4 point algorithm script 
            C = eval.CallFunction("rainflow3r", inputs2);
        }
        else if (method == "4pt-serial")
        {
            // call Beck Chen 4 point algorithm script 
            C = eval.CallFunction("rainflow4", inputs2);
        }
        else if (method == "4pt-recursive")
        {
            // call Beck Chen 4 point algorithm script 
            C = eval.CallFunction("rainflow4r", inputs2);
        }
        else
        {
            throw OML_Error(OML_ERR_OPTION, 3);
        }

        if (!C.IsMatrix())
        {
            throw OML_Error(OML_ERR_INTERNAL);
        }

        hwMatrix* CM = C.GetWritableMatrix();  // becomes an input

        hwMathStatus status = RainFlowFunc(*signal, numBins, hysteresis,
                                           CM, *H, *R, *M);

        BuiltInFuncsUtils::CheckMathStatus(eval, status);
        outputs.push_back(C);
    }

    outputs.push_back(H.release());
    outputs.push_back(R.release());
    outputs.push_back(M.release());

    return true;
}
//------------------------------------------------------------------------------
// Filters a signal with an ISO6487 filter [iso6487 command]
//------------------------------------------------------------------------------
bool OmlISO6487(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();
    if (nargin != 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_TYPE);
    }
    if (!inputs[1].IsScalar())
    {
        throw OML_Error(OML_ERR_SCALAR, 2, OML_VAR_TYPE);
    }
    if (!inputs[2].IsScalar())
    {
        throw OML_Error(OML_ERR_SCALAR, 3, OML_VAR_TYPE);
    }

    const hwMatrix* inSignal = inputs[0].ConvertToMatrix();
    double          sampFreq = inputs[1].Scalar();
    double          cfc      = inputs[2].Scalar();

    std::unique_ptr<hwMatrix> outSignal(EvaluatorInterface::allocateMatrix());

    if (inSignal->IsVector())
    {
        hwMathStatus status = ISO6487(*inSignal, sampFreq, cfc, *outSignal);
        BuiltInFuncsUtils::CheckMathStatus(eval, status);

        if (inSignal->M() != outSignal->M())
        {
            outSignal->Transpose();
        }
    }
    else
    {
        outSignal.reset(EvaluatorInterface::allocateMatrix(inSignal->M(), inSignal->N(), 0.0));
        int outM = inSignal->M();

        for (int i = 0; i < inSignal->N(); ++i)
        {
            hwMatrix inCol(outM, 1, (void*) &(*inSignal)(0,i), hwMatrix::REAL);
            hwMatrix outCol(outM, 1, (void*) &(*outSignal)(0,i), hwMatrix::REAL);

            hwMathStatus status = ISO6487(inCol, sampFreq, cfc, outCol);
            BuiltInFuncsUtils::CheckMathStatus(eval, status);
        }
    }

    outputs.push_back(outSignal.release());
    return true;
}
//------------------------------------------------------------------------------
// Filters a signal with an SAE filter [saefilter command]
//------------------------------------------------------------------------------
bool OmlSAEfilter(EvaluatorInterface           eval, 
                  const std::vector<Currency>& inputs, 
                  std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();
    if (nargin != 3 && nargin != 4)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_TYPE);
    }
    if (!inputs[1].IsScalar())
    {
        throw OML_Error(OML_ERR_SCALAR, 2, OML_VAR_TYPE);
    }
    if (!inputs[2].IsPositiveInteger())
    {
        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_TYPE);
    }

    const hwMatrix* inSignal = inputs[0].ConvertToMatrix();
    double          sampFreq = inputs[1].Scalar();
    double          cfc      = inputs[2].Scalar();
    int             fftSize  = 0;

    if (nargin > 3)
    {
        if (!inputs[3].IsPositiveInteger())
        {
           throw OML_Error(OML_ERR_POSINTEGER, 4, OML_VAR_TYPE);
        }
        fftSize = static_cast<int>(inputs[3].Scalar());
    }

    std::unique_ptr<hwMatrix> outSignal(EvaluatorInterface::allocateMatrix());

    if (inSignal->IsVector())
    {
        hwMathStatus status = SAEFilter(*inSignal, sampFreq, cfc, *outSignal, fftSize);
        BuiltInFuncsUtils::CheckMathStatus(eval, status);

        if ((inSignal->M() == 1 && outSignal->M() != 1) ||
            (inSignal->N() == 1 && outSignal->N() != 1))
        {
            outSignal->Transpose();
        }
    }
    else
    {
        int outM = (fftSize) ? fftSize : inSignal->M();

        outSignal.reset(EvaluatorInterface::allocateMatrix(outM, inSignal->N(), 0.0));

        if (inSignal->M())
        {
            for (int i = 0; i < inSignal->N(); ++i)
            {
                hwMatrix inCol(outM, 1, (void*) &(*inSignal)(0,i), hwMatrix::REAL);
                hwMatrix outCol(outM, 1, (void*) &(*outSignal)(0,i), hwMatrix::REAL);

                hwMathStatus status = SAEFilter(inCol, sampFreq, cfc, outCol, fftSize);
                BuiltInFuncsUtils::CheckMathStatus(eval, status);
            }
        }
    }

    outputs.push_back(outSignal.release());
    return true;
}
//------------------------------------------------------------------------------
// Filters a signal with an SAE class filter [saefilt95 command]
//------------------------------------------------------------------------------
bool OmlSAEfilt95(EvaluatorInterface           eval, 
                  const std::vector<Currency>& inputs, 
                  std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();
    if (nargin != 5)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_TYPE);
    }
    if (!inputs[1].IsScalar())
    {
        throw OML_Error(OML_ERR_SCALAR, 2, OML_VAR_TYPE);
    }
    if (!inputs[2].IsPositiveInteger())
    {
        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_TYPE);
    }
    if (!inputs[3].IsInteger())
    {
        throw OML_Error(OML_ERR_INTEGER, 4, OML_VAR_TYPE);
    }
    if (!inputs[4].IsPositiveInteger())
    {
        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_TYPE);
    }

    const hwMatrix* inSignal = inputs[0].ConvertToMatrix();
    double          sampFreq = inputs[1].Scalar();
    double          cfc      = inputs[2].Scalar();
    int             stdpad   = static_cast<int>(inputs[3].Scalar());
    int             dir      = static_cast<int>(inputs[4].Scalar());

    std::unique_ptr<hwMatrix> outSignal(EvaluatorInterface::allocateMatrix());

    if (inSignal->IsVector())
    {
        hwMathStatus status = SAEFilt95(*inSignal, sampFreq, cfc, stdpad, dir, 
                                      *outSignal);
        BuiltInFuncsUtils::CheckMathStatus(eval, status);

        if (inSignal->M() != outSignal->M())
        {
            outSignal->Transpose();
        }
    }
    else
    {
        outSignal.reset(EvaluatorInterface::allocateMatrix(inSignal->M(), inSignal->N(), 0.0));
        int outM = inSignal->M();

        for (int i = 0; i < inSignal->N(); ++i)
        {
            hwMatrix inCol(outM, 1, (void*) &(*inSignal)(0,i), hwMatrix::REAL);
            hwMatrix outCol(outM, 1, (void*) &(*outSignal)(0,i), hwMatrix::REAL);

            hwMathStatus status = SAEFilt95(*inSignal, sampFreq, cfc, stdpad, dir, 
                                          *outSignal);
            BuiltInFuncsUtils::CheckMathStatus(eval, status);
        }
    }

    outputs.push_back(outSignal.release());
    return true;
}
//------------------------------------------------------------------------------
// Returns toolbox version
//------------------------------------------------------------------------------
double GetToolboxVersion(EvaluatorInterface eval)
{
    return TBOXVERSION;
}
