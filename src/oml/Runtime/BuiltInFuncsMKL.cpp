/**
* @file BuiltInFuncsMKL.cpp
* @date September 2021
* Copyright (C) 2016-2021 Altair Engineering, Inc.  
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

// Begin defines/includes
#include "BuiltInFuncs.h"
#include "BuiltInFuncsMKL.h"
#include "Evaluator.h"
#include "OML_Error.h"

#include "hwComplex.h"
#include "hwMatrix.h"
#include "hwMatrixN.h"
#include "MKLutilities.h"

//------------------------------------------------------------------------------
// Helper method for elementary functions
//------------------------------------------------------------------------------
template <double    (*funcR)(double),
          hwComplex (*funcC)(const hwComplex&),
          void      (*funcR_MKL)(const MKL_INT, const double*, double*),
          void      (*funcC_MKL)(const MKL_INT, const MKL_Complex16*, MKL_Complex16*)>
static bool ElemFunc1(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& cur = inputs[0];

    if (cur.IsScalar())
    {
        double input = cur.Scalar();
        double output;
        funcR_MKL(1, &input, &output);
        outputs.push_back(output);
    }
    else if (cur.IsComplex())
    {
        outputs.push_back(funcC(cur.Complex()));
    }
    else if (cur.IsMatrix())
    {
        const hwMatrix* mtx = cur.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->IsReal());
        int size = mtx->Size();

        if (mtx->IsReal())
        {
            const double* input = mtx->GetRealData();
            double* output = result->GetRealData();

            funcR_MKL(size, input, output);
        }
        else
        {
            hwComplex* input = const_cast<hwComplex*> (mtx->GetComplexData());
            hwComplex* output = result->GetComplexData();
            MKL_Complex16* inputMKL = reinterpret_cast<MKL_Complex16*> (input);
            MKL_Complex16* outputMKL = reinterpret_cast<MKL_Complex16*> (output);

            funcC_MKL(size, inputMKL, outputMKL);
        }

        outputs.push_back(result);
    }
    else if (cur.IsNDMatrix())
    {
        const hwMatrixN* mtx = cur.MatrixN();
        const std::vector<int>& dims = mtx->Dimensions();
        hwMatrixN* result = new hwMatrixN(dims, mtx->Type());
        int size = mtx->Size();

        if (mtx->IsReal())
        {
            const double* input = mtx->GetRealData();
            double* output = result->GetRealData();

            funcR_MKL(size, input, output);
        }
        else
        {
            hwComplex* input = const_cast<hwComplex*> (mtx->GetComplexData());
            hwComplex* output = result->GetComplexData();
            MKL_Complex16* inputMKL = reinterpret_cast<MKL_Complex16*> (input);
            MKL_Complex16* outputMKL = reinterpret_cast<MKL_Complex16*> (output);

            funcC_MKL(size, inputMKL, outputMKL);
        }

        outputs.push_back(result);
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, 1, OML_VAR_TYPE);
    }

    return true;
}

//------------------------------------------------------------------------------
// Helper method for elementary functions
//------------------------------------------------------------------------------
template <double    (*funcRR)(double),
          hwComplex (*funcRC)(double),
          hwComplex (*funcCC)(const hwComplex&),
          void      (*funcR_MKL)(const MKL_INT, const double*, double*),
          void      (*funcC_MKL)(const MKL_INT, const MKL_Complex16*, MKL_Complex16*),
          bool      (*condFunc)(double)>
static bool ElemFunc2(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& cur = inputs[0];
    bool switchToComplex = false;

    if (cur.IsScalar())
    {
        double value = cur.Scalar();

        if (!IsNaN_T(value) && !condFunc(value))
        {
            outputs.push_back(funcRC(cur.Scalar()));
            return true;
        }
    }
    else if (cur.IsMatrix())
    {
        const hwMatrix* mtx = cur.Matrix();
        int size = mtx->Size();

        if (mtx->IsReal())
        {
            for (int i = 0; i < size; ++i)
            {
                if (!IsNaN_T((*mtx)(i)) && !condFunc((*mtx)(i)))
                {
                    switchToComplex = true;
                    break;
                }
            }
        }
    }
    else if (cur.IsNDMatrix())
    {
        const hwMatrixN* mtx = cur.MatrixN();
        int size = mtx->Size();

        if (mtx->IsReal())
        {
            for (int i = 0; i < size; ++i)
            {
                if (!IsNaN_T((*mtx)(i)) && !condFunc((*mtx)(i)))
                {
                    switchToComplex = true;
                    break;
                }
            }
        }
    }

    if (switchToComplex)
    {
        std::vector<Currency> inputs2;

        if (cur.IsScalar())
        {
            hwComplex* cmplx = new hwComplex(cur.Scalar(), 0.0);
            inputs2.push_back(cmplx);
        }
        else if (cur.IsMatrix())
        {
            const hwMatrix* mtx = cur.Matrix();
            hwMatrix* mtx2 = new hwMatrix;
            mtx2->PackComplex(*mtx);
            inputs2.push_back(mtx2);
        }
        else if (cur.IsNDMatrix())
        {
            const hwMatrixN* mtx = cur.MatrixN();
            hwMatrixN* mtx2 = new hwMatrixN;
            mtx2->PackComplex(*mtx);
            inputs2.push_back(mtx2);
        }

        return ElemFunc1<funcRR, funcCC, funcR_MKL, funcC_MKL>(eval, inputs2, outputs);
    }
    else
    {
        return ElemFunc1<funcRR, funcCC, funcR_MKL, funcC_MKL>(eval, inputs, outputs);
    }
}

//------------------------------------------------------------------------------
// Helper method for elementary functions
//------------------------------------------------------------------------------
template <double (*funcRR)(double),
          double (hwComplex::*funcCR)() const,
          void   (*funcRR_MKL)(const MKL_INT, const double*, double*),
          void   (*funcCR_MKL)(const MKL_INT, const MKL_Complex16*, double*)>
static bool ElemFunc3(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& cur = inputs[0];

    if (cur.IsScalar())
    {
        outputs.push_back(funcRR(cur.Scalar()));
    }
    else if (cur.IsComplex())
    {
        outputs.push_back((cur.Complex().*funcCR)());
    }
    else if (cur.IsMatrix() || cur.IsString())
    {
        const hwMatrix* mtx = cur.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), true);
        int size = mtx->Size();

        if (mtx->IsReal())
        {
            const double* input = mtx->GetRealData();
            double* output = result->GetRealData();

            funcRR_MKL(size, input, output);
        }
        else
        {
            hwComplex* input = const_cast<hwComplex*> (mtx->GetComplexData());
            MKL_Complex16* inputMKL = reinterpret_cast<MKL_Complex16*> (input);
            double* output = result->GetRealData();

            funcCR_MKL(size, inputMKL, output);
        }

        outputs.push_back(result);
    }
    else if (cur.IsNDMatrix())
    {
        const hwMatrixN* mtx = cur.MatrixN();
        const std::vector<int>& dims = mtx->Dimensions();
        hwMatrixN* result = new hwMatrixN(dims, hwMatrixN::REAL);
        int size = mtx->Size();

        if (mtx->IsReal())
        {
            const double* input = mtx->GetRealData();
            double* output = result->GetRealData();

            funcRR_MKL(size, input, output);
        }
        else
        {
            hwComplex* input = const_cast<hwComplex*> (mtx->GetComplexData());
            MKL_Complex16* inputMKL = reinterpret_cast<MKL_Complex16*> (input);
            double* output = result->GetRealData();

            funcCR_MKL(size, inputMKL, output);
        }

        outputs.push_back(result);
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, 1, OML_VAR_TYPE);
    }

    return true;
}

//------------------------------------------------------------------------------
// Helper method for elementary functions
//------------------------------------------------------------------------------
template <double (*funcR)(double),
          void   (*funcR_MKL)(const MKL_INT, const double*, double*)>
static bool ElemFunc4(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& cur = inputs[0];

    if (cur.IsScalar())
    {
        outputs.push_back(funcR(cur.Scalar()));
    }
    else if (cur.IsComplex())
    {
        hwComplex temp = cur.Complex();
        outputs.push_back(hwComplex(funcR(temp.Real()), funcR(temp.Imag())));
    }
    else if (cur.IsMatrix() || cur.IsString())
    {
        const hwMatrix* mtx = cur.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->IsReal());
        int size = mtx->Size();

        if (mtx->IsReal())
        {
            const double* input = mtx->GetRealData();
            double* output = result->GetRealData();

            funcR_MKL(size, input, output);
        }
        else
        {
            hwComplex* input = const_cast<hwComplex*> (mtx->GetComplexData());
            hwComplex* output = result->GetComplexData();
            double* dbl_in = reinterpret_cast<double*> (input);
            double* dbl_out = reinterpret_cast<double*> (output);

            funcR_MKL(2 * size, dbl_in, dbl_out);
        }

        outputs.push_back(result);
    }
    else if (cur.IsNDMatrix())
    {
        const hwMatrixN* mtx = cur.MatrixN();
        const std::vector<int>& dims = mtx->Dimensions();
        hwMatrixN* result = new hwMatrixN(dims, mtx->Type());
        int size = mtx->Size();

        if (mtx->IsReal())
        {
            const double* input = mtx->GetRealData();
            double* output = result->GetRealData();

            funcR_MKL(size, input, output);
        }
        else
        {
            hwComplex* input = const_cast<hwComplex*> (mtx->GetComplexData());
            hwComplex* output = result->GetComplexData();
            double* dbl_in = reinterpret_cast<double*> (input);
            double* dbl_out = reinterpret_cast<double*> (output);

            funcR_MKL(2 * size, dbl_in, dbl_out);
        }

        outputs.push_back(result);
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, 1, OML_VAR_TYPE);
    }

    return true;
}

//------------------------------------------------------------------------------
// Helper method for max/min functions
//------------------------------------------------------------------------------
template <double (*funcR)(double, double),
          void   (*funcR_MKL)(const MKL_INT, const double*, const double*, double*)>
static bool ElemFunc5(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&        outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];
    const Currency& input2 = inputs[1];

    if (input1.IsLogical() || input2.IsLogical())
        throw OML_Error(HW_ERROR_INPUTISLOGICAL);

    if (input1.IsScalar())
    {
        if (input2.IsScalar())
        {
            // outputs.push_back(funcR(input1.Scalar(), input2.Scalar()));
            double val1 = input1.Scalar();
            double val2 = input2.Scalar();
            double val3;
            funcR_MKL(1, &val1, &val2, &val3);  // do this for NaN support
            outputs.push_back(val3);
            return true;
        }
        else if (input2.IsMatrix())
        {
            hwMatrix* m2 = const_cast<hwMatrix*> (input2.Matrix());
            hwMatrix* temp = new hwMatrix(m2->M(), m2->N(), hwMatrix::REAL);
            temp->SetElements(input1.Scalar());

            std::vector<Currency> inputs2;
            inputs2.push_back(temp);
            m2->IncrRefCount();
            inputs2.push_back(m2);

            return ElemFunc5<funcR, funcR_MKL>(eval, inputs2, outputs);
        }
        else if (input2.IsNDMatrix())
        {
            hwMatrixN* m2 = const_cast<hwMatrixN*> (input2.MatrixN());
            hwMatrixN* temp = new hwMatrixN(m2->Dimensions(), hwMatrixN::REAL);
            temp->SetElements(input1.Scalar());

            std::vector<Currency> inputs2;
            inputs2.push_back(temp);
            m2->IncrRefCount();
            inputs2.push_back(m2);

            return ElemFunc5<funcR, funcR_MKL>(eval, inputs2, outputs);
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
        }
    }
    else if (input1.IsMatrix())
    {
        const hwMatrix* m1 = input1.Matrix();

        if (!m1->IsReal())
            throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

        if (input2.IsScalar())
        {
            hwMatrix* temp = new hwMatrix(m1->M(), m1->N(), hwMatrix::REAL);
            temp->SetElements(input2.Scalar());

            std::vector<Currency> inputs2;
            const_cast <hwMatrix*> (m1)->IncrRefCount();
            inputs2.push_back(const_cast<hwMatrix*> (m1));
            inputs2.push_back(temp);

            return ElemFunc5<funcR, funcR_MKL>(eval, inputs2, outputs);
        }
        else if (input2.IsMatrix())
        {
            const hwMatrix* m2 = input2.Matrix();

            if (!m2->IsReal())
                throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

            if (!sameSize(m1, m2))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            int size = m1->Size();
            hwMatrix* result = new hwMatrix(m1->M(), m2->N(), hwMatrix::REAL);
            funcR_MKL(size, m1->GetRealData(), m2->GetRealData(), result->GetRealData());
            outputs.push_back(result);
        }
        else if (input2.IsNDMatrix())
        {
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
        }
    }
    else if (input1.IsNDMatrix())
    {
        const hwMatrixN* m1 = input1.MatrixN();

        if (!m1->IsReal())
            throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

        if (input2.IsScalar())
        {
            hwMatrixN* temp = new hwMatrixN(m1->Dimensions(), hwMatrixN::REAL);
            temp->SetElements(input2.Scalar());

            std::vector<Currency> inputs2;
            const_cast <hwMatrixN*> (m1)->IncrRefCount();
            inputs2.push_back(const_cast<hwMatrixN*> (m1));
            inputs2.push_back(temp);

            return ElemFunc5<funcR, funcR_MKL>(eval, inputs2, outputs);
        }
        else if (input2.IsMatrix())
        {
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
        }
        else if (input2.IsNDMatrix())
        {
            const hwMatrixN* m2 = input2.MatrixN();

            if (!m2->IsReal())
                throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

            if (m1->Dimensions() != m2->Dimensions())
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            int size = m1->Size();
            hwMatrixN* result = new hwMatrixN(m1->Dimensions(), hwMatrixN::REAL);
            funcR_MKL(size, m1->GetRealData(), m2->GetRealData(), result->GetRealData());
            outputs.push_back(result);
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARMATRIX, 1, OML_VAR_TYPE);
    }

    return true;
}

//------------------------------------------------------------------------------
// Helper method to take reciprocals
//------------------------------------------------------------------------------
static void Reciprocal(const Currency& input,
                       Currency&       output)
{
    if (input.IsScalar())
    {
        output = Currency(1.0 / input.Scalar());
    }
    else if (input.IsComplex())
    {
        output = Currency(1.0 / input.Complex());
    }
    else if (input.IsMatrix())
    {
        const hwMatrix* m = input.Matrix();
        hwMatrix* r = new hwMatrix(m->M(), m->N(), m->Type());

        if (m->IsReal())
        {
            double* mdata = const_cast<double*> (m->GetRealData());
            double* rdata = r->GetRealData();
            vdInv(m->Size(), mdata, rdata);
        }
        else
        {
            // BuiltInFuncsMKL::PowerByElems(*m, -1.0, *r);
            int size = m->Size();
            hwComplex* mdata = const_cast<hwComplex*> (m->GetComplexData());
            hwComplex* rdata = r->GetComplexData();
            MKL_Complex16* mdataMKL = reinterpret_cast<MKL_Complex16*> (mdata);
            MKL_Complex16* rdataMKL = reinterpret_cast<MKL_Complex16*> (rdata);

            hwMatrix conj(size, 1, hwMatrix::COMPLEX);
            hwComplex* conjdata = conj.GetComplexData();
            MKL_Complex16* conjdataMKL = reinterpret_cast<MKL_Complex16*> (conjdata);
            vzConj(size, mdataMKL, conjdataMKL);

            hwMatrix magSq(size, 1, hwMatrix::COMPLEX);
            hwComplex* magSqdata = magSq.GetComplexData();
            MKL_Complex16* magSqdataMKL = reinterpret_cast<MKL_Complex16*> (magSqdata);
            vzMulByConj(size, mdataMKL, mdataMKL, magSqdataMKL);

            vzDiv(size, conjdataMKL, magSqdataMKL, rdataMKL);

            for (int i = 0; i < size; ++i)
            {
                if (conjdata[i] == 0.0)
                    rdata[i] = hwComplex(std::numeric_limits<double>::infinity(), 0.0);
            }
        }

        output = Currency(r);
    }
    else if (input.IsNDMatrix())
    {
        const hwMatrixN* m = input.MatrixN();
        hwMatrixN* r = new hwMatrixN(m->Dimensions(), m->Type());

        if (m->IsReal())
        {
            double* mdata = const_cast<double*> (m->GetRealData());
            double* rdata = r->GetRealData();
            vdInv(m->Size(), mdata, rdata);
        }
        else
        {
            int size = m->Size();
            hwComplex* mdata = const_cast<hwComplex*> (m->GetComplexData());
            hwComplex* rdata = r->GetComplexData();
            MKL_Complex16* mdataMKL = reinterpret_cast<MKL_Complex16*> (mdata);
            MKL_Complex16* rdataMKL = reinterpret_cast<MKL_Complex16*> (rdata);

            hwMatrix conj(size, 1, hwMatrix::COMPLEX);
            hwComplex* conjdata = conj.GetComplexData();
            MKL_Complex16* conjdataMKL = reinterpret_cast<MKL_Complex16*> (conjdata);
            vzConj(size, mdataMKL, conjdataMKL);

            hwMatrix magSq(size, 1, hwMatrix::COMPLEX);
            hwComplex* magSqdata = magSq.GetComplexData();
            MKL_Complex16* magSqdataMKL = reinterpret_cast<MKL_Complex16*> (magSqdata);
            vzMulByConj(size, mdataMKL, mdataMKL, magSqdataMKL);

            vzDiv(size, conjdataMKL, magSqdataMKL, rdataMKL);

            for (int i = 0; i < size; ++i)
            {
                if (conjdata[i] == 0.0)
                    rdata[i] = hwComplex(std::numeric_limits<double>::infinity(), 0.0);
            }
        }

        output = Currency(r);
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, 1, OML_VAR_TYPE);
    }
}

//------------------------------------------------------------------------------
// Returns cosine of input in radians [cos]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Cos(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    return ElemFunc1<cos, hwComplex::cos, vdCos, vzCos>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns sine of input in radians [sin]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Sin(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    return ElemFunc1<sin, hwComplex::sin, vdSin, vzSin>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns tangent of input in radians [tan]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Tan(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    return ElemFunc1<tan, hwComplex::tan, vdTan, vzTan>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns secant of input in radians [sec]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Sec(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    std::vector<Currency> outputs2;
    Cos(eval, inputs, outputs2);
    outputs.resize(1);
    Reciprocal(outputs2[0], outputs[0]);

    return true;
}

//------------------------------------------------------------------------------
// Returns secant of input in radians [csc]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Csc(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    std::vector<Currency> outputs2;
    Sin(eval, inputs, outputs2);
    outputs.resize(1);
    Reciprocal(outputs2[0], outputs[0]);

    return true;
}

//------------------------------------------------------------------------------
// Returns secant of input in radians [cot]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Cot(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    std::vector<Currency> outputs2;
    Tan(eval, inputs, outputs2);
    outputs.resize(1);
    Reciprocal(outputs2[0], outputs[0]);

    return true;
}

//------------------------------------------------------------------------------
// Returns cosine of input in degrees [cosd]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::CosD(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& cur = inputs[0];
    bool retv = false;

    if (cur.IsScalar())
    {
        outputs.push_back(cos(cur.Scalar() * (PI / 180.0)));
    }
    else if (cur.IsComplex())
    {
        hwComplex temp = cur.Complex() * (PI / 180.0);
        outputs.push_back(hwComplex::cos(temp));
        retv = true;
    }
    else if (cur.IsMatrix())
    {
        const hwMatrix* mtx = cur.Matrix();

        if (mtx->IsReal())
        {
            hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->IsReal());
            int size = mtx->Size();
            const double* input = mtx->GetRealData();
            double* output = result->GetRealData();

            vdCosd(size, input, output);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix* copy = new hwMatrix(*mtx);
            (*copy) *= (PI / 180.0);
            std::vector<Currency> inputs2;
            inputs2.push_back(copy);

            retv =  Cos(eval, inputs2, outputs);
        }
    }
    else if (cur.IsNDMatrix())
    {
        const hwMatrixN* mtx = cur.MatrixN();

        if (mtx->IsReal())
        {
            const std::vector<int>& dims = mtx->Dimensions();
            hwMatrixN* result = new hwMatrixN(dims, mtx->Type());
            int size = mtx->Size();
            const double* input = mtx->GetRealData();
            double* output = result->GetRealData();

            vdCosd(size, input, output);
            outputs.push_back(result);
        }
        else
        {
            hwMatrixN* copy = new hwMatrixN(*mtx);
            hwMatrix temp(copy->Size(), (void*)copy->GetComplexData(), hwMatrix::COMPLEX);
            temp *= (PI / 180.0);
            std::vector<Currency> inputs2;
            inputs2.push_back(copy);

            retv = Cos(eval, inputs2, outputs);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, 1, OML_VAR_TYPE);
    }

    if (retv)
        return retv;

    // tweak result for +/- 90
    std::vector<Currency> inputs2;
    std::vector<Currency> outputs2;

    if (cur.IsScalar())
    {
        double val = cur.Scalar() / 180.0 + 0.5;
        inputs2.push_back(val);
    }
    else if (cur.IsMatrix())
    {
        hwMatrix* mtx = new hwMatrix(*cur.Matrix());
        (*mtx) /= 180.0;
        (*mtx) += 0.5;
        inputs2.push_back(mtx);
    }
    else if (cur.IsNDMatrix())
    {
        hwMatrixN* mtxN = new hwMatrixN(*cur.MatrixN());
        int size = mtxN->Size();
        hwMatrix mtx(size, mtxN->GetRealData(), hwMatrix::REAL);
        mtx /= 180.0;
        mtx += hwComplex(0.5, 0.5);

        inputs2.push_back(mtxN);
    }

    Fix(eval, inputs2, outputs2);

    const Currency& cur1 = inputs2[0];
    Currency& cur2 = outputs2[0];

    if (cur1.IsScalar())
    {
        double val = cur1.Scalar();

        if (checkisfinite(val) && val == cur2.Scalar())
            outputs[0] = 0.0;
    }
    else if (cur1.IsMatrix())
    {
        const hwMatrix* mtx1 = cur1.Matrix();
        hwMatrix* mtx2 = cur2.GetWritableMatrix();
        hwMatrix* mtx3 = outputs[0].GetWritableMatrix();
        int size = mtx1->Size();

        if (mtx1->IsReal())
        {
            for (int i = 0; i < size; ++i)
            {
                if (checkisfinite((*mtx1)(i)) && (*mtx1)(i) == (*mtx2)(i))
                    (*mtx3)(i) = 0.0;
            }
        }
        else
        {
            for (int i = 0; i < size; ++i)
            {
                if (mtx1->z(i).IsReal() && checkisfinite(mtx1->z(i).Real()) &&
                    mtx1->z(i) == mtx2->z(i))
                    mtx3->z(i) = 0.0;
            }
        }
    }
    else if (cur1.IsNDMatrix())
    {
        const hwMatrixN* mtx1 = cur1.MatrixN();
        hwMatrixN* mtx2 = cur2.GetWritableMatrixN();
        hwMatrixN* mtx3 = outputs[0].GetWritableMatrixN();
        int size = mtx1->Size();

        if (mtx1->IsReal())
        {
            for (int i = 0; i < size; ++i)
            {
                if (checkisfinite((*mtx1)(i)) && (*mtx1)(i) == (*mtx2)(i))
                    (*mtx3)(i) = 0.0;
            }
        }
        else
        {
            for (int i = 0; i < size; ++i)
            {
                if (mtx1->z(i).IsReal() && checkisfinite(mtx1->z(i).Real()) &&
                    mtx1->z(i) == mtx2->z(i))
                    mtx3->z(i) = 0.0;
            }
        }
    }

    return true;
}

//------------------------------------------------------------------------------
// Returns sine of input in degrees [sind]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::SinD(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& cur = inputs[0];
    bool retv = false;

    if (cur.IsScalar())
    {
        outputs.push_back(sin(cur.Scalar() * (PI / 180.0)));
    }
    else if (cur.IsComplex())
    {
        hwComplex temp = cur.Complex() * (PI / 180.0);
        outputs.push_back(hwComplex::sin(temp));
        retv = true;
    }
    else if (cur.IsMatrix())
    {
        const hwMatrix* mtx = cur.Matrix();

        if (mtx->IsReal())
        {
            hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->IsReal());
            int size = mtx->Size();
            const double* input = mtx->GetRealData();
            double* output = result->GetRealData();

            vdSind(size, input, output);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix* copy = new hwMatrix(*mtx);
            (*copy) *= (PI / 180.0);
            std::vector<Currency> inputs2;
            inputs2.push_back(copy);

            retv = Sin(eval, inputs2, outputs);
        }
    }
    else if (cur.IsNDMatrix())
    {
        const hwMatrixN* mtx = cur.MatrixN();

        if (mtx->IsReal())
        {
            const std::vector<int>& dims = mtx->Dimensions();
            hwMatrixN* result = new hwMatrixN(dims, mtx->Type());
            int size = mtx->Size();
            const double* input = mtx->GetRealData();
            double* output = result->GetRealData();

            vdSind(size, input, output);
            outputs.push_back(result);
        }
        else
        {
            hwMatrixN* copy = new hwMatrixN(*mtx);
            hwMatrix temp(copy->Size(), (void*)copy->GetComplexData(), hwMatrix::COMPLEX);
            temp *= (PI / 180.0);
            std::vector<Currency> inputs2;
            inputs2.push_back(copy);

            retv = Sin(eval, inputs2, outputs);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, 1, OML_VAR_TYPE);
    }

    if (retv)
        return retv;

    // tweak result for +/- 180
    std::vector<Currency> inputs2;
    std::vector<Currency> outputs2;

    if (cur.IsScalar())
    {
        double val = cur.Scalar() / 180.0;
        inputs2.push_back(val);
    }
    else if (cur.IsMatrix())
    {
        hwMatrix* mtx = new hwMatrix(*cur.Matrix());
        (*mtx) /= 180.0;
        inputs2.push_back(mtx);
    }
    else if (cur.IsNDMatrix())
    {
        hwMatrixN* mtxN = new hwMatrixN(*cur.MatrixN());
        int size = mtxN->Size();
        hwMatrix mtx(size, mtxN->GetRealData(), hwMatrix::REAL);
        mtx /= 180.0;

        inputs2.push_back(mtxN);
    }

    Fix(eval, inputs2, outputs2);

    const Currency& cur1 = inputs2[0];
    Currency& cur2 = outputs2[0];

    if (cur1.IsScalar())
    {
        double val = cur1.Scalar();

        if (checkisfinite(val) && val == cur2.Scalar())
            outputs[0] = 0.0;
    }
    else if (cur1.IsMatrix())
    {
        const hwMatrix* mtx1 = cur1.Matrix();
        hwMatrix* mtx2 = cur2.GetWritableMatrix();
        hwMatrix* mtx3 = outputs[0].GetWritableMatrix();
        int size = mtx1->Size();

        if (mtx1->IsReal())
        {
            for (int i = 0; i < size; ++i)
            {
                if (checkisfinite((*mtx1)(i)) && (*mtx1)(i) == (*mtx2)(i))
                    (*mtx3)(i) = 0.0;
            }
        }
        else
        {
            for (int i = 0; i < size; ++i)
            {
                if (mtx1->z(i).IsReal() && checkisfinite(mtx1->z(i).Real()) &&
                    mtx1->z(i) == mtx2->z(i))
                    mtx3->z(i) = 0.0;
            }
        }
    }
    else if (cur1.IsNDMatrix())
    {
        const hwMatrixN* mtx1 = cur1.MatrixN();
        hwMatrixN* mtx2 = cur2.GetWritableMatrixN();
        hwMatrixN* mtx3 = outputs[0].GetWritableMatrixN();
        int size = mtx1->Size();

        if (mtx1->IsReal())
        {
            for (int i = 0; i < size; ++i)
            {
                if (checkisfinite((*mtx1)(i)) && (*mtx1)(i) == (*mtx2)(i))
                    (*mtx3)(i) = 0.0;
            }
        }
        else
        {
            for (int i = 0; i < size; ++i)
            {
                if (mtx1->z(i).IsReal() && checkisfinite(mtx1->z(i).Real()) &&
                    mtx1->z(i) == mtx2->z(i))
                    mtx3->z(i) = 0.0;
            }
        }
    }

    return true;
}

//------------------------------------------------------------------------------
// Returns tangent of input in degrees [tand]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::TanD(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& cur = inputs[0];
    bool retv = false;

    if (cur.IsScalar())
    {
        outputs.push_back(tan(cur.Scalar() * (PI / 180.0)));
    }
    else if (cur.IsComplex())
    {
        hwComplex temp = cur.Complex() * (PI / 180.0);
        outputs.push_back(hwComplex::tan(temp));
        retv = true;
    }
    else if (cur.IsMatrix())
    {
        const hwMatrix* mtx = cur.Matrix();

        if (mtx->IsReal())
        {
            hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->IsReal());
            int size = mtx->Size();
            const double* input = mtx->GetRealData();
            double* output = result->GetRealData();

            vdTand(size, input, output);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix* copy = new hwMatrix(*mtx);
            (*copy) *= (PI / 180.0);
            std::vector<Currency> inputs2;
            inputs2.push_back(copy);

            retv = Tan(eval, inputs2, outputs);
        }
    }
    else if (cur.IsNDMatrix())
    {
        const hwMatrixN* mtx = cur.MatrixN();

        if (mtx->IsReal())
        {
            const std::vector<int>& dims = mtx->Dimensions();
            hwMatrixN* result = new hwMatrixN(dims, mtx->Type());
            int size = mtx->Size();
            const double* input = mtx->GetRealData();
            double* output = result->GetRealData();

            vdTand(size, input, output);
            outputs.push_back(result);
        }
        else
        {
            hwMatrixN* copy = new hwMatrixN(*mtx);
            hwMatrix temp(copy->Size(), (void*)copy->GetComplexData(), hwMatrix::COMPLEX);
            temp *= (PI / 180.0);
            std::vector<Currency> inputs2;
            inputs2.push_back(copy);

            retv = Tan(eval, inputs2, outputs);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, 1, OML_VAR_TYPE);
    }

    if (retv)
        return retv;

    // tweak result for +/- 180
    std::vector<Currency> inputs2;
    std::vector<Currency> outputs2;

    if (cur.IsScalar())
    {
        double val = cur.Scalar() / 180.0;
        inputs2.push_back(val);
    }
    else if (cur.IsMatrix())
    {
        hwMatrix* mtx = new hwMatrix(*cur.Matrix());
        (*mtx) /= 180.0;
        inputs2.push_back(mtx);
    }
    else if (cur.IsNDMatrix())
    {
        hwMatrixN* mtxN = new hwMatrixN(*cur.MatrixN());
        int size = mtxN->Size();
        hwMatrix mtx(size, mtxN->GetRealData(), hwMatrix::REAL);
        mtx /= 180.0;

        inputs2.push_back(mtxN);
    }

    Fix(eval, inputs2, outputs2);

    Currency& cur1 = inputs2[0];   // cppcheck-suppress invalidContainerReference
    Currency& cur2 = outputs2[0];  // cppcheck-suppress invalidContainerReference

    if (cur1.IsScalar())
    {
        double val = cur1.Scalar();

        if (checkisfinite(val) && val == cur2.Scalar())
            outputs[0] = 0.0;
    }
    else if (cur1.IsMatrix())
    {
        const hwMatrix* mtx1 = cur1.Matrix();
        hwMatrix* mtx2 = cur2.GetWritableMatrix();
        hwMatrix* mtx3 = outputs[0].GetWritableMatrix();
        int size = mtx1->Size();

        if (mtx1->IsReal())
        {
            for (int i = 0; i < size; ++i)
            {
                if (checkisfinite((*mtx1)(i)) && (*mtx1)(i) == (*mtx2)(i))
                    (*mtx3)(i) = 0.0;
            }
        }
        else
        {
            for (int i = 0; i < size; ++i)
            {
                if (mtx1->z(i).IsReal() && checkisfinite(mtx1->z(i).Real()) &&
                    mtx1->z(i) == mtx2->z(i))
                    mtx3->z(i) = 0.0;
            }
        }
    }
    else if (cur1.IsNDMatrix())
    {
        const hwMatrixN* mtx1 = cur1.MatrixN();
        hwMatrixN* mtx2 = cur2.GetWritableMatrixN();
        hwMatrixN* mtx3 = outputs[0].GetWritableMatrixN();
        int size = mtx1->Size();

        if (mtx1->IsReal())
        {
            for (int i = 0; i < size; ++i)
            {
                if (checkisfinite((*mtx1)(i)) && (*mtx1)(i) == (*mtx2)(i))
                    (*mtx3)(i) = 0.0;
            }
        }
        else
        {
            for (int i = 0; i < size; ++i)
            {
                if (mtx1->z(i).IsReal() && checkisfinite(mtx1->z(i).Real()) &&
                    mtx1->z(i) == mtx2->z(i))
                    mtx3->z(i) = 0.0;
            }
        }
    }

    // tweak result for +/- 90
    inputs2.clear();
    outputs2.clear();

    if (cur.IsScalar())
    {
        double val = cur.Scalar() / 180.0 - 0.5;
        inputs2.push_back(val);
    }
    else if (cur.IsMatrix())
    {
        hwMatrix* mtx = new hwMatrix(*cur.Matrix());
        (*mtx) /= 180.0;
        (*mtx) -= 0.5;
        inputs2.push_back(mtx);
    }
    else if (cur.IsNDMatrix())
    {
        hwMatrixN* mtxN = new hwMatrixN(*cur.MatrixN());
        int size = mtxN->Size();
        hwMatrix mtx(size, mtxN->GetRealData(), hwMatrix::REAL);
        mtx /= 180.0;
        mtx -= 0.5;

        inputs2.push_back(mtxN);
    }

    Fix(eval, inputs2, outputs2);
    cur1 = inputs2[0];   // cppcheck-suppress invalidContainerReference
    cur2 = outputs2[0];  // cppcheck-suppress invalidContainerReference

    if (cur1.IsScalar())
    {
        double val = cur1.Scalar();

        if (checkisfinite(val) && val == cur2.Scalar())
            outputs[0] = std::numeric_limits<double>::infinity();
    }
    else if (cur1.IsMatrix())
    {
        const hwMatrix* mtx1 = cur1.Matrix();
        hwMatrix* mtx2 = cur2.GetWritableMatrix();
        hwMatrix* mtx3 = outputs[0].GetWritableMatrix();
        int size = mtx1->Size();

        if (mtx1->IsReal())
        {
            for (int i = 0; i < size; ++i)
            {
                if (checkisfinite((*mtx1)(i)) && (*mtx1)(i) == (*mtx2)(i))
                    (*mtx3)(i) = std::numeric_limits<double>::infinity();
            }
        }
        else
        {
            for (int i = 0; i < size; ++i)
            {
                if (mtx1->z(i).IsReal() && checkisfinite(mtx1->z(i).Real()) &&
                    mtx1->z(i) == mtx2->z(i))
                    mtx3->z(i) = 0.0;
            }
        }
    }
    else if (cur1.IsNDMatrix())
    {
        const hwMatrixN* mtx1 = cur1.MatrixN();
        hwMatrixN* mtx2 = cur2.GetWritableMatrixN();
        hwMatrixN* mtx3 = outputs[0].GetWritableMatrixN();
        int size = mtx1->Size();

        if (mtx1->IsReal())
        {
            for (int i = 0; i < size; ++i)
            {
                if (checkisfinite((*mtx1)(i)) && (*mtx1)(i) == (*mtx2)(i))
                    (*mtx3)(i) = std::numeric_limits<double>::infinity();
            }
        }
        else
        {
            for (int i = 0; i < size; ++i)
            {
                if (mtx1->z(i).IsReal() && checkisfinite(mtx1->z(i).Real()) &&
                    mtx1->z(i) == mtx2->z(i))
                    mtx3->z(i) = 0.0;
            }
        }
    }

    return true;
}

//------------------------------------------------------------------------------
// Returns secant of input in degrees [secd]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::SecD(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    std::vector<Currency> outputs2;
    CosD(eval, inputs, outputs2);
    outputs.resize(1);
    Reciprocal(outputs2[0], outputs[0]);

    return true;
}

//------------------------------------------------------------------------------
// Returns cosecant of input in degrees [cscd]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::CscD(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    std::vector<Currency> outputs2;
    SinD(eval, inputs, outputs2);
    outputs.resize(1);
    Reciprocal(outputs2[0], outputs[0]);

    return true;
}

//------------------------------------------------------------------------------
// Returns cotangent of input in radians [cotd]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::CotD(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    std::vector<Currency> outputs2;
    TanD(eval, inputs, outputs2);
    outputs.resize(1);
    Reciprocal(outputs2[0], outputs[0]);

    return true;
}

//------------------------------------------------------------------------------
// Returns inverse cosine of input in radians [acos]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aCos(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    return ElemFunc2<acos, hwComplex::acos_c, hwComplex::acos, vdAcos, vzAcos, absAtMostOne>
        (eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns inverse sine of input in radians [asin]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aSin(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    return ElemFunc2<asin, hwComplex::asin_c, hwComplex::asin, vdAsin, vzAsin, absAtMostOne>
        (eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns inverse tangent of input in radians [atan]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aTan(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    return ElemFunc1<atan, hwComplex::atan, vdAtan, vzAtan>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns inverse secant of input in radians [asec]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aSec(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::vector<Currency> outputs2(1);
    Reciprocal(inputs[0], outputs2[0]);

    return aCos(eval, outputs2, outputs);
}

//------------------------------------------------------------------------------
// Returns inverse cosecant of input in radians [acsc]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aCsc(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::vector<Currency> outputs2(1);
    Reciprocal(inputs[0], outputs2[0]);

    return aSin(eval, outputs2, outputs);
}

//------------------------------------------------------------------------------
// Returns inverse cotangent of input in radians [acot]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aCot(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::vector<Currency> outputs2(1);
    Reciprocal(inputs[0], outputs2[0]);

    return aTan(eval, outputs2, outputs);
}

//------------------------------------------------------------------------------
// Returns inverse cosine of input in degrees [acosd]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aCosD(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    std::vector<Currency> outputs2;

    bool retv = ElemFunc2<acos, hwComplex::acos_c, hwComplex::acos, vdAcospi, vzAcos, absAtMostOne>
        (eval, inputs, outputs2);

    Currency& cur = outputs2[0];

    if (cur.IsScalar())
    {
        outputs.push_back(cur.Scalar() * 180.0);
    }
    else if (cur.IsComplex())
    {
        retv = oml_rad2deg(eval, outputs2, outputs);
    }
    else if (cur.IsMatrix())
    {
        hwMatrix* mtx = cur.GetWritableMatrix();

        if (mtx->IsReal())
        {
            mtx->MultEquals(180.0);
            outputs.push_back(cur);
        }
        else
        {
            retv = oml_rad2deg(eval, outputs2, outputs);
        }
    }
    else if (cur.IsNDMatrix())
    {
        hwMatrixN* mtx = cur.GetWritableMatrixN();

        if (mtx->IsReal())
        {
            hwMatrix temp(mtx->Size(), (void*)mtx->GetRealData(), hwMatrix::REAL);
            temp.MultEquals(180.0);
            outputs.push_back(cur);
        }
        else
        {
            retv = oml_rad2deg(eval, outputs2, outputs);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, 1, OML_VAR_TYPE);
    }

    return retv;
}

//------------------------------------------------------------------------------
// Returns inverse sine of input in degrees [asind]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aSinD(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    std::vector<Currency> outputs2;

    bool retv = ElemFunc2<asin, hwComplex::asin_c, hwComplex::asin, vdAsinpi, vzAsin, absAtMostOne>
        (eval, inputs, outputs2);

    Currency& cur = outputs2[0];

    if (cur.IsScalar())
    {
        outputs.push_back(cur.Scalar() * 180.0);
    }
    else if (cur.IsComplex())
    {
        retv = oml_rad2deg(eval, outputs2, outputs);
    }
    else if (cur.IsMatrix())
    {
        hwMatrix* mtx = cur.GetWritableMatrix();

        if (mtx->IsReal())
        {
            mtx->MultEquals(180.0);
            outputs.push_back(cur);
        }
        else
        {
            retv = oml_rad2deg(eval, outputs2, outputs);
        }
    }
    else if (cur.IsNDMatrix())
    {
        hwMatrixN* mtx = cur.GetWritableMatrixN();

        if (mtx->IsReal())
        {
            hwMatrix temp(mtx->Size(), (void*)mtx->GetRealData(), hwMatrix::REAL);
            temp.MultEquals(180.0);
            outputs.push_back(cur);
        }
        else
        {
            retv = oml_rad2deg(eval, outputs2, outputs);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, 1, OML_VAR_TYPE);
    }

    return retv;
}

//------------------------------------------------------------------------------
// Returns inverse tangent of input in degrees [atand]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aTanD(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    std::vector<Currency> outputs2;

    bool retv = ElemFunc1<atan, hwComplex::atan, vdAtanpi, vzAtan>
        (eval, inputs, outputs2);

    Currency& cur = outputs2[0];

    if (cur.IsScalar())
    {
        outputs.push_back(cur.Scalar() * 180.0);
    }
    else if (cur.IsComplex())
    {
        retv = oml_rad2deg(eval, outputs2, outputs);
    }
    else if (cur.IsMatrix())
    {
        hwMatrix* mtx = cur.GetWritableMatrix();

        if (mtx->IsReal())
        {
            mtx->MultEquals(180.0);
            outputs.push_back(cur);
        }
        else
        {
            retv = oml_rad2deg(eval, outputs2, outputs);
        }
    }
    else if (cur.IsNDMatrix())
    {
        hwMatrixN* mtx = cur.GetWritableMatrixN();

        if (mtx->IsReal())
        {
            hwMatrix temp(mtx->Size(), (void*)mtx->GetRealData(), hwMatrix::REAL);
            temp.MultEquals(180.0);
            outputs.push_back(cur);
        }
        else
        {
            retv = oml_rad2deg(eval, outputs2, outputs);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, 1, OML_VAR_TYPE);
    }

    return retv;
}

//------------------------------------------------------------------------------
// Returns inverse secant of input in degrees [asecd]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aSecD(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::vector<Currency> outputs2(1);
    Reciprocal(inputs[0], outputs2[0]);

    return aCosD(eval, outputs2, outputs);
}

//------------------------------------------------------------------------------
// Returns inverse cosecant of input in degrees [acscd]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aCscD(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::vector<Currency> outputs2(1);
    Reciprocal(inputs[0], outputs2[0]);

    return aSinD(eval, outputs2, outputs);
}

//------------------------------------------------------------------------------
// Returns inverse cotangent of input in degrees [acotd]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aCotD(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::vector<Currency> outputs2(1);
    Reciprocal(inputs[0], outputs2[0]);

    return aTanD(eval, outputs2, outputs);
}

//------------------------------------------------------------------------------
// Returns inverse tangent of input in radians [atan2]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aTan2(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];
    const Currency& input2 = inputs[1];

    if (input1.IsLogical() || input2.IsLogical())
        throw OML_Error(HW_ERROR_INPUTISLOGICAL);

    if (input1.IsScalar())
    {
        if (input2.IsScalar())
        {
            outputs.push_back(atan2(input1.Scalar(), input2.Scalar()));
            return true;
        }
        else if (input2.IsMatrix())
        {
            hwMatrix* m2 = const_cast<hwMatrix*> (input2.Matrix());
            hwMatrix* temp = new hwMatrix(m2->M(), m2->N(), hwMatrix::REAL);
            temp->SetElements(input1.Scalar());

            std::vector<Currency> inputs2;
            inputs2.push_back(temp);
            m2->IncrRefCount();
            inputs2.push_back(m2);

            return aTan2(eval, inputs2, outputs);
        }
        else if (input2.IsNDMatrix())
        {
            hwMatrixN* m2 = const_cast<hwMatrixN*> (input2.MatrixN());
            hwMatrixN* temp = new hwMatrixN(m2->Dimensions(), hwMatrixN::REAL);
            temp->SetElements(input1.Scalar());

            std::vector<Currency> inputs2;
            inputs2.push_back(temp);
            m2->IncrRefCount();
            inputs2.push_back(m2);

            return aTan2(eval, inputs2, outputs);
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
        }
    }
    else if (input1.IsMatrix())
    {
        const hwMatrix* m1 = input1.Matrix();

        if (!m1->IsReal())
            throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

        if (input2.IsScalar())
        {
            hwMatrix* temp = new hwMatrix(m1->M(), m1->N(), hwMatrix::REAL);
            temp->SetElements(input2.Scalar());

            std::vector<Currency> inputs2;
            const_cast <hwMatrix*> (m1)->IncrRefCount();
            inputs2.push_back(const_cast<hwMatrix*> (m1));
            inputs2.push_back(temp);

            return aTan2(eval, inputs2, outputs);
        }
        else if (input2.IsMatrix())
        {
            const hwMatrix* m2 = input2.Matrix();

            if (!m2->IsReal())
                throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

            if (!sameSize(m1, m2))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            int size = m1->Size();
            hwMatrix* result = new hwMatrix(m1->M(), m2->N(), hwMatrix::REAL);
            vdAtan2(size, m1->GetRealData(), m2->GetRealData(), result->GetRealData());
            outputs.push_back(result);
        }
        else if (input2.IsNDMatrix())
        {
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
        }
    }
    else if (input1.IsNDMatrix())
    {
        const hwMatrixN* m1 = input1.MatrixN();

        if (!m1->IsReal())
            throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

        if (input2.IsScalar())
        {
            hwMatrixN* temp = new hwMatrixN(m1->Dimensions(), hwMatrixN::REAL);
            temp->SetElements(input2.Scalar());

            std::vector<Currency> inputs2;
            const_cast <hwMatrixN*> (m1)->IncrRefCount();
            inputs2.push_back(const_cast<hwMatrixN*> (m1));
            inputs2.push_back(temp);

            return aTan2(eval, inputs2, outputs);
        }
        else if (input2.IsMatrix())
        {
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
        }
        else if (input2.IsNDMatrix())
        {
            const hwMatrixN* m2 = input2.MatrixN();

            if (!m2->IsReal())
                throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

            if (m1->Dimensions() != m2->Dimensions())
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            int size = m1->Size();
            hwMatrixN* result = new hwMatrixN(m1->Dimensions(), hwMatrixN::REAL);
            vdAtan2(size, m1->GetRealData(), m2->GetRealData(), result->GetRealData());
            outputs.push_back(result);
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARMATRIX, 1, OML_VAR_TYPE);
    }

    return true;
}

//------------------------------------------------------------------------------
// Returns inverse tangent of input in degrees [atan2d]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aTan2D(EvaluatorInterface           eval,
                             const std::vector<Currency>& inputs,
                             std::vector<Currency>&       outputs)
{
    std::vector<Currency> outputs2;

    bool retv = aTan2(eval, inputs, outputs2);

    Currency& cur = outputs2[0];

    if (cur.IsScalar())
    {
        outputs.push_back(cur.Scalar() * (180.0 / PI));
    }
    else if (cur.IsComplex())
    {
        retv = oml_rad2deg(eval, outputs2, outputs);
    }
    else if (cur.IsMatrix())
    {
        hwMatrix* mtx = cur.GetWritableMatrix();

        if (mtx->IsReal())
        {
            mtx->MultEquals(180.0 / PI);
            outputs.push_back(cur);
        }
        else
        {
            retv = oml_rad2deg(eval, outputs2, outputs);
        }
    }
    else if (cur.IsNDMatrix())
    {
        hwMatrixN* mtx = cur.GetWritableMatrixN();

        if (mtx->IsReal())
        {
            hwMatrix temp(mtx->Size(), (void*)mtx->GetRealData(), hwMatrix::REAL);
            temp.MultEquals(180.0 / PI);
            outputs.push_back(cur);
        }
        else
        {
            retv = oml_rad2deg(eval, outputs2, outputs);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, 1, OML_VAR_VALUE);
    }

    return retv;
}

//------------------------------------------------------------------------------
// Returns hyperbolic cosine of input [cosh]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Cosh(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    return ElemFunc1<cosh, hwComplex::cosh, vdCosh, vzCosh>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns hyperbolic sine of input [sinh]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Sinh(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    return ElemFunc1<sinh, hwComplex::sinh, vdSinh, vzSinh>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns hyperbolic tangent of input [tanh]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Tanh(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    return ElemFunc1<tanh, hwComplex::tanh, vdTanh, vzTanh>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns inverse hyperbolic sine of input [asinh]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aSinh(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    return ElemFunc1<asinh, hwComplex::asinh, vdAsinh, vzAsinh>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns exponential of input [exp]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Exp(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    return ElemFunc1<exp, hwComplex::exp, vdExp, vzExp>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns inverse hyperbolic cosine of input [acosh]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aCosh(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    return ElemFunc2<acosh, hwComplex::acosh_c, hwComplex::acosh, vdAcosh, vzAcosh, atLeastOne>
        (eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns inverse hyperbolic tangent of input [atanh]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::aTanh(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    return ElemFunc2<atanh, hwComplex::atanh_c, hwComplex::atanh, vdAtanh, vzAtanh, absAtMostOne>
        (eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns natural logarithm of input [log]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Log(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    return ElemFunc2<log, hwComplex::log_c, hwComplex::log, vdLn, vzLn, nonnegative>
        (eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns logarithm base 2 of input [log2]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Log2(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    if (eval.GetNargoutValue() > 1)
    {
        return oml_log2(eval, inputs, outputs);
    }
        
    // MKL does not support complex for log2, so use ln(x)/ln(2.0) where needed
    bool ret = ElemFunc2<log2, hwComplex::log_c, hwComplex::log, vdLog2, vzLn, nonnegative>
        (eval, inputs, outputs);

    Currency& cur = outputs[0];

    if (cur.IsComplex())
    {
        double factor = log(2.0);
        hwComplex result = cur.Complex();
        outputs[0] = result / factor;
    }
    else if (cur.IsMatrix())
    {
        hwMatrix* mtx = cur.GetWritableMatrix();

        if (!mtx->IsReal())
        {
            double factor = log(2.0);
            mtx->DivideEquals(factor);
        }
    }
    else if (cur.IsNDMatrix())
    {
        hwMatrixN* mtx = cur.GetWritableMatrixN();

        if (!mtx->IsReal())
        {
            double factor = log(2.0);
            hwComplex* cmplx = mtx->GetComplexData();
            hwMatrix temp(mtx->Size(), static_cast<void*> (cmplx), hwMatrix::COMPLEX);
            temp.DivideEquals(factor);
        }
    }

    return ret;
}

//------------------------------------------------------------------------------
// Returns logarithm bae 10 of input [log10]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Log10(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    return ElemFunc2<log10, hwComplex::log10_c, hwComplex::log10, vdLog10, vzLog10, nonnegative>
        (eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns square root of input [sqrt]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Sqrt(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    return ElemFunc2<sqrt, hwComplex::sqrt_c, hwComplex::sqrt, vdSqrt, vzSqrt, nonnegative>
         (eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns magnitude of input [abs]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Abs(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    return ElemFunc3<fabs, &hwComplex::Mag, vdAbs, vzAbs>
        (eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Returns angle of input [arg]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Arg(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& cur = inputs[0];

    if (cur.IsScalar())
    {
        if (cur.Scalar() >= 0.0)
            outputs.push_back(0.0);
        else if (cur.Scalar() < 0.0)
            outputs.push_back(PI);
        else
            outputs.push_back(std::numeric_limits<double>::quiet_NaN());

        return true;
    }
    else if (cur.IsMatrix() || cur.IsString())
    {
        const hwMatrix* mtx = cur.Matrix();

        if (mtx->IsReal())
        {
            hwMatrix* mtx2 = new hwMatrix;
            mtx2->PackComplex(*mtx);
            std::vector<Currency> inputs2;
            inputs2.push_back(mtx2);

            return ElemFunc3<fabs, &hwComplex::Arg, vdAbs, vzArg>
                (eval, inputs2, outputs);
        }
    }
    else if (cur.IsNDMatrix())
    {
        const hwMatrixN* mtx = cur.MatrixN();

        if (mtx->IsReal())
        {
            hwMatrixN* mtx2 = new hwMatrixN;
            mtx2->PackComplex(*mtx);
            std::vector<Currency> inputs2;
            inputs2.push_back(mtx2);

            return ElemFunc3<fabs, &hwComplex::Arg, vdAbs, vzArg>
                (eval, inputs2, outputs);
        }
    }

    return ElemFunc3<fabs, &hwComplex::Arg, vdAbs, vzArg>
        (eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Rounds input to the nearest integer [round]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Round(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    return ElemFunc4<round, vdRound>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Rounds input toward inf [ceil]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Ceil(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>& outputs)
{
    return ElemFunc4<ceil, vdCeil>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Rounds input toward -inf [floor]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Floor(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>& outputs)
{
    return ElemFunc4<floor, vdFloor>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Rounds input toward zero [fix]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Fix(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    return ElemFunc4<trunc, vdTrunc>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Computes division remainder [rem]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Rem(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    return ElemFunc5<rem, vdFmod>(eval, inputs, outputs);
}

//------------------------------------------------------------------------------
// Computes modulo function [mod]
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Mod(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];
    const Currency& input2 = inputs[1];

    if (input1.IsLogical() || input2.IsLogical())
        throw OML_Error(HW_ERROR_INPUTISLOGICAL);

    if (input1.IsScalar())
    {
        if (input2.IsScalar())
        {
            outputs.push_back(mod(input1.Scalar(), input2.Scalar()));
            return true;
        }
        else if (input2.IsMatrix())
        {
            hwMatrix* m2 = const_cast<hwMatrix*> (input2.Matrix());
            hwMatrix* temp = new hwMatrix(m2->M(), m2->N(), hwMatrix::REAL);
            temp->SetElements(input1.Scalar());

            std::vector<Currency> inputs2;
            inputs2.push_back(temp);
            m2->IncrRefCount();
            inputs2.push_back(m2);

            return Mod(eval, inputs2, outputs);
        }
        else if (input2.IsNDMatrix())
        {
            hwMatrixN* m2 = const_cast<hwMatrixN*> (input2.MatrixN());
            hwMatrixN* temp = new hwMatrixN(m2->Dimensions(), hwMatrixN::REAL);
            temp->SetElements(input1.Scalar());

            std::vector<Currency> inputs2;
            inputs2.push_back(temp);
            m2->IncrRefCount();
            inputs2.push_back(m2);

            return Mod(eval, inputs2, outputs);
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
        }
    }
    else if (input1.IsMatrix())
    {
        const hwMatrix* m1 = input1.Matrix();

        if (!m1->IsReal())
            throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

        if (input2.IsScalar())
        {
            hwMatrix* temp = new hwMatrix(m1->M(), m1->N(), hwMatrix::REAL);
            temp->SetElements(input2.Scalar());

            std::vector<Currency> inputs2;
            const_cast <hwMatrix*> (m1)->IncrRefCount();
            inputs2.push_back(const_cast<hwMatrix*> (m1));
            inputs2.push_back(temp);

            return Mod(eval, inputs2, outputs);
        }
        else if (input2.IsMatrix())
        {
            const hwMatrix* m2 = input2.Matrix();

            if (!m2->IsReal())
                throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

            if (!sameSize(m1, m2))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            hwMatrix* result = new hwMatrix(m1->M(), m2->N(), hwMatrix::REAL);
            MKLutilitiesD::Mod(*m1, *m2, *result);
            outputs.push_back(result);
        }
        else if (input2.IsNDMatrix())
        {
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
        }
    }
    else if (input1.IsNDMatrix())
    {
        const hwMatrixN* m1 = input1.MatrixN();

        if (!m1->IsReal())
            throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

        if (input2.IsScalar())
        {
            hwMatrixN* temp = new hwMatrixN(m1->Dimensions(), hwMatrixN::REAL);
            temp->SetElements(input2.Scalar());

            std::vector<Currency> inputs2;
            const_cast <hwMatrixN*> (m1)->IncrRefCount();
            inputs2.push_back(const_cast<hwMatrixN*> (m1));
            inputs2.push_back(temp);

            return Mod(eval, inputs2, outputs);
        }
        else if (input2.IsMatrix())
        {
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
        }
        else if (input2.IsNDMatrix())
        {
            const hwMatrixN* m2 = input2.MatrixN();

            if (!m2->IsReal())
                throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

            if (m1->Dimensions() != m2->Dimensions())
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            int size = m1->Size();
            hwMatrixN* result = new hwMatrixN(m1->Dimensions(), hwMatrixN::REAL);
            hwMatrix temp_m1(size, const_cast<double*> (m1->GetRealData()), hwMatrix::REAL);
            hwMatrix temp_m2(size, const_cast<double*> (m2->GetRealData()), hwMatrix::REAL);
            hwMatrix temp_r(size, result->GetRealData(), hwMatrix::REAL);

            MKLutilitiesD::Mod(temp_m1, temp_m2, temp_r);
            outputs.push_back(result);
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARMATRIX, 1, OML_VAR_TYPE);
    }

    return true;
}

//------------------------------------------------------------------------------
// Returns element-wise maximum of two inputs
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Max(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
    {
        return oml_max(eval, inputs, outputs);
    }

    try
    {
        return ElemFunc5<_max, vdFmax>(eval, inputs, outputs);
    }
    catch (OML_Error& except)
    {
        if (except.ErrCode() == OML_ERR_REAL || except.ErrCode() == OML_ERR_SCALARMATRIX)
        {
            return oml_max(eval, inputs, outputs);
        }
        else
        {
            throw except;
        }
    }
}

//------------------------------------------------------------------------------
// Returns element-wise minimum of two inputs
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::Min(EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
    {
        return oml_min(eval, inputs, outputs);
    }

    try
    {
        return ElemFunc5<_min, vdFmin>(eval, inputs, outputs);
    }
    catch (OML_Error & except)
    {
        if (except.ErrCode() == OML_ERR_REAL || except.ErrCode() == OML_ERR_SCALARMATRIX)
        {
            return oml_min(eval, inputs, outputs);
        }
        else
        {
            throw except;
        }
    }
}

// Sparse matrix helper functions
//------------------------------------------------------------------------------
// Set pivot threshold for MKL sparse matrix division
//------------------------------------------------------------------------------
bool BuiltInFuncsMKL::oml_mkl_sdpt(EvaluatorInterface           eval,
                                   const std::vector<Currency>& inputs,
                                   std::vector<Currency>& outputs)
{
    if (inputs.size() == 0)
    {
        outputs.push_back(iparm9);
    }
    else if (inputs.size() == 1)
    {
        if (!inputs[0].IsPositiveInteger())
            throw OML_Error(OML_ERR_NUMARGIN);

        MKL_INT val = static_cast<MKL_INT> (inputs[0].Scalar());

        if (val < 13)
            return false;

        MKLutilitiesD::SetIparam9(val);
    }
    else
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    return true;
}
