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
#include "mkl_vml_functions.h"

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
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());
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
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);
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
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());
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
            hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());
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
            hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());
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
            hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());
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
    cur2 = outputs2[0];

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
    return ElemFunc3<abs, &hwComplex::Mag, vdAbs, vzAbs>
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

            return ElemFunc3<abs, &hwComplex::Arg, vdAbs, vzArg>
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

            return ElemFunc3<abs, &hwComplex::Arg, vdAbs, vzArg>
                (eval, inputs2, outputs);
        }
    }

    return ElemFunc3<abs, &hwComplex::Arg, vdAbs, vzArg>
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
            // vdRemainder(size, m1->GetRealData(), m2->GetRealData(), result->GetRealData());
            // vdRemainder does not follow our convention
            hwMatrix* temp1 = new hwMatrix;
            hwMatrix temp2;
            temp1->DivideByElems(*m1, *m2);
            std::vector<Currency> inputs2;
            std::vector<Currency> outputs2;
            inputs2.push_back(temp1);
            Floor(eval, inputs2, outputs2);
            temp2.MultByElems(*m2, *outputs2[0].GetWritableMatrix());
            result->Subtr(*m1, temp2);
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
            hwMatrix* temp3 = new hwMatrix;
            hwMatrix temp4;
            temp3->DivideByElems(temp_m1, temp_m2);
            std::vector<Currency> inputs2;
            std::vector<Currency> outputs2;
            inputs2.push_back(temp3);
            Floor(eval, inputs2, outputs2);
            temp4.MultByElems(temp_m2, *outputs2[0].GetWritableMatrix());
            temp_r.Subtr(temp_m1, temp4);
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

//------------------------------------------------------------------------------
// Conjugate of a complex matrix
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::Conj(const hwMatrix& A,
                           hwMatrix&       conj)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (!A.IsReal())
    {
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
        hwMathStatus status = conj.Dimension(m, n, hwMatrix::COMPLEX);
        hwComplex* pC = conj.GetComplexData();
        MKL_Complex16* pMKL_C = reinterpret_cast<MKL_Complex16*> (pC);
        vzConj(sizeA, pMKL_A, pMKL_C);
    }
    else
    {
        conj = A;
    }
}

//------------------------------------------------------------------------------
// Hypotenuse of paired matrix elements (hyp = sqrt(A^2 + B^2)
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::Hypot(const hwMatrix& A,
                            const hwMatrix& B,
                            hwMatrix&       hyp)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (A.IsReal() && B.IsReal())
    {
        const double* pA = A.GetRealData();
        const double* pB = B.GetRealData();
        hwMathStatus status = hyp.Dimension(m, n, hwMatrix::REAL);
        double* pC = hyp.GetRealData();
        vdHypot(sizeA, pA, pB, pC);
    }
}

//------------------------------------------------------------------------------
// Add two full matrices, sum = A + B
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::Add(const hwMatrix& A,
                          const hwMatrix& B,
                          hwMatrix&       sum)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (B.M() != m || B.N() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.IsReal() && B.IsReal())
    {
        const double* pA = A.GetRealData();
        const double* pB = B.GetRealData();
        hwMathStatus status = sum.Dimension(m, n, hwMatrix::REAL);
        double* pC = sum.GetRealData();
        vdAdd(sizeA, pA, pB, pC);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pB = const_cast<hwComplex*> (B.GetComplexData());
        hwMathStatus status = sum.Dimension(m, n, hwMatrix::COMPLEX);
        hwComplex* pC = sum.GetComplexData();
        MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
        MKL_Complex16* pMKL_B = reinterpret_cast<MKL_Complex16*> (pB);
        MKL_Complex16* pMKL_C = reinterpret_cast<MKL_Complex16*> (pC);
        vzAdd(sizeA, pMKL_A, pMKL_B, pMKL_C);
    }
    else if (!A.IsReal())
    {
        hwMatrix BC;
        BC.PackComplex(B);
        return Add(A, BC, sum);
    }
    else    // !B.IsReal()
    {
        hwMatrix AC;
        AC.PackComplex(A);
        return Add(AC, B, sum);
    }
}

//------------------------------------------------------------------------------
// Subtract two full matrices, diff = A - B
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::Subtr(const hwMatrix& A,
                            const hwMatrix& B,
                            hwMatrix&       diff)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (B.M() != m || B.N() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.IsReal() && B.IsReal())
    {
        const double* pA = A.GetRealData();
        const double* pB = B.GetRealData();
        hwMathStatus status = diff.Dimension(m, n, hwMatrix::REAL);
        double* pC = diff.GetRealData();
        vdSub(sizeA, pA, pB, pC);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pB = const_cast<hwComplex*> (B.GetComplexData());
        hwMathStatus status = diff.Dimension(m, n, hwMatrix::COMPLEX);
        hwComplex* pC = diff.GetComplexData();
        MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
        MKL_Complex16* pMKL_B = reinterpret_cast<MKL_Complex16*> (pB);
        MKL_Complex16* pMKL_C = reinterpret_cast<MKL_Complex16*> (pC);
        vzSub(sizeA, pMKL_A, pMKL_B, pMKL_C);
    }
    else if (!A.IsReal())
    {
        hwMatrix BC;
        BC.PackComplex(B);
        return Subtr(A, BC, diff);
    }
    else    // !B.IsReal()
    {
        hwMatrix AC;
        AC.PackComplex(A);
        return Subtr(AC, B, diff);
    }
}

//------------------------------------------------------------------------------
// Multiply two full matrices element-wise, prod = A .* B
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::MultByElems(const hwMatrix& A,
                                  const hwMatrix& B,
                                  hwMatrix&       prod)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (B.M() != m || B.N() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.IsReal() && B.IsReal())
    {
        const double* pA = A.GetRealData();
        const double* pB = B.GetRealData();
        hwMathStatus status = prod.Dimension(m, n, hwMatrix::REAL);
        double* pC = prod.GetRealData();
        vdMul(sizeA, pA, pB, pC);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pB = const_cast<hwComplex*> (B.GetComplexData());
        hwMathStatus status = prod.Dimension(m, n, hwMatrix::COMPLEX);
        hwComplex* pC = prod.GetComplexData();
        MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
        MKL_Complex16* pMKL_B = reinterpret_cast<MKL_Complex16*> (pB);
        MKL_Complex16* pMKL_C = reinterpret_cast<MKL_Complex16*> (pC);
        vzMul(sizeA, pMKL_A, pMKL_B, pMKL_C);
    }
    else if (!A.IsReal())
    {
        hwMatrix BC;
        BC.PackComplex(B);
        return MultByElems(A, BC, prod);
    }
    else    // !B.IsReal()
    {
        hwMatrix AC;
        AC.PackComplex(A);
        return MultByElems(AC, B, prod);
    }
}

//------------------------------------------------------------------------------
// Divide two full matrices element-wise, quot = A .* B
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::DivideByElems(const hwMatrix& A,
                                    const hwMatrix& B,
                                    hwMatrix& quot)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (B.M() != m || B.N() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.IsReal() && B.IsReal())
    {
        const double* pA = A.GetRealData();
        const double* pB = B.GetRealData();
        hwMathStatus status = quot.Dimension(m, n, hwMatrix::REAL);
        double* pC = quot.GetRealData();
        vdDiv(sizeA, pA, pB, pC);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pB = const_cast<hwComplex*> (B.GetComplexData());
        hwMathStatus status = quot.Dimension(m, n, hwMatrix::COMPLEX);
        hwComplex* pC = quot.GetComplexData();
        MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
        MKL_Complex16* pMKL_B = reinterpret_cast<MKL_Complex16*> (pB);
        MKL_Complex16* pMKL_C = reinterpret_cast<MKL_Complex16*> (pC);
        vzDiv(sizeA, pMKL_A, pMKL_B, pMKL_C);
    }
    else if (!A.IsReal())
    {
        hwMatrix BC;
        BC.PackComplex(B);
        return DivideByElems(A, BC, quot);
    }
    else    // !B.IsReal()
    {
        hwMatrix AC;
        AC.PackComplex(A);
        return DivideByElems(AC, B, quot);
    }
}

//------------------------------------------------------------------------------
// Power operation applied element-wise (B = A .^ P)
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::PowerByElems(const hwMatrix& A,
                                   const hwMatrix& P,
                                   hwMatrix&       B)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (P.M() != m || P.N() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.IsReal() && P.IsReal())
    {
        int count = sizeA;
        const double* pA = A.GetRealData();
        const double* pP = P.GetRealData();
        bool A_gt_zero = true;

        while (count--)
        {
            if (*pA <= 0.0)
            {
                A_gt_zero = false;
                break;
            }

            ++pA;
        }

        pA = A.GetRealData();

        if (A_gt_zero)
        {
            hwMathStatus status = B.Dimension(m, n, hwMatrix::REAL);
            double* pB = B.GetRealData();
            vdPowr(sizeA, pA, pP, pB);
        }
        else
        {
            int count = sizeA;
            bool P_is_int = true;

            while (count--)
            {
                if (!isint(*pP))
                {
                    P_is_int = false;
                    break;
                }

                ++pP;
            }

            pP = P.GetRealData();

            if (P_is_int)
            {
                pA = A.GetRealData();
                hwMathStatus status = B.Dimension(m, n, hwMatrix::REAL);
                double* pB = B.GetRealData();
                vdPow(sizeA, pA, pP, pB);
            }
            else
            {
                hwMatrix AC;
                AC.PackComplex(A);
                hwMatrix PC;
                PC.PackComplex(P);
                return PowerByElems(AC, PC, B);
            }
        }
    }
    else if (!A.IsReal() && !P.IsReal())
    {
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pP = const_cast<hwComplex*> (P.GetComplexData());
        hwMathStatus status = B.Dimension(m, n, hwMatrix::COMPLEX);
        hwComplex* pB = B.GetComplexData();
        MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
        MKL_Complex16* pMKL_P = reinterpret_cast<MKL_Complex16*> (pP);
        MKL_Complex16* pMKL_B = reinterpret_cast<MKL_Complex16*> (pB);
        vzPow(sizeA, pMKL_A, pMKL_P, pMKL_B);
    }
    else if (!A.IsReal())
    {
        hwMatrix PC;
        PC.PackComplex(P);
        return PowerByElems(A, PC, B);
    }
    else    // !P.IsReal()
    {
        hwMatrix AC;
        AC.PackComplex(A);
        return PowerByElems(AC, P, B);
    }
}

//------------------------------------------------------------------------------
// Power operation applied element-wise (B = A .^ P)
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::PowerByElems(const hwMatrix& A,
                                   double          P,
                                   hwMatrix&       B)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (A.IsReal())
    {
        const double* pA = A.GetRealData();

        if (isint(P))
        {
            hwMathStatus status = B.Dimension(m, n, hwMatrix::REAL);
            double* pB = B.GetRealData();

            if (P == 2.0)
                vdSqr(sizeA, pA, pB);
            else
                vdPowx(sizeA, pA, P, pB);
        }
        else
        {
            int count = sizeA;
            bool A_ge_zero = true;

            while (count--)
            {
                if (*pA < 0.0)
                {
                    A_ge_zero = false;
                    break;
                }

                ++pA;
            }

            pA = A.GetRealData();

            if (A_ge_zero)
            {
                hwMathStatus status = B.Dimension(m, n, hwMatrix::REAL);
                double* pB = B.GetRealData();
                vdPowx(sizeA, pA, P, pB);
            }
            else
            {
                hwMatrix AC;
                AC.PackComplex(A);
                hwComplex PC(P, 0.0);
                return PowerByElems(AC, PC, B);
            }
        }
    }
    else    // !A.IsReal()
    {
        hwComplex PC(P, 0.0);
        return PowerByElems(A, PC, B);
    }
}

//------------------------------------------------------------------------------
// Power operation applied element-wise (B = A .^ P)
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::PowerByElems(const hwMatrix&  A,
                                   const hwComplex& P,
                                   hwMatrix&        B)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (A.IsReal())
    {
        hwMatrix AC;
        AC.PackComplex(A);
        return PowerByElems(AC, P, B);
    }
    else
    {
        if (P == 2.0)
        {
            MultByElems(A, A, B);
        }
        else
        {
            hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
            hwMathStatus status = B.Dimension(m, n, hwMatrix::COMPLEX);
            hwComplex* pB = B.GetComplexData();
            MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
            const MKL_Complex16* pMKL_P = reinterpret_cast<const MKL_Complex16*> (&P);
            MKL_Complex16* pMKL_B = reinterpret_cast<MKL_Complex16*> (pB);
            vzPowx(sizeA, pMKL_A, *pMKL_P, pMKL_B);
        }
    }
}

// Sparse matrix helper functions
#include "mkl_spblas.h"
#include "mkl_pardiso.h"
//------------------------------------------------------------------------------
// Add two sparse matrices, sum = A + B
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::SparseAdd(const hwMatrixS& A,
                                const hwMatrixS& B,
                                hwMatrixS&       sum)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    hwMathStatus status;

    if (B.M() != m || B.N() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.NNZ() == 0)
    {
        sum = B;
        return;
    }
    else if (B.NNZ() == 0)
    {
        sum = A;
        return;
    }

    MKL_INT* prA = const_cast<MKL_INT*> (A.rows());
    MKL_INT* pbA = const_cast<MKL_INT*> (A.pointerB());
    MKL_INT* peA = const_cast<MKL_INT*> (A.pointerE());
    MKL_INT* prB = const_cast<MKL_INT*> (B.rows());
    MKL_INT* pbB = const_cast<MKL_INT*> (B.pointerB());
    MKL_INT* peB = const_cast<MKL_INT*> (B.pointerE());

    MKL_INT* prC = nullptr;
    MKL_INT* pbC = nullptr;
    MKL_INT* peC = nullptr;

    // MKL matrix and description
    sparse_matrix_t A_MKL = NULL;
    sparse_matrix_t B_MKL = NULL;
    sparse_matrix_t C_MKL = NULL;
    sparse_status_t mkl_status;

    if (A.IsReal() && B.IsReal())
    {
        double* pVA = const_cast<double*> (A.GetRealData());
        double* pVB = const_cast<double*> (B.GetRealData());
        double* pVC = nullptr;

        mkl_status = mkl_sparse_d_create_csr(&A_MKL, SPARSE_INDEX_BASE_ZERO,
            n, m, pbA, peA, prA, pVA);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_d_create_csr(&B_MKL, SPARSE_INDEX_BASE_ZERO,
            n, m, pbB, peB, prB, pVB);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, A_MKL, 1.0, B_MKL, &C_MKL);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

        mkl_status = mkl_sparse_d_export_csr(C_MKL, &indexing, &n, &m, &pbC, &peC, &prC, &pVC);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        sum = hwMatrixS(m, n, pbC, peC, prC, pVC);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwComplex* pVA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pVB = const_cast<hwComplex*> (B.GetComplexData());
        hwComplex* pVC = nullptr;
        struct _MKL_Complex16 one;

        one.real = 1.0;
        one.imag = 0.0;

        mkl_status = mkl_sparse_z_create_csr(&A_MKL, SPARSE_INDEX_BASE_ZERO,
                                             n, m, pbA, peA, prA,
                                             reinterpret_cast<MKL_Complex16*> (pVA));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_z_create_csr(&B_MKL, SPARSE_INDEX_BASE_ZERO,
                                             n, m, pbB, peB, prB,
                                             reinterpret_cast<MKL_Complex16*> (pVB));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, A_MKL, one, B_MKL, &C_MKL);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

        mkl_status = mkl_sparse_z_export_csr(C_MKL, &indexing, &n, &m, &pbC, &peC, &prC,
                                             reinterpret_cast<MKL_Complex16**> (&pVC));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        sum = hwMatrixS(m, n, pbC, peC, prC, pVC);
    }
    else if (!A.IsReal())
    {
        hwMatrixS BC;
        BC.PackComplex(B);
        return SparseAdd(A, BC, sum);
    }
    else    // !B.IsReal()
    {
        hwMatrixS AC;
        AC.PackComplex(A);
        return SparseAdd(AC, B, sum);
    }

    mkl_sparse_destroy(A_MKL);
    mkl_sparse_destroy(B_MKL);
    mkl_sparse_destroy(C_MKL);
}
//------------------------------------------------------------------------------
// Multiply a sparse matrix by a full matrix, prod = A * B
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::SparseMult(const hwMatrixS& A,
                                 const hwMatrix&  B,
                                 hwMatrix&        prod)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    MKL_INT k = B.N();
    hwMathStatus status;

    if (B.M() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.NNZ() == 0 || B.IsEmpty())
    {
        status = prod.Dimension(m, k, hwMatrix::REAL);
        prod.SetElements(0.0);
        return;
    }

    MKL_INT* pr = const_cast<MKL_INT*> (A.rows());
    MKL_INT* pb = const_cast<MKL_INT*> (A.pointerB());
    MKL_INT* pe = const_cast<MKL_INT*> (A.pointerE());

    // MKL matrix and description
    sparse_matrix_t A_MKL = NULL;
    struct matrix_descr DSC;
    DSC.type = SPARSE_MATRIX_TYPE_GENERAL;
    sparse_status_t mkl_status;

    if (A.IsReal() && B.IsReal())
    {
        status = prod.Dimension(m, k, hwMatrix::REAL);

        // perform MKL operations
        double* pV = const_cast<double*> (A.GetRealData());

        if (k == 1)
        {
            // populate MKL matrix
            mkl_status = mkl_sparse_d_create_csc(&A_MKL, SPARSE_INDEX_BASE_ZERO,
                                                 m, n, pb, pe, pr, pV);

            if (mkl_status != SPARSE_STATUS_SUCCESS)
                throw hwMathException(HW_MATH_ERR_UNKNOWN);

            mkl_status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,
                                         1.0, A_MKL, DSC, B.GetRealData(),
                                         0.0, prod.GetRealData());
        }
        else
        {
            // SPARSE_INDEX_BASE_ZERO is not supported in MKL for this
            // function, so convert to SPARSE_INDEX_BASE_ONE
            int nnz = A.NNZ();
            MKL_INT* pb1 = new MKL_INT[n];
            MKL_INT* pe1 = new MKL_INT[n];
            MKL_INT* pr1 = new MKL_INT[nnz];

            for (int i = 0; i < n; ++i)
            {
                pb1[i] = pb[i] + 1;
                pe1[i] = pe[i] + 1;
            }

            for (int i = 0; i < nnz; ++i)
                pr1[i] = pr[i] + 1;

            // populate MKL matrix
            mkl_status = mkl_sparse_d_create_csc(&A_MKL, SPARSE_INDEX_BASE_ONE,
                                                 m, n, pb1, pe1, pr1, pV);

            if (mkl_status != SPARSE_STATUS_SUCCESS)
                throw hwMathException(HW_MATH_ERR_UNKNOWN);

            mkl_status = mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE,
                                         1.0, A_MKL, DSC, SPARSE_LAYOUT_COLUMN_MAJOR,
                                         B.GetRealData(), k, n, 0.0,
                                         prod.GetRealData(), m);

            delete[] pb1;
            delete[] pe1;
            delete[] pr1;
        }
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        status = prod.Dimension(m, k, hwMatrix::COMPLEX);

        // perform MKL operations
        hwComplex* pV = const_cast<hwComplex*> (A.GetComplexData());
        hwTComplex<double>* pB = const_cast<hwComplex*> (B.GetComplexData());
        MKL_Complex16* pP = reinterpret_cast<MKL_Complex16*> (prod.GetComplexData());
        struct _MKL_Complex16 one;
        struct _MKL_Complex16 zero;

        one.real = 1.0;
        one.imag = 0.0;
        zero.real = 0.0;
        zero.imag = 0.0;

        if (k == 1)
        {
            // populate MKL matrix
            mkl_status = mkl_sparse_z_create_csc(&A_MKL, SPARSE_INDEX_BASE_ZERO,
                                                 m, n, pb, pe, pr,
                                                 reinterpret_cast<MKL_Complex16*> (pV));

            if (mkl_status != SPARSE_STATUS_SUCCESS)
                throw hwMathException(HW_MATH_ERR_UNKNOWN);

            mkl_status = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE,
                                         one, A_MKL, DSC,
                                         reinterpret_cast<MKL_Complex16*> (pB), zero, pP);
        }
        else
        {
            // SPARSE_INDEX_BASE_ZERO is not supported in MKL for this
            // function, so convert to SPARSE_INDEX_BASE_ONE
            int nnz = A.NNZ();
            MKL_INT* pb1 = new MKL_INT[n];
            MKL_INT* pe1 = new MKL_INT[n];
            MKL_INT* pr1 = new MKL_INT[nnz];

            for (int i = 0; i < n; ++i)
            {
                pb1[i] = pb[i] + 1;
                pe1[i] = pe[i] + 1;
            }

            for (int i = 0; i < nnz; ++i)
                pr1[i] = pr[i] + 1;

            // populate MKL matrix
            mkl_status = mkl_sparse_z_create_csc(&A_MKL, SPARSE_INDEX_BASE_ONE,
                                                 m, n, pb1, pe1, pr1,
                                                 reinterpret_cast<MKL_Complex16*> (pV));

            if (mkl_status != SPARSE_STATUS_SUCCESS)
                throw hwMathException(HW_MATH_ERR_UNKNOWN);

            mkl_status = mkl_sparse_z_mm(SPARSE_OPERATION_NON_TRANSPOSE,
                                         one, A_MKL, DSC, SPARSE_LAYOUT_COLUMN_MAJOR,
                                         reinterpret_cast<MKL_Complex16*> (pB), k, n,
                                         zero, pP, m);

            delete[] pb1;
            delete[] pe1;
            delete[] pr1;
        }
    }
    else if (!A.IsReal())
    {
        hwMatrix C;
        C.PackComplex(B);
        return SparseMult(A, C, prod);
    }
    else    // !B.IsReal()
    {
        hwMatrixS C;
        C.PackComplex(A);
        return SparseMult(C, B, prod);
    }

    if (mkl_status != SPARSE_STATUS_SUCCESS)
        throw hwMathException(HW_MATH_ERR_UNKNOWN);

    mkl_sparse_destroy(A_MKL);
}
//------------------------------------------------------------------------------
// Multiply a full matrix by a sparse matrix, prod = A * B
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::SparseMult(const hwMatrix&  A,
                                 const hwMatrixS& B,
                                 hwMatrix&        prod)
{
    // prod = A * B will be computed using trans(prod) = trans(B) * trans(A)
    // while interpreting prod and A as having row major storage

    // A: row major matrix
    MKL_INT rows_A = A.N();
    MKL_INT cols_A = A.M();

    // B: column major sparse matrix
    MKL_INT rows_B = B.M();
    MKL_INT cols_B = B.N();

    if (rows_A != rows_B)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    // prod: row major matrix
    MKL_INT rows_P = cols_B;
    MKL_INT cols_P = cols_A;

    if (B.NNZ() == 0 || A.IsEmpty())
    {
        hwMathStatus status = prod.Dimension(cols_P, rows_P, hwMatrix::REAL);
        prod.SetElements(0.0);
        return;
    }

    // compute prod
    MKL_INT* pr = const_cast<MKL_INT*> (B.rows());
    MKL_INT* pb = const_cast<MKL_INT*> (B.pointerB());
    MKL_INT* pe = const_cast<MKL_INT*> (B.pointerE());
    sparse_matrix_t B_MKL = NULL;
    struct matrix_descr DSC;
    DSC.type = SPARSE_MATRIX_TYPE_GENERAL;
    sparse_status_t mkl_status;

    MKL_INT columns = cols_P;
    MKL_INT ldx = cols_A;
    MKL_INT ldp = cols_P;

    if (A.IsReal() && B.IsReal())
    {
        hwMathStatus status = prod.Dimension(cols_P, rows_P, hwMatrix::REAL);

        // populate MKL matrix and compute
        double* pB = const_cast<double*> (B.GetRealData());

        mkl_status = mkl_sparse_d_create_csc(&B_MKL, SPARSE_INDEX_BASE_ZERO,
                                             rows_B, cols_B, pb, pe, pr, pB);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        double alpha = 1.0;
        double beta = 0.0;

        if (cols_A == 1)
        {
            mkl_status = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE,
                                         alpha, B_MKL, DSC, A.GetRealData(),
                                         beta, prod.GetRealData());
        }
        else
        {
            mkl_status = mkl_sparse_d_mm(SPARSE_OPERATION_TRANSPOSE,
                                         alpha, B_MKL, DSC, SPARSE_LAYOUT_ROW_MAJOR,
                                         A.GetRealData(), columns, ldx, beta,
                                         prod.GetRealData(), ldp);
        }
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwMathStatus status = prod.Dimension(cols_P, rows_P, hwMatrix::COMPLEX);

        // populate MKL matrix and compute
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pB = const_cast<hwComplex*> (B.GetComplexData());
        MKL_Complex16* pP = reinterpret_cast<MKL_Complex16*> (prod.GetComplexData());

        // populate MKL matrix
        mkl_status = mkl_sparse_z_create_csc(&B_MKL, SPARSE_INDEX_BASE_ZERO,
                                             rows_B, cols_B, pb, pe, pr,
                                             reinterpret_cast<MKL_Complex16*> (pB));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        struct _MKL_Complex16 alpha;
        struct _MKL_Complex16 beta;

        alpha.real = 1.0;
        alpha.imag = 0.0;
        beta.real = 0.0;
        beta.imag = 0.0;

        if (cols_A == 1)
        {
            mkl_status = mkl_sparse_z_mv(SPARSE_OPERATION_TRANSPOSE,
                                        alpha, B_MKL, DSC, reinterpret_cast<MKL_Complex16*> (pA), beta, pP);
        }
        else
        {
            mkl_status = mkl_sparse_z_mm(SPARSE_OPERATION_TRANSPOSE,
                                         alpha, B_MKL, DSC, SPARSE_LAYOUT_ROW_MAJOR,
                                         reinterpret_cast<MKL_Complex16*> (pA), columns,
                                         ldx, beta, pP, ldp);
        }
    }
    else if (!A.IsReal())
    {
        hwMatrixS BC;
        BC.PackComplex(B);
        return SparseMult(A, BC, prod);
    }
    else    // !B.IsReal()
    {
        hwMatrix AC;
        AC.PackComplex(A);
        return SparseMult(AC, B, prod);
    }

    if (mkl_status != SPARSE_STATUS_SUCCESS)
        throw hwMathException(HW_MATH_ERR_UNKNOWN);

    mkl_sparse_destroy(B_MKL);
}
//------------------------------------------------------------------------------
// Multiply a sparse matrix by a sparse matrix, prod = A * B
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::SparseMult(const hwMatrixS& A,
                                 const hwMatrixS& B,
                                 hwMatrixS& prod)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    MKL_INT k = B.N();
    hwMathStatus status;

    if (B.M() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.NNZ() == 0 || B.NNZ() == 0)
    {
        std::vector<int> ivec;
        std::vector<int> jvec;
        hwMatrix vals;
        prod = hwMatrixS(ivec, jvec, vals, m, k);
        return;
    }

    MKL_INT* prA = const_cast<MKL_INT*> (A.rows());
    MKL_INT* pbA = const_cast<MKL_INT*> (A.pointerB());
    MKL_INT* peA = const_cast<MKL_INT*> (A.pointerE());
    MKL_INT* prB = const_cast<MKL_INT*> (B.rows());
    MKL_INT* pbB = const_cast<MKL_INT*> (B.pointerB());
    MKL_INT* peB = const_cast<MKL_INT*> (B.pointerE());

    MKL_INT* prC = nullptr;
    MKL_INT* pbC = nullptr;
    MKL_INT* peC = nullptr;

    // MKL matrix and description
    sparse_matrix_t A_MKL = NULL;
    sparse_matrix_t B_MKL = NULL;
    sparse_matrix_t C_MKL = NULL;
    sparse_status_t mkl_status;

    if (A.IsReal() && B.IsReal())
    {
        double* pVA = const_cast<double*> (A.GetRealData());
        double* pVB = const_cast<double*> (B.GetRealData());
        double* pVC = nullptr;

        mkl_status = mkl_sparse_d_create_csc(&A_MKL, SPARSE_INDEX_BASE_ZERO,
                                             m, n, pbA, peA, prA, pVA);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_d_create_csc(&B_MKL, SPARSE_INDEX_BASE_ZERO,
                                             n, k, pbB, peB, prB, pVB);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, A_MKL, B_MKL, &C_MKL);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

        mkl_status = mkl_sparse_d_export_csc(C_MKL, &indexing, &n, &k, &pbC, &peC, &prC, &pVC);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        prod = hwMatrixS(n, k, pbC, peC, prC, pVC);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwComplex* pVA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pVB = const_cast<hwComplex*> (B.GetComplexData());
        hwComplex* pVC = nullptr;

        mkl_status = mkl_sparse_z_create_csc(&A_MKL, SPARSE_INDEX_BASE_ZERO,
                                             m, n, pbA, peA, prA,
                                             reinterpret_cast<MKL_Complex16*> (pVA));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_z_create_csc(&B_MKL, SPARSE_INDEX_BASE_ZERO,
                                             n, k, pbB, peB, prB, reinterpret_cast<MKL_Complex16*> (pVB));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, A_MKL, B_MKL, &C_MKL);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

        mkl_status = mkl_sparse_z_export_csc(C_MKL, &indexing, &n, &k, &pbC, &peC, &prC,
                                             reinterpret_cast<MKL_Complex16**> (&pVC));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        prod = hwMatrixS(n, k, pbC, peC, prC, pVC);
    }
    else if (!A.IsReal())
    {
        hwMatrixS BC;
        BC.PackComplex(B);
        return SparseMult(A, BC, prod);
    }
    else    // !B.IsReal()
    {
        hwMatrixS AC;
        AC.PackComplex(A);
        return SparseMult(AC, B, prod);
    }

    mkl_sparse_destroy(A_MKL);
    mkl_sparse_destroy(B_MKL);
    mkl_sparse_destroy(C_MKL);
}
//------------------------------------------------------------------------------
// Set pivot threshold for MKL sparse matrix division
//------------------------------------------------------------------------------
static int iparm9 = 13;

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

        iparm9 = val;
    }
    else
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    return true;
}
//------------------------------------------------------------------------------
// Divide a full matrix by a sparse matrix on the left side, Q = A \ B
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::SparseDivideLeft(const hwMatrixS& A,
                                       const hwMatrix&  B,
                                       hwMatrix&        Q)
{
    int nRows = A.M();
    int nCols = A.N();
    int n = B.N();

    if (nRows != nCols)
        throw hwMathException(HW_MATH_ERR_MTXNOTSQUARE, 1);

    if (nRows != B.M())
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    if (B.IsEmpty())
    {
        hwMathStatus status = Q.Dimension(nCols, n, hwMatrix::REAL);
        Q.SetElements(0.0);
        return;
    }

    if (A.NNZ() == 0)
    {
        throw hwMathException(HW_MATH_ERR_SINGMATRIX);
    }

    // Internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    void* pt[64];
    // Pardiso control parameters.
    MKL_INT iparm[64];

    // Setup Pardiso control parameters
    for (MKL_INT i = 0; i < 64; i++)
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;        // No solver default
    iparm[1] = 2;        // Fill-in reordering from METIS
    iparm[3] = 0;        // No iterative-direct algorithm
    iparm[4] = 0;        // No user fill-in reducing permutation
    iparm[5] = 0;        // Write solution into x
    iparm[6] = 0;        // Not in use
    iparm[7] = 2;        // Max numbers of iterative refinement steps
    iparm[8] = 0;        // Not in use
    iparm[9] = iparm9;   // Perturb the pivot elements with 1E-13
    iparm[10] = 1;        // Use nonsymmetric permutation and scaling MPS
    iparm[11] = 2;        // Transpose solve for systems in CSC format
    iparm[12] = 1;        // Maximum weighted matching algorithm is switched-on (default for non-symmetric)
    iparm[13] = 0;        // Output: Number of perturbed pivots
    iparm[14] = 0;        // Not in use
    iparm[15] = 0;        // Not in use
    iparm[16] = 0;        // Not in use
    iparm[17] = 0;        // Output: Number of nonzeros in the factor LU (default: -1)
    iparm[18] = 0;        // Output: Mflops for LU factorization (default: -1)
    iparm[19] = 0;        // Output: Numbers of CG Iterations
    iparm[34] = 1;        // Zero based indexing

    MKL_INT maxfct = 1;   // Maximum number of numerical factorizations.
    MKL_INT mnum = 1;   // Which factorization to use.
    MKL_INT msglvl = 0;   // Do not print statistical information in file
    MKL_INT error = 0;   // Initialize error flag
    MKL_INT phase;
    double ddum;          // Auxiliary double dummy
    MKL_INT idum;         // Auxiliary integer dummy

    // Initialize the internal solver memory pointer. This is only
    // necessary for the FIRST call of the PARDISO solver.
    for (MKL_INT i = 0; i < 64; i++)
    {
        pt[i] = 0;
    }

    // define the non-zero structure of the matrix
    int* pr = const_cast<int*> (A.rows());
    int* pe = const_cast<int*> (A.pointerE());

    std::vector<int> colCount(nCols + 1);  // TODO: rework, merging pointerB and pointerE

    for (int ii = 0; ii < nCols; ++ii)
        colCount[ii + 1] = pe[ii];

    MKL_INT* ia = colCount.data();
    MKL_INT* ja = pr;

    // choose matrix case
    if (A.IsReal() && B.IsReal())
    {
        MKL_INT mtype = 11;   // Real unsymmetric matrix

        hwMathStatus status = Q.Dimension(nCols, n, hwMatrix::REAL);

        if (!status.IsOk())
            throw hwMathException(status.GetMsgCode());

        double* a = const_cast<double*> (A.GetRealData());
        double* b = const_cast<double*> (B.GetRealData());
        double* x = Q.GetRealData();

        // Reordering and Symbolic Factorization. This step also allocates
        // all memory that is necessary for the factorization
        phase = 11;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &n, iparm, &msglvl, &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &n,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during symbolic factorization: %d", error);
        }

        // Numerical factorization
        phase = 22;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &nRows, a, ia, ja, &idum, &n, iparm, &msglvl, &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &n,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during numerical factorization: %d", error);
        }

        // Solution phase
        phase = 33;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &n, iparm, &msglvl, b, x, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &n,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during solution: %d", error);
        }

        // Termination and release of memory
        phase = -1;

        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, &ddum, ia, ja, &idum, &n,
                iparm, &msglvl, &ddum, &ddum, &error);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        MKL_INT mtype = 13;   // Complex unsymmetric matrix

        hwMathStatus status = Q.Dimension(nCols, n, hwMatrix::COMPLEX);

        if (!status.IsOk())
            throw hwMathException(status.GetMsgCode());

        hwComplex* a = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* b = const_cast<hwComplex*> (B.GetComplexData());
        hwComplex* x = Q.GetComplexData();

        // Reordering and Symbolic Factorization. This step also allocates
        // all memory that is necessary for the factorization
        phase = 11;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &n, iparm, &msglvl, &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &n,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during symbolic factorization: %d", error);
        }

        // Numerical factorization
        phase = 22;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &n, iparm, &msglvl, &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &n,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during numerical factorization: %d", error);
        }

        // Solution phase
        phase = 33;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &n, iparm, &msglvl, b, x, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &n,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during solution: %d", error);
        }

        // Termination and release of memory
        phase = -1;

        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, &ddum, ia, ja, &idum, &n,
                iparm, &msglvl, &ddum, &ddum, &error);
    }
    else if (!A.IsReal())
    {
        hwMatrix C;
        C.PackComplex(B);
        return SparseDivideLeft(A, C, Q);
    }
    else    // !B.IsReal()
    {
        hwMatrixS C;
        C.PackComplex(A);
        return SparseDivideLeft(C, B, Q);
    }
}
//------------------------------------------------------------------------------
// Divide a full matrix by a sparse matrix on the right side, Q = A / B
//------------------------------------------------------------------------------
void BuiltInFuncsMKL::SparseDivideRight(const hwMatrix&  A,
                                        const hwMatrixS& B,
                                        hwMatrix&        Q)
{
    int nRows = B.M();
    int nCols = B.N();
    int m = A.M();

    if (nRows != nCols)
        throw hwMathException(HW_MATH_ERR_MTXNOTSQUARE, 1);

    if (nCols != A.N())
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    if (A.IsEmpty())
    {
        hwMathStatus status = Q.Dimension(m, nRows, hwMatrix::REAL);
        Q.SetElements(0.0);
        return;
    }

    if (B.NNZ() == 0)
    {
        throw hwMathException(HW_MATH_ERR_SINGMATRIX);
    }

    // Internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    void* pt[64];
    // Pardiso control parameters.
    MKL_INT iparm[64];

    // Setup Pardiso control parameters
    for (MKL_INT i = 0; i < 64; i++)
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;        // No solver default
    iparm[1] = 2;        // Fill-in reordering from METIS
    iparm[3] = 0;        // No iterative-direct algorithm
    iparm[4] = 0;        // No user fill-in reducing permutation
    iparm[5] = 0;        // Write solution into x
    iparm[6] = 0;        // Not in use
    iparm[7] = 2;        // Max numbers of iterative refinement steps
    iparm[8] = 0;        // Not in use
    iparm[9] = iparm9;   // Perturb the pivot elements with 1E-13
    iparm[10] = 1;       // Use nonsymmetric permutation and scaling MPS
    iparm[11] = 0;       // Non-transpose solve for systems in CSC format
    iparm[12] = 1;       // Maximum weighted matching algorithm is switched-on (default for non-symmetric)
    iparm[13] = 0;       // Output: Number of perturbed pivots
    iparm[14] = 0;       // Not in use
    iparm[15] = 0;       // Not in use
    iparm[16] = 0;       // Not in use
    iparm[17] = 0;       // Output: Number of nonzeros in the factor LU (default: -1)
    iparm[18] = 0;       // Output: Mflops for LU factorization (default: -1)
    iparm[19] = 0;       // Output: Numbers of CG Iterations
    iparm[34] = 1;       // Zero based indexing

    MKL_INT maxfct = 1;  // Maximum number of numerical factorizations.
    MKL_INT mnum = 1;    // Which factorization to use.
    MKL_INT msglvl = 0;  // Do not print statistical information in file
    MKL_INT error = 0;   // Initialize error flag
    MKL_INT phase;
    double ddum;         // Auxiliary double dummy
    MKL_INT idum;        // Auxiliary integer dummy

    // Initialize the internal solver memory pointer. This is only
    // necessary for the FIRST call of the PARDISO solver.
    for (MKL_INT i = 0; i < 64; i++)
    {
        pt[i] = 0;
    }

    // define the non-zero structure of the matrix
    int* pr = const_cast<int*> (B.rows());
    int* pe = const_cast<int*> (B.pointerE());

    std::vector<int> colCount(nCols + 1);  // TODO: rework, merging pointerB and pointerE

    for (int ii = 0; ii < nCols; ++ii)
        colCount[ii + 1] = pe[ii];

    MKL_INT* ia = colCount.data();
    MKL_INT* ja = pr;

    hwMatrix AT;

    hwMathStatus status = AT.Transpose(A);

    if (!status.IsOk())
        throw hwMathException(status.GetMsgCode());

    // choose matrix case
    if (A.IsReal() && B.IsReal())
    {
        MKL_INT mtype = 11;   // Real unsymmetric matrix

        status = Q.Dimension(nRows, m, hwMatrix::REAL);

        if (!status.IsOk())
            throw hwMathException(status.GetMsgCode());

        double* a = const_cast<double*> (B.GetRealData());
        double* b = const_cast<double*> (AT.GetRealData());
        double* x = Q.GetRealData();

        // Reordering and Symbolic Factorization. This step also allocates
        // all memory that is necessary for the factorization
        phase = 11;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &m, iparm, &msglvl, &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &m,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during symbolic factorization: %d", error);
        }

        // Numerical factorization
        phase = 22;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &m, iparm, &msglvl, &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &m,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during numerical factorization: %d", error);
        }

        // Solution phase
        phase = 33;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &m, iparm, &msglvl, b, x, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &m,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during solution: %d", error);
        }

        status = Q.Transpose();

        // Termination and release of memory
        phase = -1;

        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, &ddum, ia, ja, &idum, &m,
                iparm, &msglvl, &ddum, &ddum, &error);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        MKL_INT mtype = 13;   // Complex unsymmetric matrix

        status = Q.Dimension(nRows, m, hwMatrix::COMPLEX);

        if (!status.IsOk())
            throw hwMathException(status.GetMsgCode());

        hwComplex* a = const_cast<hwComplex*> (B.GetComplexData());
        hwComplex* b = const_cast<hwComplex*> (AT.GetComplexData());
        hwComplex* x = Q.GetComplexData();

        // Reordering and Symbolic Factorization. This step also allocates
        // all memory that is necessary for the factorization
        phase = 11;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &m, iparm, &msglvl, &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &m,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during symbolic factorization: %d", error);
        }

        // Numerical factorization
        phase = 22;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &m, iparm, &msglvl, &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &m,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during numerical factorization: %d", error);
        }

        // Solution phase
        phase = 33;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &m, iparm, &msglvl, b, x, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &m,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during solution: %d", error);
        }

        status = Q.Transpose();

        // Termination and release of memory
        phase = -1;

        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, &ddum, ia, ja, &idum, &m,
                iparm, &msglvl, &ddum, &ddum, &error);
    }
    else if (!A.IsReal())
    {
        hwMatrixS C;
        C.PackComplex(B);
        return SparseDivideRight(A, C, Q);
    }
    else    // !B.IsReal()
    {
        hwMatrix C;
        C.PackComplex(A);
        return SparseDivideRight(C, B, Q);
    }
}
