/**
* @file PolynomialTboxFuncs.cxx
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

#include "PolynomialTboxFuncs.h"

#include "PolynomFuncs.h"

#include "BuiltInFuncs.h"
#include "BuiltInFuncsUtils.h"
#include "MatrixNUtils.h"
#include "OML_Error.h"
#include "StructData.h"
#include "hwMatrix.h"
#include "hwMatrixN.h"
#include "hwSliceArg.h"
#include "MathUtilsFuncs.h"

#include <memory>
#define POLY   "PolynomialMath"
#define TBOXVERSION 2019.0


//------------------------------------------------------------------------------
// Entry point which registers polynomial functions with oml
//------------------------------------------------------------------------------
int InitDll(EvaluatorInterface eval)
{
    eval.RegisterBuiltInFunction("roots", &OmlRoots,
        FunctionMetaData(1, 1, POLY));
    eval.RegisterBuiltInFunction("lookup", &OmlLookUp,
        FunctionMetaData(-3, 1, POLY));
    eval.RegisterBuiltInFunction("spline", &OmlSpline,
        FunctionMetaData(-3, 1, POLY));
    eval.RegisterBuiltInFunction("pchip", &OmlPchip,
        FunctionMetaData(-3, 1, POLY));
    eval.RegisterBuiltInFunction("interp1", &OmlInterp1,
        FunctionMetaData(-4, 1, POLY));
    eval.RegisterBuiltInFunction("interp2", &OmlInterp2,
        FunctionMetaData(-6, 1, POLY));
    eval.RegisterBuiltInFunction("interpn", &OmlInterpN,
        FunctionMetaData(-2, 1, POLY));
    eval.RegisterBuiltInFunction("isonormals", &OmlIsoNormals,
        FunctionMetaData(5, 1, POLY));
    eval.RegisterBuiltInFunction("deconv", &OmlDeconv,
        FunctionMetaData(2, 2, POLY));
    eval.RegisterBuiltInFunction("polyder", &OmlPolyder,
        FunctionMetaData(-2, -2, POLY));
    eval.RegisterBuiltInFunction("polyint", &OmlPolyint,
        FunctionMetaData(-2, 1, POLY));
    eval.RegisterBuiltInFunction("mkpp", &OmlMakePPoly,
        FunctionMetaData(-3, 1, POLY));
    eval.RegisterBuiltInFunction("unmkpp", &OmlUnMakePPoly,
        FunctionMetaData(3, -5, POLY));
    eval.RegisterBuiltInFunction("ppval", &OmlPPolyEval,
        FunctionMetaData(-3, 1, POLY));
    return 1;
}
//------------------------------------------------------------------------------
// Computes the roots of a polynomial and returns true
//------------------------------------------------------------------------------
bool OmlRoots(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs,
              std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);

    const hwMatrix* coefs = inputs[0].ConvertToMatrix();
    std::unique_ptr<hwMatrix> result(EvaluatorInterface::allocateMatrix());

    BuiltInFuncsUtils::CheckMathStatus(eval, PolyRoots(*coefs, *result));
    outputs.push_back(result.release());

    return true;
}
//------------------------------------------------------------------------------
// Returns true after looking up values in a sorted table
//------------------------------------------------------------------------------
bool OmlLookUp(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs,
               std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin != 2 && nargin != 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsCellArray())
    {
        HML_CELLARRAY* table = inputs[0].CellArray();

        // set up table direction
        if (!table->IsVector())
            throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_TYPE);

        int tableLength = table->Size();
        const Currency* tablePtr = table->GetRealData();
        const Currency& tableCur1 = (*table)(0);
        const Currency& tableCur2 = (*table)(tableLength - 1);
        const hwMatrix* tableMtx1 = tableCur1.Matrix();
        const hwMatrix* tableMtx2 = tableCur2.Matrix();
        int length = _min(tableMtx1->Size(), tableMtx2->Size());
        bool forward = true;

        for (int k = 0; k < length; ++k)
        {
            if ((*tableMtx1)(k) == (*tableMtx2)(k))
            {
                continue;
            }

            if ((*tableMtx1)(k) > (*tableMtx2)(k))
            {
                forward = false;
            }

            break;
        }

        // manage input options
        int dataLength = -1;
        const Currency* dataPtr = nullptr;
        double* indxPtr = nullptr;

        if (inputs[1].IsString())
        {
            dataLength = 1;
            dataPtr = &inputs[1];
            std::unique_ptr<hwMatrix> indices(EvaluatorInterface::allocateMatrix(1, 1, true));
            indxPtr = indices->GetRealData();
            outputs.push_back(indices.release());
        }
        else if (inputs[1].IsCellArray())
        {
            HML_CELLARRAY* itemCells = inputs[1].CellArray();
            dataLength = itemCells->Size();
            dataPtr = itemCells->GetRealData();
            std::unique_ptr<hwMatrix> indices(EvaluatorInterface::allocateMatrix(itemCells->M(), itemCells->N(), true));
            indxPtr = indices->GetRealData();
            outputs.push_back(indices.release());
        }
        else if (inputs[1].IsNDCellArray())
        {
            HML_ND_CELLARRAY* itemCells = inputs[1].CellArrayND();
            dataLength = itemCells->Size();
            dataPtr = itemCells->GetRealData();
            std::unique_ptr<hwMatrixN> indices(EvaluatorInterface::allocateMatrixN(itemCells->Dimensions(), true));
            indxPtr = indices->GetRealData();
            outputs.push_back(indices.release());
        }
        else
        {
            throw OML_Error(OML_ERR_STRING_STRINGCELL, 2, OML_VAR_TYPE);
        }

        // do string lookup
        for (int i = 0; i < dataLength; ++i)
        {
            const Currency& itemCur = dataPtr[i];

            if (!itemCur.IsString())
            {
                throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
            }

            const hwMatrix* itemMtx = itemCur.Matrix();
            bool insert = false;

            if (forward)
            {
                indxPtr[i] = tableLength;

                for (int j = 0; j < tableLength; ++j)
                {
                    const Currency& tableCur = (*table)(j);

                    if (!tableCur.IsString())
                    {
                        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
                    }

                    const hwMatrix* tableMtx = tableCur.Matrix();
                    int length = _min(itemMtx->Size(), tableMtx->Size());
                    bool matched = true;

                    for (int k = 0; k < length; ++k)
                    {
                        if ((*itemMtx)(k) == (*tableMtx)(k))
                        {
                            continue;
                        }

                        matched = false;

                        if ((*itemMtx)(k) < (*tableMtx)(k))
                        {
                            insert = true;
                        }

                        break;
                    }

                    if (matched)
                    {
                        if (itemMtx->Size() == tableMtx->Size())
                        {
                            indxPtr[i] = j + 1;
                        }
                        else if (itemMtx->Size() < tableMtx->Size())
                        {
                            indxPtr[i] = j;
                        }
                        else
                        {
                            continue;
                        }
                    }
                    else if (!insert)
                    {
                        continue;
                    }
                    else
                    {
                        indxPtr[i] = j;
                    }

                    break;
                }
            }
            else    // reverse
            {
                indxPtr[i] = 0;

                for (int j = tableLength - 1; j > -1; --j)
                {
                    const Currency& tableCur = (*table)(j);

                    if (!tableCur.IsString())
                    {
                        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
                    }

                    const hwMatrix* tableMtx = tableCur.Matrix();
                    int length = _min(itemMtx->Size(), tableMtx->Size());
                    bool matched = true;

                    for (int k = 0; k < length; ++k)
                    {
                        if ((*itemMtx)(k) == (*tableMtx)(k))
                        {
                            continue;
                        }

                        matched = false;

                        if ((*itemMtx)(k) < (*tableMtx)(k))
                        {
                            insert = true;
                        }

                        break;
                    }

                    if (matched)
                    {
                        if (itemMtx->Size() == tableMtx->Size())
                        {
                            indxPtr[i] = j + 1;
                        }
                        else if (itemMtx->Size() < tableMtx->Size())
                        {
                            indxPtr[i] = j + 1;
                        }
                        else
                        {
                            continue;
                        }
                    }
                    else if (!insert)
                    {
                        continue;
                    }
                    else
                    {
                        indxPtr[i] = j + 1;
                    }

                    break;
                }
            }
        }
        if (nargin == 3)
        {
            if (!inputs[2].IsString())
            {
                throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);
            }

            std::string opt = inputs[2].StringVal();

            if (opt == "m")
            {
                for (int i = 0; i < dataLength; ++i)
                {
                    if (dataPtr[i].StringVal() != tablePtr[static_cast<int>(indxPtr[i]) - 1].StringVal())
                        indxPtr[i] = 0;
                }
            }
            else if (opt == "b")
            {
                for (int i = 0; i < dataLength; ++i)
                {
                    if (dataPtr[i].StringVal() == tablePtr[static_cast<int>(indxPtr[i]) - 1].StringVal())
                        indxPtr[i] = 1;
                    else
                        indxPtr[i] = 0;
                }

                outputs[0].SetMask(Currency::MASK_LOGICAL);
            }
            else
            {
                throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_VALUE);
            }
        }
    }
    else
    {
        if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
            throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_TYPE);

        if (inputs[0].IsMatrix() && !inputs[0].Matrix()->IsEmptyOrVector())
            throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_TYPE);

        const hwMatrix* table = inputs[0].ConvertToMatrix();

        if (!table->IsReal())
            throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_TYPE);

        int tableLength = table->Size();
        const double* tablePtr = table->GetRealData();
        int dataLength = -1;
        const double* dataPtr = nullptr;
        double* indxPtr = nullptr;

        if (inputs[1].IsMatrix() || inputs[1].IsScalar())
        {
            const hwMatrix* y = inputs[1].ConvertToMatrix();

            if (!y->IsReal())
                throw OML_Error(OML_ERR_REALMATRIX, 2, OML_VAR_TYPE);

            dataLength = y->Size();
            dataPtr = y->GetRealData();
            std::unique_ptr<hwMatrix> indices(EvaluatorInterface::allocateMatrix(y->M(), y->N(), true));
            indxPtr = indices->GetRealData();
            outputs.push_back(indices.release());

            if (!tableLength)
            {
                indices->SetElements(0.0);
                return true;
            }
        }
        else if (inputs[1].IsNDMatrix())
        {
            const hwMatrixN* y = inputs[1].MatrixN();

            if (!y->IsReal())
                throw OML_Error(OML_ERR_REALMATRIX, 2, OML_VAR_TYPE);

            dataLength = y->Size();
            dataPtr = y->GetRealData();
            std::unique_ptr<hwMatrixN> indices(EvaluatorInterface::allocateMatrixN(y->Dimensions(), true));
            indxPtr = indices->GetRealData();
            outputs.push_back(indices.release());

            if (!tableLength)
            {
                indices->SetElements(0.0);
                return true;
            }
        }
        else
        {
            throw OML_Error(OML_ERR_REALMATRIX, 2, OML_VAR_TYPE);
        }

        // do numeric lookup
        if (tablePtr[0] < tablePtr[tableLength - 1])
        {
            for (int i = 0; i < dataLength; ++i)
            {
                *indxPtr++ = BinarySearch(tablePtr, tableLength, *dataPtr++) + 1;
            }
        }
        else    // reverse table
        {
            for (int i = 0; i < dataLength; ++i)
            {
                *indxPtr++ = BinarySearchR(tablePtr, tableLength, *dataPtr++) + 1;
            }
        }

        dataPtr -= dataLength;
        indxPtr -= dataLength;

        if (nargin == 3)
        {
            if (!inputs[2].IsString())
            {
                throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);
            }

            std::string opt = inputs[2].StringVal();

            if (opt == "m")
            {
                for (int i = 0; i < dataLength; ++i)
                {
                    if (dataPtr[i] != tablePtr[static_cast<int>(indxPtr[i]) - 1])
                        indxPtr[i] = 0;
                }
            }
            else if (opt == "b")
            {
                for (int i = 0; i < dataLength; ++i)
                {
                    if (dataPtr[i] == tablePtr[static_cast<int>(indxPtr[i]) - 1])
                        indxPtr[i] = 1;
                    else
                        indxPtr[i] = 0;
                }

                outputs[0].SetMask(Currency::MASK_LOGICAL);
            }
            else if (opt == "l")
            {
                for (int i = 0; i < dataLength; ++i)
                {
                    if (indxPtr[i] < 1)
                        indxPtr[i] = 1;
                }
            }
            else if (opt == "r")
            {
                for (int i = 0; i < dataLength; ++i)
                {
                    if (indxPtr[i] > tableLength - 1)
                        indxPtr[i] = tableLength - 1;
                }
            }
            else
            {
                throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_VALUE);
            }
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Interpolates (x,y) data with a cubic spline and returns true
//------------------------------------------------------------------------------
bool OmlSpline(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs,
               std::vector<Currency>&       outputs)
{
    static bool unsorted = true;
    size_t nargin = inputs.size();

    if (nargin != 2 && nargin != 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_TYPE);

    if (inputs[0].IsMatrix() && !inputs[0].Matrix()->IsVector())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_TYPE);

    if (!inputs[1].IsMatrix() && !inputs[1].IsNDMatrix() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_REALMATRIX, 2, OML_VAR_TYPE);

    const hwMatrix* x_old = inputs[0].ConvertToMatrix();

    // enforce sorted x,y lists
    if (unsorted)
    {
        std::vector<Currency> inputs2;
        inputs2.push_back(inputs[0]);           // push_back(x)
        inputs2.push_back("either");

        oml_issorted(eval, inputs2, outputs);   // issorted(x)

        if (outputs[0].Scalar() == 0)
        {
            // sort x
            outputs.clear();
            inputs2.pop_back();
            inputs2.push_back("ascend");
            oml_sort(eval, inputs2, outputs);   // sort(x)
            inputs2.clear();
            inputs2.push_back(outputs[0]);      // push_back(sorted_x)
            outputs.erase(outputs.begin());     // remove sorted_x, leaving indices

            // handle clamped spline case
            bool clamped = false;
            int  numPnts = x_old->Size();

            if (inputs[1].IsMatrix())
            {
                const hwMatrix* y_old = inputs[1].ConvertToMatrix();

                if (y_old->M() == 1)
                {
                    if (y_old->N() == numPnts + 2)
                        clamped = true;
                }
                else if (y_old->M() == numPnts + 2)
                {
                    clamped = true;
                }
            }
            else if (inputs[1].IsNDMatrix())
            {
                if (inputs[1].MatrixN()->Dimensions()[0] == numPnts + 2)
                {
                    clamped = true;
                }
            }

            if (clamped)
            {
                hwMatrix*    index = outputs[0].GetWritableMatrix();
                int          m = index->M();
                int          n = index->N();
                hwMathStatus status;

                if (m == numPnts)
                {
                    status = index->Resize(m + 2, n);
                }
                else
                {
                    status = index->Resize(m, n + 2);
                }

                BuiltInFuncsUtils::CheckMathStatus(eval, status);

                for (int i = numPnts - 1; i > -1; --i)
                {
                    (*index)(i + 1) = (*index)(i) + 1;
                }

                (*index)(0) = 1;
                (*index)(numPnts + 1) = numPnts + 2;
            }

            // reorder y (by row for spline)
            if (inputs[1].IsMatrix())
            {
                if (inputs[1].Matrix()->N() > 1)
                    outputs.insert(outputs.begin(), Currency(0.0, Currency::TYPE_COLON));
            }
            else if (inputs[1].IsNDMatrix())
            {
                const hwMatrixN* y = inputs[1].MatrixN();
                int numDims = static_cast<int>(y->Dimensions().size());

                for (int i = 1; i < numDims; ++i)
                    outputs.insert(outputs.begin(), Currency(0.0, Currency::TYPE_COLON));
            }

            Currency yOrdered = eval.VariableIndex(inputs[1], outputs);   // reordered_y = y(indices)

            // recursive call with modified arguments
            outputs.clear();
            inputs2.push_back(yOrdered);        // push_back(reordered_y)

            if (nargin == 3)
                inputs2.push_back(inputs[2]);   // push_back(xi)

            try
            {
                unsorted = false;
                bool retv = OmlSpline(eval, inputs2, outputs); // call function on sorted_x, reordered_y
                return retv;
            }
            catch (OML_Error&)
            {
                unsorted = true;    // reset
                throw;
            }
            catch (hwMathException&)
            {
                unsorted = true;    // reset
                throw;
            }
        }

        outputs.clear();
    }
    else
    {
        unsorted = true;    // reset
    }

    int n = x_old->Size();

    if (nargin == 2)
    {
        std::vector<Currency> inputs2;
        inputs2.push_back(EvaluatorInterface::allocateMatrix(x_old));
        int pieces = n - 1;

        // process yi
        if (inputs[1].IsMatrix() || inputs[1].IsScalar())
        {
            const hwMatrix* y_old = inputs[1].ConvertToMatrix();

            if (y_old->IsVector())
            {
                std::unique_ptr<hwMatrix> coefs(EvaluatorInterface::allocateMatrix());

                if (y_old->Size() == n + 2)
                {
                    // clamped spline
                    // create offset temp y vector for convenience
                    const double* y_vec = y_old->GetRealData();
                    hwMatrix y_temp(n, (void*) ++y_vec, hwMatrix::REAL);

                    BuiltInFuncsUtils::CheckMathStatus(eval, Spline(*x_old, y_temp, (*y_old)(0), (*y_old)(n + 1), *coefs));
                }
                else // (y_old->Size() == n)
                {
                    // not-a-knot spline
                    BuiltInFuncsUtils::CheckMathStatus(eval, Spline(*x_old, *y_old, *coefs));
                }

                inputs2.push_back(coefs.release());
            }
            else
            {
                int dim = y_old->M();
                std::vector<int> dims(3);

                dims[0] = dim;
                dims[1] = pieces;
                dims[2] = 4;
                hwMatrixN* coefs = eval.allocateMatrixN(dims, true);
                hwMatrix coefs2D(dim * pieces, 4, coefs->GetRealData(), hwMatrix::REAL);

                if (dim)
                {
                    std::unique_ptr<hwMatrix> yoldrow(EvaluatorInterface::allocateMatrix());
                    std::unique_ptr<hwMatrix> ppcoefs(EvaluatorInterface::allocateMatrix());
                    hwMathStatus status;

                    for (int i = 0; i < dim; ++i)
                    {
                        BuiltInFuncsUtils::CheckMathStatus(eval, y_old->ReadRow(i, *yoldrow));

                        // check length
                        if (yoldrow->Size() == n + 2)
                        {
                            // clamped spline
                            // create offset temp y vector for convenience
                            const double* y_vec = yoldrow->GetRealData();
                            hwMatrix y_temp(n, (void*) ++y_vec, hwMatrix::REAL);

                            status = Spline(*x_old, y_temp, (*y_old)(0), (*y_old)(n + 1), *ppcoefs);
                        }
                        else
                        {
                            // not-a-knot spline
                            status = Spline(*x_old, *yoldrow, *ppcoefs);
                        }

                        if (!status.IsOk())
                        {
                            if (status == HW_MATH_ERR_ARRAYSIZE)
                            {
                                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
                            }
                            else
                            {
                                BuiltInFuncsUtils::CheckMathStatus(eval, status);
                            }
                        }

                        for (int j = 0; j < pieces; ++j)
                        {
                            hwMatrix coefrow;
                            BuiltInFuncsUtils::CheckMathStatus(eval, ppcoefs->ReadRow(j, coefrow));
                            BuiltInFuncsUtils::CheckMathStatus(eval, coefs2D.WriteRow(i + j * dim, coefrow));
                        }
                    }

                    inputs2.push_back(coefs);
                    inputs2.push_back(dim);
                }
            }
        }
        else if (inputs[1].IsNDMatrix())
        {
            const hwMatrixN* y_old = inputs[1].MatrixN();
            const std::vector<int>& dims_old = y_old->Dimensions();
            int numDims = static_cast<int> (dims_old.size());
            int dim = y_old->Size() / dims_old[numDims - 1];
            hwMatrix y_old2D(dim, dims_old[numDims - 1], (void*) y_old->GetRealData(), hwMatrix::REAL);

            std::vector<int> ppdims(numDims + 1);
            hwMatrix* ppDims = eval.allocateMatrix(1, numDims - 1, true);

            for (int i = 0; i < numDims - 1; ++i)
            {
                ppdims[i] = dims_old[i];
                (*ppDims)(i) = dims_old[i];
            }

            ppdims[numDims - 1] = pieces;
            ppdims[numDims] = 4;

            hwMatrixN* coefs = eval.allocateMatrixN(ppdims, true);
            hwMatrix coefs2D(dim * pieces, 4, coefs->GetRealData(), hwMatrix::REAL);

            if (dim)
            {
                std::unique_ptr<hwMatrix> yoldrow(EvaluatorInterface::allocateMatrix());
                std::unique_ptr<hwMatrix> ppcoefs(EvaluatorInterface::allocateMatrix());
                hwMathStatus status;

                for (int i = 0; i < dim; ++i)
                {
                    BuiltInFuncsUtils::CheckMathStatus(eval, y_old2D.ReadRow(i, *yoldrow));

                    // check length
                    if (yoldrow->Size() == n + 2)
                    {
                        // clamped spline
                        // create offset temp y vector for convenience
                        const double* y_vec = yoldrow->GetRealData();
                        hwMatrix y_temp(n, (void*) ++y_vec, hwMatrix::REAL);

                        status = Spline(*x_old, y_temp, (*y_old)(0), (*y_old)(n + 1), *ppcoefs);
                    }
                    else
                    {
                        // not-a-knot spline
                        status = Spline(*x_old, *yoldrow, *ppcoefs);
                    }

                    if (!status.IsOk())
                    {
                        if (status == HW_MATH_ERR_ARRAYSIZE)
                        {
                            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
                        }
                        else
                        {
                            BuiltInFuncsUtils::CheckMathStatus(eval, status);
                        }
                    }

                    for (int j = 0; j < pieces; ++j)
                    {
                        hwMatrix coefrow;
                        BuiltInFuncsUtils::CheckMathStatus(eval, ppcoefs->ReadRow(j, coefrow));
                        BuiltInFuncsUtils::CheckMathStatus(eval, coefs2D.WriteRow(i + j * dim, coefrow));
                    }
                }
            }

            inputs2.push_back(coefs);
            inputs2.push_back(ppDims);
        }

        return OmlMakePPoly(eval, inputs2, outputs);
    }
    else // nargin == 3
    {
        // process yi when ND
        if (inputs[1].IsNDMatrix())
        {
            const hwMatrixN* y = inputs[1].MatrixN();
            int numDims = static_cast<int>(y->Dimensions().size());

            return oml_MatrixNUtil4(eval, inputs, outputs, OmlSpline, -numDims, 2);
        }

        // make sure x_new is a row when y is not a vector
        if (!inputs[2].IsMatrix() && !inputs[2].IsScalar())
            throw OML_Error(OML_ERR_REALVECTOR, 3, OML_VAR_TYPE);

        const hwMatrix* x_new = inputs[2].ConvertToMatrix();

        if (x_new->M() > 1)
        {
            if ((inputs[1].IsMatrix() && inputs[1].Matrix()->M() != 1 && inputs[1].Matrix()->N() != 1) ||
                inputs[1].IsNDMatrix())
            {
                hwMatrix* row = eval.allocateMatrix(x_new);
                row->Transpose();

                std::vector<Currency> inputs3;
                inputs3.push_back(inputs[0]);   // push_back(sorted_x)
                inputs3.push_back(inputs[1]);   // push_back(reordered_y)
                inputs3.push_back(row);         // push_back(col_xi)

                return OmlSpline(eval, inputs3, outputs);
            }
        }

        // process yi
        const hwMatrix* y_old = inputs[1].ConvertToMatrix();

        if (y_old->IsVector())
        {
            std::unique_ptr<hwMatrix> y_new(EvaluatorInterface::allocateMatrix());

            if (y_old->Size() == n + 2)
            {
                // clamped spline
                // create offset temp y vector for convenience
                const double* y_vec = y_old->GetRealData();
                hwMatrix y_temp(n, (void*) ++y_vec, hwMatrix::REAL);

                BuiltInFuncsUtils::CheckMathStatus(eval, Spline(*x_old, y_temp, (*y_old)(0), (*y_old)(n + 1), *x_new, *y_new, true));
            }
            else // (y_old->Size() == n)
            {
                // not-a-knot spline
                BuiltInFuncsUtils::CheckMathStatus(eval, Spline(*x_old, *y_old, *x_new, *y_new, true));
            }

            outputs.push_back(y_new.release());
        }
        else
        {
            if (x_old->Size() != y_old->N())
            {
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
            }

            hwMatrix* y_new = eval.allocateMatrix(y_old->M(), x_new->Size(), true);

            if (y_old->M())
            {
                std::unique_ptr<hwMatrix> yoldrow(EvaluatorInterface::allocateMatrix());
                std::unique_ptr<hwMatrix> ynewrow(EvaluatorInterface::allocateMatrix());
                hwMathStatus status;

                for (int i = 0; i < y_old->M(); ++i)
                {
                    BuiltInFuncsUtils::CheckMathStatus(eval, y_old->ReadRow(i, *yoldrow));

                    // check length
                    if (yoldrow->Size() == n + 2)
                    {
                        // clamped spline
                        // create offset temp y vector for convenience
                        const double* y_vec = yoldrow->GetRealData();
                        hwMatrix y_temp(n, (void*) ++y_vec, hwMatrix::REAL);

                        status = Spline(*x_old, y_temp, (*y_old)(0), (*y_old)(n + 1), *x_new, *y_new, true);
                    }
                    else
                    {
                        // not-a-knot spline
                        status = Spline(*x_old, *yoldrow, *x_new, *ynewrow, true);
                    }

                    if (!status.IsOk())
                    {
                        if (status == HW_MATH_ERR_ARRAYSIZE)
                        {
                            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
                        }
                        else
                        {
                            BuiltInFuncsUtils::CheckMathStatus(eval, status);
                        }
                    }

                    BuiltInFuncsUtils::CheckMathStatus(eval, y_new->WriteRow(i, *ynewrow));
                }
            }

            outputs.push_back(y_new);
        }
    }

    return true;
}
//------------------------------------------------------------------------------
//! Returns true after interpolating (x,y) data with a cubic Hermite polynomial
//------------------------------------------------------------------------------
bool OmlPchip(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs,
              std::vector<Currency>&       outputs)
{
    static bool unsorted = true;
    size_t nargin = inputs.size();

    if (nargin != 2 && nargin != 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_TYPE);

    if (inputs[0].IsMatrix() && !inputs[0].Matrix()->IsVector())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_TYPE);

    if (!inputs[1].IsMatrix() && !inputs[1].IsNDMatrix() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_REALMATRIX, 2, OML_VAR_TYPE);

    const hwMatrix* x_old = inputs[0].ConvertToMatrix();

    // enforce sorted x,y lists
    if (unsorted)
    {
        std::vector<Currency> inputs2;
        inputs2.push_back(inputs[0]);           // push_back(x)
        inputs2.push_back("either");

        oml_issorted(eval, inputs2, outputs);   // issorted(x)

        if (outputs[0].Scalar() == 0)
        {
            // sort x
            outputs.clear();
            inputs2.pop_back();
            inputs2.push_back("ascend");
            oml_sort(eval, inputs2, outputs);   // sort(x)
            inputs2.clear();
            inputs2.push_back(outputs[0]);      // push_back(sorted_x)
            outputs.erase(outputs.begin());     // remove sorted_x, leaving indices

            // reorder y (by row for spline)
            if (inputs[1].IsMatrix())
            {
                if (inputs[1].Matrix()->N() > 1)
                    outputs.insert(outputs.begin(), Currency(0.0, Currency::TYPE_COLON));
            }
            else if (inputs[1].IsNDMatrix())
            {
                const hwMatrixN* y = inputs[1].MatrixN();
                int numDims = static_cast<int>(y->Dimensions().size());

                for (int i = 1; i < numDims; ++i)
                    outputs.insert(outputs.begin(), Currency(0.0, Currency::TYPE_COLON));
            }

            Currency yOrdered = eval.VariableIndex(inputs[1], outputs);   // reordered_y = y(indices)

            // recursive call with modified arguments
            outputs.clear();
            inputs2.push_back(yOrdered);        // push_back(reordered_y)

            if (nargin == 3)
                inputs2.push_back(inputs[2]);   // push_back(xi)

            try
            {
                unsorted = false;
                bool retv = OmlPchip(eval, inputs2, outputs); // call function on sorted_x, reordered_y
                return retv;
            }
            catch (OML_Error&)
            {
                unsorted = true;    // reset
                throw;
            }
            catch (hwMathException&)
            {
                unsorted = true;    // reset
                throw;
            }
        }

        outputs.clear();
    }
    else
    {
        unsorted = true;    // reset
    }

    int n = x_old->Size();

    if (nargin == 2)
    {
        std::vector<Currency> inputs2;
        inputs2.push_back(EvaluatorInterface::allocateMatrix(x_old));
        int pieces = n - 1;

        // process yi
        if (inputs[1].IsMatrix() || inputs[1].IsScalar())
        {
            const hwMatrix* y_old = inputs[1].ConvertToMatrix();

            if (y_old->IsVector())
            {
                std::unique_ptr<hwMatrix> coefs(EvaluatorInterface::allocateMatrix());
                BuiltInFuncsUtils::CheckMathStatus(eval, PchipInterp(*x_old, *y_old, *coefs));
                inputs2.push_back(coefs.release());
            }
            else
            {
                int dim = y_old->M();
                std::vector<int> dims(3);

                dims[0] = dim;
                dims[1] = pieces;
                dims[2] = 4;
                hwMatrixN* coefs = eval.allocateMatrixN(dims, true);
                hwMatrix coefs2D(dim * pieces, 4, coefs->GetRealData(), hwMatrix::REAL);

                if (dim)
                {
                    std::unique_ptr<hwMatrix> yoldrow(EvaluatorInterface::allocateMatrix());
                    std::unique_ptr<hwMatrix> ppcoefs(EvaluatorInterface::allocateMatrix());
                    hwMathStatus status;

                    for (int i = 0; i < dim; ++i)
                    {
                        BuiltInFuncsUtils::CheckMathStatus(eval, y_old->ReadRow(i, *yoldrow));

                        // check length
                        status = PchipInterp(*x_old, *yoldrow, *ppcoefs);

                        if (!status.IsOk())
                        {
                            if (status == HW_MATH_ERR_ARRAYSIZE)
                            {
                                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
                            }
                            else
                            {
                                BuiltInFuncsUtils::CheckMathStatus(eval, status);
                            }
                        }

                        for (int j = 0; j < pieces; ++j)
                        {
                            hwMatrix coefrow;
                            BuiltInFuncsUtils::CheckMathStatus(eval, ppcoefs->ReadRow(j, coefrow));
                            BuiltInFuncsUtils::CheckMathStatus(eval, coefs2D.WriteRow(i + j * dim, coefrow));
                        }
                    }

                    inputs2.push_back(coefs);
                    inputs2.push_back(dim);
                }
            }
        }
        else if (inputs[1].IsNDMatrix())
        {
            const hwMatrixN* y_old = inputs[1].MatrixN();
            const std::vector<int>& dims_old = y_old->Dimensions();
            int numDims = static_cast<int> (dims_old.size());
            int dim = y_old->Size() / dims_old[numDims - 1];
            hwMatrix y_old2D(dim, dims_old[numDims - 1], (void*)y_old->GetRealData(), hwMatrix::REAL);

            std::vector<int> ppdims(numDims + 1);
            hwMatrix* ppDims = eval.allocateMatrix(1, numDims - 1, true);

            for (int i = 0; i < numDims - 1; ++i)
            {
                ppdims[i] = dims_old[i];
                (*ppDims)(i) = dims_old[i];
            }

            ppdims[numDims - 1] = pieces;
            ppdims[numDims] = 4;

            hwMatrixN* coefs = eval.allocateMatrixN(ppdims, true);
            hwMatrix coefs2D(dim * pieces, 4, coefs->GetRealData(), hwMatrix::REAL);

            if (dim)
            {
                std::unique_ptr<hwMatrix> yoldrow(EvaluatorInterface::allocateMatrix());
                std::unique_ptr<hwMatrix> ppcoefs(EvaluatorInterface::allocateMatrix());
                hwMathStatus status;

                for (int i = 0; i < dim; ++i)
                {
                    BuiltInFuncsUtils::CheckMathStatus(eval, y_old2D.ReadRow(i, *yoldrow));

                    status = PchipInterp(*x_old, *yoldrow, *ppcoefs);

                    if (!status.IsOk())
                    {
                        if (status == HW_MATH_ERR_ARRAYSIZE)
                        {
                            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
                        }
                        else
                        {
                            BuiltInFuncsUtils::CheckMathStatus(eval, status);
                        }
                    }

                    for (int j = 0; j < pieces; ++j)
                    {
                        hwMatrix coefrow;
                        BuiltInFuncsUtils::CheckMathStatus(eval, ppcoefs->ReadRow(j, coefrow));
                        BuiltInFuncsUtils::CheckMathStatus(eval, coefs2D.WriteRow(i + j * dim, coefrow));
                    }
                }
            }

            inputs2.push_back(coefs);
            inputs2.push_back(ppDims);
        }

        return OmlMakePPoly(eval, inputs2, outputs);
    }
    else // nargin == 3
    {
        // process yi when ND
        if (inputs[1].IsNDMatrix())
        {
            const hwMatrixN* y = inputs[1].MatrixN();
            int numDims = static_cast<int>(y->Dimensions().size());

            return oml_MatrixNUtil4(eval, inputs, outputs, OmlPchip, -numDims, 2);
        }

        // make sure x_new is a row when y is not a vector
        if (!inputs[2].IsMatrix() && !inputs[2].IsScalar())
            throw OML_Error(OML_ERR_REALVECTOR, 3, OML_VAR_TYPE);

        const hwMatrix* x_new = inputs[2].ConvertToMatrix();

        if (x_new->M() > 1)
        {
            if ((inputs[1].IsMatrix() && inputs[1].Matrix()->M() != 1 && inputs[1].Matrix()->N() != 1) ||
                inputs[1].IsNDMatrix())
            {
                hwMatrix* row = eval.allocateMatrix(x_new);
                row->Transpose();

                std::vector<Currency> inputs3;
                inputs3.push_back(inputs[0]);   // push_back(sorted_x)
                inputs3.push_back(inputs[1]);   // push_back(reordered_y)
                inputs3.push_back(row);         // push_back(col_xi)

                return OmlPchip(eval, inputs3, outputs);
            }
        }

        // process yi
        const hwMatrix* y_old = inputs[1].ConvertToMatrix();

        if (y_old->IsVector())
        {
            std::unique_ptr<hwMatrix> y_new(EvaluatorInterface::allocateMatrix());
            BuiltInFuncsUtils::CheckMathStatus(eval, PchipInterp(*x_old, *y_old, *x_new, *y_new, true));
            outputs.push_back(y_new.release());
        }
        else
        {
            if (x_old->Size() != y_old->N())
            {
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
            }

            hwMatrix* y_new = eval.allocateMatrix(y_old->M(), x_new->Size(), true);

            if (y_old->M())
            {
                std::unique_ptr<hwMatrix> yoldrow(EvaluatorInterface::allocateMatrix());
                std::unique_ptr<hwMatrix> ynewrow(EvaluatorInterface::allocateMatrix());
                hwMathStatus status;

                for (int i = 0; i < y_old->M(); ++i)
                {
                    BuiltInFuncsUtils::CheckMathStatus(eval, y_old->ReadRow(i, *yoldrow));
                    status = PchipInterp(*x_old, *yoldrow, *x_new, *ynewrow, true);

                    if (!status.IsOk())
                    {
                        if (status == HW_MATH_ERR_ARRAYSIZE)
                        {
                            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
                        }
                        else
                        {
                            BuiltInFuncsUtils::CheckMathStatus(eval, status);
                        }
                    }

                    BuiltInFuncsUtils::CheckMathStatus(eval, y_new->WriteRow(i, *ynewrow));
                }
            }

            outputs.push_back(y_new);
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Interpolates in one dimension and returns true
//------------------------------------------------------------------------------
bool OmlInterp1(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs)
{
    static bool unsorted = true;
    static bool requireUniqueX = false;     // set true if x values are unsorted
    size_t nargin = inputs.size();

    if (nargin < 3 || nargin > 5)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_TYPE);

    if (inputs[0].IsMatrix() && !inputs[0].Matrix()->IsVector())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_TYPE);

    if (!inputs[1].IsMatrix() && !inputs[1].IsNDMatrix() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_REALMATRIX, 2, OML_VAR_TYPE);

    const hwMatrix* x = inputs[0].ConvertToMatrix();

    // enforce sorted x,y lists
    if (unsorted)
    {
        std::vector<Currency> inputs2;
        inputs2.push_back(inputs[0]);           // push_back(x)
        inputs2.push_back("either");

        oml_issorted(eval, inputs2, outputs);   // issorted(x)

        if (outputs[0].Scalar() == 0)
        {
            // sort x
            outputs.clear();
            inputs2.pop_back();
            inputs2.push_back("ascend");
            oml_sort(eval, inputs2, outputs);   // sort(x)
            inputs2.clear();
            inputs2.push_back(outputs[0]);      // push_back(sorted_x)
            outputs.erase(outputs.begin());     // remove sorted_x, leaving indices

            // reorder y
            if (inputs[1].IsMatrix())
            {
                if (inputs[1].Matrix()->M() > 1)
                    outputs.push_back(Currency(0.0, Currency::TYPE_COLON));
            }
            else if (inputs[1].IsNDMatrix())
            {
                const hwMatrixN* y = inputs[1].MatrixN();
                int numDims = static_cast<int>(y->Dimensions().size());

                for (int i = 1; i < numDims; ++i)
                    outputs.push_back(Currency(0.0, Currency::TYPE_COLON));
            }

            Currency yOrdered = eval.VariableIndex(inputs[1], outputs);   // reordered_y = y(indices)

            // recursive call with modified arguments
            outputs.clear();
            inputs2.push_back(yOrdered);        // push_back(reordered_y)
            inputs2.push_back(inputs[2]);       // push_back(xi)

            if (nargin > 3)
            {
                inputs2.push_back(inputs[3]);

                if (nargin > 4)
                {
                    inputs2.push_back(inputs[4]);

                    if (nargin > 5)
                        inputs2.push_back(inputs[5]);
                }
            }

            try
            {
                unsorted       = false;
                requireUniqueX = true;
                bool retv = OmlInterp1(eval, inputs2, outputs); // call function on sorted_x, reordered_y
                return retv;
            }
            catch (OML_Error&)
            {
                unsorted       = true;    // reset
                requireUniqueX = false;
                throw;
            }
            catch (hwMathException&)
            {
                unsorted       = true;    // reset
                requireUniqueX = false;
                throw;
            }
        }

        outputs.clear();
    }
    else
    {
        unsorted = true;    // reset
    }

    // x,y lists are now sorted
    // process options
    std::string method = "linear";
    bool extrap = true;     // manage "extrap"/"noextrap" options
    int extrapOpt = 0;      // manage scalar option also
    double extrapVal = std::numeric_limits<double>::quiet_NaN();
    int xiArgOffset = 0;    // interpolated yi output

    if (nargin > 2)
    {
        if (inputs[nargin - 1].IsString() && inputs[nargin - 1].StringVal() == "pp")
        {
            xiArgOffset = -1;   // PP struct output
            --nargin;
        }
    }

    if (nargin > 3 + xiArgOffset)
    {
        bool setExtrap = true;  // extrap can be set
        bool setMethod = true;  // method can be set

        try
        {
            interpOptionsHelper(eval, inputs[3 + xiArgOffset], extrap, method, setExtrap, setMethod);

            if (!setExtrap) // extrap has been set
            {
                if (extrap)
                    extrapOpt = 1;
                else
                    extrapOpt = -1;
            }
        }
        catch (OML_Error& omlerr)
        {
            if (inputs[3 + xiArgOffset].IsScalar())
            {
                extrapVal = inputs[3 + xiArgOffset].Scalar();
            }
            else
            {
                requireUniqueX = false;
                throw OML_Error(omlerr.GetErrorMessage());
            }
        }

        if (nargin > 4 + xiArgOffset)
        {
            try
            {
                interpOptionsHelper(eval, inputs[4 + xiArgOffset], extrap, method, setExtrap, setMethod);

                if (!setExtrap) // extrap has been set
                {
                    if (extrap)
                        extrapOpt = 1;
                    else
                        extrapOpt = -1;
                }
            }
            catch (OML_Error& omlerr)
            {
                if (inputs[4 + xiArgOffset].IsScalar())
                {
                    extrapVal = inputs[4 + xiArgOffset].Scalar();
                }
                else
                {
                    requireUniqueX = false;
                    throw OML_Error(omlerr.GetErrorMessage());
                }
            }
        }
    }

    if (inputs[2].IsMatrix() || inputs[2].IsScalar())
    {
        // interpolated yi output
        const hwMatrix* xi = inputs[2].ConvertToMatrix();

        // make sure xi is a column when y is not a vector
        if (xi->N() > 1)
        {
            if ((inputs[1].IsMatrix() && inputs[1].Matrix()->M() != 1 && inputs[1].Matrix()->N() != 1) ||
                inputs[1].IsNDMatrix())
            {
                hwMatrix* col = eval.allocateMatrix(xi);
                col->Transpose();

                std::vector<Currency> inputs3;
                inputs3.push_back(inputs[0]);   // push_back(sorted_x)
                inputs3.push_back(inputs[1]);   // push_back(reordered_y)
                inputs3.push_back(col);         // push_back(col_xi)

                for (int i = 3; i < nargin; ++i)
                {
                    inputs3.push_back(inputs[i]);
                }

                return OmlInterp1(eval, inputs3, outputs);
            }
        }

        // process yi options: vector, 2D or ND
        if (inputs[1].IsMatrix() || inputs[1].IsScalar())
        {
            const hwMatrix* y = inputs[1].ConvertToMatrix();

            if (y->IsVector())
            {
                hwMatrix* yi = eval.allocateMatrix();

                if (method == "linear")
                {
                    hwMathStatus status = LinearInterp(*x, *y, *xi, *yi, requireUniqueX, extrapOpt, extrapVal);
                    requireUniqueX = false;
                    BuiltInFuncsUtils::CheckMathStatus(eval, status);
                }
                else if (method == "pchip")
                {
                    BuiltInFuncsUtils::CheckMathStatus(eval, PchipInterp(*x, *y, *xi, *yi, extrapOpt, extrapVal));
                }
                else if (method == "spline")
                {
                    BuiltInFuncsUtils::CheckMathStatus(eval, Spline(*x, *y, *xi, *yi, extrapOpt, extrapVal));
                }
                else
                {
                    requireUniqueX = false;
                    throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_NOTIMPLEMENT));
                }

                outputs.push_back(yi);
            }
            else
            {
                if (y->N() == 0 && x->Size() != y->M())
                {
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
                }

                hwMatrix* yi = eval.allocateMatrix(xi->Size(), y->N(), true);

                if (y->N())
                {
                    std::unique_ptr<hwMatrix> ycol(EvaluatorInterface::allocateMatrix());
                    std::unique_ptr<hwMatrix> yicol(EvaluatorInterface::allocateMatrix());
                    hwMathStatus status;

                    for (int i = 0; i < y->N(); ++i)
                    {
                        BuiltInFuncsUtils::CheckMathStatus(eval, y->ReadColumn(i, *ycol));

                        if (method == "linear")
                        {
                            status = LinearInterp(*x, *ycol, *xi, *yicol, requireUniqueX, extrapOpt, extrapVal);
                            requireUniqueX = false;
                        }
                        else if (method == "pchip")
                        {
                            status = PchipInterp(*x, *ycol, *xi, *yicol, extrapOpt, extrapVal);
                        }
                        else if (method == "spline")
                        {
                            status = Spline(*x, *ycol, *xi, *yicol, extrapOpt, extrapVal);
                        }
                        else
                        {
                            requireUniqueX = false;
                            throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_NOTIMPLEMENT));
                        }

                        if (!status.IsOk())
                        {
                            if (status == HW_MATH_ERR_ARRAYSIZE)
                            {
                                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
                            }
                            else
                            {
                                BuiltInFuncsUtils::CheckMathStatus(eval, status);
                            }
                        }

                        BuiltInFuncsUtils::CheckMathStatus(eval, yi->WriteColumn(i, *yicol));
                    }
                }

                outputs.push_back(yi);
            }
        }
        else // if (inputs[1].IsNDMatrix())
        {
            return oml_MatrixNUtil4(eval, inputs, outputs, OmlInterp1, -1, 2);
        }
    }
    else if (xiArgOffset == -1)
    {
        const hwMatrix* x_old = x;

        // PP struct output
        std::vector<Currency> inputs2;
        inputs2.push_back(EvaluatorInterface::allocateMatrix(x_old));
        int pieces = x_old->Size() - 1;
        hwMathStatus status;

        if (inputs[1].IsMatrix() || inputs[1].IsScalar())
        {
            const hwMatrix* y_old = inputs[1].ConvertToMatrix();

            if (y_old->IsVector())
            {
                std::unique_ptr<hwMatrix> coefs(EvaluatorInterface::allocateMatrix());

                if (method == "linear")
                {
                    status = LinearInterp(*x_old, *y_old, *coefs);
                    requireUniqueX = false;
                }
                else if (method == "pchip")
                {
                    status = PchipInterp(*x_old, *y_old, *coefs);
                }
                else if (method == "spline")
                {
                    status = Spline(*x_old, *y_old, *coefs);
                }
                else
                {
                    requireUniqueX = false;
                    throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_NOTIMPLEMENT));
                }

                if (!status.IsOk())
                {
                    if (status == HW_MATH_ERR_ARRAYSIZE)
                    {
                        throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
                    }
                    else
                    {
                        BuiltInFuncsUtils::CheckMathStatus(eval, status);
                    }
                }

                inputs2.push_back(coefs.release());
            }
            else
            {
                int dim = y_old->N();
                std::vector<int> dims(3);

                dims[0] = dim;
                dims[1] = pieces;

                if (method == "linear")
                {
                    dims[2] = 2;
                }
                else    // cubic
                {
                    dims[2] = 4;
                }

                hwMatrixN* coefs = eval.allocateMatrixN(dims, true);
                hwMatrix coefs2D(dim * pieces, dims[2], coefs->GetRealData(), hwMatrix::REAL);

                if (dim)
                {
                    std::unique_ptr<hwMatrix> yoldcol(EvaluatorInterface::allocateMatrix());
                    std::unique_ptr<hwMatrix> ppcoefs(EvaluatorInterface::allocateMatrix());

                    for (int i = 0; i < dim; ++i)
                    {
                        BuiltInFuncsUtils::CheckMathStatus(eval, y_old->ReadColumn(i, *yoldcol));

                        if (method == "linear")
                        {
                            status = LinearInterp(*x_old, *yoldcol, *ppcoefs);
                            requireUniqueX = false;
                        }
                        else if (method == "pchip")
                        {
                            status = PchipInterp(*x_old, *yoldcol, *ppcoefs);
                        }
                        else if (method == "spline")
                        {
                            status = Spline(*x_old, *yoldcol, *ppcoefs);
                        }
                        else
                        {
                            requireUniqueX = false;
                            throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_NOTIMPLEMENT));
                        }

                        if (!status.IsOk())
                        {
                            if (status == HW_MATH_ERR_ARRAYSIZE)
                            {
                                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
                            }
                            else
                            {
                                BuiltInFuncsUtils::CheckMathStatus(eval, status);
                            }
                        }

                        for (int j = 0; j < pieces; ++j)
                        {
                            hwMatrix coefrow;
                            BuiltInFuncsUtils::CheckMathStatus(eval, ppcoefs->ReadRow(j, coefrow));
                            BuiltInFuncsUtils::CheckMathStatus(eval, coefs2D.WriteRow(i + j * dim, coefrow));
                        }
                    }

                    inputs2.push_back(coefs);
                    inputs2.push_back(dim);
                }
            }
        }
        else if (inputs[1].IsNDMatrix())
        {
            // rework vector reading
            const hwMatrixN* y_old = inputs[1].MatrixN();
            const std::vector<int>& dims_old = y_old->Dimensions();
            int numDims = static_cast<int> (dims_old.size());
            int dim = y_old->Size() / dims_old[0];
            hwMatrix y_old2D(dims_old[0], dim, (void*)y_old->GetRealData(), hwMatrix::REAL);

            std::vector<int> ppdims(numDims + 1);
            hwMatrix* ppDims = eval.allocateMatrix(1, numDims - 1, true);

            for (int i = 0; i < numDims - 1; ++i)
            {
                ppdims[i] = dims_old[i + 1];
                (*ppDims)(i) = dims_old[i + 1];
            }

            ppdims[numDims - 1] = pieces;

            if (method == "linear")
            {
                ppdims[numDims] = 2;
            }
            else    // cubic
            {
                ppdims[numDims] = 4;
            }

            hwMatrixN* coefs = eval.allocateMatrixN(ppdims, true);
            hwMatrix coefs2D(dim * pieces, ppdims[numDims], coefs->GetRealData(), hwMatrix::REAL);

            if (dim)
            {
                std::unique_ptr<hwMatrix> yoldcol(EvaluatorInterface::allocateMatrix());
                std::unique_ptr<hwMatrix> ppcoefs(EvaluatorInterface::allocateMatrix());
                hwMathStatus status;

                for (int i = 0; i < dim; ++i)
                {
                    BuiltInFuncsUtils::CheckMathStatus(eval, y_old2D.ReadColumn(i, *yoldcol));

                    if (method == "linear")
                    {
                        status = LinearInterp(*x_old, *yoldcol, *ppcoefs);
                        requireUniqueX = false;
                    }
                    else if (method == "pchip")
                    {
                        status = PchipInterp(*x_old, *yoldcol, *ppcoefs);
                    }
                    else if (method == "spline")
                    {
                        status = Spline(*x_old, *yoldcol, *ppcoefs);
                    }
                    else
                    {
                        requireUniqueX = false;
                        throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_NOTIMPLEMENT));
                    }

                    if (!status.IsOk())
                    {
                        if (status == HW_MATH_ERR_ARRAYSIZE)
                        {
                            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
                        }
                        else
                        {
                            BuiltInFuncsUtils::CheckMathStatus(eval, status);
                        }
                    }

                    for (int j = 0; j < pieces; ++j)
                    {
                        hwMatrix coefrow;
                        BuiltInFuncsUtils::CheckMathStatus(eval, ppcoefs->ReadRow(j, coefrow));
                        BuiltInFuncsUtils::CheckMathStatus(eval, coefs2D.WriteRow(i + j * dim, coefrow));
                    }
                }
            }

            inputs2.push_back(coefs);
            inputs2.push_back(ppDims);
        }

        return OmlMakePPoly(eval, inputs2, outputs);
    }
    else
    {
        throw OML_Error(OML_ERR_REALVECTOR, 3, OML_VAR_TYPE);   // needs update
    }

    return true;
}
//------------------------------------------------------------------------------
// Helper method for interp1 command
//------------------------------------------------------------------------------
void interpOptionsHelper(EvaluatorInterface& eval,
                         const Currency&     input,
                         bool&               extrap,
                         std::string&        method,
                         bool&               setExtrap,
                         bool&               setMethod)
{
    // only used for interp1, could be used for interp2 also
    if (input.IsString())
    {
        std::string str = readOption(eval, input);
        if (str == "linear")
        {
            if (!setMethod)
                throw OML_Error(HW_ERROR_NOTSETMETHODMOREONCE);
            method = str;
            setMethod = false;
        }
        else if (str == "pchip")
        {
            if (!setMethod)
                throw OML_Error(HW_ERROR_NOTSETMETHODMOREONCE);
            method = str;
            setMethod = false;
        }
        else if (str == "spline")
        {
            if (!setMethod)
                throw OML_Error(HW_ERROR_NOTSETMETHODMOREONCE);
            method = str;
            setMethod = false;
        }
        else if (str == "extrap")
        {
            if (!setExtrap)
                throw OML_Error(HW_ERROR_NOTSETEXTRAPMOREONCE);
            extrap = true;
            setExtrap = false;
        }
        else if (str == "noextrap")
        {
            if (!setExtrap)
                throw OML_Error(HW_ERROR_NOTSETEXTRAPMOREONCE);
            extrap = false;
            setExtrap = false;
        }
        else
        {
            throw OML_Error(HW_ERROR_INVALIDOPTION(str));
        }
    }
    else if (input.IsScalar() && setExtrap)
    {
        throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_NOTIMPLEMENT));
    }
    else
    {
        throw OML_Error(OML_ERR_STRING);
    }
}
//------------------------------------------------------------------------------
// Interpolates in two-dimensions and returns true
//------------------------------------------------------------------------------
bool OmlInterp2(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 5 || nargin > 7)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);

    if (!inputs[1].IsMatrix())
        throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);

    if (!inputs[2].IsMatrix())
        throw OML_Error(OML_ERR_MATRIX, 3, OML_VAR_DATA);

    if (!inputs[3].IsMatrix() && !inputs[3].IsScalar())
        throw OML_Error(OML_ERR_SCALARMATRIX, 4, OML_VAR_DATA);

    if (!inputs[4].IsMatrix() && !inputs[4].IsScalar())
        throw OML_Error(OML_ERR_SCALARMATRIX, 5, OML_VAR_DATA);

    const hwMatrix* x_old = inputs[0].ConvertToMatrix();
    const hwMatrix* y_old = inputs[1].ConvertToMatrix();
    const hwMatrix* z_old = inputs[2].ConvertToMatrix();
    const hwMatrix* x_new = inputs[3].ConvertToMatrix();
    const hwMatrix* y_new = inputs[4].ConvertToMatrix();
    std::string method;
    std::string extrapStr;
    int extrapOpt = 0;
    double extrapVal = std::numeric_limits<double>::quiet_NaN();

    if (nargin > 5)
    {
        if (!inputs[5].IsString())
            throw OML_Error(OML_ERR_STRING, 6, OML_VAR_DATA);

        method = inputs[5].StringVal();

        if (nargin == 6 && method == "spline")
        {
            extrapOpt = 1;
        }
    }

    if (nargin > 6)
    {
        if (inputs[6].IsString())
        {
            extrapStr = inputs[6].StringVal();

            if (extrapStr == "extrap")
            {
                extrapOpt = 1;
            }
            else if (extrapStr == "noextrap")
            {
                extrapOpt = -1;
            }
            else
            {
                throw OML_Error(OML_ERR_BAD_STRING, 7, OML_VAR_DATA);
            }
        }
        else if (inputs[6].IsScalar())
        {
            extrapVal = inputs[6].Scalar();
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARSTRING, 7, OML_VAR_DATA);
        }
    }

    std::unique_ptr<hwMatrix> z_new(EvaluatorInterface::allocateMatrix());

    if (method.empty() || method == "linear")
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, BilinearInterp(*x_old, *y_old, *z_old,
                                           *x_new, *y_new, *z_new, extrapOpt, extrapVal));
    }
    else if (method == "spline")
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, Spline2D(*x_old, *y_old, *z_old,
                                           *x_new, *y_new, *z_new, extrapOpt, extrapVal));
    }
    else
    {
        throw OML_Error(OML_ERR_OPTIONVAL, 6, OML_VAR_STRING);
    }

    outputs.push_back(z_new.release());
    return true;
}
//------------------------------------------------------------------------------
// Interpolates in N dimensions and returns true
//------------------------------------------------------------------------------
bool OmlInterpN(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs)
{
    // check inputs
    int nargin = static_cast<int> (inputs.size());

    if (nargin < 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    int numMatrixArgs = nargin;

    if (nargin > 2)
    {
        if (inputs[nargin - 2].IsString())
        {
            numMatrixArgs -= 2;
        }
        else if (inputs[nargin - 1].IsString())
        {
            --numMatrixArgs;
        }
    }

    if (numMatrixArgs < 3)
    {
        int k = 1;

        if (numMatrixArgs == 2)
        {
            if (!inputs[1].IsPositiveInteger())
            {
                throw OML_Error(OML_ERR_POSINTEGER, 2);
            }

            k = static_cast<int> (inputs[1].Scalar());
            k = (1 << k) - 1;
            // numMatrixArgs = 1;
        }

        double inc = 1.0 / (1.0 + k);
        std::vector<Currency> inputs2;
        inputs2.push_back(inputs[0]);

        if (inputs[0].IsMatrix())
        {
            const hwMatrix* val_old = inputs[0].Matrix();
            int xv1_new_size = val_old->M() + (val_old->M() - 1) * k;
            int xv2_new_size = val_old->N() + (val_old->N() - 1) * k;
            hwMatrix* xv1_new = EvaluatorInterface::allocateMatrix(xv1_new_size, 1, true);
            hwMatrix* xv2_new = EvaluatorInterface::allocateMatrix(xv2_new_size, 1, true);

            for (int j = 0; j < xv1_new_size; ++j)
                (*xv1_new)(j) = 1 + j * inc;

            for (int j = 0; j < xv2_new_size; ++j)
                (*xv2_new)(j) = 1 + j * inc;

            inputs2.push_back(xv1_new);
            inputs2.push_back(xv2_new);
        }
        else if (inputs[0].IsNDMatrix())
        {
            const hwMatrixN* val_old = inputs[0].MatrixN();
            const std::vector<int>& dims = val_old->Dimensions();

            for (int i = 0; i < dims.size(); ++i)
            {
                int xv_new_size = dims[i] + (dims[i] - 1) * k;
                hwMatrix* xv_old = EvaluatorInterface::allocateMatrix(xv_new_size, 1, true);

                for (int j = 0; j < xv_new_size; ++j)
                    (*xv_old)(j) = 1 + j * inc;

                inputs2.push_back(xv_old);
            }
        }

        for (int i = numMatrixArgs; i < nargin; ++i)
        {
            inputs2.push_back(inputs[i]);
        }

        return OmlInterpN(eval, inputs2, outputs);
    }

    // assess whether inputs[0] is a dimension argument or V.
    // If inputs[0] is a dimension argument, determine whether the
    // dimensions are vectors or ndgrid matrices.
    bool implicitXvecs = false; // inputs[0] is not V
    bool vectorXYargs = false;  // ndgrid dimension arguments
    int numDims = (numMatrixArgs - 1) / 2;

    if (numMatrixArgs % 2 == 0)
    {
        implicitXvecs = true;
        numDims = numMatrixArgs - 1;
    }
    else
    {
        if (inputs[0].IsMatrix())
        {
            if (inputs[0].Matrix()->IsVector())
            {
                vectorXYargs = true;
            }
            else if (inputs[numDims].IsMatrix())
            {
                const hwMatrix* temp1 = inputs[0].Matrix();
                const hwMatrix* temp2 = inputs[numDims].Matrix();

                if (temp1->M() != temp2->M() || temp1->N() != temp2->N())
                {
                    implicitXvecs = true;
                    numDims = numMatrixArgs - 1;
                }
            }
            else if (inputs[numDims].IsNDMatrix())
            {
                implicitXvecs = true;
                numDims = numMatrixArgs - 1;
            }
            else
            {
                throw OML_Error(OML_ERR_REALMATRIX, numDims + 1);
            }
        }
        else if (inputs[0].IsNDMatrix())
        {
            if (inputs[numDims].IsMatrix())
            {
                implicitXvecs = true;   // inputs[0] is V
            }
            else if (inputs[numDims].IsNDMatrix())
            {
                const hwMatrixN* temp1 = inputs[0].MatrixN();
                const hwMatrixN* temp2 = inputs[numDims].MatrixN();

                if (temp1->Dimensions() != temp2->Dimensions())
                {
                    implicitXvecs = true;   // inputs[0] is V
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REALMATRIX, numDims + 1);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_REALMATRIX, 1);
        }
    }

    if (implicitXvecs)
    {
        std::vector<Currency> inputs2(numMatrixArgs - 1 + nargin);

        if (inputs[0].IsMatrix())
        {
            const hwMatrix* val_old = inputs[0].Matrix();
            hwMatrix* xv1_old = EvaluatorInterface::allocateMatrix(val_old->M(), 1, true);
            hwMatrix* xv2_old = EvaluatorInterface::allocateMatrix(val_old->N(), 1, true);

            for (int j = 0; j < val_old->M(); ++j)
                (*xv1_old)(j) = 1 + j;

            for (int j = 0; j < val_old->N(); ++j)
                (*xv2_old)(j) = 1 + j;

            inputs2[0] = xv1_old;
            inputs2[1] = xv2_old;
        }
        else if (inputs[0].IsNDMatrix())
        {
            const hwMatrixN* val_old = inputs[0].MatrixN();
            const std::vector<int>& dims = val_old->Dimensions();

            for (int i = 0; i < dims.size(); ++i)
            {
                hwMatrix* xv_old = EvaluatorInterface::allocateMatrix(dims[i], 1, true);

                for (int j = 0; j < dims[i]; ++j)
                    (*xv_old)(j) = 1 + j;

                inputs2[i] = xv_old;
            }
        }

        for (int i = 0; i < nargin; ++i)
        {
            inputs2[numMatrixArgs - 1 + i] = inputs[i];
        }

        try
        {
            OmlInterpN(eval, inputs2, outputs);
        }
        catch (OML_Error& err)
        {
            int arg1 = err.Arg1();
            err.Arg1(arg1 - (numMatrixArgs - 1));
            throw;
        }

        return true;
    }

    if (!vectorXYargs)
    {
        // matrix XY args with ndgrid format
        // convert to vector format
        std::vector<Currency> inputs2(nargin);
        int xold_stride = 1;

        if (numDims == 2)
        {
            if (!inputs[numDims].IsMatrix())
            {
                throw OML_Error(OML_ERR_REALMATRIX, numDims + 1, OML_VAR_DATA);
            }

            int val_old_m = inputs[numDims].Matrix()->M();
            int val_old_n = inputs[numDims].Matrix()->N();

            for (int i = 0; i < numDims; ++i)
            {
                if (!inputs[i].IsMatrix())
                {
                    throw OML_Error(OML_ERR_REALMATRIX, i + 1, OML_VAR_DATA);
                }

                const hwMatrix* xn_old = inputs[i].Matrix();

                if (!xn_old->IsReal())
                {
                    throw OML_Error(OML_ERR_REALMATRIX, i + 1, OML_VAR_DATA);
                }

                if (xn_old->M() != val_old_m || xn_old->N() != val_old_n)
                {
                    throw OML_Error(OML_ERR_ARRAYSIZE, i + 1, numDims + 1, OML_VAR_DATA);
                }

                int xold_size = (!i) ? val_old_m : val_old_n;

                hwMatrix* xv_old = EvaluatorInterface::allocateMatrix(xold_size, 1, true);

                inputs2[i] = xv_old;

                for (int j = 0; j < xold_size; ++j)
                {
                    (*xv_old)(j) = (*xn_old)(j * xold_stride);
                }

                xold_stride *= xold_size;
            }
        }
        else
        {
            if (!inputs[numDims].IsNDMatrix())
            {
                throw OML_Error(OML_ERR_MATRIX, numDims + 1, OML_VAR_DATA);
            }

            const std::vector<int>& val_old_dims = inputs[numDims].MatrixN()->Dimensions();

            for (int i = 0; i < numDims; ++i)
            {
                if (!inputs[i].IsNDMatrix())
                {
                    throw OML_Error(OML_ERR_MATRIX, i + 1, OML_VAR_DATA);
                }

                const hwMatrixN* xn_old = inputs[i].MatrixN();

                if (xn_old->Dimensions() != val_old_dims)
                {
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, i + 1, OML_VAR_DATA);
                }

                if (!xn_old->IsReal())
                {
                    throw OML_Error(OML_ERR_REALMATRIX, i + 1, OML_VAR_DATA);
                }

                int xold_size = xn_old->Dimensions()[i];

                hwMatrix* xv_old = EvaluatorInterface::allocateMatrix(xold_size, 1, true);

                inputs2[i] = xv_old;

                for (int j = 0; j < xold_size; ++j)
                {
                    (*xv_old)(j) = (*xn_old)(j * xold_stride);
                }

                xold_stride *= xold_size;
            }
        }

        inputs2[numDims] = inputs[numDims];

        for (int i = numDims + 1; i < nargin; ++i)
            inputs2[i] = inputs[i];

        return OmlInterpN(eval, inputs2, outputs);
    }

    bool vectorXYIargs = false;  // ndgrid dimension arguments

    if (inputs[numDims + 1].IsMatrix())
    {
        if (inputs[numDims + 1].Matrix()->IsVector())
        {
            vectorXYIargs = true;
        }
    }
    else if (inputs[numDims + 1].IsScalar())
    {
        vectorXYIargs = true;
    }

    if (!vectorXYIargs)
    {
        // matrix XY args with ndgrid format
        // convert to vector format
        std::vector<Currency> inputs2(nargin);
        int xnew_stride = 1;

        if (numDims == 2)
        {
            if (!inputs[numDims + 1].IsMatrix())
            {
                throw OML_Error(OML_ERR_REALMATRIX, numDims + 2, OML_VAR_DATA);
            }

            int val_old_m = inputs[numDims].Matrix()->M();
            int val_old_n = inputs[numDims].Matrix()->N();
            int val_new_m = inputs[numDims + 1].Matrix()->M();
            int val_new_n = inputs[numDims + 1].Matrix()->N();

            for (int i = 0; i < numDims; ++i)
            {
                if (!inputs[numDims + 1 + i].IsMatrix())
                {
                    throw OML_Error(OML_ERR_REALMATRIX, 1, numDims + 2 + i, OML_VAR_DATA);
                }

                const hwMatrix* xn_old = inputs[i].Matrix();
                const hwMatrix* xn_new = inputs[numDims + 1 + i].Matrix();

                if (i && !xn_new->IsReal())
                {
                    throw OML_Error(OML_ERR_REALMATRIX, numDims + 2 + i, OML_VAR_DATA);
                }

                if (i && (xn_new->M() != val_new_m || xn_new->N() != val_new_n))
                {
                    throw OML_Error(OML_ERR_ARRAYSIZE, numDims + 1, numDims + i + 1, OML_VAR_DATA);
                }

                int xnew_size = (!i) ? val_new_m : val_new_n;

                hwMatrix* xv_new = EvaluatorInterface::allocateMatrix(xnew_size, 1, true);

                inputs2[i] = inputs[i];
                inputs2[numDims + 1 + i] = xv_new;

                for (int j = 0; j < xnew_size; ++j)
                {
                    (*xv_new)(j) = (*xn_new)(j * xnew_stride);
                }

                xnew_stride *= xnew_size;
            }
        }
        else
        {
            if (!inputs[numDims].IsNDMatrix())
            {
                throw OML_Error(OML_ERR_MATRIX, numDims + 1, OML_VAR_DATA);
            }

            if (!inputs[numDims + 1].IsNDMatrix())
            {
                throw OML_Error(OML_ERR_MATRIX, numDims + 2, OML_VAR_DATA);
            }

            const std::vector<int>& val_old_dims = inputs[numDims].MatrixN()->Dimensions();
            const std::vector<int>& val_new_dims = inputs[numDims + 1].MatrixN()->Dimensions();

            for (int i = 0; i < numDims; ++i)
            {
                if (!inputs[numDims + 1 + i].IsNDMatrix())
                {
                    throw OML_Error(OML_ERR_MATRIX, 1, numDims + 2 + i, OML_VAR_DATA);
                }

                const hwMatrixN* xn_old = inputs[i].MatrixN();
                const hwMatrixN* xn_new = inputs[numDims + 1 + i].MatrixN();

                if (i && xn_new->Dimensions() != val_new_dims)
                {
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, numDims + 2 + i, OML_VAR_DATA);
                }

                if (i && !xn_new->IsReal())
                {
                    throw OML_Error(OML_ERR_REALMATRIX, numDims + 2 + i, OML_VAR_DATA);
                }

                int xnew_size = xn_new->Dimensions()[i];

                hwMatrix* xv_new = EvaluatorInterface::allocateMatrix(xnew_size, 1, true);

                inputs2[i] = inputs[i];
                inputs2[numDims + 1 + i] = xv_new;

                for (int j = 0; j < xnew_size; ++j)
                {
                    (*xv_new)(j) = (*xn_new)(j * xnew_stride);
                }

                xnew_stride *= xnew_size;
            }
        }

        inputs2[numDims] = inputs[numDims];

        for (int i = numMatrixArgs; i < nargin; ++i)
            inputs2[i] = inputs[i];

        return OmlInterpN(eval, inputs2, outputs);
    }

    std::string method = "linear";
    std::string argStr;
    int extrapOpt = 0;
    double extrapVal = std::numeric_limits<double>::quiet_NaN();

    if (nargin > 2)
    {
        if (inputs[nargin - 2].IsString())
        {
            method = inputs[nargin - 2].StringVal();
        }

        if (inputs[nargin - 1].IsString())
        {
            argStr = inputs[nargin - 1].StringVal();

            if (argStr == "extrap")
            {
                extrapOpt = 1;
            }
            else if (argStr == "noextrap")
            {
                extrapOpt = -1;
            }
            else if (!inputs[nargin - 2].IsString())
            {
                method = argStr;
            }
            else
            {
                throw OML_Error(OML_ERR_BAD_STRING, nargin);
            }
        }
        else if (inputs[nargin - 1].IsScalar())
        {
            extrapVal = inputs[nargin - 1].Scalar();
        }
        else if (!inputs[nargin - 1].IsMatrix())
        {
            throw OML_Error(OML_ERR_SCALARSTRING, nargin);
        }
    }

    if (method != "linear")
    {
        throw OML_Error("Error: unsupported extrapolation method; only 'linear' is currently supported");
    }

    // proceed with calling interpolation function now that all input
    // argument variations have been handled
    const hwMatrix** x_old = new const hwMatrix*[numDims];
    const hwMatrix** x_new = new const hwMatrix*[numDims];
    // std::unique_ptr<hwMatrix**> x_old(new const hwMatrix* [numDims]);
    // std::unique_ptr<hwMatrix**> x_new(new const hwMatrix* [numDims]);

    for (int i = 0; i < numDims; ++i)
    {
        // vector checks
        if (inputs[i].IsMatrix())
        {
            if (!inputs[i].Matrix()->IsVector())
            {
                delete[] x_old;
                delete[] x_new;
                throw OML_Error(OML_ERR_REALVECTOR, i + 1, OML_VAR_DATA);
            }
        }

        if (inputs[i + 1 + numDims].IsMatrix())
        {
            if (!inputs[i + 1 + numDims].Matrix()->IsVector())
            {
                delete[] x_old;
                delete[] x_new;
                throw OML_Error(OML_ERR_REALVECTOR, numDims + i + 2, OML_VAR_DATA);
            }
        }
        else if (!inputs[i + 1 + numDims].IsScalar())
        {
            delete[] x_old;
            delete[] x_new;
            throw OML_Error(OML_ERR_REALVECTOR, numDims + i + 2, OML_VAR_DATA);
        }

        x_old[i] = inputs[i].Matrix();
        x_new[i] = inputs[i + 1 + numDims].ConvertToMatrix();
    }

    std::unique_ptr<hwMatrixN> val_new(EvaluatorInterface::allocateMatrixN());

    hwMathStatus status;

    if (inputs[numDims].IsNDMatrix())
    {
        status = MultilinearInterp(x_old, *inputs[numDims].MatrixN(), x_new,
                                   *val_new, extrapOpt, extrapVal);
    }
    else if (inputs[numDims].IsMatrix())
    {
        hwMatrixN* temp = EvaluatorInterface::allocateMatrixN();
        temp->Convert2DtoND(*inputs[numDims].Matrix(), false);
        status = MultilinearInterp(x_old, *temp, x_new, *val_new, extrapOpt, extrapVal);
    }
    else
    {
        delete[] x_old;
        delete[] x_new;
        throw OML_Error(OML_ERR_REALMATRIX, numDims + 1);
    }

    delete[] x_old;
    delete[] x_new;
/*
    if (status == HW_MATH_ERR_BADRANGE)
    {
        status.SetArg1(status.GetArg1() + 1);
    }
*/
    BuiltInFuncsUtils::CheckMathStatus(eval, status);

    outputs.push_back(val_new.release());

    return true;
}
//------------------------------------------------------------------------------
// Computes gradient of trilinear interpolation
//------------------------------------------------------------------------------
bool OmlIsoNormals(EvaluatorInterface           eval,
                   const std::vector<Currency>& inputs,
                   std::vector<Currency>&       outputs)
{
    bool switchSign = true;

    if (inputs.size() == 2 || inputs.size() == 3)
    {
        if (!inputs[0].IsNDMatrix())
        {
            throw OML_Error(OML_ERR_MATRIX, 4, OML_VAR_DATA);
        }

        if (!inputs[1].IsMatrix())
        {
            throw OML_Error(OML_ERR_MATRIX, 5, OML_VAR_DATA);
        }

        if (inputs.size() == 3)
        {
            if (!inputs[2].IsString())
            {
                throw OML_Error(OML_ERR_STRING, 3, OML_VAR_DATA);
            }

            std::string str = inputs[2].StringVal();

            if (str != "negate")
                throw OML_Error(OML_ERR_OPTION, 3);

            switchSign = false;
        }

        const hwMatrixN* val_old = inputs[0].MatrixN();
        const hwMatrix*  v_new   = inputs[1].Matrix();

        const std::vector<int>& dims_old = val_old->Dimensions();

        if (dims_old.size() != 3)
            throw OML_Error(OML_ERR_ARRAYSIZE, 1);

        int nx = dims_old[1];
        int ny = dims_old[0];
        int nz = dims_old[2];

        hwMatrix x_old(nx, hwMatrix::REAL);
        hwMatrix y_old(ny, hwMatrix::REAL);
        hwMatrix z_old(nz, hwMatrix::REAL);

        for (int i = 0; i < nx; ++i)
            x_old(i) = static_cast<double>(i + 1);

        for (int i = 0; i < ny; ++i)
            y_old(i) = static_cast<double>(i + 1);

        for (int i = 0; i < nz; ++i)
            z_old(i) = static_cast<double>(i + 1);

        std::unique_ptr<hwMatrix> grad(EvaluatorInterface::allocateMatrix());

        hwMathStatus status = TrilinearInterpGrad(x_old, y_old, z_old, *val_old,
            *v_new, *grad);

        if (!status.IsOk())
        {
            int arg = status.GetArg1();
            
            if (arg > 3)
                status.SetArg1(arg - 3);

            throw OML_Error(status);
        }

        if (switchSign)
        {
            (*grad) = -(*grad);
        }

        outputs.push_back(grad.release());

        return true;
    }

    if (inputs.size() != 5 && inputs.size() != 6)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs.size() == 6)
    {
        if (!inputs[5].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 6, OML_VAR_DATA);
        }

        std::string str = inputs[5].StringVal();

        if (str != "negate")
            throw OML_Error(OML_ERR_OPTION, 6);

        switchSign = false;
    }

    if (inputs[0].IsNDMatrix() && inputs[1].IsNDMatrix() && inputs[2].IsNDMatrix())
    {
        const hwMatrixN* x_old = inputs[0].MatrixN();
        const hwMatrixN* y_old = inputs[1].MatrixN();
        const hwMatrixN* z_old = inputs[2].MatrixN();

        if (x_old->Dimensions().size() != 3)
        {
            throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
        }

        if (x_old->Dimensions() != y_old->Dimensions())
        {
            throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
        }

        if (x_old->Dimensions() != z_old->Dimensions())
        {
            throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
        }

        std::vector<Currency> inputs2;

        std::vector<int> sliceDims2D(2);
        sliceDims2D[0] = -1;
        sliceDims2D[1] = 1;

        // x_old
        std::vector<hwSliceArg> sliceArgs;   // assume permuted meshgrid inputs
        sliceArgs.push_back(0);
        sliceArgs.push_back(hwSliceArg());
        sliceArgs.push_back(0);

        hwMatrixN slice;
        hwMatrix* slice2D = new hwMatrix;
        x_old->SliceRHS(sliceArgs, slice);
        slice.Reshape(sliceDims2D);     // reshape to a column
        slice.ConvertNDto2D(*slice2D);
        inputs2.push_back(slice2D);

        // y_old
        sliceArgs.clear();
        sliceArgs.push_back(hwSliceArg());
        sliceArgs.push_back(0);
        sliceArgs.push_back(0);

        slice2D = new hwMatrix;
        y_old->SliceRHS(sliceArgs, slice);
        slice.Reshape(sliceDims2D);     // reshape to a column
        slice.ConvertNDto2D(*slice2D);
        inputs2.push_back(slice2D);

        // z_old
        sliceArgs.clear();
        sliceArgs.push_back(0);
        sliceArgs.push_back(0);
        sliceArgs.push_back(hwSliceArg());

        slice2D = new hwMatrix;
        z_old->SliceRHS(sliceArgs, slice);
        slice.Reshape(sliceDims2D);     // reshape to a column
        slice.ConvertNDto2D(*slice2D);
        inputs2.push_back(slice2D);

        for (int i = 3; i < inputs.size(); ++i)
        {
            inputs2.push_back(inputs[i]);
        }

        return OmlIsoNormals(eval, inputs2, outputs);
    }

    if (!inputs[0].IsMatrix())
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
    }

    if (!inputs[1].IsMatrix())
    {
        throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);
    }

    if (!inputs[2].IsMatrix())
    {
        throw OML_Error(OML_ERR_MATRIX, 3, OML_VAR_DATA);
    }

    if (!inputs[3].IsNDMatrix())
    {
        throw OML_Error(OML_ERR_MATRIX, 4, OML_VAR_DATA);
    }

    if (!inputs[4].IsMatrix())
    {
        throw OML_Error(OML_ERR_MATRIX, 5, OML_VAR_DATA);
    }

    const hwMatrix*  x_old   = inputs[0].Matrix();
    const hwMatrix*  y_old   = inputs[1].Matrix();
    const hwMatrix*  z_old   = inputs[2].Matrix();
    const hwMatrixN* val_old = inputs[3].MatrixN();
    const hwMatrix*  v_new   = inputs[4].Matrix();

    std::unique_ptr<hwMatrix> grad(EvaluatorInterface::allocateMatrix());

    hwMathStatus status = TrilinearInterpGrad(*x_old, *y_old, *z_old, *val_old,
                                              *v_new, *grad);

    BuiltInFuncsUtils::CheckMathStatus(eval, status);
    
    if (switchSign)
    {
        (*grad) = -(*grad);
    }

    outputs.push_back(grad.release());

    return true;
}
//------------------------------------------------------------------------------
// Computes polynomial division, or deconvolution
//------------------------------------------------------------------------------
bool OmlDeconv(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs,
               std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
    {
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);
    }

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar() && !inputs[1].IsComplex())
    {
        throw OML_Error(OML_ERR_VECTOR, 2, OML_VAR_DATA);
    }

    const hwMatrix* num = inputs[0].ConvertToMatrix();
    const hwMatrix* den = inputs[1].ConvertToMatrix();

    std::unique_ptr<hwMatrix> Q(EvaluatorInterface::allocateMatrix());
    std::unique_ptr<hwMatrix> R(EvaluatorInterface::allocateMatrix());

    BuiltInFuncsUtils::CheckMathStatus(eval, PolyDivide(*num, *den, *Q, *R));

    outputs.push_back(Q.release());
    outputs.push_back(R.release());

    return true;
}
//------------------------------------------------------------------------------
// Computes the derivative of a polynomial and returns true
//------------------------------------------------------------------------------
bool OmlPolyder(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();
    size_t nargout = eval.GetNargoutValue();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsVector() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);

    const hwMatrix* A = inputs[0].ConvertToMatrix();

    if (nargin == 1)
    {
        if (nargout > 1)
            throw OML_Error(HW_ERROR_NOMORE1OUTALLOWEDW1INP);

        std::unique_ptr<hwMatrix> D(EvaluatorInterface::allocateMatrix());
        BuiltInFuncsUtils::CheckMathStatus(eval, PolyDer(*A, *D));
        D->Transpose();
        outputs.push_back(D.release());
    }
    else // nargin == 2
    {
        if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
            throw OML_Error(OML_ERR_SCALARVECTOR, 2, OML_VAR_DATA);

        const hwMatrix* B = inputs[1].ConvertToMatrix();

        if (nargout < 2)
        {
            std::unique_ptr<hwMatrix> D(EvaluatorInterface::allocateMatrix());
            BuiltInFuncsUtils::CheckMathStatus(eval, PolyDer(*A, *B, *D));
            D->Transpose();
            outputs.push_back(D.release());
        }
        else if (nargout == 2)
        {
            std::unique_ptr<hwMatrix> P(EvaluatorInterface::allocateMatrix());
            std::unique_ptr<hwMatrix> Q(EvaluatorInterface::allocateMatrix());
            BuiltInFuncsUtils::CheckMathStatus(eval, PolyDer(*A, *B, *P, *Q));
            P->Transpose();
            Q->Transpose();
            outputs.push_back(P.release());
            outputs.push_back(Q.release());
        }
        else
        {
            throw OML_Error(OML_ERR_NUMARGOUT);
        }
    }

    return true;
}

//------------------------------------------------------------------------------
// Returns true and construct a piecewise polynomial 
//------------------------------------------------------------------------------
bool OmlMakePPoly(EvaluatorInterface           eval,
                  const std::vector<Currency>& inputs,
                  std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();
    size_t nargout = eval.GetNargoutValue();

    if (nargin < 2 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix())
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);

    const hwMatrix* input1 = inputs[0].Matrix();

    if (!input1->IsReal())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_DATA);

    if (!input1->IsVector())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_DATA);

    hwMatrix* breaks = EvaluatorInterface::allocateMatrix(input1);

    if (breaks->M() != 1)
        breaks->Transpose();

    int length = breaks->N() - 1;
    int order;

    Currency out = EvaluatorInterface::allocateStruct();
    StructData* pp = out.Struct();
    pp->SetValue(0, -1, "form", "pp");
    pp->SetValue(0, -1, "breaks", breaks);
    pp->SetValue(0, -1, "pieces", length);

    if (nargin == 2)
    {
        if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
        {
            if (inputs[1].IsNDMatrix())
                throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);
            else
                throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);
        }

        hwMatrix* coefs = EvaluatorInterface::allocateMatrix(inputs[1].ConvertToMatrix());

        if (coefs->M() != length)
        {
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DATA);
        }

        pp->SetValue(0, -1, "dim", 1);
        pp->SetValue(0, -1, "coefs", coefs);
        order = coefs->N();
    }
    else // nargin == 3
    {
        if (!inputs[1].IsNDMatrix())
        {
            if (inputs[1].IsMatrix() || inputs[1].IsScalar())
                throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);
            else
                throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);
        }

        hwMatrixN* coefs = EvaluatorInterface::allocateMatrixN(inputs[1].MatrixN());
        const std::vector<int>& dims = coefs->Dimensions();
        size_t numDims = dims.size();
        order = dims[numDims - 1];

        if (dims[numDims - 2] != length)
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DATA);

        if (numDims == 3)
        {
            if (!inputs[2].IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DATA);

            int dim = static_cast<int> (inputs[2].Scalar());

            if (dims[0] != dim)
                throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DATA);

            pp->SetValue(0, -1, "dim", dim);
        }
        else
        {
            if (!inputs[2].IsPositiveIntegralVector())
                throw OML_Error(OML_ERR_POSINTVECTOR, 3, OML_VAR_DATA);

            hwMatrix* dim = EvaluatorInterface::allocateMatrix(inputs[2].Matrix());

            if (dim->Size() != numDims - 2)
                throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DATA);

            for (int i = 0; i < dims.size() - 2; ++i)
            {
                double dim_int = (*dim)(i);

                if (dims[i] != static_cast<int> (dim_int))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DATA);
            }

            pp->SetValue(0, -1, "dim", dim);
        }

        std::vector<int> newdims(2);
        newdims[0] = -1;
        newdims[1] = order;
        coefs->Reshape(newdims);
        pp->SetValue(0, -1, "coefs", coefs);
    }

    pp->SetValue(0, -1, "order", order);

    outputs.push_back(out);

    return true;
}

//------------------------------------------------------------------------------
// Returns true and extracts details of a piecewise polynomial 
//------------------------------------------------------------------------------
bool OmlUnMakePPoly(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs)
{
    size_t nargout = eval.GetNargoutValue();

    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsStruct())
        throw OML_Error(OML_ERR_STRUCT, 1, OML_VAR_DATA);

    const StructData* pp     = inputs[0].Struct();
    const Currency& form_C   = pp->GetValue(0, -1, "form");
    const Currency& breaks_C = pp->GetValue(0, -1, "breaks");
    const Currency& coefs_C  = pp->GetValue(0, -1, "coefs");
    const Currency& pieces_C = pp->GetValue(0, -1, "pieces");
    const Currency& order_C  = pp->GetValue(0, -1, "order");
    const Currency& dim_C    = pp->GetValue(0, -1, "dim");

    if (!form_C.IsString())
    {
        throw OML_Error("Error: invalid struct; must be piecewise polynomial");
    }

    if (form_C.StringVal() != "pp")
    {
        throw OML_Error("Error: invalid struct; must be piecewise polynomial");
    }

    outputs.push_back(breaks_C);
    outputs.push_back(coefs_C);
    outputs.push_back(pieces_C);
    outputs.push_back(order_C);
    outputs.push_back(dim_C);

    return true;
}

//------------------------------------------------------------------------------
// Returns true and evaluates a piecewise polynomial 
//------------------------------------------------------------------------------
bool OmlPPolyEval(EvaluatorInterface           eval,
                  const std::vector<Currency>& inputs,
                  std::vector<Currency>&       outputs)
{
    size_t nargout = eval.GetNargoutValue();

    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsStruct())
        throw OML_Error(OML_ERR_STRUCT, 1, OML_VAR_DATA);

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar() && !inputs[1].IsNDMatrix())
        throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_DATA);

    const StructData* pp = inputs[0].Struct();
    const Currency& form_C = pp->GetValue(0, -1, "form");
    const Currency& breaks_C = pp->GetValue(0, -1, "breaks");
    const Currency& coefs_C = pp->GetValue(0, -1, "coefs");
    const Currency& pieces_C = pp->GetValue(0, -1, "pieces");
    const Currency& order_C = pp->GetValue(0, -1, "order");
    const Currency& dim_C = pp->GetValue(0, -1, "dim");

    if (!form_C.IsString())
    {
        throw OML_Error("Error: invalid struct; must be piecewise polynomial");
    }

    if (form_C.StringVal() != "pp")
    {
        throw OML_Error("Error: invalid struct; must be piecewise polynomial");
    }

    if (!breaks_C.IsMatrix())
    {
        throw OML_Error("Error: invalid struct; must be piecewise polynomial");
    }

    if (!coefs_C.IsMatrix())
    {
        throw OML_Error("Error: invalid struct; must be piecewise polynomial");
    }

    if (!pieces_C.IsPositiveInteger())
    {
        throw OML_Error("Error: invalid struct; must be piecewise polynomial");
    }

    if (!order_C.IsPositiveInteger())
    {
        throw OML_Error("Error: invalid struct; must be piecewise polynomial");
    }

    if (!dim_C.IsPositiveInteger() && !dim_C.IsPositiveIntegralVector())
    {
        throw OML_Error("Error: invalid struct; must be piecewise polynomial");
    }

    const hwMatrix* breaks = breaks_C.Matrix();
    const hwMatrix* coefs = coefs_C.Matrix();
    int pieces = static_cast<int> (pieces_C.Scalar());  // pieces per PP
    int order = static_cast<int> (order_C.Scalar());
    const hwMatrix* avDims = dim_C.ConvertToMatrix();   // av = array of vectors

    if (!breaks->IsReal())
        throw OML_Error("Error: invalid struct; pp breaks must be real");

    if (!breaks->IsVector())
        throw OML_Error("Error: invalid struct; pp breaks must be a vector");

    if (!coefs->IsReal())
        throw OML_Error("Error: invalid struct; pp coefs must be real");

    if (!avDims->IsReal())
        throw OML_Error("Error: invalid struct; pp dimensions must be integers");

    if (!avDims->IsVector())
        throw OML_Error("Error: invalid struct; pp dimensions must be a vector");

    int avSize = static_cast<int> ((*avDims)(0));

    for (int i = 1; i < avDims->Size(); ++i)
        avSize *= static_cast<int> ((*avDims)(i));

    int numPP = avSize * pieces;    // number of piecewise polynomials
    numPP = coefs->M();             // check this

    std::vector<int> yiDims;

    for (int i = 0; i < avDims->Size(); ++i)
    {
        yiDims.push_back(static_cast<int> ((*avDims)(i)));
    }

    int avNumDims = static_cast<int> (yiDims.size());
    const double* xiStart;
    int xiSize;

    if (inputs[1].IsMatrix() || inputs[1].IsScalar())
    {
        const hwMatrix* xi = inputs[1].ConvertToMatrix();

        if (!xi->IsReal())
            throw OML_Error("Error: invalid struct; must be piecewise polynomial");

        xiStart = xi->GetRealData();
        xiSize = xi->Size();

        if (xi->IsVector())
        {
            yiDims.push_back(xi->Size());
        }
        else
        {
            yiDims.push_back(xi->M());
            yiDims.push_back(xi->N());
        }
    }
    else if (inputs[1].IsNDMatrix())
    {
        const hwMatrixN* xi = inputs[1].MatrixN();

        if (!xi->IsReal())
            throw OML_Error("Error: invalid struct; must be piecewise polynomial");

        const std::vector<int>& xiDims = xi->Dimensions();
        xiStart = xi->GetRealData();
        xiSize = xi->Size();

        for (size_t i = 0; i < xiDims.size(); ++i)
        {
            yiDims.push_back(xiDims[i]);
        }
    }

    std::unique_ptr<hwMatrixN> yi(EvaluatorInterface::allocateMatrixN(yiDims, true));
    std::vector<int> yiIndex(yiDims.size());
    const double* breaks_data = breaks->GetRealData();
    int avStride = yi->Stride(avNumDims);

    for (int j = 0; j < xiSize; ++j)
    {
        // set the piece index for each xi value
        const double* coefs_start = coefs->GetRealData();
        int breakId = BinarySearch(breaks_data, pieces + 1, *xiStart);

        if (breakId < 0)
            breakId = 0;
        else if (breakId >= pieces)
            breakId = pieces - 1;

        for (int k = 0; k < avSize; ++k)
        {
            // set the matrix indices to the first index in each av slice
            int yiStart = yi->Index(yiIndex);
            double* yiData = yi->GetRealData() + yiStart;

            // perform the ppval op
            const double* coef = coefs_start + breakId * avStride;
            double x = (*xiStart) - breaks_data[breakId];
            *yiData = *coef;
            coef += coefs->M();
            int i = order - 1;

            do
            {
                *yiData = *yiData * x + *coef;
                coef += coefs->M();
            }
            while (--i);

            // advance av slice indices of yi
            for (int k = 0; k < avNumDims; ++k)
            {
                // increment index k if possible
                if (yiIndex[k] < yiDims[k] - 1)
                {
                    ++yiIndex[k];
                    break;
                }

                // index k is maxed out, so reset and continue to k+1
                yiIndex[k] = 0;
            }

            ++coefs_start;
        }

        // advance xi slice indices of yi
        for (size_t k = avNumDims; k < yiIndex.size(); ++k)
        {
            // increment index k if possible
            if (yiIndex[k] < yiDims[k] - 1)
            {
                ++yiIndex[k];
                break;
            }

            // index k is maxed out, so reset and continue to k+1
            yiIndex[k] = 0;
        }

        ++xiStart;
    }

    outputs.push_back(yi.release());

    return true;
}

//------------------------------------------------------------------------------
// Computes the integral of a polynomial and returns true
//------------------------------------------------------------------------------
bool OmlPolyint(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();
    size_t nargout = eval.GetNargoutValue();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargout > 1)
        throw OML_Error(OML_ERR_NUMARGOUT);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_SCALARVECTOR, 1, OML_VAR_DATA);

    const hwMatrix* P = inputs[0].ConvertToMatrix();
    std::unique_ptr<hwMatrix> Integral(EvaluatorInterface::allocateMatrix());

    if (nargin == 1)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, PolyInt(*P, *Integral));
        Integral->Transpose();
        outputs.push_back(Integral.release());
    }
    else // nargin == 2
    {
        if (!inputs[1].IsScalar())
            throw OML_Error(OML_ERR_SCALAR, 2, OML_VAR_VALUE);

        double k = inputs[1].Scalar();

        BuiltInFuncsUtils::CheckMathStatus(eval, PolyInt(*P, *Integral, k));
        Integral->Transpose();
        outputs.push_back(Integral.release());
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
