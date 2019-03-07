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
    eval.RegisterBuiltInFunction("spline", &OmlSpline,
        FunctionMetaData(4, 1, POLY));
    eval.RegisterBuiltInFunction("interp1", &OmlInterp1,
        FunctionMetaData(-4, 1, POLY));
    eval.RegisterBuiltInFunction("interp2", &OmlInterp2,
        FunctionMetaData(-6, 1, POLY));
    eval.RegisterBuiltInFunction("deconv", &OmlDeconv,
        FunctionMetaData(2, 2, POLY));
    eval.RegisterBuiltInFunction("polyder", &OmlPolyder,
        FunctionMetaData(-2, -2, POLY));
    eval.RegisterBuiltInFunction("polyint", &OmlPolyint,
        FunctionMetaData(-2, 1, POLY));
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
// Interpolates (x,y) data with a cubic spline and returns true
//------------------------------------------------------------------------------
bool OmlSpline(EvaluatorInterface           eval,
    const std::vector<Currency>& inputs,
    std::vector<Currency>&       outputs)
{
    static bool unsorted = true;
    size_t nargin = inputs.size();

    if (nargin != 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_TYPE);

    if (inputs[0].IsMatrix() && !inputs[0].Matrix()->IsVector())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_TYPE);

    if (!inputs[1].IsMatrix() && !inputs[1].IsNDMatrix() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_REALMATRIX, 2, OML_VAR_TYPE);

    if (!inputs[2].IsMatrix() && !inputs[2].IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 3, OML_VAR_TYPE);

    const hwMatrix* x_old = inputs[0].ConvertToMatrix();
    const hwMatrix* x_new = inputs[2].ConvertToMatrix();

    // enforce sorted x,y lists
    if (unsorted)
    {
        std::vector<Currency> inputs2;
        inputs2.push_back(inputs[0]);           // push_back(x)
        inputs2.push_back("ascending");

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
            bool            clamped = false;
            int             numPnts = x_old->Size();

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
            inputs2.push_back(inputs[2]);       // push_back(xi)

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

    // make sure x_new is a row when y is not a vector
    if (x_new->M() != 1)
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

    int n = x_old->Size();

    // process yi options: vector, 2D or ND
    if (inputs[1].IsMatrix() || inputs[1].IsScalar())
    {
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
            hwMatrix* y_new = eval.allocateMatrix(y_old->M(), x_new->Size(), hwMatrix::REAL);

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
                            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);
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
    else // if (inputs[1].IsNDMatrix())
    {
        const hwMatrixN* y = inputs[1].MatrixN();
        int numDims = static_cast<int>(y->Dimensions().size());

        return oml_MatrixNUtil4(eval, inputs, outputs, OmlSpline, -numDims, 2);
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
    size_t nargin = inputs.size();

    if (nargin < 3 || nargin > 5)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_TYPE);

    if (inputs[0].IsMatrix() && !inputs[0].Matrix()->IsVector())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_TYPE);

    if (!inputs[1].IsMatrix() && !inputs[1].IsNDMatrix() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_REALMATRIX, 2, OML_VAR_TYPE);

    if (!inputs[2].IsMatrix() && !inputs[2].IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 3, OML_VAR_TYPE);

    const hwMatrix* x = inputs[0].ConvertToMatrix();
    const hwMatrix* xi = inputs[2].ConvertToMatrix();

    // enforce sorted x,y lists
    if (unsorted)
    {
        std::vector<Currency> inputs2;
        inputs2.push_back(inputs[0]);           // push_back(x)
        inputs2.push_back("ascending");

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
                    inputs2.push_back(inputs[4]);
            }

            try
            {
                unsorted = false;
                bool retv = OmlInterp1(eval, inputs2, outputs); // call function on sorted_x, reordered_y
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

    // make sure xi is a column when y is not a vector
    if (xi->N() != 1)
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

            if (nargin > 3)
            {
                inputs3.push_back(inputs[3]);

                if (nargin > 4)
                    inputs3.push_back(inputs[4]);
            }

            return OmlInterp1(eval, inputs3, outputs);
        }
    }

    // process options
    std::string method = "linear";
    bool extrap = false;

    // work on sorted x,y lists
    if (nargin > 3)
    {
        bool setExtrap = true;
        bool setMethod = true;

        try
        {
            interpOptionsHelper(eval, inputs[3], extrap, method, setExtrap, setMethod);
        }
        catch (OML_Error& omlerr)
        {
            throw OML_Error(omlerr.GetErrorMessage());
        }

        if (nargin > 4)
        {
            try
            {
                interpOptionsHelper(eval, inputs[4], extrap, method, setExtrap, setMethod);
            }
            catch (OML_Error& omlerr)
            {
                throw OML_Error(omlerr.GetErrorMessage());
            }
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
                BuiltInFuncsUtils::CheckMathStatus(eval, LinearInterp(*x, *y, *xi, *yi, extrap));
            }
            else if (method == "pchip")
            {
                BuiltInFuncsUtils::CheckMathStatus(eval, PchipInterp(*x, *y, *xi, *yi, extrap));
            }
            else if (method == "spline")
            {
                BuiltInFuncsUtils::CheckMathStatus(eval, Spline(*x, *y, *xi, *yi, extrap));
            }
            else
            {
                throw OML_Error(HW_MATH_MSG_NOTIMPLEMENT);
            }

            outputs.push_back(yi);
        }
        else
        {
            hwMatrix* yi = eval.allocateMatrix(xi->Size(), y->N(), hwMatrix::REAL);

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
                        status = LinearInterp(*x, *ycol, *xi, *yicol, extrap);
                    }
                    else if (method == "pchip")
                    {
                        status = PchipInterp(*x, *ycol, *xi, *yicol, extrap);
                    }
                    else if (method == "spline")
                    {
                        status = Spline(*x, *ycol, *xi, *yicol, extrap);
                    }
                    else
                    {
                        throw OML_Error(HW_MATH_MSG_NOTIMPLEMENT);
                    }

                    if (!status.IsOk())
                    {
                        if (status == HW_MATH_ERR_ARRAYSIZE)
                        {
                            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);
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
    // only used for interp1, could be used for intper2 also
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
        throw OML_Error(HW_MATH_MSG_NOTIMPLEMENT);
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
    bool extrap;

    if (nargin > 5)
    {
        if (!inputs[5].IsString())
            throw OML_Error(OML_ERR_STRING, 6, OML_VAR_DATA);

        method = inputs[5].StringVal();
    }

    if (nargin > 6)
    {
        if (!inputs[6].IsString())
            throw OML_Error(OML_ERR_STRING, 7, OML_VAR_DATA);

        extrapStr = inputs[6].StringVal();
    }

    if (extrapStr == "extrap")
        extrap = true;
    else
        extrap = false;

    std::unique_ptr<hwMatrix> z_new(EvaluatorInterface::allocateMatrix());

    if (method.empty() || method == "linear")
        BuiltInFuncsUtils::CheckMathStatus(eval, BilinearInterp(*x_old, *y_old, *z_old, *x_new, *y_new, *z_new, extrap));
    else if (method == "spline")
        BuiltInFuncsUtils::CheckMathStatus(eval, Spline2D(*x_old, *y_old, *z_old, *x_new, *y_new, *z_new, extrap));
    else
        throw OML_Error(OML_ERR_OPTIONVAL, 6, OML_VAR_STRING);

    outputs.push_back(z_new.release());
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

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_SCALARVECTOR, 1, OML_VAR_DATA);

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
