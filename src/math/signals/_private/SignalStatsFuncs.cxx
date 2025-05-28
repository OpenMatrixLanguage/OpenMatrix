/**
* @file * @file SignalStatsFuncs.cxx
* @date February 2023
* Copyright (C) 2023 Altair Engineering, Inc.
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

#include "SignalStatsFuncs.h"
#include "hwMatrix.h"
#include "hwMatrixN.h"
#include "MKLutilities.h"
#include "FourierFuncs.h"

#define MKLuD MKLutilitiesD

//------------------------------------------------------------------------------
// Computes root sum squared values
//------------------------------------------------------------------------------
hwMathStatus RSSQ(const hwMatrix& data, int dim, hwMatrix& result)
{
    int m = data.M();
    int n = data.N();
    hwMathStatus status;

    if (dim == -1)
    {
        // use first non-singleton dimension
        if (data.M() == 1)
        {
            dim = 1;
        }
        else
        {
            dim = 0;
        }
    }
    else if (dim > 1)
    {
        status = result.Dimension(m, n, hwMatrix::REAL);
        result.SetElements(1.0);
        return status;
    }
    else if (dim < 0)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYDIM, 2);
    }

    std::vector<int> dims = { m, n };
    int numVecs = (dim == 0) ? n : m;
    int stride = (dim == 0) ? 1 : m;
    int row = 0;
    int col = 0;

    if (dim == 0)
    {
        status = result.Dimension(1, n, hwMatrix::REAL);
    }
    else
    {
        status = result.Dimension(m, 1, hwMatrix::REAL);
    }

    if (!status.IsOk())
        return status;

    for (int i = 0; i < numVecs; ++i)
    {
        // set the matrix indices to the first index in each vector
        int start = col * m + row;

        // perform op
        double ssq = 0.0;

        if (data.IsReal())
        {
            const double* real = data.GetRealData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                ssq += (*real) * (*real);
                real += stride;
            }
        }
        else    // complex
        {
            const hwComplex* cplx = data.GetComplexData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                ssq += cplx->MagSq();
                cplx += stride;
            }
        }

        result(row, col) = sqrt(ssq);

        // advance slice indices
        if (dim == 0)
            ++col;
        else
            ++row;
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes root sum squared values
//------------------------------------------------------------------------------
hwMathStatus RSSQ(const hwMatrixN& data, int dim, hwMatrixN& result)
{
    const std::vector<int>& dims = data.Dimensions();
    int numDim = static_cast<int> (dims.size());

    if (dim == -1)
    {
        // use first non-singleton dimension
        for (int i = 0; i < dims.size(); ++i)
        {
            if (dims[i] != 1)
            {
                dim = i;
                break;
            }
        }
    }
    else if (dim > numDim - 1)
    {
        result.Dimension(dims, hwMatrixN::REAL);
        result.SetElements(1.0);
        return hwMathStatus();
    }
    else if (dim < 0)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYDIM, 2);
    }

    int numVecs = data.Size() / dims[dim];
    int stride = data.Stride(dim);
    std::vector<int> matrixIndex(numDim);
    std::vector<int> outDims = dims;
    outDims[dim] = 1;

    result.Dimension(outDims, hwMatrixN::REAL);

    for (int i = 0; i < numVecs; ++i)
    {
        // set the matrix indices to the first index in each vector
        int start = data.Index(matrixIndex);

        // perform op
        double ssq = 0.0;

        if (data.IsReal())
        {
            const double* real = data.GetRealData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                ssq += (*real) * (*real);
                real += stride;
            }
        }
        else    // complex
        {
            const hwComplex* cplx = data.GetComplexData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                ssq += cplx->MagSq();
                cplx += stride;
            }
        }

        result(matrixIndex) = sqrt(ssq);

        // advance slice indices
        for (int j = 0; j < numDim; ++j)
        {
            if (j == dim)
                continue;

            // increment index j if possible
            if (matrixIndex[j] < static_cast<int> (dims[j]) - 1)
            {
                ++matrixIndex[j];
                break;
            }

            // index j is maxed out, so reset and continue to j+1
            matrixIndex[j] = 0;
        }
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes peak to rms values
//------------------------------------------------------------------------------
hwMathStatus Peak2RMS(const hwMatrix& data, int dim, hwMatrix& result)
{
    int m = data.M();
    int n = data.N();
    hwMathStatus status;

    if (dim == -1)
    {
        // use first non-singleton dimension
        if (data.M() == 1)
        {
            dim = 1;
        }
        else
        {
            dim = 0;
        }
    }
    else if (dim > 1)
    {
        status = result.Dimension(m, n, hwMatrix::REAL);
        result.SetElements(1.0);
        return status;
    }
    else if (dim < 0)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYDIM, 2);
    }

    std::vector<int> dims = { m, n };
    int numVecs = (dim == 0) ? n : m;
    int stride = (dim == 0) ? 1 : m;
    int row = 0;
    int col = 0;

    if (dim == 0)
    {
        status = result.Dimension(1, n, hwMatrix::REAL);
    }
    else
    {
        status = result.Dimension(m, 1, hwMatrix::REAL);
    }

    if (!status.IsOk())
        return status;

    for (int i = 0; i < numVecs; ++i)
    {
        // set the matrix indices to the first index in each vector
        int start = col * m + row;

        // perform op
        double peak = 0.0;
        double rms = 0.0;

        if (data.IsReal())
        {
            const double* real = data.GetRealData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                peak = _max(peak, fabs(*real));
                rms += (*real) * (*real);
                real += stride;
            }
        }
        else    // complex
        {
            const hwComplex* cplx = data.GetComplexData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                double magsq = cplx->MagSq();
                peak = _max(peak, magsq);
                rms += magsq;
                cplx += stride;
            }

            peak = sqrt(peak);
        }

        rms /= dims[dim];
        rms = sqrt(rms);
        result(row, col) = peak / rms;

        // advance slice indices
        if (dim == 0)
            ++col;
        else
            ++row;
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes peak to rms values
//------------------------------------------------------------------------------
hwMathStatus Peak2RMS(const hwMatrixN& data, int dim, hwMatrixN& result)
{
    const std::vector<int>& dims = data.Dimensions();
    int numDim = static_cast<int> (dims.size());

    if (dim == -1)
    {
        const std::vector<int>& dims = data.Dimensions();

        // use first non-singleton dimension
        for (int i = 0; i < dims.size(); ++i)
        {
            if (dims[i] != 1)
            {
                dim = i;
                break;
            }
        }
    }
    else if (dim > numDim - 1)
    {
        result.Dimension(dims, hwMatrixN::REAL);
        result.SetElements(1.0);
        return hwMathStatus();
    }
    else if (dim < 0)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYDIM, 2);
    }

    int numVecs = data.Size() / dims[dim];
    int stride = data.Stride(dim);
    std::vector<int> matrixIndex(numDim);
    std::vector<int> outDims = dims;
    outDims[dim] = 1;

    result.Dimension(outDims, hwMatrixN::REAL);

    for (int i = 0; i < numVecs; ++i)
    {
        // set the matrix indices to the first index in each vector
        int start = data.Index(matrixIndex);

        // perform op
        double peak = 0.0;
        double rms = 0.0;

        if (data.IsReal())
        {
            const double* real = data.GetRealData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                peak = _max(peak, fabs(*real));
                rms += (*real) * (*real);
                real += stride;
            }
        }
        else    // complex
        {
            const hwComplex* cplx = data.GetComplexData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                double magsq = cplx->MagSq();
                peak = _max(peak, magsq);
                rms += magsq;
                cplx += stride;
            }

            peak = sqrt(peak);
        }

        rms /= dims[dim];
        rms = sqrt(rms);
        result(matrixIndex) = peak / rms;

        // advance slice indices
        for (int j = 0; j < numDim; ++j)
        {
            if (j == dim)
                continue;

            // increment index j if possible
            if (matrixIndex[j] < static_cast<int> (dims[j]) - 1)
            {
                ++matrixIndex[j];
                break;
            }

            // index j is maxed out, so reset and continue to j+1
            matrixIndex[j] = 0;
        }
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes peak to peak values
//------------------------------------------------------------------------------
hwMathStatus Peak2Peak(const hwMatrix& data, int dim, hwMatrix& result)
{
    int m = data.M();
    int n = data.N();
    hwMathStatus status;

    if (dim == -1)
    {
        // use first non-singleton dimension
        if (data.M() == 1)
        {
            dim = 1;
        }
        else
        {
            dim = 0;
        }
    }
    else if (dim > 1)
    {
        status = result.Dimension(m, n, hwMatrix::REAL);
        result.SetElements(0.0);
        return status;
    }
    else if (dim < 0)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYDIM, 2);
    }

    std::vector<int> dims = { m, n };
    int numVecs = (dim == 0) ? n : m;
    int stride = (dim == 0) ? 1 : m;
    int row = 0;
    int col = 0;

    if (dim == 0)
    {
        status = result.Dimension(1, n, data.Type());
    }
    else
    {
        status = result.Dimension(m, 1, data.Type());
    }

    if (!status.IsOk())
        return status;

    for (int i = 0; i < numVecs; ++i)
    {
        // set the matrix indices to the first index in each vector
        int start = col * m + row;

        // perform op
        double min = std::numeric_limits<double>::infinity();
        double max = -min;

        if (data.IsReal())
        {
            const double* real = data.GetRealData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                max = _max(max, *real);
                min = _min(min, *real);
                real += stride;
            }

            result(row, col) = max - min;
        }
        else    // complex
        {
            const hwComplex* cplx = data.GetComplexData() + start;
            double maxph = -PI;
            double minph =  PI;
            int minidx = 0;
            int maxidx = 0;

            for (int j = 0; j < dims[dim]; ++j)
            {
                double magsq = cplx->Mag();

                if (magsq > max)
                {
                    max = magsq;
                    maxph = cplx->Arg();
                    maxidx = j;
                }
                else if (magsq == max && cplx->Arg() > maxph)
                {
                    maxph = cplx->Arg();
                    maxidx = j;
                }

                if (magsq < min)
                {
                    min = magsq;
                    minph = cplx->Arg();
                    minidx = j;
                }
                else if (magsq == min && cplx->Arg() < minph)
                {
                    minph = cplx->Arg();
                    minidx = j;
                }

                cplx += stride;
            }

            cplx = data.GetComplexData() + start;
            result.z(row, col) = cplx[maxidx * stride] - cplx[minidx * stride];
        }

        // advance slice indices
        if (dim == 0)
            ++col;
        else
            ++row;
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes peak to peak values
//------------------------------------------------------------------------------
hwMathStatus Peak2Peak(const hwMatrixN& data, int dim, hwMatrixN& result)
{
    const std::vector<int>& dims = data.Dimensions();
    int numDim = static_cast<int> (dims.size());

    if (dim == -1)
    {
        const std::vector<int>& dims = data.Dimensions();

        // use first non-singleton dimension
        for (int i = 0; i < dims.size(); ++i)
        {
            if (dims[i] != 1)
            {
                dim = i;
                break;
            }
        }
    }
    else if (dim > numDim - 1)
    {
        result.Dimension(dims, hwMatrixN::REAL);
        result.SetElements(0.0);
        return hwMathStatus();
    }
    else if (dim < 0)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYDIM, 2);
    }

    int numVecs = data.Size() / dims[dim];
    int stride = data.Stride(dim);
    std::vector<int> matrixIndex(numDim);
    std::vector<int> outDims = dims;
    outDims[dim] = 1;

    result.Dimension(outDims, data.Type());

    for (int i = 0; i < numVecs; ++i)
    {
        // set the matrix indices to the first index in each vector
        int start = data.Index(matrixIndex);

        // perform op
        double min = std::numeric_limits<double>::infinity();
        double max = -min;

        if (data.IsReal())
        {
            const double* real = data.GetRealData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                max = _max(max, *real);
                min = _min(min, *real);
                real += stride;
            }

            result(matrixIndex) = max - min;
        }
        else    // complex
        {
            const hwComplex* cplx = data.GetComplexData() + start;
            double maxph = -PI;
            double minph = PI;
            int minidx = 0;
            int maxidx = 0;

            for (int j = 0; j < dims[dim]; ++j)
            {
                double magsq = cplx->Mag();

                if (magsq > max)
                {
                    max = magsq;
                    maxph = cplx->Arg();
                    maxidx = j;
                }
                else if (magsq == max && cplx->Arg() > maxph)
                {
                    maxph = cplx->Arg();
                    maxidx = j;
                }

                if (magsq < min)
                {
                    min = magsq;
                    minph = cplx->Arg();
                    minidx = j;
                }
                else if (magsq == min && cplx->Arg() < minph)
                {
                    minph = cplx->Arg();
                    minidx = j;
                }

                cplx += stride;
            }

            cplx = data.GetComplexData() + start;
            result.z(matrixIndex) = cplx[maxidx * stride] - cplx[minidx * stride];
        }

        // advance slice indices
        for (int j = 0; j < numDim; ++j)
        {
            if (j == dim)
                continue;

            // increment index j if possible
            if (matrixIndex[j] < static_cast<int> (dims[j]) - 1)
            {
                ++matrixIndex[j];
                break;
            }

            // index j is maxed out, so reset and continue to j+1
            matrixIndex[j] = 0;
        }
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes next multiple of two
//------------------------------------------------------------------------------
int NextMultiple2(int n)
{
    int N;
    double f = std::frexp(static_cast<double> (n), &N);

    if (f == 0.5)
    {
        N -= 1;
    }

    return 1 << N;
}
//------------------------------------------------------------------------------
// Computes linear correlation of two vectors in the frequency domain
//------------------------------------------------------------------------------
hwMathStatus CorrLin(const hwMatrix& X, const hwMatrix& Y,
                     int maxlag, hwMatrix& Rxy)
{
    if (!X.IsVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);

    if (X.IsEmpty())
        return hwMathStatus(HW_MATH_ERR_EMPTYMATRIX, 1);

    if (!Y.IsVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);

    if (Y.IsEmpty())
        return hwMathStatus(HW_MATH_ERR_EMPTYMATRIX, 2);

    int xn = X.Size();
    int yn = Y.Size();
    int cn = _max(xn + yn - 1, 2 * maxlag + 1);     // size of corr
    int yzm, yzn, yzr, yzc;

    cn = NextMultiple2(cn);    // round up to avoid slow choices

    if (Y.M() == 1)
    {
        yzm = 1;
        yzn = cn;
        yzr = 0;
        yzc = _max(maxlag - yn + 1, 0);
    }
    else
    {
        yzm = cn;
        yzn = 1;
        yzr = _max(maxlag - yn + 1, 0);
        yzc = 0;
    }

    hwMathStatus status;
    hwMatrix XFT;
    status = Fft(X, XFT, cn);   // post-lagged within Fft

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    hwMatrix YReverse;
    status = YReverse.FlipVectors(Y);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    YReverse.Conjugate();

    hwMatrix YLagged(yzm, yzn, Y.Type());
    YLagged.SetElements(0.0);
    status = YLagged.WriteSubmatrix(yzr, yzc, YReverse);    // pre-lagged

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    hwMatrix YFT;
    status = Fft(YLagged, YFT, cn);

    if (YFT.M() != XFT.M())
        YFT.Transpose();

    hwMatrix RxyFT;
    MKLuD::MultByElems(XFT, YFT, RxyFT);
    status = Ifft(RxyFT, Rxy);

    if (!status.IsOk())
    {
        if (status.GetArg1() != 0)
            status.ResetArgs();
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes linear correlation of matrix columns in the frequency domain
//------------------------------------------------------------------------------
hwMathStatus CorrLin(const hwMatrix& X, int maxlag, hwMatrix& Rxy)
{
    hwMathStatus status;

    if (X.IsEmpty())
        return status(HW_MATH_ERR_EMPTYMATRIX, 1);

    if (X.IsVector())
        return status(HW_MATH_ERR_NOVECTOR, 1);

    int m = X.M();
    int n = X.N();
    int cn = _max(2 * m - 1, 2 * maxlag + 1);     // size of conv

    cn = NextMultiple2(cn);    // round up to avoid slow choices

    if (m == 0)
    {
        status = Rxy.Dimension(0, n * n, hwMatrix::REAL);
        return status;
    }

    status = Rxy.Dimension(cn, n * n, hwMatrix::REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (X.IsReal())
    {
        const double* col1;
        const double* col2;
        const double* colr;

        for (int i = 0; i < n; ++i)
        {
            col1 = &X(0, i);
            hwMatrix vec1(m, 1, (void*)col1, hwMatrix::REAL);

            for (int j = 0; j < n; ++j)
            {
                col2 = &X(0, j);
                colr = &Rxy(0, n * i + j);
                hwMatrix vec2(m, 1, (void*)col2, hwMatrix::REAL);
                hwMatrix vecr(cn, 1, (void*)colr, hwMatrix::REAL);

                status = CorrLin(vec1, vec2, maxlag, vecr);

                if (!status.IsOk())
                {
                    if (status.GetArg1() != 0)
                        status.ResetArgs();

                    return status;
                }
            }
        }
    }
    else
    {
        const hwComplex* col1;
        const hwComplex* col2;
        hwMatrix vecr;

        for (int i = 0; i < n; ++i)
        {
            col1 = &X.z(0, i);
            hwMatrix vec1(m, 1, (void*)col1, hwMatrix::COMPLEX);

            for (int j = 0; j < n; ++j)
            {
                col2 = &X.z(0, j);
                hwMatrix vec2(m, 1, (void*)col2, hwMatrix::COMPLEX);

                status = CorrLin(vec1, vec2, maxlag, vecr);

                if (!status.IsOk())
                {
                    if (status.GetArg1() != 0)
                        status.ResetArgs();

                    return status;
                }

                Rxy.WriteColumn(n * i + j, vecr);
            }
        }
    }

    return status;
}
