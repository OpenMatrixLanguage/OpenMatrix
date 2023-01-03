/**
* @file hwTMatrixSSUtil.cc
* @date December 2022
* Copyright (C) 2022 Altair Engineering, Inc.  
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

//:---------------------------------------------------------------------------
//:Description
//
//  Template sparse matrix class utility functions
//
//:---------------------------------------------------------------------------

#include "../tmpl/MKL_Tutilities.h"     // class definition

//! Multiplication operator
template < typename T1, typename T2 >
hwTMatrixS<T1, T2> operator*(const hwTMatrixS<T1, T2>& A, T1 x)
{
    hwTMatrixS<T1, T2> result(A);
    result.values() *= x;
	return result;
}

//! Multiplication operator
template < typename T1, typename T2 >
hwTMatrixS<T1, T2> operator*(const hwTMatrix<T1, T2>& A, const hwTMatrixS<T1, T2>& B)
{
    hwTMatrixS<T1, T2> result;
    SparseMult(A, B, &result);
    return result;
}

//! Multiplication operator
template < typename T1, typename T2 >
hwTMatrixS<T1, T2> operator*(const hwTMatrixS<T1, T2>& A, const hwTMatrix<T1, T2>& B)
{
    hwTMatrixS<T1, T2> result;
    SparseMult(A, B, &result);
    return result;
}

//! Division operator
template < typename T1, typename T2 >
hwTMatrixS<T1, T2> operator/(const hwTMatrixS<T1, T2>& A, T1 x)
{
    hwTMatrixS<T1, T2> result(A);
    result.values() /= x;
    return result;
}

// vector p-norm
template <typename T1, typename T2>
static T1 vectorNorm(const hwTMatrix<T1, T2>& vec, T1 p)
{
    T1 norm;

    if (p == static_cast<double> (2))
        hwMathStatus status = vec.Norm(norm, static_cast<int> (p));
    else
        hwMathStatus status = vec.Norm(norm, p);

    return norm;
}

// real vector max
template <typename T1, typename T2>
static T1 vectorMax(const hwTMatrix<T1, T2>& vec)
{
    if (!vec.IsReal())
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    if (!vec.IsVector())
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    int size = vec.Size();
    T1  max  = vec(0);

    for (int i = 1; i < size; ++i)
    {
        max = _max(max, vec(i));
    }

    return max;
}

// Find real [lambda, mu] to maximize norm(y*lambda + col*mu, p)
// subject to norm([lambda, mu], p) == 1 (p. 546).
template <typename T1, typename T2>
static void FindLambdaMU(const hwTMatrix<T1, T2>& y, const hwTMatrix<T1, T2>& col,
                         int r, T1 p, T1& lambda_max, T1& mu_max)
{
    if (p == static_cast<T1> (2))
    {
        hwTMatrix<T1, T2> W(y.M(), 2, y.Type());
        hwMathStatus status;

        status = W.WriteColumn(0, y);
        status = W.WriteColumn(1, col);

        hwTMatrix<T1, T2> S;
        hwTMatrix<T1, T2> V;
        status = W.SVD(2, nullptr, S, &V);
        lambda_max = V(0, 0);
        mu_max = V(1,0);
        return;
    }

    // find lambda and mu, p. 546.
    T1 norm_max = static_cast<T1> (0);

    for (int i = 0; i < r; i++)
    {
        T1 thetai = i * static_cast<T1> (PI) / r;
        T1 lambda = cos(thetai);
        T1 mu = sin(thetai);
        T1 nrm = pow(pow(abs(lambda), p) + pow(abs(mu), p), 1 / p);
        lambda /= nrm;
        mu /= nrm;
        T1 vnorm = vectorNorm(lambda * y + mu * col, p);

        if (vnorm > norm_max)
        {
            lambda_max = lambda;
            mu_max = mu;
            norm_max = vnorm;
        }
    }
}

// Find complex [lambda, mu] to maximize norm(y*lambda + col*mu, p)
// subject to norm([lambda, mu], p) == 1 (p. 546).
template <typename T1, typename T2>
static void FindLambdaMU(const hwTMatrix<T1, T2>& y, const hwTMatrix<T1, T2>& col,
                         int r, T1 p, T2& lambda_max, T2& mu_max)
{
    // find lambda and mu, p. 546.
    T1 norm_max = static_cast<T1> (0);

    for (int i = 0; i < r; i++)
    {
        T1 thetai = i * static_cast<T1> (PI) / r;
        T1 lambda = cos(thetai);
        T1 mu = sin(thetai);
        T1 nrm = pow(pow(abs(lambda), p) + pow(abs(mu), p), 1 / p);
        lambda /= nrm;
        mu /= nrm;
        T1 vnorm = vectorNorm(lambda * y + mu * col, p);

        if (vnorm > norm_max)
        {
            lambda_max = lambda;
            mu_max = mu;
            norm_max = vnorm;
        }
    }

    T1 lambdaMag = lambda_max.Mag();
    T2 lambdaDir;

    for (int i = 0; i < r; i++)
    {
        T1 thetai = i * static_cast<T1> (PI) / r;
        lambdaDir.Set(cos(thetai), sin(thetai));
        T1 vnorm = vectorNorm(lambdaMag * lambdaDir * y + mu_max * col, p);

        if (vnorm > norm_max)
        {
            lambda_max = lambdaMag * lambdaDir;
            norm_max = vnorm;
        }
    }
}

// vector dual
template <typename T1, typename T2>
const hwTMatrix<T1, T2> dual_p(const hwTMatrix<T1, T2>& x, T1 p, T1 q)
{
    int m = x.M();
    int n = x.N();

    hwTMatrix<T1, T2> res(m, n, x.Type());

    // compute each p-dual element
    if (x.IsReal())
    {
        for (int i = 0; i < x.Size(); i++)
        {
            if (x(i) >= 0.0)
            {
                res(i) = pow(x(i), p - static_cast<T1> (1));
            }
            else
            {
                res(i) = -pow(-x(i), p - static_cast<T1> (1));
            }
        }
    }
    else    // complex
    {
        for (int i = 0; i < x.Size(); i++)
        {
            T1 mag = x.z(i).Mag();
            hwTComplex<T1> sgn = x.z(i) / mag;

            res.z(i) = sgn * pow(mag, p - static_cast<T1> (1));
        }
    }

    return (res / vectorNorm(res, q));
}

// Higham's hybrid method from "Estimating the matrix p-norm"
template <typename T1, typename T2>
T1 HighamPnorm(const hwTMatrixS<T1, T2>& A, T1 p, T1 tol, int maxiter)
{
    int m = A.M();
    int n = A.N();
    typename hwTMatrix<T1, T2>::DataType dataType;

    if (A.IsReal())
        dataType = hwTMatrix<T1, T2>::REAL;
    else
        dataType = hwTMatrix<T1, T2>::COMPLEX;

    // find gamma and x using algorithm OSE (One Step Estimator), p. 546.
    hwTMatrix<T1, T2> x(n, 1, dataType);
    hwTMatrix<T1, T2> y(m, 1, dataType);
    hwTMatrix<T1, T2> z(m, 1, dataType);
    y.SetElements(static_cast<T1> (0));

    std::vector<hwSliceArg> sliceArg(2);
    sliceArg[0] = hwSliceArg();

    for (int k = 0; k < n; k++)
    {
        hwTMatrixS<T1, T2> col_k_sparse;
        hwTMatrix<T1, T2> col_k_full;
        sliceArg[1] = k;
        A.SliceRHS(sliceArg, col_k_sparse);
        col_k_sparse.Full(col_k_full);

        if (A.IsReal())
        {
            T1 lambda = static_cast<T1> (0);
            T1 mu = static_cast<T1> (1);

            FindLambdaMU(y, col_k_full, 4 * k, p, lambda, mu);

            for (int i = 0; i < k; i++)
                x(i) *= lambda;

            x(k) = mu;
            y = y * lambda + col_k_full * mu;
        }
        else    // complex
        {
            T2 lambda(static_cast<T1> (0), static_cast<T1> (0));
            T2 mu(static_cast<T1> (1), static_cast<T1> (0));

            FindLambdaMU(y, col_k_full, 4 * k, p, lambda, mu);

            for (int i = 0; i < k; i++)
                x.z(i) *= lambda;

            x.z(k) = mu;
            y = y * lambda + col_k_full * mu;
        }
    }

    // compute p-norm using algorithm PM (Power Method), p. 542.
    x /= vectorNorm(x, p);

    T1 q = p / (p - static_cast<T1> (1));
    T1 gamma = static_cast<T1> (0);
    T1 gamma_prev;
    int iter = 0;

    while (iter < maxiter)
    {
        MKL_Tutilities<T1, T2>::SparseMult(A, x, y);
        gamma_prev = gamma;
        gamma = vectorNorm(y, p);
        z = dual_p(y, p, q);
        hwTMatrix<T1, T2> temp;
        temp.Hermitian(z);
        MKL_Tutilities<T1, T2>::SparseMult(temp, A, z);

        if (iter > 0 && (vectorNorm(z, q) <= gamma ||
            (gamma - gamma_prev) <= tol * gamma))
            break;

        z.Hermitian();
        x = dual_p(z, q, p);
        iter++;
    }

    return gamma;
}

// Column 1-norms
template <typename T1, typename T2>
void Col1norms(const hwTMatrixS<T1, T2>& A, hwTMatrix<T1, T2>& colNorms)
{
    int n = A.N();
    hwMathStatus status;
    status = colNorms.Dimension(1, n, hwTMatrix<T1, T2>::REAL);
    colNorms.SetElements(static_cast<T1> (0));

    if (A.IsReal())
    {
        for (int col = 0; col < n; col++)
        {
            for (int idx = A.pointerB()[col]; idx < A.pointerE()[col]; idx++)
            {
                colNorms(col) += fabs(A.values()(idx));
            }
        }
    }
    else
    {
        for (int col = 0; col < n; col++)
        {
            for (int idx = A.pointerB()[col]; idx < A.pointerE()[col]; idx++)
            {
                colNorms(col) += A.values().z(idx).Mag();
            }
        }
    }
}

// Row 1-norms
template <typename T1, typename T2>
void Row1norms(const hwTMatrixS<T1, T2>& A, hwTMatrix<T1, T2>& rowNorms)
{
    int n = A.N();
    hwMathStatus status;
    status = rowNorms.Dimension(A.M(), 1, hwTMatrix<T1, T2>::REAL);
    rowNorms.SetElements(static_cast<T1> (0));

    if (A.IsReal())
    {
        for (int col = 0; col < n; col++)
        {
            for (int idx = A.pointerB()[col]; idx < A.pointerE()[col]; idx++)
            {
                rowNorms(A.rows()[idx]) += fabs(A.values()(idx));
            }
        }
    }
    else
    {
        for (int col = 0; col < n; col++)
        {
            for (int idx = A.pointerB()[col]; idx < A.pointerE()[col]; idx++)
            {
                rowNorms(A.rows()[idx]) += A.values().z(idx).Mag();
            }
        }
    }
}

//! Matrix p-norm
template<typename T1, typename T2>
T1 Norm(const hwTMatrixS<T1, T2>& A, T1 p = static_cast<T1> (2))
{
    if (A.IsVector())
    {
        return vectorNorm(A.values(), p);
    }
    else if (p == static_cast<T1> (1))
    {
        hwTMatrix<T1, T2> colNorms;
        Col1norms(A, colNorms);
        return vectorMax(colNorms);
    }
    else if (p == std::numeric_limits<double>::infinity())
    {
        hwTMatrix<T1, T2> rowNorms;
        Row1norms(A, rowNorms);
        return vectorMax(rowNorms);
    }
    else if (p > static_cast<T1> (1))
    {
        return HighamPnorm<T1, T2>(A, p, MACHEP2, 100);
    }
    else
    {
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);
    }
}

//! Matrix norm
template<typename T1, typename T2 = hwTComplex<T1>>
T1 Norm(const hwTMatrixS<T1, T2>& A, const char* type)
{
    T1 norm;

    if (A.IsVector())
    {
        hwMathStatus status = A.values().Norm(norm, type);

        if (!status.IsOk())
            throw hwMathException(status.GetMsgCode());

        if (!strcmp(type, "-inf") && A.NNZ() < A.Size())
        {
            norm = _min(norm, 0.0);
        }
        else if (!strcmp(type, "inf") && A.NNZ() < A.Size())
        {
            norm = _max(norm, 0.0);
        }
    }
    else
    {
        if (!strcmp(type, "fro"))
        {
            hwMathStatus status = A.values().Norm(norm, type);

            if (!status.IsOk())
                throw hwMathException(status.GetMsgCode());
        }
        else if (!strcmp(type, "inf"))
        {
            hwTMatrix<T1, T2> rowNorms;
            Row1norms(A, rowNorms);
            return vectorMax(rowNorms);
        }
        else
        {
            throw hwMathException(HW_MATH_ERR_INVALIDINPUT);
        }
    }

    return norm;
}
