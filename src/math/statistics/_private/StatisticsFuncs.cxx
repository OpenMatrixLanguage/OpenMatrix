/**
* @file StatisticsFuncs.cxx
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
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
#include "StatisticsFuncs.h"

#include <vector>
#include <algorithm>

#include "DistributionFuncs.h"
#include "BoxBehnken.h"
#include "FullFactorial.h"
#include "hwMatrix.h"
#include "hwNormal.h"
#include "PolynomFuncs.h"
#include "StatUtilFuncs.h"

//------------------------------------------------------------------------------
// Multiple linear regression
//------------------------------------------------------------------------------
hwMathStatus MultiRegress(const hwMatrix& y,
                          const hwMatrix& X,
                          hwMatrix&       b,
                          double          alpha,
                          hwMatrix*       bci,
                          hwMatrix*       r,
                          hwMatrix*       rci,
                          hwMatrix*       stats)
{
    if (y.N() != 1)
    {
        return hwMathStatus(HW_MATH_ERR_COLUMNVEC, 1);
    }
    if (y.M() != X.M())
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYTOOFEWROWS, 1, 2);
    }
    if (X.M() < X.N())
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYTOOFEWROWS, 2);
    }

    // find b to minimize |X*b - y|
    hwMatrix Q;
    hwMatrix R;
    hwMatrix Pinv;   // pseudo-inverse

    hwMathStatus status = X.QR(Q, R);

    if (!status.IsOk())
    {
        status.ResetArgs();
        if (!status.IsWarning())
        {
            return status;
        }
    }

    status = Q.Transpose();
    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    status = Pinv.DivideLeft(R, Q);
    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    b = Pinv * y;

    if (r)
    {
        (*r) = y - X * b;
    }

    if (r && (bci || rci || stats))
    {
        int n = X.M();
        int p = X.N() - 1;    // number of regressors
        int df = n - p - 1;
        double sse;
        double t;

        hwMathStatus status2 = r->L2NormSq(sse);
        double       varHat  = sse / df;

        status2 = T_InvCDF(1.0-alpha / 2.0, df, t);
        if (!status2.IsOk())
        {
            status2.ResetArgs();
            return status2;
        }

        if (bci)
        {
            status2 = bci->Dimension(p+1, 2, hwMatrix::REAL);
            if (!status2.IsOk())
            {
                status2.ResetArgs();
                return status2;
            }

            hwMatrix C;
            hwMatrix RT;

            status = RT.Transpose(R);
            status = C.Inverse(RT*R);  // Covariance matrix via QR decomp
            if (!status.IsOk())
            {
                if (status.GetArg1() == 1)
                {
                    status.SetArg1(2);
                }
                else
                {
                    status.ResetArgs();
                }
                if (!status.IsWarning())
                {
                    return status;
                }
            }

            for (int i = 0; i < p+1; ++i)
            {
                (*bci)(i, 0) = b(i) - t * sqrt(C(i,i) * varHat);
                (*bci)(i, 1) = b(i) + t * sqrt(C(i,i) * varHat);
            }
        }

        if (rci)
        {
            status2 = rci->Dimension(n, 2, hwMatrix::REAL);

            if (!status2.IsOk())
            {
                status2.ResetArgs();
                return status2;
            }

            double hatii;    // hat matrix diagonal element
            double width;    // confidence interval width

            df = n - p - 2;

            for (int i = 0; i < n; ++i)
            {
                hatii = 0.0;

                for (int j = 0; j < p+1; ++j)
                {
                    hatii += X(i,j) * Pinv(j,i);
                }
                // width = t * sqrt(sigmaHat2*(1-h(i,i));
                // sigmaihat2 = (sse - r(i)^ 2 / (1-h(i,i))) / df;   // omits r(i)
                width = t * sqrt((sse*(1.0-hatii) - (*r)(i)*(*r)(i))/df);
                (*rci)(i, 0)  = (*r)(i) - width;
                (*rci)(i, 1)  = (*r)(i) + width;
            }
        }

        if (stats)
        {
            status2 = stats->Dimension(1, 4, hwMatrix::REAL);

            df = n - p - 1;

            // R^2, F, p, error variance
            double sst;
            status2 = TotalSumOfSquares(y, sst);
            
            double R2 = 1.0 - sse / sst;
            double F = df * R2 / (p * (1.0 - R2));
            double pval;

            status2 = F_CDF(F, p, df, pval);

            if (!status2.IsOk())
            {
                status2.ResetArgs();
                return status2;
            }

            pval = 1.0 - pval;

            (*stats)(0) = R2;
            (*stats)(1) = F;
            (*stats)(2) = pval;
            (*stats)(3) = varHat;
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Exponential curve fit
//------------------------------------------------------------------------------
hwMathStatus ExpCurveFit(const hwMatrix& x,
                         const hwMatrix& y,
                         hwMatrix&       coef,
                         hwMatrix*       stats,
                         hwMatrix*       yEst)
{
    // exponential fit of (x,y) data using of y = ae^(bx)
    // implemented as least squares solution to log(y) = log(a) + bx
    if (!x.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!x.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    if (!y.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!y.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }
    int n = x.Size();
    if (y.Size() != n)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    hwMathStatus status = coef.Dimension(2, hwMatrix::REAL);
    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)    // renumber the argument
        {
            status.SetArg1(3);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    hwMatrix A(n, 2, hwMatrix::REAL);
    hwMatrix B(n, 1, hwMatrix::REAL);

    for (int i = 0; i < n; i++)
    {
        A(i, 0) = 1.0;
        A(i, 1) = x(i);

        if (y(i) <= 0.0)
        {
            return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
        }
        B(i, 0) = log(y(i));
    }

    status = coef.QRSolve(A, B);

    if (!status.IsOk())
    {
        if (status != HW_MATH_WARN_MATRIXDEPENDCOL && status != HW_MATH_WARN_MTXNOTFULLRANK)
        {
            status.ResetArgs();
            return status;
        }
    }

    coef(0) = exp(coef(0));     // a

    // compute statistics and fitted points
    if (stats)
    {
        double sse = 0.0;
        double sst;
        double value;

        status = stats->Dimension(2, hwMatrix::REAL);
        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
            {
                status.SetArg1(4);
            }
            else
            {
                status.ResetArgs();
            }
            return status;
        }

        if (yEst)
        {
            status = yEst->Dimension(y.M(), y.N(), hwMatrix::REAL);

            if (!status.IsOk())
            {
                if (status.GetArg1() == 0)
                {
                    status.SetArg1(5);
                }
                else
                {
                    status.ResetArgs();
                }
                return status;
            }

            for (int i = 0; i < n; i++)
            {
                (*yEst)(i) = coef(0) * exp(coef(1) * x(i));
                value = y(i) - (*yEst)(i);
                sse += value * value;
            }
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                value = y(i) - coef(0) * exp(coef(1) * x(i));
                sse += value * value;
            }
        }

        status = TotalSumOfSquares(y, sst);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 1)
            {
                status.SetArg1(2);
            }
            else
            {
                status.ResetArgs();
            }
            return status;
        }

        (*stats)(0) = sse / n;              // mse

        if (IsZero(sst, 1.0e-12))
        {
            (*stats)(1) = 1.0;
        }
        else
        {
            value = 1.0 - sse / sst;

            if (value > 0.0)
            {
                (*stats)(1) = sqrt(value);  // rho
            }
            else
            {
                (*stats)(1) = 0.0;
            }
        }
    }
    else if (yEst)
    {
        status = yEst->Dimension(n, hwMatrix::REAL);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
            {
                status.SetArg1(5);
            }
            else
            {
                status.ResetArgs();
            }
            return status;
        }

        for (int i = 0; i < n; i++)
        {
            (*yEst)(i) = coef(0) * exp(coef(1) * x(i));
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Polynomial curve fit
//------------------------------------------------------------------------------
hwMathStatus PolyCurveFit(const hwMatrix& x,
                          const hwMatrix& y,
                          int             order,
                          hwMatrix&       coef,
                          hwMatrix*       stats,
                          hwMatrix*       yEst,
                          bool            scaled,
                          double*         mu,
                          double*         sigma)
{
    // *mean and *sigma must be supplied if scaled = true
    // coefficients are in ascending order
    if (!x.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!x.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    if (!y.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!y.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }
    int n = x.Size();
    if (y.Size() != n)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }
    if (order < 0)
    {
        return hwMathStatus(HW_MATH_ERR_NONNONNEGINT, 3);
    }

    hwMathStatus status;
    if (order > n-1)
    {
        status(HW_MATH_WARN_TOOFEWPOINTS, 1, 3);
    }
    // solve Vandermonde system
    hwMatrix A(n, order+1, hwMatrix::REAL);
    hwMatrix B(n, 1, hwMatrix::REAL);

    if (!scaled)
    {
        for (int i = 0; i < n; i++)
        {
            A(i, 0) = 1.0;

            if (order > 0)
            {
                A(i, 1) = x(i);
            }
            for (int j = 2; j < order+1; j++)
            {
                A(i, j) = A(i, j-1) * x(i);
            }
            B(i, 0) = y(i);
        }
    }
    else
    {
        // normalize the input data to improve numerical stability
        if (scaled)
        {
            double factor;
            double mean;
            double stdDev;
            double value;

            Mean(x, mean);
            StdDev(x, stdDev);

            if (stdDev < MACHEP)
            {
                return status(HW_MATH_ERR_ZERORANGE, 1);
            }
            if (mu)
            {
                *mu = mean;
            }
            if (sigma)
            {
                *sigma = stdDev;
            }
            factor = 1 / stdDev;

            for (int i = 0; i < n; i++)
            {
                A(i, 0) = 1.0;
                value = (x(i) - mean) * factor;

                if (order > 0)
                {
                    A(i, 1) = value;
                }
                for (int j = 2; j < order+1; j++)
                {
                    A(i, j) = A(i, j-1) * value;
                }
                B(i, 0) = y(i);
            }
        }
    }

    hwMathStatus status2 = coef.QRSolve(A, B);

    if (!status2.IsOk())
    {
        if (status2 == HW_MATH_WARN_MATRIXDEPENDCOL || status2 == HW_MATH_WARN_MTXNOTFULLRANK)
        {
            status(HW_MATH_WARN_POORPOLYFIT);
        }
        else if (status2 != HW_MATH_WARN_NOUNIQUESOL) 
        {
            status2.ResetArgs();
            return status2;
        }
    }

    // compute statistics and fitted points
    if (stats)
    {
        double sse = 0.0;
        double sst;
        double value;
        hwMatrix yhat(n, hwMatrix::REAL);

        status2 = stats->Dimension(2, hwMatrix::REAL);

        if (!status2.IsOk())
        {
            if (status2.GetArg1() == 0)
            {
                status2.SetArg1(5);
            }
            else
            {
                status2.ResetArgs();
            }
            return status2;
        }

        if (yEst)
        {
            if (!scaled)
            {
                status2 = PolyVal(coef, x, (*yEst));

                if (!status2.IsOk())
                {
                    if (status2.GetArg1() == 1)
                    {
                        status2.SetArg1(4);
                    }
                    else if (status2.GetArg1() == 2)
                    {
                        status2.SetArg1(1);
                    }
                    else
                    {
                        status2.ResetArgs();
                    }
                    return status2;
                }
            }
            else
            {
                status2 = yEst->Dimension(y.M(), y.N(), hwMatrix::REAL);

                if (!status2.IsOk())
                {
                    status2.ResetArgs();
                    return status2;
                }

                double factor = 1 / *sigma;

                for (int i = 0; i < n; i++)
                {
                    PolyVal(coef, (x(i)-(*mu))*factor, (*yEst)(i));
                }
            }

            for (int i = 0; i < n; i++)
            {
                value = y(i) - (*yEst)(i);
                sse += value * value;
            }
        }
        else
        {
            hwMatrix yTemp;

            if (!scaled)
            {
                status2 = PolyVal(coef, x, yTemp);

                if (!status2.IsOk())
                {
                    if (status2.GetArg1() == 1)
                    {
                        status2.SetArg1(4);
                    }
                    else if (status2.GetArg1() == 2)
                    {
                        status2.SetArg1(1);
                    }
                    else
                    {
                        status2.ResetArgs();
                    }
                    return status2;
                }
            }
            else
            {
                status2 = yTemp.Dimension(y.M(), y.N(), hwMatrix::REAL);

                if (!status2.IsOk())
                {
                    status2.ResetArgs();
                    return status2;
                }

                double factor = 1 / *sigma;

                for (int i = 0; i < n; i++)
                {
                    PolyVal(coef, (x(i)-(*mu))*factor, yTemp(i));
                }
            }

            for (int i = 0; i < n; i++)
            {
                value = y(i) - yTemp(i);
                sse += value * value;
            }
        }

        (*stats)(0) = sse / n;              // mse

        status2 = TotalSumOfSquares(y, sst);

        if (!status2.IsOk())
        {
            if (status2.GetArg1() == 1)
            {
                status2.SetArg1(2);
            }
            else
            {
                status2.ResetArgs();
            }
            return status2;
        }

        if (IsZero(sst, 1.0e-12))
        {
            (*stats)(1) = 1.0;
        }
        else
        {
            value = 1.0 - sse / sst;

            (*stats)(1) = (value <= 0.0) ? 0.0 : sqrt(value);  // rho
        }
    }
    else if (yEst)
    {
        if (!scaled)
        {
            status2 = PolyVal(coef, x, (*yEst));

            if (!status2.IsOk())
            {
                if (status2.GetArg1() == 1)
                {
                    status2.SetArg1(4);
                }
                else if (status2.GetArg1() == 2)
                {
                    status2.SetArg1(1);
                }
                else if (status2.GetArg1() == 3)
                {
                    status2.SetArg1(6);
                }
                else
                {
                    status2.ResetArgs();
                }
                return status2;
            }
        }
        else
        {
            status2 = yEst->Dimension(y.M(), y.N(), hwMatrix::REAL);

            if (!status2.IsOk())
            {
                status2.ResetArgs();
                return status2;
            }

            double factor = 1 / *sigma;

            for (int i = 0; i < n; i++)
            {
                PolyVal(coef, (x(i)-(*mu))*factor, (*yEst)(i));
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Remove the mean or regression line from a set of data
//------------------------------------------------------------------------------
hwMathStatus Detrend(const hwMatrix& Y, const char* str, hwMatrix& Ymod)
{
    if (!Y.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    hwMathStatus status;
    if (!str)   // remove best fit line
    {
        if (Y.IsVector())
        {
            int size = Y.Size();
            if (size == 1)
            {
                status = Ymod.Dimension(1, 1, hwMatrix::REAL);
                if (!status.IsOk())
                {
                    if (status.GetArg1() == 0)
                    {
                        status.SetArg1(3);
                    }
                    return status;
                }

                Ymod(0) = 0.0;
                return status;
            }

            hwMatrix index(Y.M(), Y.N(), hwMatrix::REAL);
            hwMatrix coef;
            hwMatrix fit;

            for (int i = 0; i < size; ++i)
            {
                index(i) = i;
            }
            status = PolyCurveFit(index, Y, 1, coef, nullptr, &fit);

            if (!status.IsOk())
            {
                status.ResetArgs();
                return status;
            }

            Ymod = Y - fit;
        }
        else
        {
            // operate on each column vector
            int m = Y.M();
            int n = Y.N();
            hwMatrix index(m, hwMatrix::REAL);
            hwMatrix coef;
            hwMatrix fit;
            hwMatrix Ycol(m, 1, hwMatrix::REAL);

            status = Ymod.Dimension(m, n, hwMatrix::REAL);

            if (!status.IsOk())
            {
                if (status.GetArg1() == 0)
                {
                    status.SetArg1(3);
                }
                return status;
            }

            for (int i = 0; i < m; ++i)
            {
                index(i) = i;
            }
            for (int j = 0; j < n; ++j)
            {
                Y.ReadColumn(j, Ycol);

                status = PolyCurveFit(index, Ycol, 1, coef, nullptr, &fit);

                if (!status.IsOk())
                {
                    if (status.IsWarning())
                    {
                        if (status == HW_MATH_WARN_TOOFEWPOINTS)
                        {
                            status(HW_MATH_ERR_NONE);
                        }
                    }
                    {
                        status.ResetArgs();
                        return status;
                    }
                }

                for (int i = 0; i < m; ++i)
                {
                    Ymod(i,j) = Y(i,j) - fit(i);
                }
            }
        }
    }
    else if (!strcmp(str, "constant"))   // remove mean
    {
        if (Y.IsVector())
        {
            double mean;

            status = Mean(Y, mean);

            if (!status.IsOk())
            {
                status.ResetArgs();
                return status;
            }

            Ymod = Y - mean;
        }
        else
        {
            // operate on each column vector
            int m = Y.M();
            int n = Y.N();
            hwMatrix mean;

            status = Mean(Y, mean);

            if (!status.IsOk())
            {
                status.ResetArgs();
                return status;
            }

            status = Ymod.Dimension(m, n, hwMatrix::REAL);

            if (!status.IsOk())
            {
                if (status.GetArg1() == 0)
                {
                    status.SetArg1(3);
                }
                return status;
            }

            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i < m; ++i)
                {
                    Ymod(i,j) = Y(i,j) - mean(j);
                }
            }
        }
    }
    else
    {
        status(HW_MATH_ERR_INVALIDINPUT, 2);
    }

    return status;
}
//------------------------------------------------------------------------------
// Error function
//------------------------------------------------------------------------------
double ErrorFunc(double x)
{
    return hwNormal::Erf(x);
}
//------------------------------------------------------------------------------
// Compute binomial combinations
//------------------------------------------------------------------------------
hwMathStatus NchooseK(int       n,
                      int       k,
                      double&   nCk)
{
    if (n < 1)
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 1);

    if (k < 0)
        return hwMathStatus(HW_MATH_ERR_NEGATIVE, 2);

    if (n < k)
        return hwMathStatus(HW_MATH_ERR_ARG1LTARG2, 1, 2);

    k = _min(k, n - k);
    nCk = 1.0;
    double kk = static_cast<double> (k);
    double nn = static_cast<double> (n);

    while (k--)
    {
        nCk *= nn / kk;
        nn -= 1.0;
        kk -= 1.0;
    }

    nCk = floor(nCk + 0.5);

    if (nCk * 2.0 * k * MACHEP2 >= 0.5)
        return hwMathStatus(HW_MATH_WARN_DOUBLEPRECLIMIT);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Compute binomial combinations
//------------------------------------------------------------------------------
hwMathStatus NchooseK(const hwMatrix& N,
                      int             k,
                      hwMatrix&       nCk)
{
    hwMathStatus status;

    if (!N.IsEmptyOrVector())
        return status(HW_MATH_ERR_VECTOR);

    int n = N.Size();
    int m;
    double mm;

    if (k <= n)
    {
        status = NchooseK(n, k, mm);

        if (!status.IsOk())
            return status;

        m = static_cast<int> (mm);
    }
    else
    {
        m = 0;
    }

    status = nCk.Dimension(m, k, hwMatrix::REAL);

    if (!status.IsOk())
        return status;

    std::vector<int> idx(k);

    for (int col = 0; col < k; ++col)
        idx[col] = col;

    for (int row = 0; row < m; ++row)
    {
        // write a row
        for (int col = 0; col < k; ++col)
        {
            nCk(row, col) = N(idx[col]);
        }

        // update idx values
        for (int col = k - 1; col >= 0; --col)
        {
            // increment last idx[col] that is not maxed out,
            // and set values that follow it
            if (idx[col] != n - k + col)
            {
                ++idx[col];

                for (int col2 = col + 1; col2 < k; ++col2)
                {
                    idx[col2] = idx[col2 - 1] + 1;
                }

                break;
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Compute the mean of a real data vector
//------------------------------------------------------------------------------
hwMathStatus Mean(const hwMatrix& data, double& xBar)
{
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    int    n   = data.Size();
    double sum = 0.0;

    for (int i = 0; i < n; i++)
    {
        sum += data(i);
    }
    xBar = sum / (double) n;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Compute the means of the columns of a real matrix
//------------------------------------------------------------------------------
hwMathStatus Mean(const hwMatrix& A, hwMatrix& xBar)
{
    if (A.IsEmpty())
	{
		hwMathStatus status = xBar.Dimension(1, 1, hwMatrix::REAL);
		xBar(0, 0) = std::numeric_limits<double>::quiet_NaN();
		return status;
	}

    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    int m = A.M();
    int n = A.N();

    hwMathStatus status = xBar.Dimension(1, n, hwMatrix::REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    for (int j = 0; j < n; j++)
    {
        xBar(0, j) = A(0, j);

        for (int i = 1; i < m; i++)
        {
            xBar(0, j) += A(i, j);
        }
    }

    xBar /= (double) m;

    return status;
}
//------------------------------------------------------------------------------
// Compute the median of a data vector
//------------------------------------------------------------------------------
hwMathStatus Median(const hwMatrix& data, double& median)
{
    if (data.IsEmpty())
    {
		median = std::numeric_limits<double>::quiet_NaN();
		return hwMathStatus();
    }

    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    int n = data.Size();
    std::vector<double> copy(n);

    for (int i = 0; i < n; i++)
    {
        if (IsNaN_T(data(i)))
        {
            median = std::numeric_limits<double>::quiet_NaN();
            return hwMathStatus();
        }

        copy[i] = data(i);
    }

    sort(copy.begin(), copy.end());

    if (n%2)
    {
        median = copy[n/2];
    }
    else
    {
        median = 0.5 * (copy[n/2-1] + copy[n/2]);
    }
    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Compute the medians of the columns of a real matrix
//------------------------------------------------------------------------------
hwMathStatus Median(const hwMatrix& A, hwMatrix& median)
{
    hwMathStatus status;

    if (A.IsEmpty())
	{
		status = median.Dimension(1, 1, hwMatrix::REAL);
		median(0, 0) = std::numeric_limits<double>::quiet_NaN();
		return status;
	}

    if (!A.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }

    int m = A.M();
    int n = A.N();
    std::vector<double> copy(m);

    status = median.Dimension(1, n, hwMatrix::REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    for (int j = 0; j < n; j++)
    {
        bool hasNaN = false;

        for (int i = 0; i < m; i++)
        {
            if (IsNaN_T(A(i, j)))
            {
                hasNaN = true;
                break;
            }

            copy[i] = A(i, j);
        }

        if (!hasNaN)
        {
            sort(copy.begin(), copy.end());

            if (m%2)
                median(0, j) = copy[m/2];
            else
                median(0, j) = 0.5 * (copy[m/2-1] + copy[m/2]);
        }
        else
        {
            median(0, j) = std::numeric_limits<double>::quiet_NaN();
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Compute the average absolute deviation of a data vector
//------------------------------------------------------------------------------
hwMathStatus AvgDev(const hwMatrix& data, double& avgDev)
{
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    double mean;
    hwMathStatus status = Mean(data, mean);
    if (!status.IsOk())
    {
        return status;
    }

    int    n   = data.Size();
    double sum = 0.0;

    for (int i = 0; i < n; i++)
    {
        sum += fabs(data(i) - mean);
    }

    avgDev = sum / static_cast<double>(n);

    return status;
}
//------------------------------------------------------------------------------
// Compute the median absolute deviation of a data vector
//------------------------------------------------------------------------------
hwMathStatus MedianDev(const hwMatrix& data, double& medianDev)
{
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    double median;
    hwMathStatus status = Median(data, median);
    if (!status.IsOk())
    {
        return status;
    }

    hwMatrix devs;

    status = devs.Abs(data - median);

    if (!status.IsOk())
    {
        return status;
    }

    status = Median(devs, medianDev);

    return status;
}
//------------------------------------------------------------------------------
// Compute the standard deviation of a data vector
//------------------------------------------------------------------------------
hwMathStatus StdDev(const hwMatrix& data, double& stdDev, bool sampleStd)
{
    double       variance;
    hwMathStatus status = Variance(data, variance, sampleStd);
    if (!status.IsOk())
    {
        return status;
    }

    stdDev = sqrt(variance);
    return status;
}
//------------------------------------------------------------------------------
// Compute the variance of a data vector
//------------------------------------------------------------------------------
hwMathStatus Variance(const hwMatrix& data, double& variance, bool sampleVar)
{
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    int n = data.Size();
    double data_zero;
    double sum = 0.0;
    double sumSq = 0.0;
    double value;

    if (n)
    {
        data_zero = data(0);
    }

    for (int i = 1; i < n; i++)
    {
        value = data(i) - data_zero;        // shift mean to avoid overflow
        sum += value;
        sumSq += value * value;
    }

    if (sampleVar)
    {
        if (n != 1)
        {
            variance = ((double) n * sumSq - sum * sum) / ((double) n * ((double) (n - 1)));
        }
        else // (n == 1)
        {
            variance = 0.0;
        }
    }
    else // population variance
    {
        variance = ((double) n * sumSq - sum * sum) / ((double) (n * n));
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Compute the variances of the columns of a real matrix
//------------------------------------------------------------------------------
hwMathStatus Variance(const hwMatrix& A, hwMatrix& variance, bool sampleVar)
{
    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    int m = A.M();
    int n = A.N();
    hwMathStatus status = (n) ? 
                          variance.Dimension(1, n, hwMatrix::REAL):
                          variance.Dimension(0, 0, hwMatrix::REAL);
    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    double sum;
    double sumSq;
    double value;
    double data_zero;

    for (int j = 0; j < n; j++)
    {
        sum   = 0.0;
        sumSq = 0.0;

        if (m > 0)
        {
            data_zero = A(0, j);

            for (int i = 0; i < m; i++)
            {
                value = A(i, j) - data_zero;
                sum += value;
                sumSq += value * value;
            }
        }

        if (sampleVar)
        {
            if (m != 1)
            {
                variance(0, j) = ((double) m * sumSq - sum * sum) / ((double) m * ((double) (m - 1)));
            }
            else // (m == 1)
            {
                variance(0, j) = 0.0;
            }
        }
        else // population variance
        {
            variance(0, j) = ((double) m * sumSq - sum * sum) / ((double) (m * m));
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Compute the skewness of a data vector
//------------------------------------------------------------------------------
hwMathStatus Skewness(const hwMatrix& data, double& skewness, bool correctBias)
{
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    int n = data.Size();
    double x;
    double value;
    double xSum = 0.0;
    double xSqSum = 0.0;
    double xCuSum = 0.0;

    if (correctBias && n < 3)
    {
        skewness = std::numeric_limits<double>::quiet_NaN();
        return hwMathStatus();
    }

    double data_zero = data(0);

    for (int i = 1; i < n; i++)
    {
        x = data(i);
        x -= data_zero;         // shift data
        xSum += x;
        value = x * x;
        xSqSum += value;
        value *= x;
        xCuSum += value;
    }

    double nn = static_cast<double> (n);
    double mean = xSum / nn;
    double variance = xSqSum / nn - mean * mean;    // population, not sample
    double stdDev = sqrt(variance);

    skewness = xCuSum - 3.0 * mean * xSqSum + 2.0 * nn * mean * mean * mean;
    skewness /= nn * variance * stdDev;

    if (correctBias)
    {
        skewness *= sqrt(nn * (nn - 1.0)) / (nn - 2.0);
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Compute the kurtosis of a data vector
//------------------------------------------------------------------------------
hwMathStatus Kurtosis(const hwMatrix& data, double& kurtosis, bool correctBias)
{
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    int n = data.Size();
    double x;
    double value;
    double xSum = 0.0;
    double xSqSum = 0.0;
    double xCuSum = 0.0;
    double xQdSum = 0.0;

    if (correctBias && n < 4)
    {
        kurtosis = std::numeric_limits<double>::quiet_NaN();
        return hwMathStatus();
    }

    double data_zero = data(0);

    for (int i = 1; i < n; i++)
    {
        x = data(i);
        x -= data_zero;         // shift data
        xSum += x;
        value = x * x;
        xSqSum += value;
        value *= x;
        xCuSum += value;
        value *= x;
        xQdSum += value;
    }

    double nn = static_cast<double> (n);
    double mean = xSum / nn;  // shifted mean
    double mean2 = mean * mean;
    double mean4 = mean2 * mean2;
    double variance = xSqSum / nn - mean * mean;    // population, not sample
    double stdDev = sqrt(variance);

    kurtosis = xQdSum - 4.0 * mean * xCuSum + 6.0 * mean2 * xSqSum
             - 3.0 * nn * mean4;
    kurtosis /= nn * variance * variance;

    if (correctBias)
    {
        kurtosis = (nn + 1.0) * kurtosis - 3.0 * (nn - 1.0);
        kurtosis = (nn - 1.0) / ((nn - 2.0) * (nn - 3.0)) * kurtosis + 3.0;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Compute the root mean square value of a data vector
//------------------------------------------------------------------------------
hwMathStatus RMS(const hwMatrix& data, double& rms)
{
    if (!data.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    int n = data.Size();
    double sumSq = 0.0;

    if (data.IsReal())
    {
        for (int i = 0; i < n; i++)
        {
            sumSq += data(i) * data(i);
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            sumSq += data.z(i).MagSq();
        }
    }

    rms = sqrt(sumSq / static_cast<double>(n));

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Compute covariance matrix for the columns of a matrix
//------------------------------------------------------------------------------
hwMathStatus Cov(const hwMatrix& X, hwMatrix& cov, bool sampleVar)
{
    if (!X.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    int n = X.N();
    hwMathStatus status = cov.Dimension(n, n, hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(2);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    int m = X.M();
    int inc_i = 0;
    int inc_j;
    double covar;
    const double* x_start = X.GetRealData();

    for (int i = 0; i < n; i++)
    {
        inc_j = 0;

        for (int j = 0; j < i; j++)
        {
            covar = Covariance(x_start + inc_i, x_start + inc_j, m, sampleVar);
            cov(i, j) = covar;
            cov(j, i) = covar;
            inc_j += m;
        }

        cov(i, i) = Covariance(x_start + inc_i, x_start + inc_i, m, sampleVar);
        inc_i += m;
    }

    return status;
}
//------------------------------------------------------------------------------
// Compute covariance matrix for the columns of two matrices
//------------------------------------------------------------------------------
hwMathStatus Cov(const hwMatrix& X,
                 const hwMatrix& Y,
                 hwMatrix&       cov,
                 bool            sampleCovar)
{
    if (!X.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!Y.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }

    int m = X.M();
    if (Y.M() != m)
    {
        return hwMathStatus(HW_MATH_ERR_COLUMNDIM, 1, 2);
    }

    int xn = X.N();
    int yn = Y.N();

    hwMathStatus status = cov.Dimension(xn, yn, hwMatrix::REAL);

    if (!status.IsOk())
    {
        status.SetArg1(3);
        return status;
    }

    int inc_i = 0;
    int inc_j;
    double covar;
    const double* x_start = X.GetRealData();
    const double* y_start = Y.GetRealData();

    for (int i = 0; i < xn; i++)
    {
        inc_j = 0;

        for (int j = 0; j < yn; j++)
        {
            covar = Covariance(x_start + inc_i, y_start + inc_j, m, sampleCovar);
            cov(i, j) = covar;
            inc_j += m;
        }

        inc_i += m;
    }

    return status;
}
//------------------------------------------------------------------------------
// Compute correlation matrix for the columns of a matrix
//------------------------------------------------------------------------------
hwMathStatus Corr(const hwMatrix& X, hwMatrix& corr)
{
    hwMatrix var;
    hwMathStatus status = Cov(X, corr, true);
    if (!status.IsOk())
    {
        return status;
    }

    status = Variance(X, var, true);
    if (!status.IsOk())
    {
        return status;
    }

    int n = X.N();

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            corr(i,j) /= sqrt(var(i) * var(j));
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Compute correlation matrix for the columns of two matrices
//------------------------------------------------------------------------------
hwMathStatus Corr(const hwMatrix& X, const hwMatrix& Y, hwMatrix& corr)
{
    hwMatrix varX;
    hwMatrix varY;

    hwMathStatus status = Cov(X, Y, corr, true);
    if (!status.IsOk())
    {
        return status;
    }

    status = Variance(X, varX, true);
    if (!status.IsOk())
    {
        return status;
    }

    status = Variance(Y, varY, true);
    if (!status.IsOk())
    {
        return status;
    }

    int xn = X.N();
    int yn = Y.N();

    for (int i = 0; i < xn; ++i)
    {
        for (int j = 0; j < yn; ++j)
        {
            corr(i,j) /= sqrt(varX(i) * varY(j));
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Compute the geomtric mean of a data vector
//------------------------------------------------------------------------------
hwMathStatus GeoMean(const hwMatrix& data, double& geoMean)
{
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    bool zeroVal = false;
    int n = data.Size();
    double temp = 0.0;
    double value;

    for (int i = 0; i < n; i++)
    {
        value = data(i);

        if (value <= 0.0)
        {
            if (value < 0.0)
            {
                return hwMathStatus(HW_MATH_ERR_NEGATIVE, 1);
            }
            zeroVal = true;
            break;
        }

        temp += log(value);
    }

    geoMean = (zeroVal) ? 0.0 : exp(temp / (double) n);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generate the histogram bin count for a data sample
//------------------------------------------------------------------------------
hwMathStatus PDF(const hwMatrix& data, hwMatrixI& bin)
{
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    if (!bin.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!bin.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }

    int n = data.Size();
    if (n < 1)
    {
        bin.SetElements(0);
        return hwMathStatus();
    }

    int numBins = bin.Size();
    double temp_max = data(0);
    double temp_min = data(0);

    for (int i = 1; i < n; i++)
    {
        if (data(i) > temp_max)
        {
            temp_max = data(i);
        }
        else if (data(i) < temp_min)
        {
            temp_min = data(i);
        }
    }

    int index;
    double bin_size = (temp_max - temp_min) / numBins;

    if (IsZero(bin_size, 1.0e-12))
    {
        if (numBins == 1)
        {
            bin(0) = n;
            return hwMathStatus();
        }
        return hwMathStatus(HW_MATH_ERR_ZERORANGE, 1);
    }

    bin.SetElements(0);

    for (int i = 0; i < n; i++)
    {
        index = (int) floor((data(i) - temp_min) / bin_size);

        if (index == numBins)
        {
            index--;
        }
        bin(index)++;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Construct a full factorial design matrix
//------------------------------------------------------------------------------
hwMathStatus FullDoe(const hwMatrixI& runs, hwMatrixI& matrix)
{
    return FullFact(runs, matrix);
}
//------------------------------------------------------------------------------
// Construct a Box-Behnken design matrix
//------------------------------------------------------------------------------
hwMathStatus BBDoe(int n, hwMatrixI& matrix)
{
    return BoxBehnken(n, matrix);
}