/**
* @file StatUtilFuncs.cxx
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
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
#include "StatUtilFuncs.h"

#include "hwMatrix.h"

//------------------------------------------------------------------------------
// Returns the status and the covariance of two vectors
//------------------------------------------------------------------------------
hwMathStatus Covariance(const hwMatrix& X, 
                        const hwMatrix& Y, 
                        double&         covar, 
                        bool            sampleCov)
{
    // vectors are stored as hwMatrix objects
    if (!X.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!X.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    if (!Y.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!Y.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }
    if (X.Size() != Y.Size())
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }
    int n = X.Size();
    const double* x = X.GetRealData();
    const double* y = Y.GetRealData();

    covar = Covariance(x, y, n, sampleCov);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Returns the covariance of two vectors
//------------------------------------------------------------------------------
double Covariance(const double* x,
                  const double* y,
                  int           n,
                  bool          sampleCov)
{
    // vectors are stored as double* arrays
    double numer = 0.0;
    double denom = 0.0;

    if (x && y)
    {
        if (n == 1)
        {
            return 0.0;
        }
        double xx;
        double yy;
        double xSum = 0.0;
        double ySum = 0.0;
        double xySum = 0.0;

        for (int i = 0; i < n; i++)
        {
            xx = x[i];
            yy = y[i];
            xSum += xx;
            ySum += yy;
            xySum += xx * yy;
        }

        numer = n * xySum - xSum * ySum;

        denom = (sampleCov) ? n * (n - 1.0) : n * n;
    }

    return numer / denom;
}
//------------------------------------------------------------------------------
// Returns the status and the total sum of squares of a data vector
//------------------------------------------------------------------------------
hwMathStatus TotalSumOfSquares(const hwMatrix& data, double& sst)
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
    double sum = 0.0;
    double sumSq = 0.0;
    double value;
    double data_zero;

    if (n)
    {
        data_zero = data(0);
    }
    for (int i = 1; i < n; i++)
    {
        value = data(i) - data_zero;
        sum += value;
        sumSq += value * value;
    }

    sst = sumSq - sum * sum / (double) n;

    return hwMathStatus();
}
