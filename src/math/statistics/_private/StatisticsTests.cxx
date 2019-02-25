/**
* @file StatisticsTests.cxx
* @date June 2012
* Copyright (C) 2012-2018 Altair Engineering, Inc.  
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
#include "StatisticsTests.h"

#include <math.h>

#include "hwMatrix.h"
#include "DistributionFuncs.h"
#include "StatisticsFuncs.h"

//------------------------------------------------------------------------------
// z Test
//------------------------------------------------------------------------------
hwMathStatus ZTest(const hwMatrix& data,
                   double          mu,
                   double          sigma,
                   bool&           reject,
                   double&         p,
                   hwMatrix&       CI,
                   double          alpha)
{
    int n = data.Size();
    double xBar;
    double testStat;

    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    if (sigma <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }
    if (alpha <= 0.0 || alpha >= 0.5)
    {
        return hwMathStatus(HW_MATH_ERR_STATTESTALPHA, 7);
    }

    hwMathStatus status = CI.Dimension(1, 2, hwMatrix::REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    // compute test statistic
    double s = sigma / sqrt((double) n);
    Mean(data, xBar);
    testStat = (xBar - mu) / (sigma / sqrt((double) n));

    // compute P value
    NormCDF(-fabs(testStat), 0.0, 1.0, p);
    p *= 2.0;

    // compute confidence interval and perform test
    NormInvCDF(0.5 * alpha, xBar, s, CI(0));
    NormInvCDF(1.0 - 0.5 * alpha, xBar, s, CI(1));

    // strictly, the test applies to the test statistic, but this is equivalent
    reject = (mu < CI(0) || mu > CI(1)) ? true: false;

    return status;
}
//------------------------------------------------------------------------------
// One sample t Test
//------------------------------------------------------------------------------
hwMathStatus TTest(const hwMatrix& data,
                   double          mu,
                   bool&           reject,
                   double&         p,
                   hwMatrix&       CI,
                   double          alpha)
{
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!data.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    if (alpha <= 0.0 || alpha >= 0.5)
    {
        return hwMathStatus(HW_MATH_ERR_STATTESTALPHA, 6);
    }

    hwMathStatus status = CI.Dimension(1, 2, hwMatrix::REAL);
    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    int n = data.Size();
    double xBar;
    double s;
    double testStat;

    // compute test statistic
    double testSigma;
    Mean(data, xBar);
    StdDev(data, s);

    testSigma = s / sqrt((double) n);
    testStat = (xBar - mu) / testSigma;

    // compute P value
    T_CDF(-fabs(testStat), n-1, p);
    p *= 2.0;

    // compute confidence interval and perform test
    T_InvCDF(0.5 * alpha, n-1, CI(0));
    CI(1) = -CI(0);

    // CI = CI * s + xBar; 
    CI(0) = CI(0) * testSigma + xBar; 
    CI(1) = CI(1) * testSigma + xBar; 

    // strictly, the test applies to the test statistic, but this is equivalent
    reject = (mu < CI(0) || mu > CI(1)) ? true : false;

    return status;
}
//------------------------------------------------------------------------------
// Two sample t Test
//------------------------------------------------------------------------------
hwMathStatus TTest2(const hwMatrix& data1,
                    const hwMatrix& data2,
                    bool&           reject,
                    double&         p,
                    hwMatrix&       CI,
                    double          alpha)
{
    if (!data1.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!data1.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (!data2.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }

    if (!data2.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }

    if (alpha <= 0.0 || alpha >= 0.5)
    {
        return hwMathStatus(HW_MATH_ERR_STATTESTALPHA, 6);
    }

    hwMathStatus status = CI.Dimension(1, 2, hwMatrix::REAL);
    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    int n = data1.Size();
    int m = data2.Size();
    double xBar1, xBar2;
    double var1, var2;
    double testStat;

    // compute test statistic
    double Sp;  // pooled standard deviation
    double testSigma;

    Mean(data1, xBar1);
    Mean(data2, xBar2);
    Variance(data1, var1);
    Variance(data2, var2);

    Sp = sqrt(((n-1)*var1 + (m-1)*var2) / (double)(n+m-2));
    testSigma = Sp * sqrt(1.0/n + 1.0/m);
    testStat = (xBar1 - xBar2) / testSigma;

    // compute P value
    T_CDF(-fabs(testStat), n+m-2, p);
    p *= 2.0;

    // compute confidence interval and perform test
    T_InvCDF(0.5 * alpha, n+m-2, CI(0));
    CI(1) = -CI(0);

    CI(0) = CI(0) * testSigma + (xBar1 - xBar2); 
    CI(1) = CI(1) * testSigma + (xBar1 - xBar2); 

    // strictly, the test applies to the test statistic, but this is equivalent
    reject = (0 < CI(0) || 0 > CI(1)) ? true : false;

    return status;
}
//------------------------------------------------------------------------------
// Chi-squared Test
//------------------------------------------------------------------------------
hwMathStatus ChiSqTest(const hwMatrix& data,
                       double          var,
                       bool&           reject,
                       double&         p,
                       hwMatrix&       CI,
                       double          alpha)
{
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    if (var <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    if (alpha <= 0.0 || alpha >= 0.5)
    {
        return hwMathStatus(HW_MATH_ERR_STATTESTALPHA, 6);
    }

    hwMathStatus status = CI.Dimension(1, 2, hwMatrix::REAL);
    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    int n = data.Size();
    double s2;
    double testStat;

    // compute test statistic
    Variance(data, s2);
    testStat = (n-1) * s2 / var;

    // compute P value
    Chi2CDF(testStat, n-1, p);

    if (p > 0.5)
    {
       p = 1.0 - p;
    }
    p *= 2.0;

    // compute confidence interval and perform test
    Chi2InvCDF(0.5 * alpha, n-1, CI(1));
    Chi2InvCDF(1.0 - 0.5 * alpha, n-1, CI(0));

    CI(0) = (n-1) * s2 / CI(0);
    CI(1) = (n-1) * s2 / CI(1);

    // strictly, the test applies to the test statistic, but this is equivalent
    reject = (var < CI(0) || var > CI(1)) ? true : false;
    return status;
}
//------------------------------------------------------------------------------
// F Test
//------------------------------------------------------------------------------
hwMathStatus FTest(const hwMatrix& data1,
                   const hwMatrix& data2,
                   bool&           reject,
                   double&         p,
                   hwMatrix&       CI,
                   double          alpha)
{
    if (!data1.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data1.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    if (!data2.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!data2.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }
    if (alpha <= 0.0 || alpha >= 0.5)
    {
        return hwMathStatus(HW_MATH_ERR_STATTESTALPHA, 6);
    }

    hwMathStatus status = CI.Dimension(1, 2, hwMatrix::REAL);
    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    int n = data1.Size();
    int m = data2.Size();
    double var1;
    double var2;
    double testStat;

    // compute test statistic
    Variance(data1, var1);
    Variance(data2, var2);

    testStat = var2 / var1;

    // compute P value
    F_CDF(testStat, m-1, n-1, p);

    if (p > 0.5)
    {
       p = 1.0 - p;
    }
    p *= 2.0;

    // compute confidence interval and perform test
    F_InvCDF(0.5 * alpha, m-1, n-1, CI(0));
    F_InvCDF(1.0 - 0.5 * alpha, m-1, n-1, CI(1));

    CI /= testStat;

    // strictly, the test applies to the test statistic, but this is equivalent
    reject = (1 < CI(0) || 1 > CI(1)) ? true : false;

    return status;
}