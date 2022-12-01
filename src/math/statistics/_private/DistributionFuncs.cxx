/**
* @file DistributionFuncs.cxx
* @date June 2009
* Copyright (C) 2009-2018 Altair Engineering, Inc.  
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
#include "hwMatrix.h"
#include "DistributionFuncs.h"

#include <sys/timeb.h>
#include <unordered_map>

#include "StatisticsFuncs.h"
#include "hwMersenneTwister.h"
#include "hwUniform.h"
#include "hwNormal.h"
#include "hwBeta.h"
#include "hwGamma.h"
#include "hwExponential.h"
#include "hwChiSquared.h"
#include "hwStudent_t.h"
#include "hw_F.h"
#include "hwPearson.h"
#include "hwLogNormal.h"
#include "hwWeibull.h"
#include "hwPoisson.h"
#include "StatUtilFuncs.h"
#include "SpecialFuncs.h"
#include "hwOptimizationFuncs.h"
#include "MathUtilsFuncs.h"
#include "hwMathStatus.h"

#ifndef OS_WIN
#if defined(OS_UNIX) || defined(_DARWIN)
#define _ftime ftime
#define _timeb timeb
#endif // _DARWIN
static void ftime_s_wrapper(struct _timeb *timeptr)
{
    _ftime(timeptr);
    return;
}
#define _ftime_s ftime_s_wrapper
#endif

double nearZero = 1.0e-12;

//------------------------------------------------------------------------------
// Generates seed for random number generator using the system clock
//------------------------------------------------------------------------------
unsigned long GetSeedFromClock()
{
    // seed static generator if first time
    static bool firstTime = true;
    static hwMersenneTwisterState MTstate;

    if (firstTime)
    {
        unsigned long seed;
        struct _timeb timebuffer;
        _ftime_s( &timebuffer );

        seed  = (unsigned long) timebuffer.time;
        seed += (unsigned long) timebuffer.millitm;

        MTstate.Initialize(seed);
        firstTime = false;
    }

    // cycle generator
    hwUniform uniDist(&MTstate);

    uniDist.GetDeviate();
    int index = MTstate.mti - 1;

    return MTstate.mt[index];
}
//------------------------------------------------------------------------------
// Uniform distribution probability density function
//------------------------------------------------------------------------------
hwMathStatus UnifPDF(double x, double a, double b, double& density)
{
    double width = b - a;
    if (width < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINTERVAL, 2, 3);
    }

    density = (x < a || x > b) ? 0.0 : 1.0 / width;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Uniform distribution cumulative density function
//------------------------------------------------------------------------------
hwMathStatus UnifCDF(double x, double a, double b, double& prob)
{
    double width = b - a;
    if (width < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINTERVAL, 2, 3);
    }

    if (x < a)
    {
        prob = 0.0;
    }
    else if (x > b)
    {
        prob = 1.0;
    }
    else
    {
        prob = (x - a) / width;
    }
    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Uniform distribution inverse cumulative density function
//------------------------------------------------------------------------------
hwMathStatus UnifInvCDF(double prob, double a, double b, double& x)
{
    if (prob <= 0.0 || prob >= 1.0)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDPROB, 1);
    }

    double width = b - a;
    if (width < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINTERVAL, 2, 3);
    }

    if (prob < 0.0)
    {
        x = a;
    }
    else if (prob > 1.0)
    {
        x = b;
    }
    else
    {
        x = a + prob * width;
    }
    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Uniform distribution parameter estimation function
//------------------------------------------------------------------------------
hwMathStatus UnifFit(const hwMatrix& data, 
                     double&         aHat, 
                     double&         bHat,
                     hwMatrix*       aCI, 
                     hwMatrix*       bCI)
{
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    hwMathStatus status;

    // allocate output
    if (aCI)
    {
        status = aCI->Dimension(2, hwMatrix::REAL);
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
    }

    if (bCI)
    {
        status = bCI->Dimension(2, hwMatrix::REAL);
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
    }

    int numPnts = data.Size();

    // compute MLE
    aHat = (numPnts == 0) ? 
           std::numeric_limits<double>::quiet_NaN() : data(0);
    bHat = aHat;

    double value;
    for (int i = 1; i < numPnts; ++i)
    {
        value = data(i);

        if (value < aHat)
        {
            aHat = value;
        }
        if (value > bHat)
        {
            bHat = value;
        }
    }

    // compute confidence intervals
    value = (bHat - aHat) / pow(0.05, 1.0 / (double) numPnts);  // Beta(1,numPnts)

    if (aCI)
    {
        (*aCI)(0) = bHat - value;
        (*aCI)(1) = aHat;
    }

    if (bCI)
    {
        (*bCI)(0) = bHat;
        (*bCI)(1) = aHat + value;
    }

    return status;
}
//------------------------------------------------------------------------------
// Generates uniform distribution random number, a single value on (0,1)
//------------------------------------------------------------------------------
hwMathStatus UnifRnd(hwMersenneTwisterState* MTstate,
                     unsigned long*          seed, 
                     double&                 value)
{
    // prepare random number generator
    bool createState = false;

    if (!MTstate)
    {
        MTstate     = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values from U(0,1)
    hwUniform uniDist(MTstate);

    value = uniDist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates uniform distribution random number, a single value on (A,B)
//------------------------------------------------------------------------------
hwMathStatus UnifRnd(double                  A, 
                     double                  B, 
                     hwMersenneTwisterState* MTstate,
                     unsigned long*          seed,
                     double&                 value)
{
    if (B - A < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINTERVAL, 1, 2);
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate     = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    hwUniform uniDist(MTstate);

    value = A + (B - A) * uniDist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates uniform distribution random numbers, a matrix of values on (A,B)
//------------------------------------------------------------------------------
hwMathStatus UnifRnd(double                  A, 
                     double                  B, 
                     hwMersenneTwisterState* MTstate,
                     unsigned long*          seed, 
                     hwMatrix&               matrix)
{
    if (B - A < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINTERVAL, 1, 2);
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate     = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int m = matrix.M();
    int n = matrix.N();
    int size = m * n;
    hwUniform uniDist(MTstate);

    for (int i = 0; i < size; ++i)
    {
        matrix(i) = A + (B - A) * uniDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates uniform distribution random numbers, a matrix of values on (A,B)
//------------------------------------------------------------------------------
hwMathStatus UnifRnd(const hwMatrix&         A, 
                     const hwMatrix&         B,
                     hwMersenneTwisterState* MTstate, 
                     unsigned long*          seed,
                     hwMatrix&               matrix)
{
    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!B.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }

    int m = A.M();
    int n = A.N();
    if (B.M() != m || B.N() != n)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    hwMathStatus status = matrix.Dimension(m, n, hwMatrix::REAL);
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

    if (matrix.IsEmpty())
    {
        return status;
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    double width;
    int size = m * n;
    hwUniform uniDist(MTstate);

    for (int i = 0; i < size; ++i)
    {
        width = B(i) - A(i);
        if (width < 0.0)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return status(HW_MATH_ERR_INVALIDINTERVAL, 1, 2);
        }
        matrix(i) = A(i) + width * uniDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return status;
}
//------------------------------------------------------------------------------
// Normal distribution probability density function
//------------------------------------------------------------------------------
hwMathStatus NormPDF(double x, double mu, double sigma, double& density)
{
    if (sigma < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }

    hwNormal normDist;
    double xx = (x - mu) / sigma;

    density = normDist.Pdf(xx) / sigma;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Normal distribution cumulative density function
//------------------------------------------------------------------------------
hwMathStatus NormCDF(double x, double mu, double sigma, double& prob)
{
    if (sigma < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }

    hwNormal normDist;
    double xx = (x - mu) / sigma;

    prob = normDist.Cdf(xx);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Normal distribution inverse cumulative density function
//------------------------------------------------------------------------------
hwMathStatus NormInvCDF(double prob, double mu, double sigma, double& x)
{
    if (prob < MACHEP || prob > 1.0 - MACHEP)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDPROB, 1);
    }
    if (sigma < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }

    hwNormal normDist;
    x = normDist.CdfInv(prob);
    x = mu + x * sigma;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Normal distribution parameter estimation function
//------------------------------------------------------------------------------
hwMathStatus NormFit(const hwMatrix& data, 
                     double&         muHat, 
                     double&         sigmaHat,
                     hwMatrix*       muCI,
                     hwMatrix*       sigmaCI)
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
    hwMathStatus status;

    // allocate output
    if (muCI)
    {
        status = muCI->Dimension(2, hwMatrix::REAL);
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
    }

    if (sigmaCI)
    {
        status = sigmaCI->Dimension(2, hwMatrix::REAL);
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
    }

    // compute MLE
    double data0 = (n == 0) ? 
                   std::numeric_limits<double>::quiet_NaN() : data(0);

    double nRec = 1.0 / (double) n;
    double sum   = 0.0;
    double sumSq = 0.0;
    double value;

    for (int i = 1; i < n; ++i)
    {
        value = data(i) - data0;
        sum += value;
        sumSq += value * value;
    }

    muHat = sum * nRec + data0;

    if (n != 1)
    {
        sigmaHat = sqrt((sumSq - sum * sum / (double) n) / (double) (n-1));
    }
    else
    {
        sigmaHat = 0.0;
    }

    // compute confidence intervals
    if (muCI)
    {
        if (n > 1)
        {
            hwStudent_t stdTDist(n - 1);

            value = stdTDist.CdfInv(0.975) / sqrt((double) n);
            (*muCI)(0) = muHat - value * sigmaHat;
            (*muCI)(1) = muHat + value * sigmaHat;
        }
        else
        {
            (*muCI)(0) = muHat;
            (*muCI)(1) = muHat;
        }
    }

    if (sigmaCI)
    {
        if (n > 1)
        {
            hwChiSquared chiSqDist(n - 1);

            value = (n - 1) / chiSqDist.CdfInv(0.975);
            (*sigmaCI)(0) = sigmaHat * sqrt(value);

            value = (n - 1) / chiSqDist.CdfInv(0.025);
            (*sigmaCI)(1) = sigmaHat * sqrt(value);
        }
        else
        {
            (*sigmaCI)(0) = sigmaHat;
            (*sigmaCI)(1) = sigmaHat;
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Generates Normal distribution random number, a single value from N(0,1)
//------------------------------------------------------------------------------
hwMathStatus NormRnd(hwMersenneTwisterState* MTstate,
                     unsigned long*          seed, 
                     double&                 value)
{
    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values from N(0,1)
    hwNormal normDist(MTstate);

    value = normDist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Normal distribution random number, a single value from 
// N(mu,sigma^2)
//------------------------------------------------------------------------------
hwMathStatus NormRnd(double                  mu, 
                     double                  sigma,
                     hwMersenneTwisterState* MTstate,
                     unsigned long*          seed, 
                     double&                 value)
{
    if (sigma < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    hwNormal normDist(MTstate);

    value = mu + sigma * normDist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Normal distribution random numbers, a matrix of values
// from N(mu,sigma^2)
//------------------------------------------------------------------------------
hwMathStatus NormRnd(double                  mu, 
                     double                  sigma, 
                     hwMersenneTwisterState* MTstate,
                     unsigned long*          seed, 
                     hwMatrix&               matrix)
{
    if (sigma < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int m = matrix.M();
    int n = matrix.N();
    int size = m * n;
    hwNormal normDist(MTstate);

    for (int i = 0; i < size; ++i)
    {
        matrix(i) = mu + sigma * normDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Normal distribution random numbers, a matrix of values from N(mu,sigma^2)
//------------------------------------------------------------------------------
hwMathStatus NormRnd(const hwMatrix&         mu,
                     const hwMatrix&         sigma,
                     hwMersenneTwisterState* MTstate, 
                     unsigned long*          seed, 
                     hwMatrix&               matrix)
{
    if (!mu.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!sigma.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }

    int m = mu.M();
    int n = mu.N();
    if (sigma.M() != m || sigma.N() != n)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    hwMathStatus status = matrix.Dimension(m, n, hwMatrix::REAL);

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

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    double sigma_temp;
    int size = m * n;
    hwNormal normDist(MTstate);

    for (int i = 0; i < size; ++i)
    {
        sigma_temp = sigma(i);
        if (sigma_temp < nearZero)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return status(HW_MATH_ERR_NONPOSITIVE, 2);
        }
        matrix(i) = mu(i) + sigma_temp * normDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return status;
}
//------------------------------------------------------------------------------
// Beta distribution probability density function
//------------------------------------------------------------------------------
hwMathStatus BetaPDF(double x, double a, double b, double& density)
{
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    if (b < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }

    hwBeta betaDist(a, b);
    density = betaDist.Pdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Beta distribution cumulative density function
//------------------------------------------------------------------------------
hwMathStatus BetaCDF(double x, double a, double b, double& prob)
{
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    if (b < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }

    hwBeta betaDist(a, b);
    prob = betaDist.Cdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Beta distribution inverse cumulative density function
//------------------------------------------------------------------------------
hwMathStatus BetaInvCDF(double prob, double a, double b, double& x)
{
    if (prob < MACHEP || prob > 1.0 - MACHEP)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDPROB, 1);
    }
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    if (b < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }

    hwBeta betaDist(a, b);
    x = betaDist.CdfInv(prob);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Beta distribution log-likelihood function
//------------------------------------------------------------------------------
static hwMathStatus BetaLogLike(const hwMatrix& P, 
                                const hwMatrix* userData, 
                                double&         logL)
{
    double a   = P(0) * P(0);
    double b   = P(1) * P(1);
    double n   = (*userData)(0);
    double ld  = (*userData)(1);
    double l1d = (*userData)(2);

    logL = n * BetaLog(a,b) + (1-a) * ld + (1-b) * l1d;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Beta distribution parameter estimation function
//------------------------------------------------------------------------------
hwMathStatus BetaFit(const hwMatrix& data, 
                     double&         aHat, 
                     double&         bHat,
                     hwMatrix*       aCI, 
                     hwMatrix*       bCI)
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

    // allocate output
    hwMathStatus status;
    if (aCI)
    {
        status = aCI->Dimension(2, hwMatrix::REAL);
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
    }

    if (bCI)
    {
        status = bCI->Dimension(2, hwMatrix::REAL);
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
    }

    // compute MLE
    if (n == 0)
    {
        aHat = std::numeric_limits<double>::quiet_NaN();
        bHat = std::numeric_limits<double>::quiet_NaN();

        if (aCI)
        {
            aCI->SetElements(aHat);
        }
        if (bCI)
        {
            bCI->SetElements(bHat);
        }
        return status;
    }

    if (n == 1)
    {
        return status(HW_MATH_ERR_DISTINCTVALS3, 1);
    }

    double max = data(0);
    double min = data(0);

    for (int i = 1; i < n; i++)
    {
        if (data(i) < min)
        {
            min = data(i);
        }
        if (data(i) > max)
        {
            max = data(i);
        }
    }

    if (min <= 0.0)
    {
        return status(HW_MATH_ERR_BADRANGE, 1);
    }

    if (max >= 1.0)
    {
        return status(HW_MATH_ERR_BADRANGE, 1);
    }
    if (min == max)
    {
        return status(HW_MATH_ERR_ZERORANGE, 1);
    }

    // check for enough distinct points
    int i;
    int index;
    double value;

    for (i = 1; i < n; i++)
    {
        if (IsZero(data(i) - data(0), 1.0e-12))
        {
            continue;
        }
        index = i;
        value = data(i);
        break;
    }

    if (i == n)
    {
        return status(HW_MATH_ERR_DISTINCTVALS3, 1);
    }

    for (i = 1; i < n; i++)
    {
        if (i == index)
        {
            continue;
        }

        if (!IsZero(data(i) - value, 1.0e-12))
        {
            if (!IsZero(data(i) - data(0), 1.0e-12))
            {
                break;
            }
        }
    }

    if (i == n)
    {
        return status(HW_MATH_ERR_DISTINCTVALS3, 1);
    }

    // compute initial estimates
    double nRec   = 1.0 / (double) n;
    double value1 = pow(1.0 - data(0), nRec);
    double value2 = pow(data(0), nRec);
    double value3;
    hwMatrix param(2, hwMatrix::REAL);

    for (i = 1; i < n; i++)
    {
        value1 *= pow(1.0 - data(i), nRec);
        value2 *= pow(data(i), nRec);
    }

    value3 = 1.0 - value1 - value2;
    param(0) = 0.5 * (1.0 - value1) / value3;   // aHat
    param(1) = 0.5 * (1.0 - value2) / value3;   // bHat

    param(0) = sqrt(param(0));
    param(1) = sqrt(param(1));

    // optimization step
    hwMatrix userData(3, hwMatrix::REAL);

    value1 = log(data(0));
    value2 = log(1.0 - data(0));

    for (i = 1; i < n; i++)
    {
        value1 += log(data(i));
        value2 += log(1.0 - data(i));
    }

    userData(0)     = (double) n;
    userData(1)     = value1;
    userData(2)     = value2;
    int maxIter     = 200;
    int maxFuncEval = 1000;

    status = FMinUncon(BetaLogLike, nullptr, param, min, maxIter, maxFuncEval,
                       1.0e-8, 1.0e-6, nullptr, nullptr, &userData);

    if (!status.IsOk() && !status.IsInfoMsg() && !status.IsWarning())
    {
        status.ResetArgs();
        return status(HW_MATH_ERR_NOTCONVERGE);
    }

    status = hwMathStatus();

    aHat = param(0) * param(0);
    bHat = param(1) * param(1);

    // compute confidence intervals
    if (aCI || bCI)
    {
        double eps = pow(MachPrecision(1.0), 1.0/3.0);
        hwMatrix J(n, 2, hwMatrix::REAL);
        hwMatrix JT;
        hwMatrix H(2, 2, hwMatrix::REAL);
        hwMatrix Cov;

        double delta1 = _max(eps * fabs(aHat), 1.0e-10);
        double temp = aHat + delta1;
        delta1 = temp - aHat;

        double delta2 = _max(eps * fabs(bHat), 1.0e-10);
        temp = bHat + delta2;
        delta2 = temp - bHat;

        double deriv1 = (BetaLog(aHat + delta1, bHat)
                      - BetaLog(aHat - delta1, bHat)) / (2.0 * delta1);
        double deriv2 = (BetaLog(aHat, bHat + delta2)
                      - BetaLog(aHat, bHat - delta2)) / (2.0 * delta2);

        for (i = 0; i < n; i++)
        {
            J(i, 0) = log(data(i)) - deriv1;
            J(i, 1) = log(1.0 - data(i)) - deriv2;
        }

        JT.Transpose(J);
        H = JT * J;     // modify to use QR later

        hwMathStatus status2 = Cov.Inverse(H);

        if (!status2.IsOk())
        {
            status2.ResetArgs();
            return status2;
        }

        if (aCI)
        {
            NormInvCDF(0.025, aHat, sqrt(Cov(0,0)), (*aCI)(0));
            NormInvCDF(0.975, aHat, sqrt(Cov(0,0)), (*aCI)(1));
        }

        if (bCI)
        {
            NormInvCDF(0.025, bHat, sqrt(Cov(1,1)), (*bCI)(0));
            NormInvCDF(0.975, bHat, sqrt(Cov(1,1)), (*bCI)(1));
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Generates Beta distribution random number, a single value from Beta(a,b)
//------------------------------------------------------------------------------
hwMathStatus BetaRnd(double                  a, 
                     double                  b, 
                     hwMersenneTwisterState* MTstate,
                     unsigned long*          seed, 
                     double&                 value)
{
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }
    if (b < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    // prepare random number generator
    bool createState =  false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    hwBeta betaDist(a, b, MTstate);

    value = betaDist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Beta distribution random numbers, a matrix of values from Beta(a,b)
//------------------------------------------------------------------------------
hwMathStatus BetaRnd(double                  a, 
                     double                  b, 
                     hwMersenneTwisterState* MTstate,
                     unsigned long*          seed, 
                     hwMatrix&               matrix)
{
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }
    if (b < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int m = matrix.M();
    int n = matrix.N();
    int size = m * n;
    hwBeta betaDist(a, b, MTstate);

    for (int i = 0; i < size; i++)
    {
        matrix(i) = betaDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Beta distribution random numbers, a matrix of values from Beta(a,b)
//------------------------------------------------------------------------------
hwMathStatus BetaRnd(const hwMatrix&         A, 
                     const hwMatrix&         B,
                     hwMersenneTwisterState* MTstate, 
                     unsigned long*          seed, 
                     hwMatrix&               matrix)
{
    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!B.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }

    int m = A.M();
    int n = A.N();

    if (B.M() != m || B.N() != n)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    hwMathStatus status = matrix.Dimension(m, n, hwMatrix::REAL);

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

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    double A_temp;
    double B_temp;
    int size = m * n;

    for (int i = 0; i < size; i++)
    {
        A_temp = A(i);
        B_temp = B(i);

        if (A_temp < nearZero)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
        }
        if (B_temp < nearZero)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
        }
        hwBeta betaDist(A_temp, B_temp, MTstate);

        matrix(i) = betaDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return status;
}
//------------------------------------------------------------------------------
// Gamma distribution probability density function
//------------------------------------------------------------------------------
hwMathStatus GammaPDF(double x, double a, double b, double& density)
{
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    if (b <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }

    hwGamma gammaDist(a);

    density = gammaDist.Pdf(x / b) / b;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Gamma distribution cumulative density function
//------------------------------------------------------------------------------
hwMathStatus GammaCDF(double x, double a, double b, double& prob)
{
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    if (b <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }
    hwGamma gammaDist(a);

    prob = gammaDist.Cdf(x / b);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Gamma distribution inverse cumulative density function
//------------------------------------------------------------------------------
hwMathStatus GammaInvCDF(double prob, double a, double b, double& x)
{
    if (prob < MACHEP || prob > 1.0 - MACHEP)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDPROB, 1);
    }

    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    if (b <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }

    hwGamma gammaDist(a);

    x = b * gammaDist.CdfInv(prob);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Gamma distribution log-likelihood function
//------------------------------------------------------------------------------
static hwMathStatus GammaLogLikeUtil(const hwMatrix& P, 
                                     const hwMatrix& X,
                                     const hwMatrix* userData,
                                     hwMatrix&       result)
{
    double a  = P(0);
    double n  = (*userData)(0);
    double ld = (*userData)(1);
    double sd = (*userData)(2);

    result(0) = log(a) - Digamma(a) - log(sd / n) + ld / n;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Gamma distribution parameter estimation function
//------------------------------------------------------------------------------
hwMathStatus GammaFit(const hwMatrix& data, 
                      double&         aHat, 
                      double&         bHat,
                      hwMatrix*       aCI, 
                      hwMatrix*       bCI)
{
    hwMathStatus status;

    if (!data.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 1);
    }

    int n = data.Size();

    // allocate output
    if (aCI)
    {
        status = aCI->Dimension(2, hwMatrix::REAL);

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
    }

    if (bCI)
    {
        status = bCI->Dimension(2, hwMatrix::REAL);

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
    }

    // compute MLE
    if (n == 0)
    {
        aHat = std::numeric_limits<double>::quiet_NaN();
        bHat = std::numeric_limits<double>::quiet_NaN();

        if (aCI)
        {
            aCI->SetElements(aHat);
        }
        if (bCI)
        {
            bCI->SetElements(bHat);
        }
        return status;
    }

    if (n == 1)
    {
        return status(HW_MATH_ERR_TOOFEWPOINTS, 1);
    }

    double min;
    status = Min(data, &min);
    if (!status.IsOk())
    {
        if (status != HW_MATH_ERR_EMPTYMATRIX)
        {
            return status;
        }
        status = hwMathStatus();    // reset and continue
    }
    else if (min <= 0.0)
    {
        return status(HW_MATH_ERR_NONPOSITIVE, 1);
    }

    double max;
    status = Max(data, &max);

    if (!status.IsOk())
    {
        if (status != HW_MATH_ERR_EMPTYMATRIX)
        {
            return status;
        }
        status = hwMathStatus();    // reset and continue
    }
    else if (min == max)
    {
        return status(HW_MATH_ERR_ZERORANGE, 1);
    }

    // compute initial estimates
    double ld = 0.0;
    double sd = 0.0;
    double nRec = 1.0 / (double) n;
    double mean;
    double value;
    hwMatrix param(1, hwMatrix::REAL);

    Mean(data, mean);

    for (int i = 0; i < n; i++)
    {
        value = data(i) / mean;     // scale by mean
        ld += log(value);
        sd += value;
    }

    double s = log(sd * nRec) - ld * nRec;
    value = 3.0 - s;

    param(0) = (value + sqrt(value*value + 24.0*s)) / (12.0*s);     // aHat

    // optimization step
    hwMatrix userData(3, hwMatrix::REAL);

    userData(0)     = (double) n;
    userData(1)     = ld;
    userData(2)     = sd;
    int maxIter     = 200;
    int maxFuncEval = 1000;

    status = NLSolve(GammaLogLikeUtil, nullptr, param, min, maxIter, maxFuncEval,
                     1, 1.0e-6, 1.0e-6, nullptr, nullptr, &userData);

    if (!status.IsOk() && !status.IsInfoMsg() && !status.IsWarning())
    {
        status.ResetArgs();
        return status(HW_MATH_ERR_NOTCONVERGE);
    }

    status = hwMathStatus();

    aHat = param(0);
    bHat = sd * nRec / aHat;

    // compute confidence intervals
    if (aCI || bCI)
    {
        hwMatrix H(2, 2, hwMatrix::REAL);
        hwMatrix Cov;

        value = bHat * bHat;

        H(0, 0) = n * Trigamma(aHat);
        H(0, 1) = n / bHat;
        H(1, 0) = H(0, 1);
        H(1, 1) = (2.0 * sd / bHat - n * aHat) / value;

        hwMathStatus status2 = Cov.Inverse(H);

        if (!status2.IsOk())
        {
            status2.ResetArgs();
            return status2;
        }

        if (aCI)
        {
            NormInvCDF(0.025, log(aHat), sqrt(Cov(0,0)) / aHat, (*aCI)(0));
            NormInvCDF(0.975, log(aHat), sqrt(Cov(0,0)) / aHat, (*aCI)(1));
            (*aCI)(0) = exp((*aCI)(0));
            (*aCI)(1) = exp((*aCI)(1));
        }

        if (bCI)
        {
            NormInvCDF(0.025, log(bHat), sqrt(Cov(1,1)) / bHat, (*bCI)(0));
            NormInvCDF(0.975, log(bHat), sqrt(Cov(1,1)) / bHat, (*bCI)(1));
            (*bCI)(0) = exp((*bCI)(0));
            (*bCI)(1) = exp((*bCI)(1));
            (*bCI) *= mean;
        }
    }

    bHat *= mean;

    return status;
}
//------------------------------------------------------------------------------
// Generates Gamma distribution random number, a single value from Gamma(a,b)
//------------------------------------------------------------------------------
hwMathStatus GammaRnd(double                  a, 
                      double                  b, 
                      hwMersenneTwisterState* MTstate,
                      unsigned long*          seed, 
                      double&                 value)
{
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }
    if (b <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    hwGamma gammaDist(a, MTstate);

    value = b * gammaDist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Gamma distribution random numbers, a matrix of values from Gamma(a,b)
//------------------------------------------------------------------------------
hwMathStatus GammaRnd(double                  a, 
                      double                  b, 
                      hwMersenneTwisterState* MTstate,
                      unsigned long*          seed, 
                      hwMatrix&               matrix)
{
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }
    if (b <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    { 
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int m = matrix.M();
    int n = matrix.N();
    int size = m * n;
    hwGamma gammaDist(a, MTstate);

    for (int i = 0; i < size; i++)
    {
        matrix(i) = b * gammaDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Gamma distribution random numbers, a matrix of values from Gamma(a,b)
//------------------------------------------------------------------------------
hwMathStatus GammaRnd(const hwMatrix&         A,
                      const hwMatrix&         B, 
                      hwMersenneTwisterState* MTstate,
                      unsigned long*          seed, 
                      hwMatrix&               matrix)
{
    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!B.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }

    int m = A.M();
    int n = A.N();
    if (B.M() != m || B.N() != n)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    hwMathStatus status = matrix.Dimension(m, n, hwMatrix::REAL);
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

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else  if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    double A_temp;
    double B_temp;
    int size = m * n;

    for (int i = 0; i < size; i++)
    {
        A_temp = A(i);
        B_temp = B(i);

        if (A_temp < nearZero)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
        }
        if (B_temp <= 0.0)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
        }
        hwGamma gammaDist(A_temp, MTstate);

        matrix(i) = B_temp * gammaDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return status;
}
//------------------------------------------------------------------------------
// Exponential distribution probability density function
//------------------------------------------------------------------------------
hwMathStatus ExpPDF(double x, double mu, double& density)
{
    if (mu <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    hwExponential expDist(mu);

    density = expDist.Pdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Exponential distribution cumulative density function
//------------------------------------------------------------------------------
hwMathStatus ExpCDF(double x, double mu, double& prob)
{
    if (mu <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    hwExponential expDist(mu);

    prob = expDist.Cdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Exponential distribution inverse cumulative density function
//------------------------------------------------------------------------------
hwMathStatus ExpInvCDF(double prob, double mu, double& x)
{
    if (prob < MACHEP || prob > 1.0 - MACHEP)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDPROB, 1);
    }
    if (mu <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    hwExponential expDist(mu);

    x = expDist.CdfInv(prob);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Exponential distribution parameter estimation function
//------------------------------------------------------------------------------
hwMathStatus ExpFit(const hwMatrix& data, double& muHat, hwMatrix* muCI)
{
    hwMathStatus status;

    if (!data.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 1);
    }
    int n = data.Size();

    // allocate output
    if (muCI)
    {
        status = muCI->Dimension(2, hwMatrix::REAL);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
            {
                status.SetArg1(3);
            }
            else
            {
                status.ResetArgs();
            }
            return status;
        }
    }

    // compute MLE
    if (n == 0)
    {
        muHat = std::numeric_limits<double>::quiet_NaN();
        if (muCI)
        {
            muCI->SetElements(muHat);
        }
        return status;
    }

    if (n == 1)
    {
        return status(HW_MATH_ERR_TOOFEWPOINTS, 1);
    }

    double min;
    status = Min(data, &min);

    if (!status.IsOk())
    {
        if (status != HW_MATH_ERR_EMPTYMATRIX)
        {
            return status;
        }
        status = hwMathStatus();    // reset and continue
    }
    else if (min <= 0.0)
    {
        return status(HW_MATH_ERR_NONPOSITIVE, 1);
    }

    Mean(data, muHat);

    // compute confidence intervals
    if (muCI)
    {
        hwGamma gammaDist(static_cast<double>(n));

        double value = gammaDist.CdfInv(0.975);
        (*muCI)(0) = (double) n * muHat / value;  //\todo: Check value for near zero

        value = gammaDist.CdfInv(0.025);
        (*muCI)(1) = (double) n * muHat / value;
    }

    return status;
}
//------------------------------------------------------------------------------
// Generates exponential distribution random number, a single value with mean mu
//------------------------------------------------------------------------------
hwMathStatus ExpRnd(double                  mu, 
                    hwMersenneTwisterState* MTstate,
                    unsigned long*          seed, 
                    double&                 value)
{
    if (mu <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {     
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    hwExponential expDist(mu, MTstate);

    value = expDist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates exponential distribution random numbers, a matrix of values
// with mean mu
//------------------------------------------------------------------------------
hwMathStatus ExpRnd(double                  mu, 
                    hwMersenneTwisterState* MTstate,
                    unsigned long*          seed, 
                    hwMatrix&               matrix)
{
    if (mu <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int m = matrix.M();
    int n = matrix.N();
    int size = m * n;
    hwExponential expDist(mu, MTstate);

    for (int i = 0; i < size; i++)
    {
        matrix(i) = expDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates exponential distribution random numbers, a matrix of values
// with mean mu
//------------------------------------------------------------------------------
hwMathStatus ExpRnd(const hwMatrix&         Mu, 
                    hwMersenneTwisterState* MTstate,
                    unsigned long*          seed, 
                    hwMatrix&               matrix)
{
    if (!Mu.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    int m = Mu.M();
    int n = Mu.N();

    hwMathStatus status = matrix.Dimension(m, n, hwMatrix::REAL);

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

    if (matrix.IsEmpty())
    {
        return status;
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int size = m * n;
    double mu;

    for (int i = 0; i < size; i++)
    {
        mu = Mu(i);

        if (mu <= 0.0)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
        }
        hwExponential expDist(mu, MTstate);
        matrix(i) = expDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return status;
}
//------------------------------------------------------------------------------
// Chi-squared distribution probability density function
//------------------------------------------------------------------------------
hwMathStatus Chi2PDF(double x, int ndof, double& density)
{
    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);
    }
    hwChiSquared chi2Dist(ndof);

    density = chi2Dist.Pdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Chi-squared distribution cumulative density function
//------------------------------------------------------------------------------
hwMathStatus Chi2CDF(double x, int ndof, double& prob)
{
    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);
    }
    hwChiSquared chi2Dist(ndof);

    prob = chi2Dist.Cdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Chi-squared distribution inverse cumulative density function
//------------------------------------------------------------------------------
hwMathStatus Chi2InvCDF(double prob, int ndof, double& x)
{
    if (prob < MACHEP || prob > 1.0 - MACHEP)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDPROB, 1);
    }
    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);
    }
    hwChiSquared chi2Dist(ndof);

    x = chi2Dist.CdfInv(prob);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Chi-squared distribution random number, a single value with ndof
// degrees of freedom
//------------------------------------------------------------------------------
hwMathStatus Chi2Rnd(int ndof, hwMersenneTwisterState* MTstate,
                     unsigned long* seed, double& value)
{
    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 1);
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    hwChiSquared chi2Dist(ndof, MTstate);

    value = chi2Dist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Chi-squared distribution random numbers, a matrix of values
// with ndof degrees of freedom
//------------------------------------------------------------------------------
hwMathStatus Chi2Rnd(int                     ndof, 
                     hwMersenneTwisterState* MTstate,
                     unsigned long*          seed, 
                     hwMatrix&               matrix)
{
    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 1);
    }
    // prepare random number generator
    bool createState = false;

    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int m = matrix.M();
    int n = matrix.N();
    int size = m * n;
    hwChiSquared chi2Dist(ndof, MTstate);

    for (int i = 0; i < size; i++)
    {
        matrix(i) = chi2Dist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Chi-squared distribution random numbers, a matrix of values
// with ndof degrees of freedom
//------------------------------------------------------------------------------
hwMathStatus Chi2Rnd(const hwMatrixI&        Ndof, 
                     hwMersenneTwisterState* MTstate,
                     unsigned long*          seed, 
                     hwMatrix&               matrix)
{
    if (!Ndof.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    int m = Ndof.M();
    int n = Ndof.N();

    hwMathStatus status = matrix.Dimension(m, n, hwMatrix::REAL);

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

    if (matrix.IsEmpty())
    {
        return status;
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int size = m * n;
    int ndof;

    for (int i = 0; i < size; i++)
    {
        ndof = Ndof(i);
        if (ndof < 1)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return status(HW_MATH_ERR_NONPOSINT, 1);
        }
        hwChiSquared chi2Dist(ndof, MTstate);

        matrix(i) = chi2Dist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return status;
}
//------------------------------------------------------------------------------
// Student T distribution probability density function
//------------------------------------------------------------------------------
hwMathStatus T_PDF(double x, int ndof, double& density)
{
    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);
    }
    hwStudent_t tDist(ndof);

    density = tDist.Pdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Student T distribution cumulative density function
//------------------------------------------------------------------------------
hwMathStatus T_CDF(double x, int ndof, double& prob)
{
    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);
    }
    hwStudent_t tDist(ndof);

    prob = tDist.Cdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Student T distribution inverse cumulative density function
//------------------------------------------------------------------------------
hwMathStatus T_InvCDF(double prob, int ndof, double& x)
{
    if (prob < MACHEP || prob > 1.0 - MACHEP)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDPROB, 1);
    }
    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);
    }
    hwStudent_t tDist(ndof);

    x = tDist.CdfInv(prob);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Student T distribution random number, a single value with ndof
// degrees of freedom
//------------------------------------------------------------------------------
hwMathStatus T_Rnd(int                     ndof, 
                   hwMersenneTwisterState* MTstate,
                   unsigned long*          seed, 
                   double&                 value)
{
    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 1);
    }
    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }
    // generate random values
    hwStudent_t tDist(ndof, MTstate);

    value = tDist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Student T distribution random numbers, a matrix of valuex with
// ndof degrees of freedom
//------------------------------------------------------------------------------
hwMathStatus T_Rnd(int                     ndof, 
                   hwMersenneTwisterState* MTstate,
                   unsigned long*          seed, 
                   hwMatrix&               matrix)
{
    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 1);
    }
    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int m = matrix.M();
    int n = matrix.N();
    int size = m * n;
    hwStudent_t tDist(ndof, MTstate);

    for (int i = 0; i < size; i++)
    {
        matrix(i) = tDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Student T distribution random numbers, a matrix of values
// with Ndof degrees of freedom
//------------------------------------------------------------------------------
hwMathStatus T_Rnd(const hwMatrixI&        Ndof, 
                   hwMersenneTwisterState* MTstate,
                   unsigned long*          seed, 
                   hwMatrix&               matrix)
{
    if (!Ndof.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    int m = Ndof.M();
    int n = Ndof.N();

    hwMathStatus status = matrix.Dimension(m, n, hwMatrix::REAL);

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

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int size = m * n;
    int ndof;

    for (int i = 0; i < size; i++)
    {
        ndof = Ndof(i);

        if (ndof < 1)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return status(HW_MATH_ERR_NONPOSINT, 1);
        }

        hwStudent_t tDist(ndof, MTstate);
        matrix(i) = tDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return status;
}
//------------------------------------------------------------------------------
// F distribution probability density function
//------------------------------------------------------------------------------
hwMathStatus F_PDF(double x, int mdof, int ndof, double& density)
{
    if (mdof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);
    }

    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 3);
    }
    hw_F F_dist(mdof, ndof);

    density = F_dist.Pdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// F distribution cumulative density function
//------------------------------------------------------------------------------
hwMathStatus F_CDF(double x, int mdof, int ndof, double& prob)
{
    if (mdof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);
    }
    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 3);
    }
    hw_F F_dist(mdof, ndof);

    prob = F_dist.Cdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// F distribution inverse cumulative density function
//------------------------------------------------------------------------------
hwMathStatus F_InvCDF(double prob, int mdof, int ndof, double& x)
{
    if (prob < MACHEP || prob > 1.0 - MACHEP)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDPROB, 1);
    }
    if (mdof < 1)
    {
        x = 0.0;
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);
    }

    if (ndof < 1)
    {
        x = 0.0;
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 3);
    }

    hw_F F_dist(mdof, ndof);

    x = F_dist.CdfInv(prob);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates F distribution random number, a single value with mdof and ndof
// degrees of freedom
//------------------------------------------------------------------------------
hwMathStatus F_Rnd(int                     mdof, 
                   int                     ndof, 
                   hwMersenneTwisterState* MTstate,
                   unsigned long*          seed, 
                   double&                 value)
{
    if (mdof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 1);
    }
    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    hw_F F_dist(mdof, ndof, MTstate);

    value = F_dist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates F distribution random numbers, a matrix of values with mdof and
// ndof degrees of freedom
//------------------------------------------------------------------------------
hwMathStatus F_Rnd(int                     mdof, 
                   int                     ndof, 
                   hwMersenneTwisterState* MTstate,
                   unsigned long*          seed, 
                   hwMatrix&               matrix)
{
    if (mdof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 1);
    }
    if (ndof < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);
    }
    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int m = matrix.M();
    int n = matrix.N();
    int size = m * n;
    hw_F F_dist(mdof, ndof, MTstate);

    for (int i = 0; i < size; i++)
    {
        matrix(i) = F_dist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates F distribution random numbers, a matrix of values with Mdof
// and Ndof degrees of freedom
//------------------------------------------------------------------------------
hwMathStatus F_Rnd(const hwMatrixI&        Mdof, 
                   const hwMatrixI&        Ndof,
                   hwMersenneTwisterState* MTstate,
                   unsigned long*          seed, 
                   hwMatrix&               matrix)
{
    if (!Mdof.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!Ndof.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }

    int m = Ndof.M();
    int n = Ndof.N();

    if (Mdof.M() != m)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }
    if (Mdof.N() != n)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }
    hwMathStatus status = matrix.Dimension(m, n, hwMatrix::REAL);

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

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int size = m * n;
    int mdof;
    int ndof;

    for (int i = 0; i < size; i++)
    {
        mdof = Mdof(i);
        if (mdof < 1)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return status(HW_MATH_ERR_NONPOSINT, 1);
        }
        ndof = Ndof(i);
        if (ndof < 1)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return status(HW_MATH_ERR_NONPOSINT, 1);
        }
        hw_F F_dist(mdof, ndof, MTstate);

        matrix(i) = F_dist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return status;
}
//------------------------------------------------------------------------------
// Compute four moments needed for Pearson distributions
//------------------------------------------------------------------------------
hwMathStatus FourMoments(const hwMatrix& data, hwMatrix& moment)
{
    // compute the mean, variance, skewness, and kurtosis
    // skewness and kurtosis are normalized
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    int size = data.Size();
    if (size < 2)
    {
        return hwMathStatus(HW_MATH_ERR_TOOFEWPOINTS, 1);
    }

    hwMathStatus status = moment.Dimension(4, hwMatrix::REAL);
    if (!status.IsOk())
    {
        status.SetArg1(2);
        return status;
    }

    double x;
    double value;
    double xSum = 0.0;
    double xSqSum = 0.0;
    double xCuSum = 0.0;
    double xQdSum = 0.0;

    double data_zero = data(0);

    for (int i = 1; i < size; i++)
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

    double mean = xSum / size;  // shifted mean

    double variance = (xSqSum - xSum * xSum / size) / (size - 1);
    double stdDev = sqrt(variance);

    double mean2 = mean * mean;
    value = (size - 1) * variance * stdDev;
    double skewness = xCuSum - 3.0 * mean * xSqSum + 2.0 * size * mean2 * mean;
    skewness /= value;

    double mean4 = mean2 * mean2;
    value = (size - 1) * variance * variance;
    double kurtosis = xQdSum - 4.0 * mean * xCuSum + 6.0 * mean2 * xSqSum
        - 3.0 * size * mean4;
    kurtosis /= value;
    kurtosis -= 3.0;

    mean += data_zero;          // reverse mean shift

    moment(0) = mean;
    moment(1) = variance;
    moment(2) = skewness;
    moment(3) = kurtosis;

    return status;
}
//------------------------------------------------------------------------------
// Gets the Pearson distribution family type and parameters
//------------------------------------------------------------------------------
hwMathStatus PearsonInfo(const hwMatrix& moment, int& type, hwMatrix& param)
{
    hwPearson pearsonDist(moment);

    type = pearsonDist.Type();
    pearsonDist.GetParams(param);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Pearson distribution probability density function
//------------------------------------------------------------------------------
hwMathStatus PearsonPDF(double x, const hwMatrix& moment, double& density)
{
    hwPearson pearsonDist(moment);
    if (pearsonDist.Type() == -1)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT);
    }
    density = pearsonDist.Pdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Pearson distribution cumulative density function
//------------------------------------------------------------------------------
hwMathStatus PearsonCDF(double x, const hwMatrix& moment, double& prob)
{
    hwPearson pearsonDist(moment);
    if (pearsonDist.Type() == -1)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT);
    }
    prob = pearsonDist.Cdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Pearson distribution inverse cumulative density function
//------------------------------------------------------------------------------
hwMathStatus PearsonInvCDF(double prob, const hwMatrix& moment, double& x)
{
    if (prob < MACHEP || prob > 1.0 - MACHEP)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDPROB, 1);
    }

    hwPearson pearsonDist(moment);
    if (pearsonDist.Type() == -1)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT);
    }
    x = pearsonDist.CdfInv(prob);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Pearson distribution median
//------------------------------------------------------------------------------
hwMathStatus PearsonMedian(const hwMatrix& moment, double& median)
{
    hwPearson pearsonDist(moment);
    if (pearsonDist.Type() == -1)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT);
    }

    median = pearsonDist.Median();

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Pearson distribution mode
//------------------------------------------------------------------------------
hwMathStatus PearsonMode(const hwMatrix& moment, double& mode)
{
    hwPearson pearsonDist(moment);
    if (pearsonDist.Type() == -1)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT);
    }

    mode = pearsonDist.Mode();

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Pearson distribution random number, a single value
//------------------------------------------------------------------------------
hwMathStatus PearsonRnd(const hwMatrix&         moment, 
                        hwMersenneTwisterState* MTstate,
                        unsigned long*          seed, 
                        double&                 value)
{
    hwPearson pearsonDist(moment, MTstate);
    if (pearsonDist.Type() == -1)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT);
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    value = pearsonDist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Pearson distribution random numbers, a matrix of values
//------------------------------------------------------------------------------
hwMathStatus PearsonRnd(const hwMatrix&         moment, 
                        hwMersenneTwisterState* MTstate,
                        unsigned long*          seed, 
                        hwMatrix&               matrix)
{
    int m = matrix.M();
    int n = matrix.N();

    hwPearson pearsonDist(moment, MTstate);
    if (pearsonDist.Type() == -1)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT);
    }

    // prepare random number generator
    bool createState =  false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int size = m * n;

    for (int i = 0; i < size; i++)
    {
        matrix(i) = pearsonDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Lognormal distribution probability density function
//------------------------------------------------------------------------------
hwMathStatus LogNormPDF(double x, double mu, double sigma, double& density)
{
    if (x <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }
    if (sigma < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }
    if (sigma < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }

    hwLogNormal logNormDist(mu, sigma);

    density = logNormDist.Pdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Lognormal distribution cumulative density function
//------------------------------------------------------------------------------
hwMathStatus LogNormCDF(double x, double mu, double sigma, double& prob)
{
    if (x <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }
    if (sigma < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }
    hwLogNormal logNormDist(mu, sigma);

    prob = logNormDist.Cdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Lognormal distribution inverse cumulative density function
//------------------------------------------------------------------------------
hwMathStatus LogNormInvCDF(double prob, double mu, double sigma, double& x)
{
    if (prob < MACHEP || prob > 1.0 - MACHEP)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDPROB, 1);
    }
    if (sigma < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }
    hwLogNormal logNormDist(mu, sigma);

    x = logNormDist.CdfInv(prob);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Lognormal distribution parameter estimation function
//------------------------------------------------------------------------------
hwMathStatus LogNormFit(const hwMatrix& data, 
                        double&         muHat, 
                        double&         sigmaHat,
                        hwMatrix*       muCI, 
                        hwMatrix*       sigmaCI)
{
    hwMathStatus status;

    if (!data.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 1);
    }
    int n = data.Size();

    // allocate output
    if (muCI)
    {
        status = muCI->Dimension(2, hwMatrix::REAL);

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
    }

    if (sigmaCI)
    {
        status = sigmaCI->Dimension(2, hwMatrix::REAL);

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
    }

    // compute MLE
    double min;
    double lndata0 = (n == 0) ? 
                     std::numeric_limits<double>::quiet_NaN() : log(data(0));
    status = Min(data, &min);

    if (!status.IsOk())
    {
        if (status != HW_MATH_ERR_EMPTYMATRIX)
        {
            return status;
        }
        status = hwMathStatus();    // reset and continue
    }
    else if (min <= 0.0)
    {
        return status(HW_MATH_ERR_NONPOSITIVE, 1);
    }

    double nRec = 1.0 / (double) n;
    double sum = 0.0;
    double sumSq = 0.0;
    double value;

    for (int i = 1; i < n; i++)
    {
        value = log(data(i)) - lndata0;
        sum += value;
        sumSq += value * value;
    }

    muHat    = sum * nRec + lndata0;
    sigmaHat = (n != 1) ? 
               sqrt((sumSq - sum * sum * nRec) / (double) (n-1)) : 0.0;

    // compute confidence intervals
    if (muCI)
    {
        if (n > 1)
        {
            hwStudent_t stdTDist(n - 1);

            value = stdTDist.CdfInv(0.975) * sqrt(nRec);
            (*muCI)(0) = muHat - value * sigmaHat;
            (*muCI)(1) = muHat + value * sigmaHat;
        }
        else
        {
            (*muCI)(0) = muHat;
            (*muCI)(1) = muHat;
        }
    }

    if (sigmaCI)
    {
        if (n > 1)
        {
            hwChiSquared chiSqDist(n - 1);

            value = (n - 1) / chiSqDist.CdfInv(0.975);
            (*sigmaCI)(0) = sigmaHat * sqrt(value);

            value = (n - 1) / chiSqDist.CdfInv(0.025);
            (*sigmaCI)(1) = sigmaHat * sqrt(value);
        }
        else
        {
            (*sigmaCI)(0) = sigmaHat;
            (*sigmaCI)(1) = sigmaHat;
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Generates Lognormal distribution random number, a single value from N(mu,sigma^2)
//------------------------------------------------------------------------------
hwMathStatus LogNormRnd(double                  mu, 
                        double                  sigma, 
                        hwMersenneTwisterState* MTstate,
                        unsigned long*          seed,
                        double&                 value)
{
    if (sigma < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    hwLogNormal logNormDist(mu, sigma, MTstate);
    value = logNormDist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Lognormal distribution random numbers, a matrix of values
// from N(mu,sigma^2)
//------------------------------------------------------------------------------
hwMathStatus LogNormRnd(double                  mu, 
                        double                  sigma, 
                        hwMersenneTwisterState* MTstate,
                        unsigned long*          seed, 
                        hwMatrix&               matrix)
{
    if (sigma < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    // prepare random number generator
    bool createState =  false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int m = matrix.M();
    int n = matrix.N();
    int size = m * n;
    hwLogNormal logNormDist(mu, sigma, MTstate);

    for (int i = 0; i < size; i++)
    {
        matrix(i) = logNormDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Lognormal distribution random numbers, a matrix of values
// from N(mu,sigma^2)
//------------------------------------------------------------------------------
hwMathStatus LogNormRnd(const hwMatrix&         mu, 
                        const hwMatrix&         sigma,
                        hwMersenneTwisterState* MTstate, 
                        unsigned long*          seed, 
                        hwMatrix&               matrix)
{
    if (!mu.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!sigma.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    int m = mu.M();
    int n = mu.N();

    if (sigma.M() != m || sigma.N() != n)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }
    hwMathStatus status = matrix.Dimension(m, n, hwMatrix::REAL);

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

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    double sigma_temp;
    int size = m * n;
    hwNormal normDist(MTstate);

    for (int i = 0; i < size; i++)
    {
        sigma_temp = sigma(i);
        if (sigma_temp < nearZero)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return status(HW_MATH_ERR_NONPOSITIVE, 2);
        }
        matrix(i) = exp(mu(i) + sigma_temp * normDist.GetDeviate());
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return status;
}
//------------------------------------------------------------------------------
// Weibull distribution probability density function
//------------------------------------------------------------------------------
hwMathStatus WeibullPDF(double x, double a, double b, double& density)
{
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    if (b < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }
    hwWeibull weibullDist(a, b);

    density = weibullDist.Pdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Weibull distribution cumulative density function
//------------------------------------------------------------------------------
hwMathStatus WeibullCDF(double x, double a, double b, double& prob)
{
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    if (b < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }
    hwWeibull weibullDist(a, b);

    prob = weibullDist.Cdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Weibull distribution inverse cumulative density function
//------------------------------------------------------------------------------
hwMathStatus WeibullInvCDF(double prob, double a, double b, double& x)
{
    if (prob < MACHEP || prob > 1.0 - MACHEP)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDPROB, 1);
    }
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    if (b < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }
    hwWeibull weibullDist(a, b);

    x = weibullDist.CdfInv(prob);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Weibull utility function for method of moments estimation
//------------------------------------------------------------------------------
static hwMathStatus WeibullFitUtil_MM(const hwMatrix& b, 
                                      const hwMatrix& X,
                                      const hwMatrix* userData, 
                                      hwMatrix&       res)
{
    double cv2 = (*userData)(0);    // coefficient of variation^2
	double gn = GammaFunc(1.0 + 2.0 / b(0));
	double gd = GammaFunc(1.0 + 1.0 / b(0));

	res(0) = gn / (gd * gd) - 1.0 - cv2;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Weibull distribution log-likelihood function (actually its negative)
//------------------------------------------------------------------------------
static hwMathStatus WeibullLogLike(const hwMatrix& P, 
                                   const hwMatrix* userData, 
                                   double&         neglogL)
{
	int n = userData->Size();
	double a = P(0);
	double b = P(1);
    double x;
	double term1 = 0.0;
	double term2 = 0.0;
	double term3;

    for (int i = 0; i < n; i++)
    {
        x = (*userData)(i);
        term1 += pow(x, b);
        term2 += log(x);
    }

    term1 /= pow(a, b);
    term2 *= (b - 1.0);
    term3 = n * log(b) - n * b * log(a);
	neglogL = term1 - term2 - term3;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Gradient of WeibullLogLike
//------------------------------------------------------------------------------
static hwMathStatus WeibullLLGrad(const hwMatrix& P, 
                                  const hwMatrix* userData, 
                                  hwMatrix&       grad)
{
	int n = userData->Size();
	double a = P(0);
	double b = P(1);
    double x;
	double term1 = 0.0;
	double term2 = 0.0;
	double term3 = 0.0;

    for (int i = 0; i < n; i++)
    {
        x = (*userData)(i);
        term1 += pow(x, b);
        term2 += log(x);
        term3 += (log(x) - log(a)) * pow(x, b);
    }

    term1 /= pow(a, b);
    term3 /= pow(a, b);

    grad(0) = b / a * (n - term1);                  // pLpa
    grad(1) = -n / b + n * log(a) - term2 + term3;  // pLpb

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Weibull distribution parameter estimation function
//------------------------------------------------------------------------------
hwMathStatus WeibullFit(const hwMatrix& data, 
                        double&         aHat, 
                        double&         bHat,
                        hwMatrix*       aCI, 
                        hwMatrix*       bCI)
{
    hwMathStatus status;

    if (!data.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 1);
    }
    int n = data.Size();

    // allocate output
    if (aCI)
    {
        status = aCI->Dimension(2, hwMatrix::REAL);

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
    }

    if (bCI)
    {
        status = bCI->Dimension(2, hwMatrix::REAL);

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
    }

    // compute MLE
    if (n == 0)
    {
        aHat = std::numeric_limits<double>::quiet_NaN();
        bHat = std::numeric_limits<double>::quiet_NaN();

        if (aCI)
        {
            aCI->SetElements(aHat);
        }
        if (bCI)
        {
            bCI->SetElements(bHat);
        }
        return status;
    }

    if (n == 1)
    {
        return status(HW_MATH_ERR_TOOFEWPOINTS, 1);
    }
    double min;
    double max;

    status = Min(data, &min);

    if (!status.IsOk())
    {
        if (status != HW_MATH_ERR_EMPTYMATRIX)
        {
            return status;
        }
        status = hwMathStatus();    // reset and continue
    }
    else if (min <= 0.0)
    {
        return status(HW_MATH_ERR_NONPOSITIVE, 1);
    }

    status = Max(data, &max);

    if (!status.IsOk())
    {
        if (status != HW_MATH_ERR_EMPTYMATRIX)
        {
            return status;
        }
        status = hwMathStatus();    // reset and continue
    }
    else if (min == max)
    {
        return status(HW_MATH_ERR_ZERORANGE, 1);
    }

    // first fit with method of moments (may switch to regression fit)
    double mean;
    double var;
    hwMatrix b(1, hwMatrix::REAL);
    hwMatrix cv2(1, hwMatrix::REAL);    // coefficient of variation^2

    status = Mean(data, mean);

    if (!status.IsOk())
        return status;

    status = Variance(data, var);
    if (!status.IsOk())
    {
        return status;
    }

    cv2(0)          = var / (mean*mean);
    b(0)            = mean / sqrt(var);  // = 1/sqrt(cv2(0)) - initial estimate
    int maxIter     = 200;
    int maxFuncEval = 1000;

    status = NLSolve(WeibullFitUtil_MM, nullptr, b, min, maxIter, maxFuncEval,
                     1, 1.0e-6, 1.0e-6, nullptr, nullptr, &cv2);

    if (!status.IsOk() && !status.IsInfoMsg() && !status.IsWarning())
    {
        status.ResetArgs();
        return status(HW_MATH_ERR_NOTCONVERGE);
    }

    status = hwMathStatus();

    aHat = mean / GammaFunc(1.0 + 1.0 / b(0));
    bHat = b(0);

    // next fit with maximum likelihood estimation (MLE)
    hwMatrix param(2, hwMatrix::REAL);

    param(0) = aHat;
    param(1) = bHat;
    maxIter     = 200;
    maxFuncEval = 1000;

    status = FMinUncon(WeibullLogLike, WeibullLLGrad, param, min, maxIter,
                       maxFuncEval, 1.0e-6, 1.0e-6, nullptr, nullptr, &data);

    if (!status.IsOk() && !status.IsInfoMsg() && !status.IsWarning())
    {
        status.ResetArgs();
        return status(HW_MATH_ERR_NOTCONVERGE);
    }

    status = hwMathStatus();

    aHat = param(0);
    bHat = param(1);

    // allocate output
    if (!aCI && !bCI)
    {
        return status;
    }
    // compute MLE confidence intervals
    double temp1;
    double temp2;
    double term1 = 0.0;
    double term2 = 0.0;
    double term3 = 0.0;
    double x;
    double z;
    hwMatrix H(2, 2, hwMatrix::REAL);
    hwMatrix cov;

    for (int i = 0; i < n; i++)
    {
        x = data(i);
        temp1 = pow(x, bHat);
        temp2 = log(x) - log(aHat);
        term1 += temp1;
        term2 += temp1 * temp2;
        term3 += temp1 * temp2 * temp2;
    }

    temp1 = pow(aHat, bHat);
    term1 /= temp1;
    term2 /= temp1;
    term3 /= temp1;

    H(0, 0) = bHat / (aHat*aHat) * (n - (bHat+1) * term1);
    H(0, 1) = -1 / aHat * (n - term1 - bHat * term2);
    H(1, 0) = H(0, 1);
    H(1, 1) = -n / (bHat*bHat) - term3;
    hwMathStatus status2 = cov.Inverse(-H);				    // Fisher information matrix

    if (!status2.IsOk())
    {
        status2.ResetArgs();
        return status2;
    }

    status2 = NormInvCDF(0.025, 0.0, 1.0, z);    // 95% confidence interval

    if (!status2.IsOk())
    {
        status2.ResetArgs();
        return status2;
    }

    if (aCI)
    {
        (*aCI)(0) = aHat * exp(z/aHat*sqrt(cov(0,0)));
        (*aCI)(1) = aHat * exp(-z/aHat*sqrt(cov(0,0)));
    }

    if (bCI)
    {
        (*bCI)(0) = bHat * exp(z/bHat*sqrt(cov(1,1)));
        (*bCI)(1) = bHat * exp(-z/bHat*sqrt(cov(1,1)));
    }

    return status;
}
//------------------------------------------------------------------------------
// Generates Weibull distribution random number, a single value from Weibull(a,b)
//------------------------------------------------------------------------------
hwMathStatus WeibullRnd(double                  a, 
                        double                  b, 
                        hwMersenneTwisterState* MTstate,
                        unsigned long*          seed, 
                        double&                 value)
{
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }
    if (b < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    hwWeibull weibullDist(a, b, MTstate);

    value = weibullDist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Weibull distribution random numbers, a matrix of values
// from Weibull(a,b)
//------------------------------------------------------------------------------
hwMathStatus WeibullRnd(double                  a, 
                        double                  b, 
                        hwMersenneTwisterState* MTstate,
                        unsigned long*          seed, 
                        hwMatrix&               matrix)
{
    if (a < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }
    if (b < nearZero)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int m = matrix.M();
    int n = matrix.N();
    int size = m * n;
    hwWeibull weibullDist(a, b, MTstate);

    for (int i = 0; i < size; i++)
    {
        matrix(i) = weibullDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Weibull distribution random numbers, a matrix of values
// from Weibull(a,b)
//------------------------------------------------------------------------------
hwMathStatus WeibullRnd(const hwMatrix&         A, 
                        const hwMatrix&         B,
                        hwMersenneTwisterState* MTstate, 
                        unsigned long*          seed, 
                        hwMatrix&               matrix)
{
    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!B.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    int m = A.M();
    int n = A.N();

    if (B.M() != m || B.N() != n)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }
    hwMathStatus status = matrix.Dimension(m, n, hwMatrix::REAL);

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

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int size = m * n;
    double A_temp;
    double B_temp;

    for (int i = 0; i < size; i++)
    {
        A_temp = A(i);
        B_temp = B(i);

        if (A_temp < nearZero)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
        }

        if (B_temp < nearZero)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
        }
        hwWeibull weibullDist(A_temp, B_temp, MTstate);

        matrix(i) = weibullDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return status;
}
//------------------------------------------------------------------------------
// Poisson distribution probability density function
//------------------------------------------------------------------------------
hwMathStatus PoissonPDF(int x, double lambda, double& density)
{
    if (x < 0)
    {
        return hwMathStatus(HW_MATH_ERR_NONNONNEGINT, 1);
    }
    if (lambda <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    hwPoisson poissonDist(lambda);
    return poissonDist.Pdf(x, density);
}
//------------------------------------------------------------------------------
// Poisson distribution cumulative density function
//------------------------------------------------------------------------------
hwMathStatus PoissonCDF(double x, double lambda, double& prob)
{
    if (x < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONNONNEGINT, 1);
    }
    if (lambda <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    hwPoisson poissonDist(lambda);

    prob = poissonDist.Cdf(x);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Poisson distribution inverse cumulative density function
//------------------------------------------------------------------------------
hwMathStatus PoissonInvCDF(double prob, double lambda, int& x)
{
    if (prob < MACHEP || prob > 1.0 - MACHEP)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDPROB, 1);
    }
    if (lambda <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONNONNEGINT, 2);
    }
    hwPoisson poissonDist(lambda);

    return poissonDist.CdfInv(prob, x);
}
//------------------------------------------------------------------------------
// Poisson distribution parameter estimation function
//------------------------------------------------------------------------------
hwMathStatus PoissonFit(const hwMatrix& data, 
                        double&         lambdaHat, 
                        hwMatrix*       lambdaCI)
{
    hwMathStatus status;

    if (!data.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsEmptyOrVector())
    {   
        return status(HW_MATH_ERR_VECTOR, 1);
    }
    int n = data.Size();

    // allocate output
    if (lambdaCI)
    {
        lambdaCI->Dimension(2, hwMatrix::REAL);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
            {
                status.SetArg1(3);
            }
            else
            {
                status.ResetArgs();
            }
            return status;
        }
    }

    // compute MLE
    if (n == 0)
    {
        lambdaHat = std::numeric_limits<double>::quiet_NaN();
        if (lambdaCI)
        {
            lambdaCI->SetElements(lambdaHat);
        }
        return status;
    }

    if (n == 1)
    {
        return status(HW_MATH_ERR_TOOFEWPOINTS, 1);
    }
    double min;

    status = Min(data, &min);

    if (!status.IsOk())
    {
        if (status != HW_MATH_ERR_EMPTYMATRIX)
        {
            return status;
        }
        status = hwMathStatus();    // reset and continue
    }
    else if (min < 0.0)
    {
        return status(HW_MATH_ERR_NONNONNEGINT, 1);
    }

    Mean(data, lambdaHat);

    // compute confidence intervals
    if (lambdaCI)
    {
        int size = data.Size();
        double sum = lambdaHat * size;

        if (size < 100)
        {
            status = Chi2InvCDF(0.025, 2 * (int) sum, (*lambdaCI)(0));
            status = Chi2InvCDF(0.975, 2 * ((int) sum + 1), (*lambdaCI)(1));
            (*lambdaCI) *= 0.5;
        }
        else
        {
            // normal approximation to ChiSquared
            status = NormInvCDF(0.025, sum, sqrt(sum), (*lambdaCI)(0));
            status = NormInvCDF(0.975, sum, sqrt(sum), (*lambdaCI)(1));
        }

        (*lambdaCI) /= (double) size;
    }

    for (int i = 0; i < n; i++)
    {
        if (IsInteger(data(i), 1.0e-12) == HW_MATH_ERR_NONINTEGER)
        {
            status(HW_MATH_WARN_NONNONNEGINT, 1);
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Generates Poisson distribution random number, a single value with parameter
// lambda
//------------------------------------------------------------------------------
hwMathStatus PoissonRnd(double                  lambda, 
                        hwMersenneTwisterState* MTstate,
                        unsigned long*          seed, 
                        double&                 value)
{
    if (lambda <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    hwPoisson poissonDist(lambda, MTstate);

    value = poissonDist.GetDeviate();

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Poisson distribution random numbers, a matrix of values with
// parameter lambda
//------------------------------------------------------------------------------
hwMathStatus PoissonRnd(double                  lambda, 
                        hwMersenneTwisterState* MTstate,
                        unsigned long*          seed, 
                        hwMatrixI&              matrix)
{
    if (lambda <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int m = matrix.M();
    int n = matrix.N();
    int size = m * n;
    hwPoisson poissonDist(lambda, MTstate);

    for (int i = 0; i < size; i++)
    {
        matrix(i) = poissonDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Generates Poisson distribution random numbers, a matrix of values with
// parameters Lambda
//------------------------------------------------------------------------------
hwMathStatus PoissonRnd(const hwMatrix&         Lambda, 
                        hwMersenneTwisterState* MTstate,
                        unsigned long*          seed, 
                        hwMatrixI&              matrix)
{
    if (!Lambda.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    int m = Lambda.M();
    int n = Lambda.N();

    hwMathStatus status = matrix.Dimension(m, n, hwMatrixI::REAL);

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

    // prepare random number generator
    bool createState = false;
    if (!MTstate)
    {
        MTstate = new hwMersenneTwisterState;
        createState = true;
    }

    if (seed)
    {
        MTstate->Initialize(*seed);
    }
    else if (!MTstate->Initialized())
    {
        MTstate->Initialize(GetSeedFromClock());
    }

    // generate random values
    int size = m * n;
    double lambda;

    for (int i = 0; i < size; i++)
    {
        lambda = Lambda(i);
        if (lambda <= 0.0)
        {
            if (createState)
            {
                delete MTstate;
                MTstate = nullptr;
            }
            return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
        }
        hwPoisson poissonDist(lambda, MTstate);

        matrix(i) = poissonDist.GetDeviate();
    }

    if (createState)
    {
        delete MTstate;
        MTstate = nullptr;
    }

    return status;
}
//------------------------------------------------------------------------------
// Generate a random permutation vector on [1:max]
//------------------------------------------------------------------------------
hwMathStatus RandPerm(int                     max,
                      int                     numPts,
                      hwMersenneTwisterState* pMTState,
                      hwMatrixI&              permVec)
{
    // check inputs
    if (max < 0)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 1);
    }

    if (numPts < 0)
    {
        return hwMathStatus(HW_MATH_ERR_NONNONNEGINT, 2);
    }

    if (numPts > max)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 2);
    }

    if (!pMTState || !pMTState->Initialized())
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 4);
    }

    // Fisher-Yates shuffle
    hwUniform uniDist(pMTState);
    hwMathStatus status = permVec.Dimension(1, numPts, hwMatrixI::REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (numPts == max)
    {
        for (int i = 1; i <= numPts; ++i)
            permVec(i - 1) = i;

        for (int i = 0; i < numPts - 1; ++i)
        {
            // get random integer on [i, numPts-1] and reorder
            int indx      = static_cast<int> (floor(uniDist.GetDeviate() * (numPts - i))) + i;
            int temp      = permVec(indx);
            permVec(indx) = permVec(i);
            permVec(i)    = temp;
        }
    }
    else if (numPts >= max / 5)
    {
        hwMatrixI index(max, hwMatrixI::REAL);

        for (int i = 1; i <= max; ++i)
            index(i - 1) = i;

        for (int i = 0; i < numPts; ++i)
        {
            // get random integer on [i, max-1] and reorder
            int indx    = static_cast<int> (floor(uniDist.GetDeviate() * (max - i))) + i;
            int temp    = index(indx);
            index(indx) = index(i);
            permVec(i)  = temp;
        }
    }
    else
    {
        // Entries > numPnts are tracked with a map
        std::unordered_map<int, int> map(numPts);
        std::vector<int> idx(numPts);

        for (int i = 0; i < numPts; i++)
            idx[i] = i;

        for (int i = 0; i < numPts; i++)
        {
            int k = static_cast<int> (floor(uniDist.GetDeviate() * (max - i))) + i;

            if (k < numPts)
            {
                std::swap(idx[i], idx[k]);
            }
            else
            {
                if (map.find(k) == map.end())
                    map[k] = k;

                std::swap(idx[i], map[k]);
            }
        }

        for (int i = 0; i < numPts; i++)
            permVec(i) = idx[i] + 1;
    }

    return status;
}
//------------------------------------------------------------------------------
// Generate a random permutation vector on [1:max]
//------------------------------------------------------------------------------
hwMathStatus RandPerm(int64_t                 max,
                      int                     numPts,
                      hwMersenneTwisterState* pMTState,
                      hwMatrixI64&            permVec)
{
    // check inputs
    if (max < 0)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 1);
    }

    if (numPts < 0)
    {
        return hwMathStatus(HW_MATH_ERR_NONNONNEGINT, 2);
    }

    if (numPts > max)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 2);
    }

    if (!pMTState || !pMTState->Initialized())
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 3);
    }

    // Fisher-Yates shuffle
    hwUniform uniDist(pMTState);
    hwMathStatus status = permVec.Dimension(1, numPts, hwMatrixI64::REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    // Entries > numPnts are tracked with a map
    std::unordered_map<int64_t, int64_t> map(numPts);
    std::vector<int64_t> idx(numPts);

    for (int i = 0; i < numPts; i++)
        idx[i] = i;

    for (int i = 0; i < numPts; i++)
    {
        int64_t k = static_cast<int64_t> (floor(uniDist.GetDeviate() * (max - i))) + i;

        if (k < numPts)
        {
            std::swap(idx[i], idx[k]);
        }
        else
        {
            if (map.find(k) == map.end())
                map[k] = k;

            std::swap(idx[i], map[k]);
        }
    }

    for (int i = 0; i < numPts; i++)
        permVec(i) = idx[i] + 1;

    return status;
}
//------------------------------------------------------------------------------
// D'Agostino-Pearson omnibus normality test
//------------------------------------------------------------------------------
hwMathStatus NormalityTestDP(const hwMatrix& data, double sigLevel, bool& reject)
{
    hwMatrix moment;
    hwMathStatus status = FourMoments(data, moment);

    if (!status.IsOk())
    {
        return status;
    }
    int n = data.Size();
    if (n < 20)
    {
        return status(HW_MATH_ERR_NORMTESTPNTS, 1);
    }
    if (sigLevel <= 0.0 || sigLevel >= 0.5)
    {
        status(HW_MATH_WARN_STATTESTALPHA, 2);
    }

    // compute skewness transform
    double factor = n / (n - 1.0);
    double skew   = moment(2);
    double beta1  = skew * skew * factor;

    double Y = sqrt(beta1 * (n+1.0)* (n+3.0) / (6.0 * (n-2.0)));
    double b2 = 3.0 * (n*n+27.0*n-70.0)/((n-2.0)*(n+5.0))*(n+1.0)/(n+7.0)*(n+3.0)/(n+9.0);
    double Wsq = -1.0 + sqrt(2.0 * (b2 - 1.0));
    double delta = 1.0 / sqrt(0.5*log(Wsq));
    double alpha = sqrt(2.0 / (Wsq - 1.0));
    double value = Y / alpha;
    double Z1 = delta * log(sqrt(value*value + 1.0) + value);

    // compute kurtosis transform
    double kurt = moment(3);
    double beta2 = (kurt + 3.0) * factor;
    double ev_beta2 = 3.0 * (n-1.0) / (n+1.0);
    double var_beta2 = 24.0 * n / (n+1.0) * (n-2.0) / (n+1.0) * (n-3.0) / (n+3.0) / (n+5.0);
    double x = (beta2 - ev_beta2) / sqrt(var_beta2);
    double b1sqrt = 6.0 * (n*n-5.0*n+2.0) / ((n+7.0)*(n+9.0)) * sqrt(6.0*(n+3.0)/n*(n+5.0)/(n-2.0)/(n-3.0));
    double A = 6.0 + 8.0 / b1sqrt * (2.0 / b1sqrt + sqrt(4.0 / b1sqrt + 1.0));

    value = (1.0 - 2.0 / A) / (1.0 + x * sqrt(2.0 / (A - 4.0)));
    double Z2 = (1.0 - 2.0 / A - pow(value,1.0/3.0)) / sqrt(2.0 / (9.0 * A));

    // run test
    double testStat = Z1 * Z1 + Z2 * Z2;   // chiSq with 2 dof is exponential
    double chiSqValue = -2.0 * log(0.5 * sigLevel);

    if (testStat > chiSqValue)
    {
        reject = true;
    }
    else
    {
        reject = false;
    }
    return status;
}
