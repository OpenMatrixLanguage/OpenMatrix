/**
* @file hwGammaInv.cxx
* @date May 2009
* Copyright (C) 2009-2018 Altair Engineering, Inc.  
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
#include "hwGamma.h"

#include <math.h>

#include "GeneralFuncs.h"
#include "hwNormal.h"
#include "SpecialFuncs.h"
#include "hwUniform.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwGamma::hwGamma(double alpha, hwMersenneTwisterState* pMTState)
    : m_pUniform (nullptr)
    , m_alpha    (alpha)
{
    if (alpha <= 0.0)
    {
        return;
    }

    if (pMTState)
    {
        m_pUniform = new hwUniform(pMTState);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwGamma::~hwGamma()
{
    delete m_pUniform;
}
//------------------------------------------------------------------------------
// Returns probability density function
//------------------------------------------------------------------------------
double hwGamma::Pdf(double x)
{
    if (x < 0.0)
    {
        return 0.0;
    }

    double gamma;

    if (m_alpha < MAXGAM)
    {
        gamma = GammaFunc(m_alpha);
    }
    else
    {
        gamma = GammaLog(m_alpha);
        gamma = (gamma < MINLOG) ? 0.0 : exp(gamma);
    }

    double density  = 0.0;
    double nearZero = 1.0e-12;
    if (fabs(gamma) > nearZero)
    {
        density = pow(x, m_alpha - 1.0) * exp(-x) / gamma;
    }

    return density;
}
//------------------------------------------------------------------------------
// Returns cumulative density function
//------------------------------------------------------------------------------
double hwGamma::Cdf(double x)
{
    return IncGamma(x);
}
//------------------------------------------------------------------------------
// Returns inverse cumulative density function
//------------------------------------------------------------------------------
double hwGamma::CdfInv(double prob)
{
    return IncGammaCInv(1.0 - prob);
}
//------------------------------------------------------------------------------
// Returns random value
//------------------------------------------------------------------------------
double hwGamma::GetDeviate()
{
    int    maxLoop  = 500;
    double largeNum = 1.0e+30;
    double nearZero = 1.0e-12;
    double x        = 0.0;

    if (m_alpha < nearZero)
    {
    }
    else if (m_alpha < 1.0)
    {
        double c = 1 + m_alpha * 0.3678794411714;	// 1/e = 0.3678794411714
        double oneOverAlpha = 1.0 / m_alpha;

        double u1;
        double u2;
        double t;

        for (int i = 0; i < maxLoop; i++)
        {
            u1 = m_pUniform->GetDeviate();
            u2 = m_pUniform->GetDeviate();

            t = c * u1;

            if (t <= 1.0)
            {
                x = pow(t, oneOverAlpha);

                if (u2 <= exp(-x)) 
                    break;
            }
            else
            {
                x = -log((c - t) * oneOverAlpha);

                if (u2 <= pow(x, m_alpha - 1)) 
                    break;
            }
        }
    }
    else if (fabs(m_alpha - 1.0) < nearZero)
    {
        double u1 = m_pUniform->GetDeviate();
        if (u1 >= nearZero)
        {
            x = -log(u1);
        }
    }
    else // if (m_alpha > 1.0)
    {
        // This section uses Fishman's method. For more information see
        // Communications of the ACM Vol. 19, Number 7 (1976)
        double logU1;
        double u1;
        double u2;
        double t;

        for (int i = 0; i < maxLoop; i++)
        {
            u1    = m_pUniform->GetDeviate();
            logU1 = largeNum;

            if (u1 >= nearZero)
            {
                logU1 = -log(u1);
            }

            t = logU1 / exp(logU1 - 1);
            u2 = m_pUniform->GetDeviate();

            if (u2 <= pow(t, m_alpha - 1.0)) 
            {
                x = m_alpha * logU1;
                break;
            }
        }
    }

    return x;
}
//------------------------------------------------------------------------------
// Incomplete Gamma Integral
// The function is defined by
//
//                           x
//                            -
//                   1       | |  -t  a-1
//  igam(a,x)  =   -----     |   e   t   dt.
//                  -      | |
//                 | (a)    -
//                           0
//
// In this implementation both arguments must be positive. The integral is 
// evaluated by either a power series or continued fraction expansion, depending 
// on the relative values of a and x.
//
// left tail of incomplete gamma function:
//
//          inf.      k
//   a  -x   -       x
//  x  e     >   ----------
//           -     -
//          k=0   | (a+k+1)
//
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0,30       200000       3.6e-14     2.9e-15
//    IEEE      0,100      300000       9.9e-14     1.5e-14
//
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1995, 2000 by Stephen L. Moshier
//The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------
double hwGamma::IncGamma(double x)
{
    if (x <= 0 || m_alpha <= 0)
    {
        return 0.0;
    }
    if (x > 1.0 && x > m_alpha)
    {
        return (1.0 - IncGammaC(x));
    }

    // Compute  x**m_alpha * exp(-x) / gamma(m_alpha)
    double ax = m_alpha * log(x) - x - GammaLog(m_alpha);
    ax = exp(ax);

    // power series
    double r   = m_alpha;
    double c   = 1.0;
    double ans = 1.0;

    do
    {
        r += 1.0;
        c *= x / r;
        ans += c;
    }
    while (c / ans > MACHEP);

    return ans * ax / m_alpha;
}
//------------------------------------------------------------------------------
// Complementary Incomplete Gamma Integral
// The function is defined by
//
//                           x
//                            -
//                   1       | |  -t  a-1
//  igam(a,x)  =   -----     |   e   t   dt.
//                  -      | |
//                 | (a)    -
//                           0
//
// In this implementation both arguments must be positive. The integral is 
// evaluated by either a power series or continued fraction expansion, depending 
// on the relative values of a and x.
// ACCURACY:
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0,30       200000       3.6e-14     2.9e-15
//    IEEE      0,100      300000       9.9e-14     1.5e-14
//
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1995, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------
double hwGamma::IncGammaC(double x)
{
    if (x <= 0 || m_alpha <= 0)
    {
        return 1.0;
    }
    if (x < 1.0 || x < m_alpha)
    {
        return 1.0 - IncGamma(x);
    }
    if (IsInf_T(x))
    {
        return 0.0;
    }
    double yc;
    double r;
    double t;
    double pk;
    double qk;
    double big    = 4.503599627370496e+15;
    double biginv = 2.22044604925031308085e-16;


    double ax = m_alpha * log(x) - x - GammaLog(m_alpha);
    ax = exp(ax);

    // continued fraction
    double y    = 1.0 - m_alpha;
    double z    = x + y + 1.0;
    double c    = 0.0;
    double pkm2 = 1.0;
    double qkm2 = x;
    double pkm1 = x + 1.0;
    double qkm1 = z * x;
    double ans  = pkm1 / qkm1;

    do
    {
        c += 1.0;
        y += 1.0;
        z += 2.0;
        yc = y * c;
        pk = pkm1 * z - pkm2 * yc;
        qk = qkm1 * z - qkm2 * yc;

        if (qk != 0)
        {
            r = pk / qk;
            t = fabs((ans - r) / r);
            ans = r;
        }
        else
        {
            t = 1.0;
        }

        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        if (fabs(pk) > big)
        {
            pkm2 *= biginv;
            pkm1 *= biginv;
            qkm2 *= biginv;
            qkm1 *= biginv;
        }
    }
    while (t > MACHEP);

    if (ax)     // avoid limit issues with ans
    {
        return ans * ax;
    }

    return 0.0;
}
//------------------------------------------------------------------------------
// Inverse Complementary Incomplete Gamma Integral
// Given p, the function finds x such that
//  igamc( a, x ) = p.
// Starting with the approximate value
//         3
//  x = a t
// where t = 1 - d - ndtri(p) sqrt(d) and d = 1/9a,
// the routine performs up to 10 Newton iterations to find the root of 
// igamc(a,x) - p = 0.
// ACCURACY:
// Tested at random a, p in the intervals indicated.
//
//                a        p                      Relative error:
// arithmetic   domain   domain     # trials      peak         rms
//    IEEE     0.5,100   0,0.5       100000       1.0e-14     1.7e-15
//    IEEE     0.01,0.5  0,0.5       100000       9.0e-14     3.4e-15
//    IEEE    0.5,10000  0,0.5        20000       2.3e-13     3.8e-14
//
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1995, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------
double hwGamma::IncGammaCInv(double y0)
{
    // bound the solution
    double x0       = MAXNUM;
    double yl       = 0;
    double x1       = 0;
    double yh       = 1.0;
    double dithresh = 5.0 * MACHEP;

    // approximation to inverse function
    hwNormal normDist;

    double d = 1.0 / (9.0 * m_alpha);
    double y = (1.0 - d - normDist.CdfInv(y0) * sqrt(d));
    double x = m_alpha * y * y * y;

    double lgm = GammaLog(m_alpha);

    int dir;

    for (int i = 0; i < 10; i++)
    {
        if (x > x0 || x < x1)
        {
            goto ihalve;
        }
        y = IncGammaC(x);

        if (y < yl || y > yh)
            goto ihalve;

        if (y < y0)
        {
            x0 = x;
            yl = y;
        }
        else
        {
            x1 = x;
            yh = y;
        }

        // compute the derivative of the function at this point
        d = (m_alpha - 1.0) * log(x) - x - lgm;

        if (d < -MAXLOG)
        {
            goto ihalve;
        }
        d = -exp(d);

        // compute the step to the next approximation of x
        d = (y - y0) / d;

        if (fabs(d / x) < MACHEP)
        {
            goto done;
        }
        x -= d;
    }

    // Resort to interval halving if Newton iteration did not converge.
ihalve:

    d = 0.0625;

    if (x0 == MAXNUM)
    {
        if (x <= 0.0)
        {
            x = 1.0;
        }

        while (x0 == MAXNUM)
        {
            x = (1.0 + d) * x;
            y = IncGammaC(x);

            if (y < y0)
            {
                x0 = x;
                yl = y;
                break;
            }

            d += d;
        }
    }

    d   = 0.5;
    dir = 0;

    for (int i = 0; i < 400; i++)
    {
        x = x1  +  d * (x0 - x1);
        y = IncGammaC(x);
        lgm = (x0 - x1) / (x1 + x0);

        if (fabs(lgm) < dithresh)
        {
            break;
        }
        lgm = (y - y0) / y0;

        if (fabs(lgm) < dithresh)
        {
            break;
        }
        if (x <= 0.0)
        {
            break;
        }
        if (y >= y0)
        {
            x1 = x;
            yh = y;

            if (dir < 0)
            {
                dir = 0;
                d = 0.5;
            }
            else if (dir > 1)
            {
                d = 0.5 * d + 0.5; 
            }
            else
            {
                d = (y0 - yl) / (yh - yl);
            }
            dir += 1;
        }
        else
        {
            x0 = x;
            yl = y;

            if (dir > 0)
            {
                dir = 0;
                d = 0.5;
            }
            else if (dir < -1)
            {
                d = 0.5 * d;
            }
            else
            {
                d = (y0 - yl) / (yh - yl);
            }
            dir -= 1;
        }
    }

done:
    return x;
}
