/**
* @file hwBeta.cxx
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
#include "hwBeta.h"

#include <math.h>

#include "Globals.h"
#include "hwGamma.h"
#include "hwNormal.h"
#include "SpecialFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwBeta::hwBeta(double                  alpha,
               double                  beta,
               hwMersenneTwisterState* pMTState)
    : m_pGamma1 (nullptr)
    , m_pGamma2 (nullptr)
    , m_alpha   (alpha)
    , m_beta    (beta)
{
    if (alpha <= 0.0 || beta <= 0.0)
    {
        return;
    }

    if (pMTState)
    {
        m_pGamma1 = new hwGamma(m_alpha, pMTState);
        m_pGamma2 = new hwGamma(m_beta,  pMTState);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwBeta::~hwBeta()
{
    delete m_pGamma1;
    delete m_pGamma2;
}
//------------------------------------------------------------------------------
// Returns probability density function
//------------------------------------------------------------------------------
double hwBeta::Pdf(double x)
{
    if (x < 0 || x > 1) 
    {
        return 0.0;
    }

    double density  = 0.0;
    double nearZero = 1.0e-12;
    double beta     = BetaFunc(m_alpha, m_beta);
    if (fabs(beta) > nearZero)
    {
        density = pow(x, m_alpha - 1.0) * pow(1.0 - x, m_beta - 1.0) / beta;
    }

    return density;
}
//------------------------------------------------------------------------------
// Returns cumulative density function
//------------------------------------------------------------------------------
double hwBeta::Cdf(double x)
{
    return IncBeta(x);
}
//------------------------------------------------------------------------------
// Returns inverse cumulative density function
//------------------------------------------------------------------------------
double hwBeta::CdfInv(double prob)
{
    return IncBetaInv(prob);
}
//------------------------------------------------------------------------------
// Returns distribution mean
//------------------------------------------------------------------------------
double hwBeta::Mean()
{
    double value = m_alpha + m_beta;

    return m_alpha / value;
}
//------------------------------------------------------------------------------
// Compute the variance of the distribution
//------------------------------------------------------------------------------
double hwBeta::Variance()
{
    double value = m_alpha + m_beta;

    return m_alpha * m_beta / (value * value * (1 + value));
}
//------------------------------------------------------------------------------
// Returns distribution variance
//------------------------------------------------------------------------------
double hwBeta::Mode()
{
    double mode = (m_alpha - 1.0) / (m_alpha + m_beta - 2.0);

    if (mode < 0.0)
    {
        mode = 0.0;
    }
    else if (mode > 1.0)
    {
        mode = 1.0;
    }
    return mode;
}
//------------------------------------------------------------------------------
// Compute a random deviate from the distribution
//------------------------------------------------------------------------------
double hwBeta::GetDeviate()
{
    double nearZero = 1.0e-12;
    double g1       = m_pGamma1->GetDeviate();
    double g2       = m_pGamma2->GetDeviate();
    double b        = g1 + g2;	

    if (b > nearZero)
    {
        return g1 / b;
    }

    return 0.0;
}
//------------------------------------------------------------------------------
// Compute natural log of Beta(a, b)
//------------------------------------------------------------------------------
double hwBeta::LogBeta(double a, double b)
{
    return GammaLog(a) + GammaLog(b) - GammaLog(a + b);
}
//------------------------------------------------------------------------------
// Incomplete Beta Integral
// Returns incomplete beta integral of the arguments, evaluated from zero to x.  
// The function is defined as
//                  x
//     -            -
//    | (a+b)      | |  a-1     b-1
//  -----------    |   t   (1-t)   dt.
//   -     -     | |
//  | (a) | (b)   -
//                 0
//
// The domain of definition is 0 <= x <= 1. In this implementation a and b are 
// restricted to positive values. The integral from x to 1 may be obtained by 
// the symmetry relation
//    1 - incbet( a, b, x )  =  incbet( b, a, 1-x ).
//
// The integral is evaluated by a continued fraction expansion or, when b*x is 
// small, by a power series.
// ACCURACY:
// Tested at uniformly distributed random points (a,b,x) with a and b in 
// "domain" and x between 0 and 1.
//                                        Relative error
// arithmetic   domain     # trials      peak         rms
//    IEEE      0,5         10000       6.9e-15     4.5e-16
//    IEEE      0,85       250000       2.2e-13     1.7e-14
//    IEEE      0,1000      30000       5.3e-12     6.3e-13
//    IEEE      0,10000    250000       9.3e-11     7.1e-12
//    IEEE      0,100000    10000       8.7e-10     4.8e-11
// Outputs smaller than the IEEE gradual underflow threshold were excluded from 
// these statistics.
//
// ERROR MESSAGES:
//   message         condition      value returned
// incbet domain      x<0, x>1          0.0
// incbet underflow                     0.0
//
//Cephes Math Library, Release 2.8:  June, 2000
//Copyright 1984, 1995, 2000 by Stephen L. Moshier
//The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Incomplete beta function (Cephes library function)
//------------------------------------------------------------------------------
double hwBeta::IncBeta(double xx)
{
    return IncBeta(m_alpha, m_beta, xx);
}
//------------------------------------------------------------------------------
// Incomplete beta function (Cephes library function)
//------------------------------------------------------------------------------
double hwBeta::IncBeta(double aa, double bb, double xx)
{
    if (aa <= 0.0 || bb <= 0.0)
    {
        return 0.0;
    }
    if (xx <= 0.0)
    {
        return 0.0;
    }
    if (xx >= 1.0)
    {
        return 1.0;
    }

    int    flag = 0;
    double a;
    double b;
    double t;
    double x;
    double xc;
    double w;
    double y;

    if (bb * xx <= 1.0 && xx <= 0.95)
    {
        t = Pseries(aa, bb, xx);
    }
    else
    {
        w = 1.0 - xx;

        // Reverse a and b if x is greater than the mean
        if (xx > aa / (aa+bb))
        {
            flag = 1;
            a    = bb;
            b    = aa;
            xc   = xx;
            x    = w;
        }
        else
        {
            a  = aa;
            b  = bb;
            xc = w;
            x  = xx;
        }

        if (flag == 1 && b * x <= 1.0 && x <= 0.95)
        {
            t = Pseries(a, b, x);
        }
        else
        {
            // Choose expansion for better convergence
            y = x * (a+b-2.0) - (a-1.0);

            w = (y < 0.0) ? Incbcf(a, b, x) : Incbd(a, b, x) / xc;

            // Multiply w by the factor
            //	 a      b   _             _     _
            //	x  (1-x)   | (a+b) / ( a | (a) | (b) )

            y = a * log(x);
            t = b * log(xc);

            if (a + b < MAXGAM && fabs(y) < MAXLOG && fabs(t) < MAXLOG)
            {
                t = pow(xc,b);
                t *= pow(x,a);
                t /= a;
                t *= w;
                t *= GammaFunc(a+b) / (GammaFunc(a) * GammaFunc(b));
            }
            else
            {
                // Resort to logarithms
                y += t + GammaLog(a+b) - GammaLog(a) - GammaLog(b);
                y += log(w/a);

                t = (y < MINLOG) ? 0.0 : exp(y);
            }
        }
    }

    if (flag == 1)
    {
        t = (t <= MACHEP) ? 1.0 - MACHEP : 1.0 - t;
    }

    return t;
}
//------------------------------------------------------------------------------
// Inverse Incomplete Beta Integral
// Given y, the function finds x such that
//  incbet( a, b, x ) = y .
//
// The routine performs interval halving or Newton iterations to find the
// root of incbet(a,b,x) - y = 0.
//
// ACCURACY:
//                      Relative error:
//                x     a,b
// arithmetic   domain  domain  # trials    peak       rms
//    IEEE      0,1    .5,10000   50000    5.8e-12   1.3e-13
//    IEEE      0,1   .25,100    100000    1.8e-13   3.9e-15
//    IEEE      0,1     0,5       50000    1.1e-12   5.5e-15
//    VAX       0,1    .5,100     25000    3.5e-14   1.1e-15
// With a and b constrained to half-integer or integer values:
//    IEEE      0,1    .5,10000   50000    5.8e-12   1.1e-13
//    IEEE      0,1    .5,100    100000    1.7e-14   7.9e-16
// With a = .5, b constrained to half-integer or integer values:
//    IEEE      0,1    .5,10000   10000    8.3e-11   1.0e-11
//
//Cephes Math Library Release 2.8:  June, 2000
//Copyright 1984, 1996, 2000 by Stephen L. Moshier
//The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Returns inverse incomplete beta function value
//------------------------------------------------------------------------------
double hwBeta::IncBetaInv(double yy0)
{
    if (yy0 <= 0)
    {
        return 0.0;
    }
    if (yy0 >= 1.0)
    {
        return 1.0;
    }

    double a;
    double b;
    double y0;
    double x;
    double y;
    double dithresh = 1.0e-4;
    double di;
    double xt;
    double x0 = 0.0;
    double yl = 0.0;
    double x1 = 1.0;
    double yh = 1.0;

    int rflg;
    int dir;
    int nflg = 0;

    hwNormal normDist;

    double d;
    double lgm;
    double yp;

    if (m_alpha <= 1.0 || m_beta <= 1.0)
    {
        dithresh = 1.0e-6;
        rflg     = 0;
        a        = m_alpha;
        b        = m_beta;
        y0       = yy0;
        x        = a/(a+b);
        y        = IncBeta( a, b, x );
        goto ihalve;
    }

    // approximation to inverse function
    yp = -normDist.CdfInv(yy0);

    if (yy0 > 0.5)
    {
        rflg = 1;
        a    = m_beta;
        b    = m_alpha;
        y0   = 1.0 - yy0;
        yp   = -yp;
    }
    else
    {
        rflg = 0;
        a    = m_alpha;
        b    = m_beta;
        y0   = yy0;
    }

    lgm = (yp * yp - 3.0)/6.0;
    x = 2.0 / (1.0 / (2.0*a-1.0) + 1.0 / (2.0*b-1.0));
    d = yp * sqrt( x + lgm ) / x
        - (1.0 / (2.0*b-1.0) - 1.0 / (2.0*a-1.0))
        * (lgm + 5.0 / 6.0 - 2.0 / (3.0*x));
    d += d;

    if (d < MINLOG)
    {
        x = 1.0;
        goto under;
    }

    x = a / (a + b * exp(d));
    y = IncBeta(a, b, x);
    yp = (y - y0) / y0;

    if (fabs(yp) < 0.2)
    {
        goto newt;
    }

    // Resort to interval halving if not close enough
ihalve:

    dir = 0;
    di  = 0.5;

    for (int i = 0; i < 100; i++)
    {
        if (i != 0)
        {
            x = x0 + di * (x1 - x0);

            if (x == 1.0)
            {
                x = 1.0 - MACHEP;
            }
            if (x == 0.0)
            {
                di = 0.5;
                x = x0  +  di * (x1 - x0);

                if (x == 0.0)
                {
                    goto under;
                }
            }

            y  = IncBeta(a, b, x);
            yp = (x1 - x0) / (x1 + x0);

            if (fabs(yp) < dithresh)
            {
                goto newt;
            }
            yp = (y - y0) / y0;

            if (fabs(yp) < dithresh)
            {
                goto newt;
            }
        }

        if (y < y0)
        {
            x0 = x;
            yl = y;

            if (dir < 0)
            {
                dir = 0;
                di = 0.5;
            }
            else if (dir > 3)
            {
                di = 1.0 - (1.0 - di) * (1.0 - di);
            }
            else if (dir > 1)
            {
                di = 0.5 * di + 0.5; 
            }
            else
            {
                di = (y0 - y) / (yh - yl);
            }
            dir += 1;

            if (x0 > 0.75)
            {
                if (rflg == 1)
                {
                    rflg = 0;
                    a    = m_alpha;
                    b    = m_beta;
                    y0   = yy0;
                }
                else
                {
                    rflg = 1;
                    a    = m_beta;
                    b    = m_alpha;
                    y0   = 1.0 - yy0;
                }

                x  = 1.0 - x;
                y  = IncBeta( a, b, x );
                x0 = 0.0;
                yl = 0.0;
                x1 = 1.0;
                yh = 1.0;
                goto ihalve;
            }
        }
        else
        {
            x1 = x;

            if (rflg == 1 && x1 < MACHEP)
            {
                x = 0.0;
                goto done;
            }

            yh = y;

            if(dir > 0)
            {
                dir = 0;
                di = 0.5;
            }
            else if (dir < -3)
            {
                di = di * di;
            }
            else if (dir < -1)
            {
                di = 0.5 * di;
            }
            else
            {
                di = (y - y0) / (yh - yl);
            }
            dir -= 1;
        }
    }
    if (x0 >= 1.0)
    {
        x = 1.0 - MACHEP;
        goto done;
    }

    if (x <= 0.0)
    {

under:
        x = 0.0;
        goto done;
    }

newt:
    if (nflg)
    {
        goto done;
    }
    nflg = 1;
    lgm  = GammaLog(a+b) - GammaLog(a) - GammaLog(b);

    for (int i = 0; i < 8; i++)
    {
        // Compute the function at this point
        if (i != 0)
        {
            y = IncBeta(a, b, x);
        }
        if (y < yl)
        {
            x = x0;
            y = yl;
        }
        else if (y > yh)
        {
            x = x1;
            y = yh;
        }
        else if (y < y0)
        {
            x0 = x;
            yl = y;
        }
        else
        {
            x1 = x;
            yh = y;
        }

        if (x == 1.0 || x == 0.0)
            break;

        // Compute the derivative of the function at this point
        d = (a - 1.0) * log(x) + (b - 1.0) * log(1.0-x) + lgm;

        if (d < MINLOG)
        {
            goto done;
        }
        if (d > MAXLOG)
        {
            break;
        }
        d = exp(d);

        // Compute the step to the next approximation of x
        d = (y - y0) / d;
        xt = x - d;

        if (xt <= x0)
        {
            y = (x - x0) / (x1 - x0);
            xt = x0 + 0.5 * y * (x - x0);

            if (xt <= 0.0)
            {
                break;
            }
        }

        if (xt >= x1)
        {
            y = (x1 - x) / (x1 - x0);
            xt = x1 - 0.5 * y * (x1 - x);

            if (xt >= 1.0)
            {
                break;
            }
        }

        x = xt;
        if (fabs(d / x) < 128.0 * MACHEP)
        {
            goto done;
        }
    }

    // Did not converge
    dithresh = 256.0 * MACHEP;
    goto ihalve;

done:

    if (rflg)
    {
        x = (x <= MACHEP) ? 1.0 - MACHEP : 1.0 - x;
    }

    return x;
}
//------------------------------------------------------------------------------
// Incomplete Beta Integral (continued fraction #1)
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1996, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//
// Incomplete Beta function using a continued fraction
//------------------------------------------------------------------------------
double hwBeta::Incbcf(double a, double b, double x)
{
    double big    = 4.503599627370496e+15;
    double biginv = 2.22044604925031308085e-16;
    double xk;
    double pk;
    double qk;
    double t;

    double k1 = a;
    double k2 = a + b;
    double k3 = a;
    double k4 = a + 1.0;
    double k5 = 1.0;
    double k6 = b - 1.0;
    double k7 = k4;
    double k8 = a + 2.0;

    double pkm2   = 0.0;
    double qkm2   = 1.0;
    double pkm1   = 1.0;
    double qkm1   = 1.0;
    double ans    = 1.0;
    double r      = 1.0;
    double thresh = 3.0 * MACHEP;

    int n = 0;

    do
    {
        xk = -(x * k1 * k2) / (k3 * k4);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        xk = (x * k5 * k6) / (k7 * k8);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        if (qk != 0)
        {
            r = pk / qk;
        }
        if (r != 0)
        {
            t = fabs((ans - r)/r);
            ans = r;
        }
        else
        {
            t = 1.0;
        }

        if (t < thresh)
        {
            break;
        }
        k1 += 1.0;
        k2 += 1.0;
        k3 += 2.0;
        k4 += 2.0;
        k5 += 1.0;
        k6 -= 1.0;
        k7 += 2.0;
        k8 += 2.0;

        if (fabs(qk) + fabs(pk) > big)
        {
            pkm2 *= biginv;
            pkm1 *= biginv;
            qkm2 *= biginv;
            qkm1 *= biginv;
        }

        if (fabs(qk) < biginv || fabs(pk) < biginv)
        {
            pkm2 *= big;
            pkm1 *= big;
            qkm2 *= big;
            qkm1 *= big;
        }
    }
    while (++n < 300);

    return ans;
}
//------------------------------------------------------------------------------
// Incomplete Beta Integral (continued fraction #2)
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1996, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//
// Incomplete Beta function using a continued fraction
//------------------------------------------------------------------------------
double hwBeta::Incbd(double a, double b, double x)
{
    double big    = 4.503599627370496e+15;
    double biginv = 2.22044604925031308085e-16;
    double xk;
    double pk;
    double qk;
    double t;
    int    n = 0;

    double k1 = a;
    double k2 = b - 1.0;
    double k3 = a;
    double k4 = a + 1.0;
    double k5 = 1.0;
    double k6 = a + b;
    double k7 = a + 1.0;;
    double k8 = a + 2.0;

    double pkm2   = 0.0;
    double qkm2   = 1.0;
    double pkm1   = 1.0;
    double qkm1   = 1.0;
    double z      = x / (1.0-x);
    double ans    = 1.0;
    double r      = 1.0;
    double thresh = 3.0 * MACHEP;

    do
    {
        xk = -(z * k1 * k2) / (k3 * k4);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        xk = (z * k5 * k6) / (k7 * k8);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        if (qk != 0)
        {
            r = pk/qk;
        }
        if (r != 0)
        {
            t = fabs((ans - r) / r);
            ans = r;
        }
        else
        {
            t = 1.0;
        }
        if (t < thresh)
        {
            break;
        }
        k1 += 1.0;
        k2 -= 1.0;
        k3 += 2.0;
        k4 += 2.0;
        k5 += 1.0;
        k6 += 1.0;
        k7 += 2.0;
        k8 += 2.0;

        if (fabs(qk) + fabs(pk) > big)
        {
            pkm2 *= biginv;
            pkm1 *= biginv;
            qkm2 *= biginv;
            qkm1 *= biginv;
        }

        if (fabs(qk) < biginv || fabs(pk) < biginv)
        {
            pkm2 *= big;
            pkm1 *= big;
            qkm2 *= big;
            qkm1 *= big;
        }
    }
    while (++n < 300);

    return ans;
}
//------------------------------------------------------------------------------
// Incomplete Beta Integral (Power Series)
// Power series for incomplete beta integral.
// Use when b*x is small and x not too close to 1.
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1996, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
// Incomplete Beta function using a power series
//------------------------------------------------------------------------------
double hwBeta::Pseries(double a, double b, double x)
{
    double ai = 1.0 / a;
    double u  = (1.0 - b) * x;
    double v  = u / (a + 1.0);
    double t1 = v;
    double t  = u;
    double n  = 2.0;
    double s  = 0.0;
    double z  = MACHEP * ai;

    while (fabs(v) > z)
    {
        u = (n - b) * x / n;
        t *= u;
        v = t / (a + n);
        s += v; 
        n += 1.0;
    }

    s += t1;
    s += ai;

    u = a * log(x);

    if (a + b < MAXGAM && fabs(u) < MAXLOG)
    {
        t = GammaFunc(a+b) / (GammaFunc(a) * GammaFunc(b));
        s = s * t * pow(x,a);
    }
    else
    {
        t = GammaLog(a+b) - GammaLog(a) - GammaLog(b) + u + log(s);
        s = (t < MINLOG) ? 0.0 : exp(t);
    }

    return s;
}
