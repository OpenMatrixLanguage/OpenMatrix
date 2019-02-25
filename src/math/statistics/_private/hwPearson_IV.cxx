/**
* @file hwPearson_IV.cxx
* @date October 2014
* Copyright (C) 2014-2018 Altair Engineering, Inc.  
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
#include "hwPearson_IV.h"

#include <math.h>

#include "hwMatrix.h"
#include "hwUniform.h"
#include "SpecialFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwPearson_IV::hwPearson_IV(double                  m,
                           double                  nu,
                           hwMersenneTwisterState* pMTState)
    : m_m        (m)
    , m_nu       (nu)
    , m_pUniform (nullptr)
{
    m_mean     = -m_nu / (2.0*m_m - 2.0);
    m_variance = (1.0 + m_mean*m_mean) / (2.0*m_m-3.0);
    m_k        = type4norm();

    if (pMTState)
    {
        m_pUniform = new hwUniform(pMTState);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwPearson_IV::~hwPearson_IV()
{
    delete m_pUniform;
}
//------------------------------------------------------------------------------
// Probability density function
//------------------------------------------------------------------------------
double hwPearson_IV::Pdf(double x)
{
    return m_k * pow((1.0+x*x), -m_m) * exp(-m_nu*atan(x));
}
//------------------------------------------------------------------------------
// Gauss-Legender quadrature reproduced from /hwcommon/math/calculus
//------------------------------------------------------------------------------
double hwPearson_IV::GaussQuad(double          a, 
                               double          b, 
                               const hwMatrix& X,
                               const hwMatrix& W)
{
    int    n    = W.M();
    double c    = 0.5 * (b + a);
    double d    = 0.5 * (b - a);
    double area = 0.0;

    for (int i = 0; i < n; i++)
    {
        area += W(i) * Pdf(c + d * X(i));
    }
    area *= d;

    return area;
}
//------------------------------------------------------------------------------
// Returns cumulative density function
//------------------------------------------------------------------------------
double hwPearson_IV::Cdf(double x)
{
    // from Stavroyiannis
    double tail = 1.1;
    // note : Stavroyiannis indicates tail=sqrt(3), but
    // emperical results show that 1.0 works. Using 1.1
    // here because convergence is very slow near 1.0.

    if (x <= -tail)
    {
        int n = 0;
        double pdf;
        double yn = -(1.0+x*x)/x;           // y[n]
        double ynp1 = -0.5*yn*m_nu/(m_m*x); // y[n+1]
        double ynp2;                        // y[n+2]
        double omega = yn + ynp1;

        do
        {
            ynp2 = -(m_nu*ynp1/x + (n+1)*yn/(x*x)) / (2.0*m_m+(n+1));
            omega += ynp2;
            yn = ynp1;
            ynp1 = ynp2;
            ++n;
        } while (fabs(ynp2) > 1.0e-14);

        pdf = Pdf(x);

        return omega*pdf / (2.0*m_m-1.0);
    }
    else if (x >= tail)
    {
        int n = 0;
        double pdf;
        double yn = (1.0+x*x)/x;            // y[n]
        double ynp1 = -0.5*yn*m_nu/(m_m*x); // y[n+1]
        double ynp2;                        // y[n+2]
        double omega = yn + ynp1;

        do
        {
            ynp2 = -(m_nu*ynp1/x + (n+1)*yn/(x*x)) / (2.0*m_m+(n+1));
            omega += ynp2;
            yn = ynp1;
            ynp1 = ynp2;
            ++n;
        } while (fabs(ynp2) > 1.0e-14);

        pdf = Pdf(x);

        return 1.0 - omega*pdf / (2.0*m_m-1.0);
    }

    // use Gauss-Legendre quadrature to integrate
    // when -tail < x < tail
    double area;
    int n = 10;
    hwMatrix X(n, hwMatrix::REAL);   // locations
    hwMatrix W(n, hwMatrix::REAL);   // weights

    X(0)= -0.973906528517172;
    X(1)= -0.865063366688985;
    X(2)= -0.679409568299024;
    X(3)= -0.433395394129247;
    X(4)= -0.148874338981631;
    X(5)=  0.148874338981631;
    X(6)=  0.433395394129247;
    X(7)=  0.679409568299024;
    X(8)=  0.865063366688985;
    X(9)=  0.973906528517172;
    W(0) = 0.066671344308688;
    W(1) = 0.149451349150581;
    W(2) = 0.219086362515982;
    W(3) = 0.269266719309996;
    W(4) = 0.295524224714753;
    W(5) = 0.295524224714753;
    W(6) = 0.269266719309996;
    W(7) = 0.219086362515982;
    W(8) = 0.149451349150581;
    W(9) = 0.066671344308688;

    if (x < 0.0)
    {
        area = Cdf(-tail);
        area += GaussQuad(-tail, x, X, W);
    }
    else
    {
        area = Cdf(tail);
        area -= GaussQuad(x, tail, X, W);
    }

    return area;
}
//------------------------------------------------------------------------------
// Returns inverse cumulative density function
//------------------------------------------------------------------------------
double hwPearson_IV::CdfInv(double prob)
{
    double mode     = Mode();
    double delta    = sqrt((4.0*m_m*m_m+m_nu*m_nu)/(2.0*m_m+1)) / (2.0*m_m);
    double inflect1 = mode - delta;
    double inflect2 = mode + delta;
    double modeCDF  = Cdf(mode);
    double x        = (prob < modeCDF) ? inflect1 :  // use Newton's method
                                         inflect2;
    double f;
    double fp;

    for (int i = 0; i < 30; i++)
    {
        f = Cdf(x) - prob;
        if (fabs(f) < 1.0e-10)
        {
            break;          // function convergence
        }
        fp = Pdf(x);

        if (fabs(f) < 1.0e-10 * fabs(fp))
        {
            break;          // step convergence
        }
        x -= f / fp;
    }

    return x;
}
//------------------------------------------------------------------------------
// Compute the mode of the distribution
//------------------------------------------------------------------------------
double hwPearson_IV::Mode()
{
    return -0.5 * m_nu / m_m;
}
//------------------------------------------------------------------------------
// Compute a random deviate from the distribution
//------------------------------------------------------------------------------
double hwPearson_IV::GetDeviate()
{
    // code from Heinrich, theory is from Luc Devroye,
    // "Non-Uniform Random Variate Generation", Chapter 7, p. 308
    // if X is Pearson IV, then arctan(x) is log-concave with density
    // g(x) = c * [cos(x)]^[2*(m-1)] * exp(-nu*x)
    const double b = 2 * m_m-2;
    const double M = atan(-m_nu / b);
    const double cosM = b / sqrt(b*b + m_nu*m_nu);
    const double r = b * log(cosM) - m_nu * M;
    const double rc = exp(-r) / m_k;
    double x;
    double z;

    do
    {
        int s = 0;
        z = 0;

        if( (x=4.0*m_pUniform->GetDeviate()) > 2.0 )
        {
            x -= 2;
            s = 1;
        }

        if (x > 1)
        {
            x = 1 - (z=log(x-1));
        }
        x = (s) ? M + rc*x : M - rc*x;
    }
    while (fabs(x) >= PI / 2.0 ||
           z + log(m_pUniform->GetDeviate()) > b*log(cos(x)) - m_nu*x - r);

    return tan(x);
}
//------------------------------------------------------------------------------
// Compute abs(gamma(x+iy)/gamma(x))^2
//------------------------------------------------------------------------------
double hwPearson_IV::gammar2(double x, double y)
{
    const double y2 = y*y, xmin = (2*y2 > 10.0) ? 2*y2 : 10.0;
    double r = 1;
    double s = 1;
    double p = 1;
    double f = 0;

    while (x < xmin)
    {
        const double t = y/x++;
        r *= 1 + t*t;
    }

    while (p > s*MACHEP2)
    {
        p *= y2 + f*f;
        p /= x++ * ++f;
        s += p;
    }

    return 1.0/(r*s);
}
//------------------------------------------------------------------------------
// Compute pdf normalization constant
//------------------------------------------------------------------------------
double hwPearson_IV::type4norm()
{
    return 1.0/sqrt(PI)*gammar2(m_m,0.5*m_nu)*exp(GammaLog(m_m)-GammaLog(m_m-0.5));
}
