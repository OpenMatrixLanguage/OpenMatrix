/**
* @file EllipticFuncs.cxx
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
#include <EllipticFuncs.h>

#include <GeneralFuncs.h>
#include <math.h>

//------------------------------------------------------------------------------
// Complete elliptic integral of the first kind
// DESCRIPTION:
// Approximates the integral
//            pi/2
//             -
//            | |
//            |           dt
// K(m)  =    |    ------------------
//            |                   2
//          | |    sqrt( 1 - m sin t )
//           -
//            0
// where m = 1 - m1, using the approximation
//     P(x)  -  log x Q(x).
// The argument m1 is used rather than m so that the logarithmic singularity at
// m = 1 will be shifted to the origin; this preserves maximum accuracy.
// K(0) = pi/2.
// ACCURACY:
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    DEC        0,1        16000       3.5e-17     1.1e-17
//    IEEE       0,1        30000       2.5e-16     6.8e-17
// ERROR MESSAGES:
//   message         condition      value returned
// ellpk domain       x<0, x>1           0.0
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1987, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------
double ellpk(double x)
{
    double P[] =
    {
        1.37982864606273237150E-4,
        2.28025724005875567385E-3,
        7.97404013220415179367E-3,
        9.85821379021226008714E-3,
        6.87489687449949877925E-3,
        6.18901033637687613229E-3,
        8.79078273952743772254E-3,
        1.49380448916805252718E-2,
        3.08851465246711995998E-2,
        9.65735902811690126535E-2,
        1.38629436111989062502E0
    };
    double Q[] =
    {
        2.94078955048598507511E-5,
        9.14184723865917226571E-4,
        5.94058303753167793257E-3,
        1.54850516649762399335E-2,
        2.39089602715924892727E-2,
        3.01204715227604046988E-2,
        3.73774314173823228969E-2,
        4.88280347570998239232E-2,
        7.03124996963957469739E-2,
        1.24999999999870820058E-1,
        4.99999999999999999821E-1
    };

    double C1 = 1.3862943611198906188;            // log(4)

    if (x < 0.0 || x > 1.0)
    {
        return 0.0;
    }

    if (x > MACHEP)
    {
        return (polynomFunc1(P, 10, x) - log(x) * polynomFunc1(Q, 10, x));
    }

    if (x == 0.0)
    {
        return MAXNUM;
    }
   
    double val = C1 - 0.5 * log(x);
    return val;
}
//------------------------------------------------------------------------------
// Complete elliptic integral of the second kind
// DESCRIPTION:
// Approximates the integral
//            pi/2
//             -
//            | |                 2
// E(m)  =    |    sqrt( 1 - m sin t ) dt
//          | |    
//           -
//            0
//
// Where m = 1 - m1, using the approximation
//
//      P(x)  -  x log x Q(x).
//
// Though there are no singularities, the argument m1 is used rather than m for 
// compatibility with ellpk().
// E(1) = 1; E(0) = pi/2.
// ACCURACY:
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    DEC        0, 1       13000       3.1e-17     9.4e-18
//    IEEE       0, 1       10000       2.1e-16     7.3e-17
//
// ERROR MESSAGES:
//   message         condition      value returned
// ellpe domain      x<0, x>1            0.0
// Cephes Math Library, Release 2.8: June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------
double ellpe(double x)
{
    double P[] = 
    {
        1.53552577301013293365E-4,
        2.50888492163602060990E-3,
        8.68786816565889628429E-3,
        1.07350949056076193403E-2,
        7.77395492516787092951E-3,
        7.58395289413514708519E-3,
        1.15688436810574127319E-2,
        2.18317996015557253103E-2,
        5.68051945617860553470E-2,
        4.43147180560990850618E-1,
        1.00000000000000000299E0
    };
    double Q[] = 
    {
        3.27954898576485872656E-5,
        1.00962792679356715133E-3,
        6.50609489976927491433E-3,
        1.68862163993311317300E-2,
        2.61769742454493659583E-2,
        3.34833904888224918614E-2,
        4.27180926518931511717E-2,
        5.85936634471101055642E-2,
        9.37499997197644278445E-2,
        2.49999999999888314361E-1
    };

    if (x <= 0.0 || x > 1.0)
    {
        if (x == 0.0)
        {
            return 1.0;
        }
        return 0.0;
    }

    double val = polynomFunc1(P, 10, x) - log(x) * (x * polynomFunc1(Q, 9, x));
    return val;
}
//------------------------------------------------------------------------------
// Returns the incomplete elliptic integral of the first kind
// DESCRIPTION:
// Approximates the integral
//                phi
//                 -
//                | |
//                |           dt
// F(phi_\m)  =    |    ------------------
//                |                   2
//              | |    sqrt( 1 - m sin t )
//               -
//                0
// of amplitude phi and modulus m, using the arithmetic - geometric mean 
// algorithm.
// ACCURACY:
// Tested at random points with m in [0, 1] and phi as indicated.
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE     -10,10       200000      7.4e-16     1.0e-16
//
// Cephes Math Library, Release 2.8: June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------
double ellik(double phi,
             double m)
{
    if (m == 0.0)
    {
        return phi;
    }

    double a = 1.0 - m;
    if (a == 0.0)
    {
        if (fabs(phi) >= 0.5 * PI)
        {
            // ASSERT(false);
            return MAXNUM;
        }

        return log(tan((0.5 * PI + phi) / 2.0));
    }

    int npio2 = static_cast<int>(floor(2.0 * phi / PI));
    if (npio2 & 1)
    {
        npio2 += 1;
    }

    double K = 0.0;
    if (npio2)
    {
        K = ellpk(a);
        phi = phi - 0.5 * npio2 * PI;
    }

    int sign = 0;
    if (phi < 0.0)
    {
        phi = -phi;
        sign = -1;
    }

    double b = sqrt(a);
    double t = tan(phi);
    double e;
    double temp;
    double c;

    int d   = 1;
    int mod = 0;
    if (fabs(t) > 10.0)
    {
        // Transform the amplitude
        e = 1.0/(b*t);

        // ... but avoid multiple recursions
        if (fabs(e) < 10.0)
        {
            e = atan(e);

            if (npio2 == 0)
                K = ellpk(a);

            temp = K - ellik(e, m);
        }
    }
    else
    {
        a = 1.0;
        c = sqrt(m);

        while (fabs(c / a) > MACHEP)
        {
            temp = b/a;
            phi  = phi + atan(t * temp) + mod * PI;
            mod  = static_cast<int>((phi + 0.5 * PI) / PI);
            t    = t * (1.0 + temp) / (1.0 - temp * t * t);
            c    = (a - b) / 2.0;
            temp = sqrt(a * b);
            a    = (a + b) / 2.0;
            b    = temp;
            d   += d;
        }

        temp = (atan(t) + mod * PI) / (d * a);
    }

    if (sign < 0)
    {
        temp = -temp;
    }

    temp += npio2 * K;

    return temp;
}
//------------------------------------------------------------------------------
// Jacobian elliptic functions
// DESCRIPTION: Evaluates the Jacobian elliptic functions sn(u|m), cn(u|m),
// and dn(u|m) of parameter m between 0 and 1, and real argument u.
// These functions are periodic, with quarter-period on the real axis equal to 
// the complete elliptic integral ellpk(1.0-m).
// Relation to incomplete elliptic integral:
// If u = ellik(phi,m), then sn(u|m) = sin(phi), and cn(u|m) = cos(phi).  
// Phi is called the amplitude of u.
// Computation is by means of the arithmetic-geometric mean algorithm, except 
// when m is within 1e-9 of 0 or 1. In the latter case with m close to 1, the
// approximation applies only for phi < pi/2.
// ACCURACY:
// Tested at random points with u between 0 and 10, m between 0 and 1.
//            Absolute error (* = relative error):
// arithmetic   function   # trials      peak         rms
//    DEC       sn           1800       4.5e-16     8.7e-17
//    IEEE      phi         10000       9.2e-16*    1.4e-16*
//    IEEE      sn          50000       4.1e-15     4.6e-16
//    IEEE      cn          40000       3.6e-15     4.4e-16
//    IEEE      dn          10000       1.3e-12     1.8e-14
// Peak error observed in consistency check using addition theorem for sn(u+v) 
// was 4e-16 (absolute). Also tested by the above relation to the incomplete
// elliptic integral. Accuracy deteriorates when u is large.
// Cephes Math Library, Release 2.8: June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------
int ellpj(double  u, 
          double  m, 
          double& sn, 
          double& cn, 
          double& dn, 
          double& ph)
{
    // Check for special cases
    if (m < 0.0 || m > 1.0)
    {
        // ASSERT(false);
        sn = 0.0;
        cn = 0.0;
        ph = 0.0;
        dn = 0.0;

        return -1;
    }

    double ai;
    double b;
    double t;
    if (m < 1.0e-9)
    {
        t = sin(u);
        b = cos(u);
        ai = 0.25 * m * (u - t*b);
        sn = t - ai * b;
        cn = b + ai * t;
        ph = u - ai;
        dn = 1.0 - 0.5 * m * t * t;

        return 0;
    }

    double phi;
    double twon;
    if (m >= 0.9999999999)
    {
        ai = 0.25 * (1.0-m);
        b = cosh(u);
        t = tanh(u);
        phi = 1.0 / b;
        twon = b * sinh(u);
        sn = t + ai * (twon - u) / (b*b);
        ph = 2.0 * atan(exp(u)) - 0.5 * PI + ai*(twon - u) / b;
        ai *= t * phi;
        cn = phi - ai * (twon - u);
        dn = phi + ai * (twon + u);

        return 0;
    }

    double a[9];
    double c[9];

    //	A. G. M. scale
    a[0] = 1.0;
    b    = sqrt(1.0 - m);
    c[0] = sqrt(m);
    twon = 1.0;
    
    int i = 0;
    while (fabs(c[i] / a[i]) > MACHEP)
    {
        if( i > 7 )
        {
            // ASSERT(false);
            break;
        }

        ai = a[i];
        i++;

        c[i]  = (ai - b) / 2.0;
        t     = sqrt(ai * b);
        a[i]  = (ai + b) / 2.0;
        b     = t;
        twon *= 2.0;
    }

    // backward recurrence
    phi = twon * a[i] * u;

    do
    {
        t   = c[i] * sin(phi) / a[i];
        b   = phi;
        phi = (asin(t) + phi) / 2.0;
    }
    while (--i);

    sn = sin(phi);
    t  = cos(phi);
    cn = t;
    dn = t / cos(phi - b);
    ph = phi;

    return 0;
}
//------------------------------------------------------------------------------
// Find parameter corresponding to given nome by expansion in theta functions:
// AMS55 #16.38.5, 16.38.7
//       1/2
// ( 2K )                   4     9
// ( -- )     =  1 + 2q + 2q  + 2q  + ...  =  Theta (0,q)
// ( pi )                                          3
//       1/2
// ( 2K )     1/4       1/4        2    6    12    20
// ( -- )    m     =  2q    ( 1 + q  + q  + q   + q   + ...) = Theta (0,q)
// ( pi )                                                           2
//
// The nome q(m) = exp( - pi K(1-m)/K(m) ).
//
//                                1/2
// Given q, this program returns m   .
//
// Cephes Math Library, Release 2.8: June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------
double cay(double q)
{
    double a = 1.0;
    double b = 1.0;
    double r = 1.0;
    double p = q;

    double t1;
    double t2;

    do
    {
        r *= p;
        a += 2.0 * r;
        t1 = fabs(r / a);

        r *= p;
        b += r;
        p *= q;
        t2 = fabs(r / b);
        if (t2 > t1)
        {
            t1 = t2;
        }
    }
    while (t1 > MACHEP);

    a = b / a;
    a = 4.0 * sqrt(q) * a * a;      // see above formulas, solved for m

    return a;
}
