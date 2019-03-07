/**
* @file hwNormal.cxx
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
#include "hwNormal.h"

#ifndef OS_WIN
#   include <math.h>
#endif

#include "GeneralFuncs.h"
#include "hwUniform.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwNormal::hwNormal(hwMersenneTwisterState* pMTState)
    : m_pUniform      (nullptr)
    , m_bExtraDeviate (false)
{
    if (pMTState)
    {
        m_pUniform = new hwUniform(pMTState);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwNormal::~hwNormal()
{
    delete m_pUniform;
}
//------------------------------------------------------------------------------
// Returns probability density function
//------------------------------------------------------------------------------
double hwNormal::Pdf(double x)
{
    double pi = acos(-1.0);

    return 1.0 / sqrt(2.0 * pi) * exp(-0.5 * x * x);
}
//------------------------------------------------------------------------------
// Standard Normal Distribution CDF
// Returns the area under the Gaussian probability density function, integrated 
// from minus infinity to x:
//                            x
//                             -
//                   1        | |          2
//    ndtr(x)  = ---------    |    exp( - t /2 ) dt
//               sqrt(2pi)  | |
//                           -
//                          -inf.
//
//             =  ( 1 + erf(z) ) / 2
//             =  erfc(z) / 2
// where z = x/sqrt(2). Computation is via the functions erf and erfc with care 
// to avoid error amplification in computing exp(-x^2).
// ACCURACY:
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE     -13,0        30000       1.3e-15     2.2e-16
//
// ERROR MESSAGES:
//
//   message         condition         value returned
// erfc underflow    x > 37.519379347       0.0
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1995, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------
double hwNormal::Cdf(double a)
{
    double y;
    double SQRTH = 0.707106781186547524401;

    double x = a * SQRTH;
    double z = fabs(x);

    if (z < 1.0)
    {
        y = 0.5 + 0.5 * Erf(x);
    }
    else
    {
        y = 0.5 * Erfc(z);
        if (x > 0)
        {
            y = 1.0 - y;
        }
    }

    return y;
}
//------------------------------------------------------------------------------
// Error Function
// The integral is
//
//                           x 
//                            -
//                 2         | |          2
//   erf(x)  =  --------     |    exp( - t  ) dt.
//              sqrt(pi)   | |
//                          -
//                           0
//
// The magnitude of x is limited to 9.231948545 for DEC arithmetic; 1 or -1 is 
// returned outside this range.
// For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise
// erf(x) = 1 - erfc(x).
// ACCURACY:
//
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    DEC       0,1         14000       4.7e-17     1.5e-17
//    IEEE      0,1         30000       3.7e-16     1.0e-16
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1995, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------
double hwNormal::Erf(double x)
{
    double T[] = {
        9.60497373987051638749E0,
        9.00260197203842689217E1,
        2.23200534594684319226E3,
        7.00332514112805075473E3,
        5.55923013010394962768E4
    };
    double U[] = {
        // 1.00000000000000000000E0,
        3.35617141647503099647E1,
        5.21357949780152679795E2,
        4.59432382970980127987E3,
        2.26290000613890934246E4,
        4.92673942608635921086E4
    };

    if (fabs(x) > 1.0)
    {
        return 1.0 - Erfc(x);
    }
    double z = x * x;
    double y = x * polynomFunc1(T, 4, z) / polynomFunc2(U, 5, z);

    return y;
}
//------------------------------------------------------------------------------
// Complementary Error Function
//  1 - erf(x) =
//
//                           inf. 
//                             -
//                  2         | |          2
//   erfc(x)  =  --------     |    exp( - t  ) dt
//               sqrt(pi)   | |
//                           -
//                            x
//
// For small x, erfc(x) = 1 - erf(x); otherwise rational approximations are computed.
//
// A special function expx2.c is used to suppress error amplification
// in computing exp(-x^2).
// ACCURACY:
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0,26.6417   30000       1.3e-15     2.2e-16
//
// ERROR MESSAGES:
//   message         condition              value returned
// erfc underflow    x > 9.231948545 (DEC)       0.0
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1995, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------
double hwNormal::Erfc(double a)
{
    double x = (a < 0.0) ? -a : a;
    if (x < 1.0)
    {
        return 1.0 - Erf(a);
    }

    double P[] = {
        2.46196981473530512524E-10,
        5.64189564831068821977E-1,
        7.46321056442269912687E0,
        4.86371970985681366614E1,
        1.96520832956077098242E2,
        5.26445194995477358631E2,
        9.34528527171957607540E2,
        1.02755188689515710272E3,
        5.57535335369399327526E2
    };
    double Q[] = {
        // 1.00000000000000000000E0,
        1.32281951154744992508E1,
        8.67072140885989742329E1,
        3.54937778887819891062E2,
        9.75708501743205489753E2,
        1.82390916687909736289E3,
        2.24633760818710981792E3,
        1.65666309194161350182E3,
        5.57535340817727675546E2
    };
    double R[] = {
        5.64189583547755073984E-1,
        1.27536670759978104416E0,
        5.01905042251180477414E0,
        6.16021097993053585195E0,
        7.40974269950448939160E0,
        2.97886665372100240670E0
    };
    double S[] = {
        // 1.00000000000000000000E0,
        2.26052863220117276590E0,
        9.39603524938001434673E0,
        1.20489539808096656605E1,
        1.70814450747565897222E1,
        9.60896809063285878198E0,
        3.36907645100081516050E0
    };

    double p;
    double q;

    double z = -a * a;
    z = exp(z);

    if (x < 8.0)
    {
        p = polynomFunc1(P, 8, x);
        q = polynomFunc2(Q, 8, x);
    }
    else
    {
        p = polynomFunc1(R, 5, x);
        q = polynomFunc2(S, 6, x);
    }

    double y = 0.0;
    if (z)      // avoid limit issues with p,q
    {
        y = (z * p)/q;
    }

    if (a < 0)
    {
        y = 2.0 - y;
    }
    return y;
}

//////////////////////////////////////////////////////////////////////
// Standard Normal Distribution Inverse CDF

/*
*
* DESCRIPTION:
*
* Returns the argument, x, for which the area under the
* Gaussian probability density function (integrated from
* minus infinity to x) is equal to y.
*
* For small arguments 0 < y < exp(-2), the program computes
* z = sqrt( -2.0 * log(y) );  then the approximation is
* x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
* There are two rational functions P/Q, one for 0 < y < exp(-32)
* and the other for y up to exp(-2).  For larger arguments,
* w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).
*
* ACCURACY:
*
*                      Relative error:
* arithmetic   domain        # trials      peak         rms
*    DEC      0.125, 1         5500       9.5e-17     2.1e-17
*    DEC      6e-39, 0.135     3500       5.7e-17     1.3e-17
*    IEEE     0.125, 1        20000       7.2e-16     1.3e-16
*    IEEE     3e-308, 0.135   50000       4.6e-16     9.8e-17
*
*
* ERROR MESSAGES:
*
*   message         condition    value returned
* ndtri domain       x <= 0        -MAXNUM
* ndtri domain       x >= 1         MAXNUM
*

Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1995, 2000 by Stephen L. Moshier

The Cephes Library is public domain, available from http://www.netlib.org/cephes.
*/
//------------------------------------------------------------------------------
// Inverse cummulative density function
//------------------------------------------------------------------------------
double hwNormal::CdfInv(double prob)
{
    int code;
    double x, y, z, y2, x0, x1;

    // sqrt(2pi)
    static double s2pi = 2.50662827463100050242E0;

    // approximation for 0 <= |y - 0.5| <= 3/8
    static double P0[5] = {
        -5.99633501014107895267E1,
        9.80010754185999661536E1,
        -5.66762857469070293439E1,
        1.39312609387279679503E1,
        -1.23916583867381258016E0,
    };
    static double Q0[8] = {
        // 1.00000000000000000000E0,
        1.95448858338141759834E0,
        4.67627912898881538453E0,
        8.63602421390890590575E1,
        -2.25462687854119370527E2,
        2.00260212380060660359E2,
        -8.20372256168333339912E1,
        1.59056225126211695515E1,
        -1.18331621121330003142E0,
    };

    // Approximation for interval z = sqrt(-2 log y ) between 2 and 8
    // i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
    static double P1[9] = {
        4.05544892305962419923E0,
        3.15251094599893866154E1,
        5.71628192246421288162E1,
        4.40805073893200834700E1,
        1.46849561928858024014E1,
        2.18663306850790267539E0,
        -1.40256079171354495875E-1,
        -3.50424626827848203418E-2,
        -8.57456785154685413611E-4,
    };
    static double Q1[8] = {
        // 1.00000000000000000000E0,
        1.57799883256466749731E1,
        4.53907635128879210584E1,
        4.13172038254672030440E1,
        1.50425385692907503408E1,
        2.50464946208309415979E0,
        -1.42182922854787788574E-1,
        -3.80806407691578277194E-2,
        -9.33259480895457427372E-4,
    };

    // Approximation for interval z = sqrt(-2 log y ) between 8 and 64
    // i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
    static double P2[9] = {
        3.23774891776946035970E0,
        6.91522889068984211695E0,
        3.93881025292474443415E0,
        1.33303460815807542389E0,
        2.01485389549179081538E-1,
        1.23716634817820021358E-2,
        3.01581553508235416007E-4,
        2.65806974686737550832E-6,
        6.23974539184983293730E-9,
    };
    static double Q2[8] = {
        // 1.00000000000000000000E0,
        6.02427039364742014255E0,
        3.67983563856160859403E0,
        1.37702099489081330271E0,
        2.16236993594496635890E-1,
        1.34204006088543189037E-2,
        3.28014464682127739104E-4,
        2.89247864745380683936E-6,
        6.79019408009981274425E-9,
    };

    if (prob <= 0.0)
        return -99999.0;

    if (prob >= 1.0)
        return 99999.0;

    code = 1;
    y = prob;

    if (y > 1.0 - 0.13533528323661269189) // 0.135... = exp(-2)
    {
        y = 1.0 - y;
        code = 0;
    }

    if (y > 0.13533528323661269189)
    {
        y -= 0.5;
        y2 = y * y;
        x = y + y * (y2 * polynomFunc1(P0, 4, y2) / polynomFunc2(Q0, 8, y2));
        x *= s2pi;

        return x;
    }

    x = sqrt(-2.0 * log(y));
    x0 = x - log(x) / x;
    z = 1.0 / x;

    if (x < 8.0) // y > exp(-32) = 1.2664165549e-14
        x1 = z * polynomFunc1(P1, 8, z) / polynomFunc2(Q1, 8, z);
    else
        x1 = z * polynomFunc1(P2, 8, z) / polynomFunc2(Q2, 8, z);

    x = x0 - x1;

    if (code != 0)
        x = -x;

    return x;
}
//------------------------------------------------------------------------------
// Compute a random deviate from the distribution
//------------------------------------------------------------------------------
double hwNormal::GetDeviate()
{
    if (m_bExtraDeviate)
    {
        m_bExtraDeviate = false;
        return m_extraDeviate;
    }

    double radiusSq;
    double mag;
    double u1, u2;

    do
    {
        u1 = 2.0 * m_pUniform->GetDeviate() - 1.0;
        u2 = 2.0 * m_pUniform->GetDeviate() - 1.0;
        radiusSq = u1 * u1 + u2 * u2;
    }
    while (radiusSq >= 1.0 || radiusSq == 0.0);

    mag = sqrt(-2.0 * log(radiusSq) / radiusSq);
    m_extraDeviate = mag * u2;		// deviates generated in pairs, save the extra one
    m_bExtraDeviate = true;

    return mag * u1;
}
