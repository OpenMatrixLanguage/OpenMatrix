/**
* @file SpecialFuncs.cxx
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
#include <SpecialFuncs.h>

#include <math.h>
#include <GeneralFuncs.h>
#include <Globals.h>

//------------------------------------------------------------------------------
// Returns gamma function of the argument. The result is correctly signed, and 
// the sign (+1 or -1) is also returned in a global (extern) variable named 
// sgngam. This variable is also filled in by the logarithmic gamma function 
// lgam().
// Arguments |x| <= 34 are reduced by recurrence and the function approximated 
// by a rational function of degree 6/7 in the interval (2,3). Large arguments 
// are handled by Stirling's formula. Large negative arguments are made positive
// using a reflection formula.  
//
// ACCURACY:
//                      Relative error:
// arithmetic   domain     # trials      peak         rms
//    DEC      -34, 34      10000       1.3e-16     2.5e-17
//    IEEE    -170,-33      20000       2.3e-15     3.3e-16
//    IEEE     -33,  33     20000       9.4e-16     2.2e-16
//    IEEE      33, 171.6   20000       2.3e-15     3.2e-16
//
// Error for arguments outside the test range will be larger owing to error 
// amplification by the exponential function.
// Cephes Math Library, Release 2.8: June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from 
// http://www.netlib.org/cephes.
//------------------------------------------------------------------------------
double GammaFunc(double x)
{
    if (IsNaN_T(x) || IsNegInf_T(x))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (IsInf_T(x))
    {
        return std::numeric_limits<double>::infinity();
    }

    double P[] = {
        1.60119522476751861407E-4,
        1.19135147006586384913E-3,
        1.04213797561761569935E-2,
        4.76367800457137231464E-2,
        2.07448227648435975150E-1,
        4.94214826801497100753E-1,
        9.99999999999999996796E-1
    };
    double Q[] = {
       -2.31581873324120129819E-5,
        5.39605580493303397842E-4,
       -4.45641913851797240494E-3,
        1.18139785222060435552E-2,
        3.58236398605498653373E-2,
       -2.34591795718243348568E-1,
        7.14304917030273074085E-2,
        1.00000000000000000320E0
    };
		
	double q = fabs(x);
	double z = 1.0;
    double p;

	while (x >= 3.0)
	{
		x -= 1.0;
		z *= x;
	}
	
	while (x < 0.0)
	{
		if (x > -1.0e-9)
        {
        	return z / ((1.0 + 0.5772156649015329 * x) * x);
        }
		z /= x;
		x += 1.0;
	}
	
	while (x < 2.0)
	{
		if (x < 1.0e-9)
        {
        	return z / ((1.0 + 0.5772156649015329 * x) * x);
        }
		z /= x;
		x += 1.0;
	}
	
	if (x == 2.0)
    {
		return z;
    }

	x -= 2.0;
	p  = polynomFunc1(P, 6, x);
	q  = polynomFunc1(Q, 7,x);

	return z * p / q;
}
//------------------------------------------------------------------------------
// Returns the base e (2.718...) logarithm of the absolute value of the gamma 
// function of the argument. The sign (+1 or -1) of the gamma function is 
// returned in a global (extern) variable named sgngam.
// For arguments greater than 13, the logarithm of the gamma function is 
// approximated by the logarithmic version of Stirling's formula using a 
// polynomial approximation of degree 4. Arguments between -33 and +33 are 
// reduced by recurrence to the interval [2,3] of a rational approximation. The 
// cosecant reflection formula is employed for arguments less than -33.
// Arguments greater than MAXLGM return MAXNUM and an error message. 
// MAXLGM = 2.035093e36 for DEC arithmetic or 2.556348e305 for IEEE arithmetic.
// ACCURACY:
// arithmetic      domain        # trials     peak         rms
//    DEC     0, 3                  7000     5.2e-17     1.3e-17
//    DEC     2.718, 2.035e36       5000     3.9e-17     9.9e-18
//    IEEE    0, 3                 28000     5.4e-16     1.1e-16
//    IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
// The error criterion was relative when the function magnitude was greater than
// one but absolute when it was less than one.
// The following test used the relative error criterion, though
// at certain points the relative error could be much higher than
// indicated.
//    IEEE    -200, -4             10000     4.8e-16     1.3e-16

// Cephes Math Library, Release 2.8: June, 2000
// Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// The Cephes Library is public domain, available from http://www.netlib.org/cephes.
//------------------------------------------------------------------------------
double GammaLog(double x)
{
    if (IsNaN_T(x) || IsNegInf_T(x))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (IsInf_T(x))
    {
        return std::numeric_limits<double>::infinity();
    }

    double pi    = acos(-1.0);
	double LS2PI = log(sqrt(2.0*pi));
	
    double A[] = 
    {
        8.11614167470508450300E-4,
       -5.95061904284301438324E-4,
        7.93650340457716943945E-4,
       -2.77777777730099687205E-3,
        8.33333333333331927722E-2
    };

    double B[] =
    {
        -1.37825152569120859100E3,
        -3.88016315134637840924E4,
        -3.31612992738871184744E5,
        -1.16237097492762307383E6,
        -1.72173700820839662146E6,
        -8.53555664245765465627E5
    };

    double C[] = 
    {
        // 1.00000000000000000000E0,
        -3.51815701436523470549E2,
        -1.70642106651881159223E4,
        -2.20528590553854454839E5,
        -1.13933444367982507207E6,
        -2.53252307177582951285E6,
        -2.01889141433532773231E6
    };

    double p;
    double q;
    double u;
    double z;

	if (x < 13.0)
	{
		z = 1.0;
		p = 0.0;
		u = x;

		while (u >= 3.0)
		{
			p -= 1.0;
			u  = x + p;
			z *= u;
		}

		while (u < 2.0)
		{
			if (u == 0.0)
            {
				// return 0.0;
                return std::numeric_limits<double>::infinity();
            }

			z /= u;
			p += 1.0;
			u  = x + p;
		}

		if (z < 0.0)
        {
			z = -z;
        }

		if (u == 2.0)
        {
			return log(z);
        }

		p -= 2.0;
		x += p;
		p  = x * polynomFunc1(B, 5, x) / polynomFunc2(C, 6, x);

		return log(z) + p;
	}
	
	q = (x - 0.5) * log(x) - x + LS2PI;

	if (x > 1.0e8)
    {
		return q;
    }

	p = 1.0 / (x * x);

	if (x >= 1000.0)
    {
		q += ((7.9365079365079365079365e-4    * p
		       - 2.7777777777777777777778e-3) * p
		       + 0.0833333333333333333333) / x;
    }
	else
    {
		q += polynomFunc1(A, 4, p) / x;
    }
	return q;
}
//------------------------------------------------------------------------------
// Psi (digamma) function
// SYNOPSIS:
// double x, y, psi();
// y = psi( x );
// DESCRIPTION:
//
//              d      -
//   psi(x)  =  -- ln | (x)
//              dx
// is the logarithmic derivative of the gamma function.
// For integer x,
//                   n-1
//                    -
// psi(n) = -EUL  +   >  1/k.
//                    -
//                   k=1
//
// This formula is used for 0 < n <= 10.  If x is negative, it is transformed to 
// a positive argument by the reflection formula  psi(1-x) = psi(x) + pi cot(pi x).
// For general positive x, the argument is made greater than 10 using the 
// recurrence  psi(x+1) = psi(x) + 1/x. Then the following asymptotic expansion 
// is applied:
//                           inf.   B
//                            -      2k
// psi(x) = log(x) - 1/2x -   >   -------
//                            -        2k
//                           k=1   2k x
// where the B2k are Bernoulli numbers.
// ACCURACY:
//    Relative error (except absolute when |psi| < 1):
// arithmetic   domain     # trials      peak         rms
//    DEC       0,30         2500       1.7e-16     2.0e-17
//    IEEE      0,30        30000       1.3e-15     1.4e-16
//    IEEE      -30,0       40000       1.5e-15     2.2e-16
// ERROR MESSAGES:
//     message         condition      value returned
// psi singularity    x integer <=0      MAXNUM
// Cephes Math Library Release 2.8:  June, 2000
// Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
//------------------------------------------------------------------------------
double Digamma(double x)
{
    // first derivative of the gamma function
    static double A[] =
    {
         8.33333333333333333333E-2,
        -2.10927960927960927961E-2,
         7.57575757575757575758E-3,
        -4.16666666666666666667E-3,
         3.96825396825396825397E-3,
        -8.33333333333333333333E-3,
         8.33333333333333333333E-2
    };

    double p;
    double q;
    double s;
    double euler = 0.57721566490153286061;
    double nz    = 0.0;

    int negative = 0;

    if ( x <= 0.0 )
    {
        negative = 1;
        q        = x;
        p        = floor(q);

        if (p == q)
        {
            return MAXNUM;
        }

        // Remove the zeros of tan(PI x)
        // by subtracting the nearest integer from x
        nz = q - p;

        if (nz != 0.5)
        {
            if (nz > 0.5)
            {
                p += 1.0;
                nz = q - p;
            }

            nz = PI / tan(PI*nz);
        }
        else
        {
            nz = 0.0;
        }

        x = 1.0 - x;
    }

    int    n;
    double w;
    double y;

    // check for positive integer up to 10
    if (x <= 10.0 && x == floor(x))
    {
        y = 0.0;
        n = static_cast<int>(x);

        for (int i = 1; i < n; i++)
        {
            w = i;
            y += 1.0 / w;
        }

        y -= euler;
    }
    else
    {
        s = x;
        w = 0.0;

        while (s < 10.0)
        {
            w += 1.0/s;
            s += 1.0;
        }

        double z;
        if (s < 1.0e+17)
        {
            z = 1.0 / (s * s);
            y = z * polynomFunc1(A, 6, z);
        }
        else
        {
            y = 0.0;
        }

        y = log(s) - (0.5/s) - y - w;
    }

    if (negative)
    {
        y -= nz;
    }

    return y;
}
//------------------------------------------------------------------------------
// Trigamma function
//------------------------------------------------------------------------------
double Trigamma(double x)
{
    if (x <= 0.0)
    {
        return 0.0;
    }

    // second derivative of the gamma function
    int    count = 0;
    double sum   = 0.0;
    double value;

    while (true)
    {
        value = 1.0 / (x + (double) count);

        if (value < 1.0e-7)
        {
            break;
        }

        sum += value * value;
        count++;
    }

    return sum;
}
//------------------------------------------------------------------------------
// Beta function
//------------------------------------------------------------------------------
double BetaFunc(double x,
                double y)
{
    double beta;

	if (x + y < MAXGAM)
	{
		beta = GammaFunc(x) * GammaFunc(y) / GammaFunc(x + y);
	}
	else
	{
		beta = GammaLog(x) + GammaLog(y) - GammaLog(x + y);
					
		if (beta < MINLOG)
        {
			beta = 0.0;
        }
		else
        {
			beta = exp(beta);
        }
	}
				
	return beta;
}
//------------------------------------------------------------------------------
// Logarithm of the Beta function
//------------------------------------------------------------------------------
double BetaLog(double x,
               double y)
{
    double betaLog;
				
	if (x + y < MAXGAM)
	{
		betaLog = GammaFunc(x) * GammaFunc(y) / GammaFunc(x + y);

		if (betaLog > 0.0)
        {
            betaLog = log(betaLog);
        }
		else
        {
			betaLog = 0.0;
        }
	}
	else
	{
		betaLog = GammaLog(x) + GammaLog(y) - GammaLog(x + y);
	}
				
	return betaLog;
}
//------------------------------------------------------------------------------
// Zero order modified Bessel function of the first kind (Io)
//------------------------------------------------------------------------------
double Bessel(double alpha)
{
    // alpha is the domain variable, not the function order
    if (alpha < 0.0)
    {
        return 0.0;
    }

    double fact = 1.0;
    double sum  = 1.0;
    double value;

    for (int j = 1; j < 33; j++)
    {
        fact *= j;
        value = pow(0.5 * alpha, j) / fact;
        sum  += value * value;
    }

    return sum;
}
