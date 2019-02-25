/**
* @file GeneralFuncs.cxx
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

//:---------------------------------------------------------------------------
//:Description
//
//  General functions
//
//:---------------------------------------------------------------------------

#include <math.h>
#include <limits>
#include <GeneralFuncs.h>

//////////////////////////////////////////////////////////////////////
// Polynomial (Cephes)

/*
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 *                     2          N
 * y  =  C  + C x + C x  +...+ C x
 *        0    1     2          N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C  , ..., coef[N] = C  .
 *            N                   0
 *
 * SPEED:
 *
 * In the interest of speed, there are no checks for out
 * of bounds arithmetic.  This routine is used by most of
 * the functions in the library.  Depending on available
 * equipment features, the user may wish to rewrite the
 * program in microcode or assembly language.
 *

 Cephes Math Library, Release 2.8: June, 2000
 Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier

 The Cephes Library is public domain, available from http://www.netlib.org/cephes.
 This function was extracted and reformatted by Brian Shier.
*/

// -------------------------------------------------------------------
//*******************************************************************
//                     Polynomial functions
//*******************************************************************
//! Evaluate a polynomial of degree n at point x. 
//! The coefficients are stored in descending order.
double polynomFunc1(double* coef, int n, double x)
{
	double ans;
	int i;
	
	ans = *coef++;
	i = n;
	
	do
	    ans = ans * x  +  *coef++;
	while (--i);
	
	return ans;
}

//! Evaluate a polynomial of degree n at point x. 
//! The coefficients are stored in descending order.
//! The first coefficient is assumed to equal 1.0 and is omitted.
double polynomFunc2(double* coef, int n, double x)
{
	double ans;
	int i;
	
	ans = x + *coef++;
	i = n - 1;
	
	do
	    ans = ans * x  + *coef++;
	while (--i);
	
	return ans;
}

//! Multiply polynomial 'a' times polynomial 'b' 
void polyNomMultAccum(double* a, int n, double* b, int m)
{
    // The coefficients are stored in descending order.
    // n is the order of a on input, not the size of a
    // m is the order of b on input, not the size of b
    // n+m is the order of a on output
    int i, j;

    for (i = n + m; i > n; --i)
    {
        a[i] = a[n] * b[i-n];

        for (j = i - n + 1; j <= m; ++j)
            a[i] += a[i-j] * b[j];
    }

    for (i = n; i >= m; --i)
    {
        a[i] *= b[0];

        for (j = 1; j <= m; ++j)
            a[i] += a[i-j] * b[j];
    }

    for (i = m - 1; i >= 0; --i)
    {
        a[i] *= b[0];

        for (j = 1; j < i + 1; ++j)
            a[i] += a[i-j] * b[j];
    }
}

//! Find the real roots of ax^2 + bx + c 
bool quadraticRoots(double a, double b, double c,
					double& zero1, double& zero2)
{
	double discriminant = b * b - 4.0 * a * c;
	
	if (discriminant < 0.0)
		return false;
	
	double value = sqrt(discriminant);
	double q;
	
	if (b < 0.0)
		q = -0.5 * (b - value);
	else
		q = -0.5 * (b + value);
	
	zero1 = q / a;
	zero2 = c / q;
	
	return true;
}

//*******************************************************************
//                    General utility functions
//*******************************************************************
//! Check to see if a value is an integer to within a tolerance
hwMathStatus IsInteger(double value, double tol)
{
    if (tol < 0.0)
        return hwMathStatus(HW_MATH_ERR_NEGATIVE, 2);

    if (value >= 0.0)
    {
        if (value > std::numeric_limits<int>::max())
            return hwMathStatus(HW_MATH_ERR_BADRANGE, 1);

        if (fabs(value - floor(value + 0.5)) > tol)
            return hwMathStatus(HW_MATH_ERR_NONINTEGER, 1);
    }
    else
    {
        if (value < std::numeric_limits<int>::min())
            return hwMathStatus(HW_MATH_ERR_BADRANGE, 1);

        if (fabs(value - ceil(value - 0.5)) > tol)
            return hwMathStatus(HW_MATH_ERR_NONINTEGER, 1);
    }
    
    return hwMathStatus();
}

