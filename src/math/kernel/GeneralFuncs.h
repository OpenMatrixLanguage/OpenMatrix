/**
* @file GeneralFuncs.h
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
#ifndef _GeneralFuncs_h
#define _GeneralFuncs_h

#include <limits>
#include <MathCoreExports.h>
#include <hwMathStatus.h>

#ifndef OS_WIN
  #include <stdlib.h>
  #include <cmath>
#endif

//*******************************************************************
//                     Polynomial functions
//*******************************************************************
MATHCORE_DECLS double polynomFunc1(double* coef, int n, double x);
MATHCORE_DECLS double polynomFunc2(double* coef, int n, double x);
MATHCORE_DECLS void polyNomMultAccum(double* a, int n, double* b, int m);
MATHCORE_DECLS bool quadraticRoots(double a, double b, double c,
								   double& zero1, double& zero2);


//*******************************************************************
//                    General utility functions
//*******************************************************************
template<typename T>
bool IsNaN_T(T x)
{
    return (x != x);
}

template<typename T>
bool IsInf_T(T x)
{
    return std::numeric_limits<T>::has_infinity &&
           (x == std::numeric_limits<T>::infinity());
}

template<typename T>
bool IsNegInf_T(T x)
{
    return std::numeric_limits<T>::has_infinity &&
           ((-x) == std::numeric_limits<T>::infinity());
}

template<typename T>
bool IsFinite_T(T x)
{
    return !(IsNaN_T(x) || IsInf_T(x) || IsNegInf_T(x));
}

MATHCORE_DECLS hwMathStatus IsInteger(double value, double tol = 0.0);

//! Return the nearest integer value of type T
template < typename T >
static T RoundT(T value)
{
    if (value < (T) 0)
        return ceil(value - (T) 0.5);
    else
        return floor(value + (T) 0.5);
}

//! Fix function (round toward zero)
template < typename T >
static T FixT(T value)
{
    if (value < (T) 0)
        return ceil(value);
    else
        return floor(value);
}

static inline int absT(int value)
{
    return abs(value);
}

static inline float absT(float value)
{
    return fabsf(value);
}

static inline double absT(double value)
{
    return fabs(value);
}

template<typename T>
bool IsZero(T value, double tol = 0.0)
{
    T mag = absT(value);

    if (IsNaN_T(mag) || IsNaN_T(tol))
        return false;

    if (mag > tol)
        return false;
    else
        return true;
}

template<typename T>
bool AreEqual(T value1, T value2, double tol = 0.0)
{
    if (IsNaN_T(value1) || IsNaN_T(value2) || IsNaN_T(tol))
        return false;

    if (tol)
    {
        if (value1 > value2 + tol || value2 > value1 + tol)
            return false;
        else
            return true;
    }
    else
    {
        return (value1 == value2);
    }
}

template< typename T >
T _min(T a, T b)
{
	if (a < b)
		return a;
	else
		return b;
}

template< typename T >
T _max(T a, T b)
{
	if (a > b)
		return a;
	else
		return b;
}

//! Return a, but with the sign of b
template< typename T >
T sign(T a, T b)
{
	if (a < 0)
    {
        if (b < 0)
            return a;
        else
            return -a;
    }
	else
    {
        if (b < 0)
            return -a;
        else
            return a;
    }
}

//*******************************************************************
//                    Machine related functions
//*******************************************************************
//! Return smallest non-zero number relative to reference
template< typename T >
T MachPrecision(T ref) 
{
    if (!IsFinite_T(ref))
        return std::numeric_limits<T>::quiet_NaN();

    if (ref < (T) 0)
        ref = -ref;

    T eps;
    T temp = (T) 1;

    while (ref + temp != ref)
    {
        eps = temp;
        temp *= (T) 0.5;
    }

    return eps;
}

//*******************************************************************
//            Custom functions (overriding stdc++ lib)
//*******************************************************************
//! Return the power, allowing for a NaN input
template < typename T >
T CustomPow(T base, T exp)
{
	if (base != base)
		return base;     // base is a NaN

	if (exp != exp)
		return exp;     // exponent is a NaN

	return pow(base, exp);
}

#endif // _GeneralFuncs_h
