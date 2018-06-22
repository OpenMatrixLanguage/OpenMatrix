/**
* @file hwTComplex.cc
* @date March 2012
* Copyright (C) 2012-2018 Altair Engineering, Inc.  
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

//:---------------------------------------------------------------------------
//:Description
//
//  Template complex number class implementation file
//
//:---------------------------------------------------------------------------

#include <GeneralFuncs.h>
#include <hwMathStatus.h>

#ifndef OS_WIN
  #define _hypot hypot
  #define _hypotf hypotf
#endif

//*******************************************************************
//            hwTComplex constructor/destructor functions
//*******************************************************************
//! Constructor for complex number class
template < typename T >
hwTComplex<T>::hwTComplex()
    : x((T) 0), y((T) 0)
{
}

//! Constructor
template < typename T >
hwTComplex<T>::hwTComplex(T real, T imag)
    : x(real), y(imag)
{
}

//! Copy constructor
template < typename T >
hwTComplex<T>::hwTComplex(const hwTComplex<T>& z)
    : x(z.x), y(z.y)
{
}

//! Destructor
template < typename T >
hwTComplex<T>::~hwTComplex()
{
}

//*******************************************************************
//                 hwTComplex public functions
//*******************************************************************
//! Assignment operator
template < typename T >
hwTComplex<T>& hwTComplex<T>::operator=(const hwTComplex<T>& z)
{
    x = z.x;
    y = z.y;
    return *this;
}

//! Assignment operator
template < typename T >
hwTComplex<T>& hwTComplex<T>::operator=(T real)
{
    x = real;
    y = (T) 0;
    return *this;
}

//! Compute the magnitude
template <>
inline double hwTComplex<double>::Mag() const
{
    return _hypot(x, y);
}

//! Compute the magnitude
template <>
inline float hwTComplex<float>::Mag() const
{
    return _hypotf(x, y);
}

//! Compute the phase
template < typename T >
inline T hwTComplex<T>::Arg() const
{
    return atan2(y, x);
}

//! Determine if the complex number is real
template < typename T >
inline bool hwTComplex<T>::IsReal(T tol) const
{
    return (abs(y) <= tol ? true : false);
}

//! Determine if the complex number is real
template <>
inline bool hwTComplex<float>::IsReal(float tol) const
{
    return (fabsf(y) <= tol ? true : false);
}

//! Determine if the complex number is real
template <>
inline bool hwTComplex<double>::IsReal(double tol) const
{
    return (fabs(y) <= tol ? true : false);
}

//! Determine if the complex number is imaginary (and non-zero)
template < typename T >
inline bool hwTComplex<T>::IsImag(T tol) const
{
    if (y == (T)0)
        return false;

    return (abs(x) <= tol ? true : false);
}

//! Determine if the complex number is imaginary (and non-zero)
template <>
inline bool hwTComplex<float>::IsImag(float tol) const
{
    if (y == 0.0f)
        return false;

    return (fabsf(x) <= tol ? true : false);
}

//! Determine if the complex number is imaginary (and non-zero)
template <>
inline bool hwTComplex<double>::IsImag(double tol) const
{
    if (y == 0.0)
        return false;

    return (fabs(x) <= tol ? true : false);
}

//! Check equality with a complex number
template < typename T >
inline bool hwTComplex<T>::IsEqual(const hwTComplex<T>& z, T tol) const
{
    if (IsNaN_T(x) || IsNaN_T(y) || IsNaN_T(z.x) || IsNaN_T(z.y) || IsNaN_T(tol))
        return false;

    if (!AreEqual(x, z.x, tol))
        return false;

    if (!AreEqual(y, z.y, tol))
        return false;

    return true;
}

//! Check equality with a real number
template < typename T >
inline bool hwTComplex<T>::IsEqual(T a, T tol) const
{
    if (IsNaN_T(x) || IsNaN_T(y) || IsNaN_T(a) || IsNaN_T(tol))
        return false;

    if (!AreEqual(x, a, tol))
        return false;

    if (!AreEqual(y, 0.0, tol))
        return false;

    return true;
}

//*******************************************************************
//                 hwTComplex operator functions
//*******************************************************************
//! Equality operator with complex right hand side
template < typename T >
inline bool hwTComplex<T>::operator==(const hwTComplex<T>& z) const
{
    if (x != z.x || y != z.y)
        return false;

    return true;
}

//! Equality operator with real right hand side
template < typename T >
inline bool hwTComplex<T>::operator==(T a) const
{
    if (x != a || y != 0.0)
        return false;

    return true;
}

//! Inequality operator with complex right hand side
template < typename T >
inline bool hwTComplex<T>::operator!=(const hwTComplex<T>& z) const
{
    if (x != z.x || y != z.y)
        return true;

    return false;
}

//! Inequality operator with real right hand side
template < typename T >
inline bool hwTComplex<T>::operator!=(T a) const
{
    if (x != a || y != 0.0)
        return true;

    return false;
}

//! Addition assignment operator
template < typename T >
inline hwTComplex<T>& hwTComplex<T>::operator+=(const hwTComplex<T>& z)
{
    x += z.x;
    y += z.y;
    return *this;
}

//! Addition assignment operator
template < typename T >
inline hwTComplex<T>& hwTComplex<T>::operator+=(T a)
{
    x += a;
    return *this;
}

//! Addition operator
template < typename T >
inline hwTComplex<T> hwTComplex<T>::operator+(const hwTComplex<T>& z) const
{
    hwTComplex<T> c;

    c.x = x + z.x;
    c.y = y + z.y;

    return c;
}

//! Addition operator
template < typename T >
inline hwTComplex<T> hwTComplex<T>::operator+(T a) const
{
    hwTComplex<T> c;

    c.x = x + a;
    c.y = y;

    return c;
}

//! Subtraction assignment operator
template < typename T >
inline hwTComplex<T>& hwTComplex<T>::operator-=(const hwTComplex<T>& z)
{
    x -= z.x;
    y -= z.y;
    return *this;
}

//! Subtraction assignment operator
template < typename T >
inline hwTComplex<T>& hwTComplex<T>::operator-=(T a)
{
    x -= a;
    return *this;
}

//! Subtraction operator
template < typename T >
inline hwTComplex<T> hwTComplex<T>::operator-(const hwTComplex<T>& z) const
{
    hwTComplex<T> c;

    c.x = x - z.x;
    c.y = y - z.y;

    return c;
}

//! Subtraction operator
template < typename T >
inline hwTComplex<T> hwTComplex<T>::operator-(T a) const
{
    hwTComplex<T> c;

    c.x = x - a;
    c.y = y;

    return c;
}

//! Negation operator (urnary prefix)
template < typename T >
inline hwTComplex<T> hwTComplex<T>::operator-() const
{
    hwTComplex<T> c(-x, -y);

    return c;
}

//! Multiplication operator
template < typename T >
inline hwTComplex<T> hwTComplex<T>::operator*(const hwTComplex<T>& z) const
{
    hwTComplex<T> c;

    c.x = x * z.x - y * z.y;
    c.y = x * z.y + y * z.x;

    return c;
}

//! Multiplication operator
template < typename T >
inline hwTComplex<T> hwTComplex<T>::operator*(T a) const
{
    hwTComplex<T> c;

    c.x = x * a;
    c.y = y * a;

    return c;
}

//! Multiplication assignment operator
template < typename T >
inline hwTComplex<T>& hwTComplex<T>::operator*=(const hwTComplex<T>& z)
{
    *this = (*this) * z;
    return *this;
}

//! Multiplication assignment operator
template < typename T >
inline hwTComplex<T>& hwTComplex<T>::operator*=(T a)
{
    x *= a;
    y *= a;
    return *this;
}

//! Division operator
template < typename T >
inline hwTComplex<T> hwTComplex<T>::operator/(T a) const
{
    hwTComplex<T> c;

    c.x = x / a;
    c.y = y / a;

    return c;
}

//! Division assignment operator
template < typename T >
inline hwTComplex<T>& hwTComplex<T>::operator/=(T a)
{
    x /= a;
    y /= a;
    return *this;
}

//! Division operator
template < typename T >
inline hwTComplex<T> hwTComplex<T>::operator/(const hwTComplex<T>& z) const
{
    T magSq = z.MagSq();
    hwTComplex<T> c;

    //c = (*this) * z.Conjugate();

    // if (magSq > 0.0)
    //c /= magSq;

    // ieee 754 complex division with 0, NaN
    if ((z.x == 0.) && (z.y == 0))
    {
        if (((x == 0.) && (y == 0.)) || ((x != x) || (y != y)))
        {
            // NaN
            double d = 0;
            c.x = 0. / d;
            c.y = 0.;

        }
        else
        {
            // Inf
            double d = 0;
            c.x = 1. / d;
            c.y = 0.;
        }
    }
    else
    {
        if (fabs(z.x) > fabs(z.y))
        {
            T den = z.x + z.y*(z.y/z.x);
            T realNum = x + y*(z.y/z.x);
            T imagNum = y - x*(z.y/z.x);

            c.x = realNum / den;
            c.y = imagNum / den;
        }
        else
        {
            T den = z.x*(z.x/z.y) + z.y;
            T realNum = x*(z.x/z.y) + y;
            T imagNum = y*(z.x/z.y) - x;

            c.x = realNum / den;
            c.y = imagNum / den;
        }
    }
    return c;
}

//! Division assignment operator
template < typename T >
inline hwTComplex<T>& hwTComplex<T>::operator/=(const hwTComplex<T>& z)
{
    *this = (*this) / z;
    return *this;
}

//! Complex sine
template < typename T >
inline hwTComplex<T> hwTComplex<T>::sin(const hwTComplex<T>& z)
{
    return hwTComplex<T>(::sin(z.x) * ::cosh(z.y), ::cos(z.x) * ::sinh(z.y));
}

//! Complex cosine
template < typename T >
inline hwTComplex<T> hwTComplex<T>::cos(const hwTComplex<T>& z)
{
    return hwTComplex<T>(::cos(z.x) * ::cosh(z.y), -::sin(z.x) * ::sinh(z.y));
}

//! Complex tangent
template < typename T >
inline hwTComplex<T> hwTComplex<T>::tan(const hwTComplex<T>& z)
{
    T value1 = ::tan(z.x);
    T value2 = ::tanh(z.y);
    hwTComplex<T> denom((T)1, -value1 * value2);

    return hwTComplex<T>(value1, value2) / denom;
}

//! Complex inverse sine
template < typename T >
inline hwTComplex<T> hwTComplex<T>::asin(const hwTComplex<T>& z)
{
    hwTComplex<T> i((T)0, (T)1);

    return -i * log(i*z + sqrt(-z*z + (T)1));
}

//! Complex inverse cosine
template < typename T >
inline hwTComplex<T> hwTComplex<T>::acos(const hwTComplex<T>& z)
{
    hwTComplex<T> i((T)0, (T)1);

    return i * log(i*z + sqrt(-z*z + 1.0)) + (T)0.5 * (T)PI;
}

//! Complex inverse tangent
template < typename T >
inline hwTComplex<T> hwTComplex<T>::atan(const hwTComplex<T>& z)
{
    hwTComplex<T> i((T)0, (T)1);

    return i * (T)0.5 * (log(-i*z + (T)1) - log(i*z + (T)1));
}

//! Complex hyperbolic sine
template < typename T >
inline hwTComplex<T> hwTComplex<T>::sinh(const hwTComplex<T>& z)
{
    return hwTComplex<T>(::sinh(z.x) * ::cos(z.y), ::cosh(z.x) * ::sin(z.y));
}

//! Complex hyperbolic cosine
template < typename T >
inline hwTComplex<T> hwTComplex<T>::cosh(const hwTComplex<T>& z)
{
    return hwTComplex<T>(::cosh(z.x) * ::cos(z.y), ::sinh(z.x) * ::sin(z.y));
}

//! Complex hyperbolic tangent
template < typename T >
inline hwTComplex<T> hwTComplex<T>::tanh(const hwTComplex<T>& z)
{
    T value1 = ::tanh(z.x);
    T value2 = ::tan(z.y);
    hwTComplex<T> denom((T)1, value1 * value2);

    return hwTComplex<T>(value1, value2) / denom;
}

//! Complex inverse hyperbolic sine
template < typename T >
inline hwTComplex<T> hwTComplex<T>::asinh(const hwTComplex<T>& z)
{
    return log(z + sqrt(z*z + (T)1));
}

//! Complex inverse hyperbolic cosine
template < typename T >
inline hwTComplex<T> hwTComplex<T>::acosh(const hwTComplex<T>& z)
{
    return log(z + sqrt(z + (T)1) * sqrt(z - (T)1));
}

//! Complex inverse hyperbolic tangent
template < typename T >
inline hwTComplex<T> hwTComplex<T>::atanh(const hwTComplex<T>& z)
{
    return (log(z + (T)1) - log(-z + (T)1)) * (T)0.5;
}

//! Complex exponential
template < typename T >
inline hwTComplex<T> hwTComplex<T>::exp(const hwTComplex<T>& z)
{
    T value = ::exp(z.x);

    return hwTComplex<T>(value * ::cos(z.y), value * ::sin(z.y));
}

//! Complex logarithm (base e)
template < typename T >
inline hwTComplex<T> hwTComplex<T>::log(const hwTComplex<T>& z)
{
    // the principal value is on (-pi, pi]
    T mag = (T) z.Mag();
    T arg = (T) z.Arg();

    return hwTComplex<T>(::log(mag), arg);
}

//! Complex logarithm (base 2)
template < typename T >
inline hwTComplex<T> hwTComplex<T>::log2(const hwTComplex<T>& z)
{
    return log(z) / ::log((T)2);
}

//! Complex logarithm (base 10)
template < typename T >
inline hwTComplex<T> hwTComplex<T>::log10(const hwTComplex<T>& z)
{
    return log(z) / ::log((T)10);
}

//! Complex square root
template < typename T >
inline hwTComplex<T> hwTComplex<T>::sqrt(const hwTComplex<T>& z)
{
    // the principal value is on (-pi, pi]
    // this method bypasses the trig using the half angle
    // identities to improve accuracy
    T mag = (T) z.Mag();

    if (z.y < 0)
        return hwTComplex<T>(::sqrt((T)0.5 * (mag + z.x)), -::sqrt((T)0.5 * (mag - z.x)));
    else
        return hwTComplex<T>(::sqrt((T)0.5 * (mag + z.x)), ::sqrt((T)0.5 * (mag - z.x)));
}

//! Power function for complex base to a real exponent
template < typename T >
inline hwTComplex<T> hwTComplex<T>::pow(const hwTComplex<T>& base, T power)
{
    T mag = base.Mag();
    mag = CustomPow(mag, power);
    T arg = base.Arg();
    arg *= power;

    if (arg != (T)0)
    {
        // check for pure real or imaginary solutions
        T factor = (T)2 * arg / PI; // = arg / (PI/2)

        if (IsInteger(factor, 1.0e-15).IsOk())      // needs modification if T != double
        {
            int numQuads = ((int) RoundT(factor)) %4;

            if (numQuads == 0)   // arg = 0
                return hwTComplex<T>(mag, (T)0);
            if (numQuads == 1)   // arg = PI/2
                return hwTComplex<T>((T)0, mag);
            if (numQuads == 2)   // arg = PI
                return hwTComplex<T>(-mag, (T)0);
            if (numQuads == 3)   // arg = -PI/2
                return hwTComplex<T>((T)0, -mag);
        }

        return hwTComplex<T>(mag * ::cos(arg), mag * ::sin(arg));
    }

    return hwTComplex<T>(mag, (T)0);
}

//! Power function for real base to a complex exponent
template < typename T >
inline hwTComplex<T> hwTComplex<T>::pow(T base, const hwTComplex<T>& power)
{
	if (base != base)     // base is a NaN
        return hwTComplex<T>(std::numeric_limits<T>::quiet_NaN(), 0);

	if (power.x != power.x || power.y != power.y) // power is a NaN
        return hwTComplex<T>(std::numeric_limits<T>::quiet_NaN(), 0);

    if (power.x == (T)0 && power.y == (T)0)
        return hwTComplex<T>((T)1, (T)0);
    if (base == (T)0)
        return hwTComplex<T>((T)0, (T)0);
    if (base < (T)0)
    {
        if (power.y == (T) 0)
            return pow_c(base, power.x);

        return exp(power * log_c(base));
    }
    // base > (T)0
    return exp(power * ::log(base));
}

//! Power function for complex base to a complex exponent
template < typename T >
inline hwTComplex<T> hwTComplex<T>::pow(const hwTComplex<T>& base, const hwTComplex<T>& power)
{
    if (power.y == (T)0)
        return pow(base, power.x);

    if (base.y == (T)0)
        return pow(base.x, power);

	if (base.x != base.x || base.y != base.y) // base is a NaN
        return hwTComplex<T>(std::numeric_limits<T>::quiet_NaN(), 0);

	if (power.x != power.x || power.y != power.y) // power is a NaN
        return hwTComplex<T>(std::numeric_limits<T>::quiet_NaN(), 0);

    return exp(power * log(base));
}

//! Inverse sine of a real argument
template < typename T >
inline hwTComplex<T> hwTComplex<T>::asin_c(T x)
{
    // for |x| > 1
    hwTComplex<T> i((T) 0, (T) 1);

    return -i * log(i*x + sqrt_c(1.0 - x*x));
}

//! Inverse cosine of a real argument
template < typename T >
inline hwTComplex<T> hwTComplex<T>::acos_c(T x)
{
    // for |x| > 1
    hwTComplex<T> i((T)0, (T)1);

    return i * log(i*x + sqrt_c((T)1 - x*x)) + (T)0.5 * (T)PI;
}

//! Real inverse hyperbolic cosine producing complex output
template < typename T >
inline hwTComplex<T> hwTComplex<T>::acosh_c(T x)
{
    // for x < 1
    if (x < (T)-1)
        // return log(sqrt_c(x + 1.0) * sqrt_c(x - 1.0) + x);
        return hwTComplex<T>(::log(::sqrt(x*x - (T)1) - x), (T)PI); // simplified

    else
        // return log(sqrt_c(x - 1.0) * sqrt(x + 1.0) + x);
        return hwTComplex<T>(0, atan2(::sqrt((T)1 - x*x), x));    // simplified
}

//! Real inverse hyperbolic tangent producing complex output
template < typename T >
inline hwTComplex<T> hwTComplex<T>::atanh_c(T x)
{
    // for |x| > 1
    if (x > 1)
        // return (log(x + 1.0) - log_c(1.0 - x)) * 0.5;
        return hwTComplex<T>(::log(x+(T)1) - ::log(x-(T)1), (T)PI) * (T)0.5;    // simplified
    else
        // return (log_c(x + 1.0) - log(1.0 - x)) * 0.5;
        return hwTComplex<T>(::log(-x-(T)1) - ::log((T)1-x), (T)PI) * (T)0.5;   // simplified
}

//! Logarithm base e of a real argument producing complex value
template < typename T >
inline hwTComplex<T> hwTComplex<T>::log_c(T x)
{
    // for x < 0
    return hwTComplex<T>(::log(-x), (T)PI);
}

//! Logarithm base 2 of a real argument producing complex value
template < typename T >
inline hwTComplex<T> hwTComplex<T>::log2_c(T x)
{
    // for x < 0
    return log_c(x) / ::log((T) 2);
}

//! Logarithm base 10 of a real argument producing complex value
template < typename T >
inline hwTComplex<T> hwTComplex<T>::log10_c(T x)
{
    // for x < 0
    return log_c(x) / ::log((T) 10);
}

//! Square root of a real argument producing complex output
template < typename T >
inline hwTComplex<T> hwTComplex<T>::sqrt_c(T x)
{
    // for x < 0
    return hwTComplex<T>((T)0, ::sqrt(-x));
}

//! Power function for real base and exponent, producing complex result 
template < typename T >
inline hwTComplex<T> hwTComplex<T>::pow_c(T base, T power)
{
    T mag = CustomPow(fabs(base), power);

    if (base < (T) 0)
    {
        // check for pure real or imaginary solutions
        T factor = (T)2 * power; // = arg / (PI/2)

        if (IsInteger(factor, 1.0e-15).IsOk())      // needs modification of T != double
        {
            int numQuads = ((int) RoundT(factor)) %4;

            if (numQuads == 0)   // arg = 0
                return hwTComplex<T>(mag, (T)0);
            if (numQuads == 1)   // arg = PI/2
                return hwTComplex<T>((T)0, mag);
            if (numQuads == 2)   // arg = PI
                return hwTComplex<T>(-mag, (T)0);
            if (numQuads == 3)   // arg = -PI/2
                return hwTComplex<T>((T)0, -mag);
        }

        return hwTComplex<T>(mag * ::cos(PI*power), mag * ::sin(PI*power));
    }
    else
        return hwTComplex<T>(mag, (T)0);
}

//! Write the value to the output stream
template < typename T >
inline void hwTComplex<T>::emitVarVal (std::ostringstream &os)
{ 
    os << x;

    if (y)
        os << " + " << y << "i";
}
