/**
* @file hwTComplex.h
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
//  Template complex number class
//
//:---------------------------------------------------------------------------
#ifndef _hwTComplex_h
#define _hwTComplex_h

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

class hwMathStatus; // forward declaration

//*****************************************************************************
//********************************* WARNING ***********************************
//*****************************************************************************
//
// This class has been designed for use with HyperMath. It is **CRITICAL** that
// no virtual functions or other data members be added to it. Please contact the
// HyperMath team if you have a need to modify the class.
//
//*****************************************************************************
//********************************* WARNING ***********************************
//*****************************************************************************
template < typename T >
class hwTComplex
{
public:
    // The following functions are implemented below the class definition
    //! Default constructor
    hwTComplex();
    //! Constructor specific object with specified value
    hwTComplex(T real, T imag);
    //! Copy constructor
    hwTComplex(const hwTComplex<T>& z);
    //! Destructor
    ~hwTComplex();

protected:
    //! Real component
    T x;
    //! Imaginary component
    T y;

public:
    //! Set the real and imaginary components
    void Set(T real, T imag) { x = real; y = imag; }
    //! Return the real component
    inline T& Real() { return x; }
    //! Return a reference to the real component
    inline const T& Real() const { return x; }
    //! Set the real component
    inline void Real(T a) { x = a; }
    //! Return the imaginary component
    inline T& Imag() { return y; }
    //! Return a reference to the imaginary component
    inline const T& Imag() const { return y; }
    //! Set the real component
    inline void Imag(T a) { y = a; }
    //! Implement the = operator
    hwTComplex<T>& operator=(const hwTComplex<T>& z);
    //! Implement the = operator
    hwTComplex<T>& operator=(T x);
    //! Compute the magnitude
    T Mag() const;
    //! Compute the squared magnitude
    inline T MagSq() const { return x*x + y*y; }
    //! Compute the phase
    T Arg() const;
    //! Return the conjugate
    hwTComplex<T> Conjugate() const { return hwTComplex<T>(x, -y); }
    //! Determine if the complex number is real
    bool IsReal(T tol = (T) 0) const;
    //! Determine if the complex number is imaginary (and non-zero)
    bool IsImag(T tol = (T) 0) const;
    //! Determine if the complex number equals a complex value
    bool IsEqual(const hwTComplex<T>& z, T tol = (T) 0) const;
    //! Determine if the complex number equals a real value
    bool IsEqual(T a, T tol = (T) 0) const;
    //! Implement the == operator
    bool operator==(const hwTComplex<T>& z) const;
    //! Implement the == operator
    bool operator==(T a) const;
    //! Implement the != operator
    bool operator!=(const hwTComplex<T>& z) const;
    //! Implement the != operator
    bool operator!=(T a) const;
    //! Implement the += operator
    hwTComplex<T>& operator+=(const hwTComplex<T>& z);
    //! Implement the += operator
    hwTComplex<T>& operator+=(T a);
    //! Implement the + operator
    hwTComplex<T> operator+(const hwTComplex<T>& z) const;
    //! Implement the + operator
    hwTComplex<T> operator+(T a) const;
    //! Implement the -= operator
    hwTComplex<T>& operator-=(const hwTComplex<T>& z);
    //! Implement the -= operator
    hwTComplex<T>& operator-=(T a);
    //! Implement the - operator
    hwTComplex<T> operator-(const hwTComplex<T>& z) const;
    //! Implement the - operator
    hwTComplex<T> operator-(T a) const;
    //! Implement the - urnary prefix operator
    hwTComplex<T> operator-() const;
    //! Implement the * operator
    hwTComplex<T> operator*(const hwTComplex<T>& z) const;
    //! Implement the * operator
    hwTComplex<T> operator*(T a) const;
    //! Implement the *= operator
    hwTComplex<T>& operator*=(const hwTComplex<T>& z);
    //! Implement the *= operator
    hwTComplex<T>& operator*=(T a);
    //! Implement the / operator
    hwTComplex<T> operator/(T a) const;
    //! Implement the /= operator
    hwTComplex<T>& operator/=(T a);
    //! Implement the / operator
    hwTComplex<T> operator/(const hwTComplex<T>& z) const;
    //! Implement the /= operator
    hwTComplex<T>& operator/=(const hwTComplex<T>& z);

    //! Sine
    static hwTComplex<T> sin(const hwTComplex<T>& z);
    //! Cosine
    static hwTComplex<T> cos(const hwTComplex<T>& z);
    //! Tangent
    static hwTComplex<T> tan(const hwTComplex<T>& z);
    //! Inverse sine
    static hwTComplex<T> asin(const hwTComplex<T>& z);
    //! Inverse cosine
    static hwTComplex<T> acos(const hwTComplex<T>& z);
    //! Inverse tangent
    static hwTComplex<T> atan(const hwTComplex<T>& z);
    //! Hyperbolic sine
    static hwTComplex<T> sinh(const hwTComplex<T>& z);
    //! Hyperbolic cosine
    static hwTComplex<T> cosh(const hwTComplex<T>& z);
    //! Hyperbolic tangent
    static hwTComplex<T> tanh(const hwTComplex<T>& z);
    //! Inverse hyperbolic sine
    static hwTComplex<T> asinh(const hwTComplex<T>& z);
    //! Inverse hyperbolic cosine
    static hwTComplex<T> acosh(const hwTComplex<T>& z);
    //! Inverse hyperbolic tangent
    static hwTComplex<T> atanh(const hwTComplex<T>& z);
    //! Exponential
    static hwTComplex<T> exp(const hwTComplex<T>& z);
    //! Logarithm base e
    static hwTComplex<T> log(const hwTComplex<T>& z);
    //! Logarithm base 2
    static hwTComplex<T> log2(const hwTComplex<T>& z);
    //! Logarithm base 10
    static hwTComplex<T> log10(const hwTComplex<T>& z);
    //! Square root
    static hwTComplex<T> sqrt(const hwTComplex<T>& z);
    //! Power function for complex base to a real exponent
    static hwTComplex<T> pow(const hwTComplex<T>& base, T power);
    //! Power function for real base to a complex exponent
    static hwTComplex<T> pow(T base, const hwTComplex<T>& power);
    //! Power function for complex base to a complex exponent
    static hwTComplex<T> pow(const hwTComplex<T>& base, const hwTComplex<T>& power);
    //! Inverse sine of a real argument
    static hwTComplex<T> asin_c(T x);
    //! Inverse cosine of a real argument
    static hwTComplex<T> acos_c(T x);
    //! Inverse hyperbolic cosine of a real argument
    static hwTComplex<T> acosh_c(T x);
    //! Inverse hyperbolic tangent of a real argument
    static hwTComplex<T> atanh_c(T x);
    //! Logarithm base e of a real argument
    static hwTComplex<T> log_c(T x);
    //! Logarithm base 2 of a real argument
    static hwTComplex<T> log2_c(T x);
    //! Logarithm base 10 of a real argument
    static hwTComplex<T> log10_c(T x);
    //! Square root of a real argument
    static hwTComplex<T> sqrt_c(T x);
    //! Power function for real base and exponent, producing complex result 
    static hwTComplex<T> pow_c(T base, T power);
    //! Write the value to the output stream
    void emitVarVal(std::ostringstream &os);
};

//! template implementation file
#include <tmpl/hwTComplex.cc>

//! template utility function file
#include <utl/hwTComplexUtil.cc>

#endif // _hwTComplex_h

