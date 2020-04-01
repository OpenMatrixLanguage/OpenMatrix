/**
* @file DisplayFormatVars.cpp
* @date June, 2019
* Copyright (C) 2019 Altair Engineering, Inc.
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

#include "DisplayFormatVars.h"

#include <cassert>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <string>
#include <limits.h>

#include "BuiltInFuncsUtils.h"
#include "OutputFormat.h"

#include "math/kernel/GeneralFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
DisplayFormatVars::DisplayFormatVars()
{
    Reset();
}
//------------------------------------------------------------------------------
// Gets format for given double
//------------------------------------------------------------------------------
DisplayFormatVars::Type DisplayFormatVars::GetFormatType(double val) 
{
    assert(_type != TypeScientific); // Should never get into this situation
    if (_type == TypeScientific)   
    {
        return TypeScientific; // If format is scientific, we are done
    }

    if (IsInf_T(val) || IsNegInf_T(val) || IsNaN_T(val))
    {
        return TypeInt;
    }

    double aval = fabs(val);
    BuiltInFuncsUtils utils;
    if (utils.IsInt(val)) // Integer, compare by values
    {
        // Check for very small floats
        if (IsZero(aval) && !(aval < 1e-16))
        {
            return TypeScientific;
        }

        // Max digits for integers is 12 (short format), 15 (long format)
        if ((_precision == OutputFormat::PRECISION_SHORT  && !(aval < _maxint))    ||
            (_precision == OutputFormat::PRECISION_LONG   && !(aval < _maxdigits)) ||
            (_precision == OutputFormat::PRECISION_SCALAR && !(aval < _maxdigits)))
        {
            return TypeScientific;
        }

        // Check for large integers
        if (aval > LLONG_MAX) 
        {
            return TypeScientific;
        }

        // Additional check if we have large integers
        if (!_haslargeint && _precision == OutputFormat::PRECISION_SHORT &&
            !(aval < _maxfloat))
        {
            _haslargeint = true;
        }
        else if (!_haslargeint && _precision == OutputFormat::PRECISION_SCALAR &&
            !(aval < _maxfloat))
        {
            _haslargeint = true;
        }
        std::string tmp(RealToString(val, TypeInt));
        if (tmp.size() > _totaldigits)
        {
            return TypeScientific;
        }

        if (!_haslargeint || _type == TypeInt)
        {
            return TypeInt;
        }
    }

    if (_haslargeint)
    {
        return TypeScientific;
    }

    // Additional check if we have very small floats
    if (!(aval > _minfloat))
    {
        return TypeScientific;
    }

    std::string tmp(RealToString(val, TypeFloat));
    if (tmp.size() > _totaldigits)
    {
        return TypeScientific;
    }

    return TypeFloat;
}
//------------------------------------------------------------------------------
// Returns format for complex numbers
//------------------------------------------------------------------------------
DisplayFormatVars::Type DisplayFormatVars::GetFormatType(const hwComplex& val)
{
    Type realfmt = GetFormatType(val.Real());
    if (realfmt == TypeScientific)
    {
        return TypeScientific;
    }

    Type imagfmt = GetFormatType(val.Imag());
    if (imagfmt == TypeScientific)
    {
        return TypeScientific;
    }

    return std::max(realfmt, imagfmt);
}
//------------------------------------------------------------------------------
// Returns string to display if value is Inf/Nan/NegInf
//------------------------------------------------------------------------------
std::string DisplayFormatVars::GetOutOfRangeOutput(double val)
{
    if (IsInf_T(val))
    {
        return "Inf";
    }
    else if (IsNaN_T(val))
    {
        return "NaN";
    }
    else if (IsNegInf_T(val))
    {
        return "-Inf";
    }
    return "";
}
//------------------------------------------------------------------------------
// Converts real value to formatted string
//------------------------------------------------------------------------------
std::string DisplayFormatVars::RealToString(double val, Type fmt) const
{
    std::string output(GetOutOfRangeOutput(val));
    if (!output.empty())
    {
        return output;
    }
    else if (fmt == TypeInt)
    {
        std::ostringstream os;
        os.setf(static_cast<std::ios_base::fmtflags>(0), std::ios::floatfield);
        os << std::setprecision(static_cast<std::streamsize>(_totaldigits));

        if (std::signbit(static_cast<long double>(val))) // Detect -0
        {
            os << "-" << fabs(val);
        }
        else
        {
            os << val;
        }
        return std::string(os.str());
    }

    char tmp[1024];

    if (fmt == TypeFloat && _precision == OutputFormat::PRECISION_SHORT)
    {
        sprintf(tmp, "%.5f", val);
        return std::string(tmp);
    }
    else if (fmt == TypeFloat && _precision == OutputFormat::PRECISION_LONG)
        sprintf(tmp, "%.8f", val);

    else if (fmt == TypeFloat && _precision == OutputFormat::PRECISION_SCALAR)
    {
        std::ostringstream os;
        os << "%" << _formatinteger << "." << _formatdecimal << "f";
        sprintf(tmp, std::string(os.str()).c_str(), val);
    }

    else if (fmt == TypeScientific && _precision == OutputFormat::PRECISION_SHORT)
    {
        _uppercase ? sprintf(tmp, "%.5E", val) : sprintf(tmp, "%.5e", val);
    }
    else if (fmt == TypeScientific && _precision == OutputFormat::PRECISION_LONG)
    {
        _uppercase ? sprintf(tmp, "%.8E", val) : sprintf(tmp, "%.8e", val);
    }
    else
    {
        std::ostringstream os;
        os << "%" << _formatinteger << "." << _formatdecimal;
        if (_uppercase)
            os << "E";
        else
            os << "e";
        sprintf(tmp, std::string(os.str()).c_str(), val);
    }
    return std::string(tmp);
}
//------------------------------------------------------------------------------
// Utility which returns real value to string
//------------------------------------------------------------------------------
std::string DisplayFormatVars::RealToString(double val, const std::string& precision)
{
    std::string output(GetOutOfRangeOutput(val));
    if (!output.empty())
    {
        return output;
    }

    std::string fmt(precision);
    if (fmt.empty())
    {
        fmt = "%lf";
    }

    bool isint = BuiltInFuncsUtils::IsInt(val);
    if (!isint)
    {
        size_t pos = fmt.find("d");
        if (pos != std::string::npos)
        {
            fmt.replace(pos, 2, "lf"); // For sprintf invalid formats
        }
    }

    // Handle double precision issues
    if (!isint)
    {
        if (fmt.find("lf") == std::string::npos)
        {
            size_t pos = fmt.find("g");
            if (pos == std::string::npos)
                pos = fmt.find("f");
            if (pos != std::string::npos)
                fmt.replace(pos, 2, "lf"); // For sprintf invalid formats
        }
    }
    else
    {
        size_t pos1 = fmt.find("d");
        if (pos1 != std::string::npos)
        {
            size_t pos = fmt.find(".");
            if (pos != std::string::npos)
            {
                fmt.replace(pos1, 2, "lf"); // For sprintf invalid formats
            }
        }
    }
    char tmp[1028];
    sprintf(tmp, fmt.c_str(), val);
    std::string out(tmp);
    return tmp;
}
//------------------------------------------------------------------------------
// Initializes variables based on output format
//------------------------------------------------------------------------------
void DisplayFormatVars::Initialize(const OutputFormat* fmt)
{
    if (!fmt)
    {
        return;
    }
    _formatinteger = std::min(18, fmt->GetIntegerPart());
    _formatdecimal = std::min(18, fmt->GetDecimalPart());

    if (_formatinteger > 0 && _formatdecimal > 0)
    {   // Initialization for user specified format
        _totaldigits = _formatinteger + _formatdecimal;
        _maxdigits = static_cast<long long>(std::pow(10.0, 
                                            static_cast<double>(_totaldigits)));
        _maxint = _maxdigits;
        _maxfloat = std::pow(10.0, _formatinteger);
        _minfloat = std::pow(10.0, -(_formatdecimal + 1));
        _precision = OutputFormat::PRECISION_SCALAR;
    }
    else if (fmt->GetPrecision() == OutputFormat::PRECISION_LONG)
    {   // Initialization for format long
        _formatinteger = 0;
        _formatdecimal = 0;
        _totaldigits = 24;
        _maxdigits = static_cast<long long>(1e+15);
        _minfloat = 1e-9;
        _maxint = static_cast<long long>(1e+12);
        _precision = OutputFormat::PRECISION_LONG;
    }
    else
    {   // Initialization for short format
        InitializeForShortFormat();
    }

    if (fmt->GetFlags() == std::ios::scientific)
    {
        _type = TypeScientific;
        return;
    }

    else if (fmt->GetFlags() == (std::ios::scientific | std::ios::uppercase))
    {
        _type = TypeScientific;
        _uppercase = true;
    }
}
//------------------------------------------------------------------------------
// Initialize for short format
//------------------------------------------------------------------------------
void DisplayFormatVars::InitializeForShortFormat()
{
    _formatinteger = 0;
    _formatdecimal = 0;
    _totaldigits   = 12;
    _maxdigits     = static_cast<long long>(1e+15);
    _maxint        = static_cast<long long>(1e+12);
    _maxfloat      = 1e+6;
    _minfloat      = 1e-6;
    _precision     = OutputFormat::PRECISION_SHORT;
}
//------------------------------------------------------------------------------
// Resets variables
//------------------------------------------------------------------------------
void DisplayFormatVars::Reset()
{
    _type        = TypeInt;
    _haslargeint = false;
    _uppercase   = false;
    InitializeForShortFormat();
}
