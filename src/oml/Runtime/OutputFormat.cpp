/**
* @file OutputFormat.cpp
* @date February 2016
* Copyright (C) 2016-2018 Altair Engineering, Inc.  
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

// Begin defines/includes

#include "OutputFormat.h"

#include <algorithm>

std::streamsize OutputFormat::PRECISION_SCALAR = 9;
std::streamsize OutputFormat::PRECISION_SHORT  = 5;  // %6.5f;
std::streamsize OutputFormat::PRECISION_LONG   = 8;  // %15.8f;

// End defines/includes

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
OutputFormat::OutputFormat(std::streamsize prec, std::ios_base::fmtflags flags) 
    : _precision   (prec)
    , _flags       (flags)
    , _integerpart (0)
    , _decimalpart (0)
{
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
OutputFormat::OutputFormat(int intpart, int decpart) 
    : _precision   (PRECISION_SCALAR)
    , _flags       (static_cast<std::ios_base::fmtflags>(0))
{
    // Limit the number of digits displayed to 18
    _integerpart = (std::min)(18, intpart);
    _decimalpart = (std::min)(18, decpart);
}
//------------------------------------------------------------------------------
// Resets format
//------------------------------------------------------------------------------
void OutputFormat::Reset()
{
    _precision   = PRECISION_SCALAR;
    _flags       = static_cast<std::ios_base::fmtflags>(0);
    _integerpart = 0;
    _decimalpart = 0;
}
//------------------------------------------------------------------------------
// Sets long format
//------------------------------------------------------------------------------
void OutputFormat::SetLongFormat()
{
    Reset();
    _precision = PRECISION_LONG;       
}
//------------------------------------------------------------------------------
// Sets short format
//------------------------------------------------------------------------------
void OutputFormat::SetShortFormat()
{
    Reset();
    _precision = PRECISION_SHORT;       
}
//------------------------------------------------------------------------------
// Sets custom format, specifying integer/decimal digits
//------------------------------------------------------------------------------
void OutputFormat::SetCustomFormat(int intpart, int decpart)
{
    Reset();
    // Limit the number of digits displayed to 18
    _integerpart = std::min(18, intpart);
    _decimalpart = std::min(18, decpart);
}
//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
OutputFormat::OutputFormat(const OutputFormat& src)
    : _precision   (src._precision)
    , _flags       (src._flags)
    , _integerpart (src._integerpart)
    , _decimalpart (src._decimalpart)
{
} 
//------------------------------------------------------------------------------
// Assignment op
//------------------------------------------------------------------------------
OutputFormat& OutputFormat::operator=(const OutputFormat& src)
{
    _precision   = src._precision;
    _flags       = src._flags;
    _integerpart = src._integerpart;
    _decimalpart = src._decimalpart;
    return *this;
} 
