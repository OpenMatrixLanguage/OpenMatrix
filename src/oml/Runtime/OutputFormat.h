/**
* @file OutputFormat.h
* @date February 2016
* Copyright (C) 2016-2018 Altair Engineering, Inc.  
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

// Begin defines/includes
#ifndef __OUTPUTFORMAT_H__
#define __OUTPUTFORMAT_H__

// Begin defines/includes
#include "Hml2Dll.h"

#include <iostream>
#include <string>

// End defines/includes

//------------------------------------------------------------------------------
//!
//! \class OutputFormat
//! \brief Specifies output format
//!
//------------------------------------------------------------------------------
class HML2DLL_DECLS OutputFormat
{
public:
    //!
    //! Default constructor
    //!
    OutputFormat() { Reset(); }
    //!
    //! Constructor
    //! \param prec  Precision
    //! \param flags Format flags
    //!
    OutputFormat( std::streamsize          prec, 
                  std::ios_base::fmtflags flags);
    //!
    //! Constructor
    //! \param intpart Max number of digits in integer part <= 18
    //! \param decpart Max number of digits in decimal part <= 18
    //!
    OutputFormat( int intpart,
                  int decpart); 
    //!
    //! Destructor
    //!
    ~OutputFormat() {}
        
    //!
    //! Copy constructor
    //! \param src Source format
    //!
	OutputFormat(const OutputFormat& src);
    //!
    //! Assignment op
    //! \param src Source format
    //!
	OutputFormat& operator=(const OutputFormat& src);

    static std::streamsize PRECISION_SCALAR;  //!< Scalar precision = 9
    static std::streamsize PRECISION_SHORT;   //!< Scalar precision = 5
    static std::streamsize PRECISION_LONG;    //!< Scalar precision = 8

    //!
    //! Gets precision
    //!
    std::streamsize GetPrecision() const { return _precision; }
    //!
    //! Gets format flags
    //!
    std::ios_base::fmtflags GetFlags() const { return _flags; }
    //!
    //! Sets format flags
    //! \param val Value to set, 0 is free format
    //!
    void SetFlags( std::ios_base::fmtflags val) { _flags = val; }

    //!
    //! Gets integer part for format, if applicable
    //!
    int GetIntegerPart() const { return _integerpart; }
    //!
    //! Gets decimal part for format, if applicable
    //!
    int GetDecimalPart() const { return _decimalpart; }

    //!
    //! Resets format
    //!
    void Reset();
    //!
    //! Sets long format
    //!
    void SetLongFormat();
    //!
    //! Sets short format
    //!
    void SetShortFormat();
    //! 
    //! Sets custom format, specifying integer/decimal digits
    //! \param intpart Maximum digits in integer part, <= 18
    //! \param decpart Maximum digits in decimal part, <= 18
    //!
    void SetCustomFormat( int intpart,
                          int decpart);
private:
    std::streamsize          _precision;    //!< Precision
    std::ios_base::fmtflags  _flags;        //!< Format flags
    int                      _integerpart;  //!< Int precision,     user defined
    int                      _decimalpart;  //!< Decimal precision, user defined
};

#endif

