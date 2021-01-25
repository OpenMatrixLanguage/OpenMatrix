/**
* @file DisplayFormatVars.h
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
#ifndef __DISPLAYFORMATVARS_H__
#define __DISPLAYFORMATVARS_H__

// Begin defines/includes
#include "OMLDll.h"

#include <string>

#include "hwComplex.h"

class OutputFormat;
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//! \class DisplayFormatVars
//! \brief Contains variables needed for formatting displays
//------------------------------------------------------------------------------
class OMLDLL_DECLS DisplayFormatVars
{
public:

    //!
    //! Constructor
    //!
    DisplayFormatVars();
    //!
    //! Destructor
    //!
    ~DisplayFormatVars() {}

    //!
    //! \enum Format type
    //!
    enum Type
    {
        TypeInt,              //!< Integer format
        TypeFloat,            //!< Float format
        TypeScientific,       //!< Scientific format
    };

    //!
    //! Initializes variables based on output format
    //! \param fmt Output format
    //!
    void Initialize(const OutputFormat* fmt);
    //!
    //! Resets variables
    //!
    void Reset();

    //!
    //! Gets format type for scalars
    //! \param val Given value
    //!
    Type GetFormatType(double val);
    //!
    //! Gets format type for complex numbers
    //! \param val Given value
    //!
    Type GetFormatType(const hwComplex& val);
    //!
    //! Gets the default format type
    //!
    Type GetFormatType() const { return _type;  }
    //!
    //! Sets the format type
    //! \param val Value to set
    //!
    void SetFormatType(Type val) { _type = val;  }

    //!
    //! Converts real value to formatted string
    //! \param val Value
    //!
    std::string RealToString(double val) const { return RealToString(val, _type); }
    //!
    //! Converts real value to formatted string
    //! \param val  Value
    //! \param type Type
    //!
    std::string RealToString(double val,
                             Type   fmt) const;
    //!
    //! Utility which returns real value to string
    //! \param val       Value
    //! \param precision Custom precision
    //!
    static std::string RealToString(double             val,
                                    const std::string& precision);
    //!
    //! Returns string to display if value is Inf/Nan/NegInf
    //! \param val
    //!
    static std::string GetOutOfRangeOutput(double val);

private:
    bool             _haslargeint;   //!< True if data has large ints
    bool             _uppercase;     //!< True if scientific uppercase
    int              _formatinteger; //!< Integer part for format, if applicable
    int              _formatdecimal; //!< Decimal part for format, if applicable
    Type             _type;          //!< Display type 
    std::streamsize  _precision;     //!< Precision

    long long   _maxdigits;          //!< Max digits
    long double _maxfloat;           //!< Max possible float
    long long   _maxint;             //!< Max possible int
    long double _minfloat;           //!< Min possible float
    size_t      _totaldigits;        //!< Total digits

    // Stubbed out default, copy constructors and assignment operator
    DisplayFormatVars(const DisplayFormatVars& src);
    DisplayFormatVars& operator=(const DisplayFormatVars& src);

    //!
    //! Initializes for short format
    //!
    void InitializeForShortFormat();
};
#endif
// End of file:

