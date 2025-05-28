/**
* @file omljson.h
* @date September 2023
* Copyright (C) 2023 Altair Engineering, Inc.
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
#ifndef __OMLJSON_H__
#define __OMLJSON_H__

#include "omljsondefs.h"

#include <map>

#include "EvaluatorInt.h"

//!
//! Entry point which registers library with oml
//! \param Evaluator interface
//!
extern "C" OMLJSON_DECLS int InitDll(EvaluatorInterface);
//!
//! Returns library version
//! \param Evaluator interface
//!
extern "C" OMLJSON_DECLS double GetToolboxVersion(EvaluatorInterface);

//------------------------------------------------------------------------------
//! \class OmlJSON
//! \brief Utility class for JSON related functions in oml
//------------------------------------------------------------------------------
class OMLJSON_DECLS OmlJSON
{
public:
    //!
    //! Destructor
    //!
    ~OmlJSON() {}

    //!
    //! Decodes JSON text to oml data and returns true [jsondecode]
    //! \param Evaluator
    //! \param Inputs
    //! \param Outputs
    //!
    static bool Decode(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Encodes oml data to JSON text and returns true [jsonencode]
    //! \param Evaluator
    //! \param Inputs
    //! \param Outputs
    //!
    static bool Encode(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);

private:

    //!
    //! \enum JSON decode replacement style for invalid characters
    //!
    enum REPLACESTYLE
    {
        REPLACESTYLE_DELETE,     //! Delete invalid chars
        REPLACESTYLE_HEX,        //! Replace invalid chars with hex notation
        REPLACESTYLE_UNDERSCORE  //! Replace invalid chars with underscore
    };

    bool         _convert;   //!< Encode: Convert invalid scalars to null
    bool         _format;    //!< Encode: Output has indentation and newlines
    std::string  _newline;   //!< Encode: New line char for formatting

    bool         _checkname; //!< Decode: Check for valid field names
    std::string  _prefix;    //!< Decode: Prefix added for numeric field names  
    REPLACESTYLE _style;     //!< Decode: Handling invalid characters in struct fields

    //!
    //! Constructor
    //! 
    OmlJSON();

    OmlJSON(const OmlJSON&);              // Stubbed out 
    OmlJSON& operator=(const OmlJSON&);   // Stubbed out

    // Utilities
   
    //!
    //! Gets comma string, if needed
    //! \param Index
    //! \param Max
    //! 
    static std::string Comma(int, int);
    //!
    //! Gets value between open brackets/parenthesis or quotes
    //! \param Input string
    //!
    static std::string ContainerValue(std::string&);
    //!
    //! Returns true if this is a null string
    //! \param Input
    //!
    static bool IsNull(const std::string&);
    //!
    //! Returns true if this is a a valid prefix
    //! \param Input
    //!
    static bool IsValidPrefix(const std::string&);
    //!
    //! Trims space and new line characters from beginning and end
    //! \param Input string
    //! 
    static void Trim(std::string&);

    //!
    //! Gets the indent string
    //! \param Level of the currency hierarchy
    //!
    std::string Indent(int) const;
    //!
    //! Increments current indent, if formatting
    //! \param Current indent
    //! 
    std::string IncrIndent(const std::string&) const;
    //!
    //! Splits a row into substrings with comma delimiter
    //! \param Input
    //! 
    std::vector<std::string> SplitRow(const std::string&) const;

    // Conversion functions

    //!
    //! Converts currency to JSON format string
    //! \param Input
    //! \param Level of the currency hierarchy
    //! 
    std::string Currency2JSON(const Currency&, int = 0) const;
    //!
    //! Converts JSON format string to Currency
    //! \param Input
    //! 
    Currency JSON2Currency(const std::string&) const;

    //!
    //! Converts cell to JSON format string and returns it
    //! \param Input
    //! \param Level of the currency hierarchy
    //! 
    std::string Cell2JSON(const Currency&, int) const;
    //!
    //! Converts JSON format string to a cell and returns it
    //! \param Input
    //! 
    Currency JSON2Cell(const std::string&) const;

    //!
    //! Returns true if successful in converting from string to matrix
    //! \param Input
    //! 
    bool JSON2Matrix(const std::string&, Currency&) const;
    //!
    //! Returns JSON format string after converting a matrix
    //! \param Input
    //! \param Level of the currency hierarchy
    //! 
    std::string Matrix2JSON(const Currency&, int) const;

    //! 
    //! Converts complex to JSON string
    //! \param Input
    //! \param Level of indent
    //! 
    std::string Complex2JSON(const Currency&, int) const;
    //!
    //! Returns true if conversion to scalar is valid
    //! \param Input
    //! \param Scalar currency
    //! 
    bool JSON2Scalar(const std::string&, Currency&) const;
    //!
    //! Converts scalar to JSON format string with no indent
    //! \param Input
    //! \param True if input is logical
    //! 
    std::string Scalar2JSON(double, bool = false) const;

    //!
    //! Converts multiline strings to JSON format string
    //! \param Input
    //! \param Level of the currency, 0 if root level
    //! 
    std::string StringND2JSON(const Currency&, int) const;
    //!
    //! Converts single dimension string to JSON format string
    //! \param Input
    //! 
    std::string String2JSON(const std::string&) const;
    //!
    //! Converts JSON string to string currency
    //! \param Input
    //! 
    Currency JSON2String(const std::string&) const;


    //!
    //! Creates the field name
    //! \param Input
    //! 
    std::string FieldName(const std::string&) const;
    //!
    //! Gets the value, given the field
    //! \param Field name
    //! 
    std::string GetFieldValueString(std::string&) const;
    //!
    //! Gets name-value pairs for a struct
    //! \param Input string
    //! 
    void ParseStructData(const std::string&, std::map<std::string, std::vector<Currency> >&) const;
    //!
    //! Converts struct array to JSON format string
    //! \param Input
    //! \param Output currency
    //! 
    bool JSON2StructArray(const std::string&, Currency&) const;
    //!
    //! Utility to convert to JSON format string to struct
    //! \param Input
    //! 
    Currency JSON2Struct(const std::string&) const;
    //!
    //! Converts struct to JSON format string
    //! \param Input
    //! \param Level of the currency, 0 if root level
    //! 
    std::string Struct2JSON(const Currency&, int) const;
};
#endif

