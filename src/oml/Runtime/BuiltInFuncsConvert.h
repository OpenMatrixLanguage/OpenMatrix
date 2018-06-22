/**
* @file BuiltInFuncsConvert.h
* @date November 2015
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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

#ifndef __BUILTINFUNCSCONVERT__
#define __BUILTINFUNCSCONVERT__

#include "EvaluatorInt.h"

//------------------------------------------------------------------------------
//!
//! \brief Utility class for built-in functions related to datatype conversions
//!
//------------------------------------------------------------------------------
class HML2DLL_DECLS BuiltInFuncsConvert
{
public:
    //!
    //! Destructor
    //!
    ~BuiltInFuncsConvert() {}
    //!
    //! Returns true after converting decimal to binary [dec2bin command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool hml_dec2bin( EvaluatorInterface           eval,
                             const std::vector<Currency>& inputs,
                             std::vector<Currency>&       outputs);
    //!
    //! Returns true after converting binary string to decimal [bin2dec command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool hml_bin2dec( EvaluatorInterface           eval,
                             const std::vector<Currency>& inputs,
                             std::vector<Currency>&       outputs);
    //!
    //! Returns true after converting decimal to hexdecimal [dec2hex command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool hml_dec2hex( EvaluatorInterface           eval,
                             const std::vector<Currency>& inputs,
                             std::vector<Currency>&       outputs);
    //!
    //! Returns true after converting hexdecimal to decimal [hex2dec command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool hml_hex2dec( EvaluatorInterface           eval,
                             const std::vector<Currency>& inputs,
                             std::vector<Currency>&       outputs);

private:
    //!
    //! Constructor
    //!
    BuiltInFuncsConvert() {}
    //!
    //! Returns true if given value is a non-negative integer
    //! \param val Given value
    //!
    static bool IsValidNaturalNumber( double val);

    //!
    // Helper funtions for dec2hex
    //! Returns hex string given value and minimum width to pad element with
    //! \param value Value
    //! \param width Minimum width to pad strings with
    //!
    static std::string Dec2HexHelper( double value, 
                                      int    width);
    //!
    //! Returns vector of hex strings, given matrix and minimum width
    //! \param mtx   Given matrix
    //! \param width Minimum width to pad strings with
    //!
    static std::vector<std::string> Dec2HexHelper( const hwMatrix* mtx, 
                                                   int             width);
    //!
    //! Returns vector of hex strings, given cell array and minimum width
    //! \param cell  Given cell
    //! \param width Minimum width to pad strings with
    //!
    static std::vector<std::string> Dec2HexHelper( HML_CELLARRAY* cell, 
                                                   int            width);
    //!
    // Helper funtions for hex2dec
    //! Returns currency given hex string
    //! \param hexstr Hexadecimal string
    //!
    static Currency Hex2DecHelper( const std::string& hexstr); 
    //!
    //! Returns currency (column matrix) for given input
    //! \param eval Evaluator interface
    //! \param cur  Given currency
    //!
    static Currency Hex2DecHelper( EvaluatorInterface& eval,
                                   const Currency&     cur); 
    //!
    // Helper funtions for dec2bin
    //! Returns binary string given value and minimum width to pad element with
    //! \param val   Value
    //! \param width Minimum width to pad strings with
    //!
    static std::string Dec2BinHelper( double val, 
                                      int    width);
    //!
    //! Returns vector of binary strings, given matrix and minimum width
    //! \param mtx   Given matrix
    //! \param width Minimum width to pad strings with
    //!
    static std::vector<std::string> Dec2BinHelper( const hwMatrix* mtx, 
                                                   int             width);
    //!
    //! Returns vector of binary strings, given cell array and minimum width
    //! \param cell  Given cell
    //! \param width Minimum width to pad strings with
    //!
    static std::vector<std::string> Dec2BinHelper( HML_CELLARRAY* cell, 
                                                   int            width);
    //!
    // Helper funtions for bin2dec - returns currency given hex string
    //! \param hexstr Hexadecimal string
    //!
    static Currency Bin2DecHelper( const std::string& hexstr); 
    //!
    //! Returns currency (column matrix) for given input
    //! \param eval Evaluator interface
    //! \param cur  Given currency
    //!
    static Currency Bin2DecHelper( EvaluatorInterface& eval,
                                   const Currency&     cur); 



};
#endif // __BUILTINFUNCSCONVERT__
