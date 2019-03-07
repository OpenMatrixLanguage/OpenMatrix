/**
* @file BuiltInFuncsString.h
* @date November 2015
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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

#ifndef __BUILTINFUNCSSTRING__
#define __BUILTINFUNCSSTRING__

#include "EvaluatorInt.h"

#include <regex>

class OutputFormat;

//------------------------------------------------------------------------------
//!
//! \brief Utility class for built-in functions implementing string commands
//!
//------------------------------------------------------------------------------
class HML2DLL_DECLS BuiltInFuncsString
{
public:
    //!
    //! Constructor
    //!
    BuiltInFuncsString() {}
    //!
    //! Destructor
    //!
    ~BuiltInFuncsString() {}

    //!
    //! Returns true and returns a string of blanks
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool hml_blanks( EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs);
    //!
    //! Returns true and tokenizes given string
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool hml_strtok( EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs);
    //!
    //! Returns true and vertically concatenates given strings
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool hml_strvcat( EvaluatorInterface           eval,
                             const std::vector<Currency>& inputs,
                             std::vector<Currency>&       outputs);
    //!
    //! Returns true and reads formatted input from a string
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Sscanf( EvaluatorInterface           eval,
                        const std::vector<Currency>& inputs,
                        std::vector<Currency>&       outputs);
    //!
    //! Returns true and trims leading/trailing spaces from strings/cells
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Strtrim( EvaluatorInterface           eval,
                        const std::vector<Currency>& inputs,
                        std::vector<Currency>&       outputs);
    //!
    //! Returns true and converts input string to ascii row matrix
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool ToAscii( EvaluatorInterface           eval,
                         const std::vector<Currency>& inputs,
                         std::vector<Currency>&       outputs);
    //!
    //! Returns true and writes contents of a matrix to a string
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Mat2Str( EvaluatorInterface           eval,
                         const std::vector<Currency>& inputs,
                         std::vector<Currency>&       outputs);
    //!
    //! Returns true after matching/replacing  substrings
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Regexprep( EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs);
    //!
    //! Utility to throw regex error
    //! \param code Error code
    //!
    static void ThrowRegexError( std::regex_constants::error_type code);
    //!
    //! Returns true after converting a numerical value to a string
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Num2Str( EvaluatorInterface           eval,
                         const std::vector<Currency>& inputs,
                         std::vector<Currency>&       outputs);
    //!
    //! Returns true and creates single matrix from string inputs [str2mat
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Str2mat(EvaluatorInterface           eval,
                        const std::vector<Currency>& inputs,
                        std::vector<Currency>&       outputs);

private:
    //!
    //! Helper function for strvcat
    //! \param eval  Evaluator interface
    //! \param in    Input currency
    //! \param index Index of input currency being processed
    //! \param out   Vector of output strings
    //!
    void StrvcatHelper( EvaluatorInterface        eval,
                        const Currency&           input,
                        int                       index,
                        std::vector<std::string>& out);
    //!
    //! Processess cell array for strvcat
    //! \param eval  Evaluator interface
    //! \param cell  Cell array
    //! \param index Index of input currency being processed
    //! \param out   Vector of output strings
    //!
    void StrvcatHelperCellArray( EvaluatorInterface        eval,
                                 HML_CELLARRAY*            cell,
                                 int                       index,
                                 std::vector<std::string>& out);
    //!
    //! Helper function for strvcat, returns string of outputs
    //! \param eval  Evaluator interface
    //! \param in    Input currency
    //! \param index Index of input currency being processed
    //! \param out   Vector of output strings
    //!
    void StrvcatHelperMatrix( EvaluatorInterface eval,
                              const Currency&    in,
                              int                index,
                              std::vector<std::string>& out);
    //!
    //! Returns string from processing a scalar value for strvcat
    //! \param eval  Evaluator interface
    //! \param in    Input currency
    //! \param index Index of input currency being processed
    //!
    std::string StrvcatHelperScalar( EvaluatorInterface eval,
                                     const Currency&    in,
                                     int                index);
    //!
    //! Helper function for strtrim
    //! \param eval Evaluator interface
    //! \param cur  Input currency
    //!
    static Currency StrtrimHelper( EvaluatorInterface& eval, 
                                   const Currency&     cur);
    //!
    //! Helper function to get vector of strings from currency
    //! \param cur Input currency
    //! \param idx Input index for error messages
    //!
    std::vector<std::string> Currency2StringVec( const Currency& cur,
                                                 int             idx);
    //!
    //! Converts double to string
    //! \param cur    Input currency
    //! \param fmtstr String format, if available
    //! \param fmt    Output format
    //!
    std::string Dbl2Str( const Currency&     cur,
                         const std::string&  fmtstr,
                         const OutputFormat* fmt) const;
    //!
    //! Converts double to string
    //! \param cur         Input currency
    //! \param fmtstr      String format, if available
    //! \param fmt         Output format
    //! \param totaldigits Total digits, if available
    //!
    std::string Dbl2Str( const Currency&     cur,
                         const std::string&  fmtstr,
                         const OutputFormat* fmt,
                         int                 totaldigits) const;
    //!
    //! Converts complex to string
    //! \param cur    Input currency
    //! \param fmtstr String format, if available
    //! \param fmt    Output format
    //! \param totaldigits Total digits, if available
    //!
    std::string Complex2Str( const Currency&     cur,
                             const std::string&  fmtstr,
                             const OutputFormat* fmt,
                             int                 totaldigits) const;
};
#endif

// End of file:


