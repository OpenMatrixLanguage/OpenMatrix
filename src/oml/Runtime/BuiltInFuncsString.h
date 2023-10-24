/**
* @file BuiltInFuncsString.h
* @date November 2015
* Copyright (C) 2015-2023 Altair Engineering, Inc.  
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
class OMLDLL_DECLS BuiltInFuncsString
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
    //! Returns true after converting a numerical value to a string
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Num2Str( EvaluatorInterface           eval,
                         const std::vector<Currency>& inputs,
                         std::vector<Currency>&       outputs);
	//!
	//! Returns true and creates a cell array from an input string array [cellstr]
	//! \param eval    Evaluator interface
	//! \param inputs  Vector of inputs
	//! \param outputs Vector of outputs
	//!
	static bool CellStr(EvaluatorInterface           eval,
		const std::vector<Currency>& inputs,
		std::vector<Currency>& outputs);
	//! Returns true and creates a cell array from an input string array [cellstr2]
	//! \param eval    Evaluator interface
	//! \param inputs  Vector of inputs
	//! \param outputs Vector of outputs
	//!
	static bool CellStr2(EvaluatorInterface           eval,
		const std::vector<Currency>& inputs,
		std::vector<Currency>& outputs);
	//!
    //! Returns true and creates single matrix from string inputs [str2mat]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Str2mat(EvaluatorInterface           eval,
                        const std::vector<Currency>& inputs,
                        std::vector<Currency>&       outputs);
    //!
    //! Returns true and converts string to double without eval [str2double]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Str2Double(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs);
    //!
    //! Returns true if given pattern is found in the input [contains]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Contains(EvaluatorInterface           eval,
                         const std::vector<Currency>& inputs,
                         std::vector<Currency>&       outputs);
    //!
    //! Returns true and strips leading/trailing characters from input [strip]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Strip(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns true successful in converting string to scalar/complex
    //! \param in       Input string
    //! \param rval     Real part
    //! \param ival     Imaginary part, if it exists
    //! \param isscalar True if value is scalar
    //!
    bool Str2Num(const std::string& in,
                 double& rval,
                 double& ival,
                 bool& isscalar);
    //!
    //! Returns true if successful in converting string to scalar/complex
    //! \param Input string
    //! \param Output currency
    //!
    static bool IsNumber(const std::string&,
                         Currency&);
    //!
    //! Returns true after doing a case sensitive comparison of strings [strcmp]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //!
    static bool Strcmp(EvaluatorInterface,
                       const std::vector<Currency>&,
                       std::vector<Currency>&);
    //!
    //! Returns true after doing a case insensitive comparison of strings [strcmpi]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //!
    static bool Strcmpi(EvaluatorInterface,
                        const std::vector<Currency>&,
                        std::vector<Currency>&);
    //!
    //! Returns a matrix indicating which elements are printable [isprint]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //!
    static bool IsPrint(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Helper method to split an input string into substrings
    //! \param Input string
    //! \param Delimiter
    //! 
    static std::vector<std::string> Split(const std::string&, const std::string&);

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
    //!
    //! Returns true if strtod conversion is successful
    //! \param in  Input string
    //! \param end End string
    //! \param val Value converted by strtod
    //!
    bool IsValidStrtodResult(const std::string& in,
                             const std::string& end,
                             double             val);
    //!
    //! Helper method to left trim input
    //! \param Input string
    //! \param Pattern(s) to trim
    //!
    void LeftTrim(std::string&, 
                  const std::vector<std::string>&);
    //!
    //! Helper method to right trim input
    //! \param Input string
    //! \param Pattern(s) to trim
    //!
    void RightTrim(std::string&       in,
                   const std::vector<std::string>& trim);

    //!
    //! Utility function to convert matrix to a lower case string
    //! \param Input matrix
    //!
    static Currency ToLowerString(const hwMatrix*);
    //!
    //! Returns true if this is a valid numeric format
    //! \param Format for sprintf
    //! 
    static bool IsValidNumericFormat(const std::string&);

    static Currency IsPrintImpl(const Currency&);
    //!
    //! Parse input string to get formats
    //! \param Input
    //! \param Base format
    //! \param Full format
    //! \param Skip format flag
    //! \param True if warning needs to be shown
    //! 
    static bool ParseFormat(const std::string&, std::vector<std::string>&, std::vector<std::string>&, std::vector<bool>&, bool&);
    //!
    //! Returns true after scanning double data
    //! \param Input string
    //! \param Format
    //! \param True if value is skipped
    //! \param Outputs
    //! \param String read
    //! \param Number of chars read
    static bool SscanfFloat(const std::string&,
        const std::string& fmt,
        bool usefmt,
        std::vector<Currency>& outvals,
        std::string& stringread,
        int& numread);
    //!
    //! sscanf helper function, returns true after reading formatted string
    //! \param Input string
    //! \param Format
    //! \param Base format
    //! \param True if value is skipped
    //! \param Output values
    //! \param True if warning about invalid formats needs to be shown
    //! 
    static bool Sscanf(std::string&, const std::string&, const std::string&, bool, std::vector<Currency>&);
};
#endif

// End of file:


