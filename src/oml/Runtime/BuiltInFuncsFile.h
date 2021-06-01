/**
* @file BuiltInFuncsFile.h
* @date March 2016
* Copyright (C) 2016-2021 Altair Engineering, Inc.  
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

#ifndef __BUILTINFUNCSFILE__
#define __BUILTINFUNCSFILE__

#include "OMLDll.h"
#include "EvaluatorInt.h"

//------------------------------------------------------------------------------
//!
//! \class BuiltInFuncsFile
//! \brief Utility class for built-in functions related to file I/O commands
//!
//------------------------------------------------------------------------------
class OMLDLL_DECLS BuiltInFuncsFile
{
public:
    //!
    //! Destructor
    //!
    ~BuiltInFuncsFile() {}
    //!
    //! Returns true after copying files/directories (copyfile command)
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //! 
    static bool Copyfile(EvaluatorInterface           eval,
	                     const std::vector<Currency>& inputs,
		                 std::vector<Currency>&       outputs);
    //!
    //! Returns true after writing a matrix to a file (dlmwrite command)
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //! 
    static bool Dlmwrite(EvaluatorInterface           eval,
	                     const std::vector<Currency>& inputs,
		                 std::vector<Currency>&       outputs);
    //!
    //! Returns true after importing data from files (importdata command)
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //! 
    static bool Importdata(EvaluatorInterface           eval,
	                       const std::vector<Currency>& inputs,
		                   std::vector<Currency>&       outputs);
    //!
    //! Returns true after reading a text file (textread command)
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //! 
    static bool Textread(EvaluatorInterface           eval,
	                     const std::vector<Currency>& inputs, 
		                 std::vector<Currency>&       outputs);
    //!
    //! Returns true after reading a text file/string (textscan command)
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //! 
    static bool Textscan(EvaluatorInterface           eval,
	                     const std::vector<Currency>& inputs, 
		                 std::vector<Currency>&       outputs);
    //!
    //! Returns true after reading formatted input from file/string (fscanf)
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //! 
    static bool Fscanf(EvaluatorInterface           eval,
	                   const std::vector<Currency>& inputs, 
		               std::vector<Currency>&       outputs);
    //!
    //! Returns true after renaming a file (rename)
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //! 
    static bool Rename(EvaluatorInterface           eval,
	                   const std::vector<Currency>& inputs, 
		               std::vector<Currency>&       outputs);
    //!
    //! Returns true if given file extension is Excel compatible
    //! \param ext File extension
    //! 
    static bool IsExtXlsCompatible(const std::string& ext);
    //!
    //! Returns true if given filename is an Ascii file
    //! \param name File name
    //! 
    static bool IsAsciiFile(const std::string& name);
    //!
    //! Returns true and reads float/integer data from a file (dlmread)
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //! 
    static bool Dlmread(EvaluatorInterface           eval,
	                    const std::vector<Currency>& inputs, 
		                std::vector<Currency>&       outputs);
    //!
    //! Returns true and opens a file for read/write [fopen]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //! 
    static bool Fopen(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns true and reads a comma separated file [csvread]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Csvread(EvaluatorInterface,
                        const std::vector<Currency>&,
                        std::vector<Currency>&);
    //!
    //! Returns true and writes a comma separated file [csvwrite]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Csvwrite(EvaluatorInterface,
                         const std::vector<Currency>&,
                         std::vector<Currency>&);

private:
    //!
    //! Constructor
    //!
    BuiltInFuncsFile() {}
    //!
    //! Stubbed out copy constructor
    //!
	BuiltInFuncsFile(const BuiltInFuncsFile&);
    //!
    //! Stubbed out assignment operator
    //!
    BuiltInFuncsFile& operator=(const BuiltInFuncsFile&);

    //!
    //! Returns string value of a currency, if applicable
    //! \param inputs Vector of input currencies
    //! \param nargin Number of inputs
    //! \param index  Updated currency index
    //!
    std::string GetStringValue(const std::vector<Currency>& inputs, 
                               int                          nargin,
                               int&                         index) const;
    //!
    //! Returns integer value of a currency, if applicable
    //! \param inputs Vector of input currencies
    //! \param nargin Number of inputs
    //! \param index  Updated currency index
    //!
    int GetIntegerValue(const std::vector<Currency>& inputs, 
                        int                          nargin,
                        int&                         index) const;
    //!
    //! Returns result and scans number from a file
    //! \param file       File
    //! \param fmt        Format
    //! \param iscomplex  True if this is a complex number
    //! \param realval    Real value
    //! \param complexval Complex value, if applicable
    //!
    int ReadDouble(FILE*              file, 
                   const std::string& fmt,
                   bool               isdouble,
                   bool&              iscomplex,
                   double&            realval,
                   double&            complexval);
    //!
    //! Returns true if this is a complex number along with the value
    //! \param in  Input string
    //! \param val Value
    //!
    bool IsComplex(const std::string& in,
                   double&            val);
    //!
    //! Utility returning true if string represents complex value
    //! \param in  Input string
    //! \param rval Real value
    //! \param ival Imaginary value
    //!
    static bool IsComplex(const std::string& in,
                          double&            rval,
                          double&            ival);

    //!
    //! Returns a vector range, given a string range in excel format
    //! \param strRange Given range in the form A1..B5 or A1:B5
    //! \param idx      Index of input argument for error 
    //!
    std::vector<int> String2VectorRange(const std::string& strRange,
                                        int                idx);
    //!
    //! Parses given string to get row/column information
    //! \param str Given string in the form A1
    //! \param err Error to throw if index is incorrect
    //! \param row Row - 0 based
    //! \param col Column - 0 based
    //!
    void GetRowColInfo(const std::string& str,
                       const std::string& err,
                       int&               row,
                       int&               col);
    //!
    //! Utility to get format
    //! \param Input
    //! \param Vector of base formats
    //! \param Vector of raw  formats
    //! \param True if specific format can be used
    //! \param Error
    //! \param Prefix
    //!
    bool GetFormats(const std::string&,
                    std::vector<std::string>&,
                    std::vector<std::string>&,
                    std::vector<bool>&,
                    std::string&,
                    std::string&);
    //!
    //! Sets default delimiter for certain file types
    //! \param File name
    //! \param Delimiter
    //!
    static void SetDefaultDelimiter(const std::string&,
                                    std::string&);
    //!
    //! Escapes delimiters
    //! \param Current delimiter(s)
    //! \param True if there is a newline delimiter
    //!
    static void EscapeDelimiters(std::string&,
                                 bool&);
    //!
    //! Returns a line read from the given file, with the given end of line character
    //! \param File pointer
    //! \param End of line character
    //!
   static std::string GetLine(std::FILE*,
                              const std::string&);
}; 

#endif  // __BUILTINFUNCSFILE__


