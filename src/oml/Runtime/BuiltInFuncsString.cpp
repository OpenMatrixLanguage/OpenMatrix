/**
* @file BuiltInFuncsString.cpp
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

// Begin defines/includes

#include "BuiltInFuncsString.h"

#include <cassert>
#include <iomanip>
#include <memory>     // For std::unique_ptr
#include <sstream>
#include <string>

#include "BuiltInFuncsUtils.h"
#include "CurrencyDisplay.h"
#include "ErrorInfo.h"
#include "MatrixDisplay.h"
#include "OML_Error.h"
#include "OutputFormat.h"

#include "hwComplex.h"
#include "hwMatrix.h"

// End defines/includes

//------------------------------------------------------------------------------
// Returns true and returns a string of blanks
//------------------------------------------------------------------------------
bool BuiltInFuncsString::hml_strvcat(EvaluatorInterface           eval,
                                     const std::vector<Currency>& inputs,
                                     std::vector<Currency>&       outputs)
{
    if (inputs.empty()) throw OML_Error(OML_ERR_NUMARGIN);

    BuiltInFuncsString funcs;

    // Creates character array from one or more matrices, cell arrays or strings
    // with one row by concatinating across all columns
    std::vector<std::string> out;
    std::vector<Currency>::const_iterator itr = inputs.begin();
    for (int index = 0; itr != inputs.end(); ++itr, ++index)
        funcs.StrvcatHelper(eval, *itr, index, out);

    Currency result = BuiltInFuncsUtils::FormatOutput(
                      out, false, false, false, ' ');
    outputs.push_back(result);

    return true;
}
//------------------------------------------------------------------------------
// Helper function for strvcat
//------------------------------------------------------------------------------
void BuiltInFuncsString::StrvcatHelper(EvaluatorInterface        eval,
                                       const Currency&           in,
                                       int                       index,
                                       std::vector<std::string>& out)
{
    if (in.IsNothing()) return;

    if (in.IsComplex()) throw OML_Error(OML_ERR_NOCOMPLEX, index+1, OML_VAR_VALUE);

    else if (in.IsScalar())
        out.push_back(StrvcatHelperScalar(eval, in, index));
        
    else if (in.IsString())
    {
        std::string val (in.StringVal());
        if (!val.empty())
            out.push_back(in.StringVal());       
    }

    else if (in.IsMatrix())
        StrvcatHelperMatrix(eval, in, index, out);

    else if (in.IsCellArray())
        StrvcatHelperCellArray(eval, in.CellArray(), index, out);
}
//------------------------------------------------------------------------------
// Returns string from processing a scalar value for strvcat
//------------------------------------------------------------------------------
std::string BuiltInFuncsString::StrvcatHelperScalar(EvaluatorInterface eval, 
                                                    const Currency&    in,
                                                    int                index)
{
    if (in.IsPositiveInteger())
    {
        int val = static_cast<int>(in.Scalar());
        std::ostringstream os;
        os << static_cast<char>(val);
        std::string result (os.str());
        return result;
    }

    if (!in.IsInteger()) return "";  // Just return an empty string

    if (in.Scalar() < 0)
    {
        std::ostringstream os;
        os << "Warning: invalid input in argument " << index + 1
           << "; input must be a nonnegative integer";
        BuiltInFuncsUtils::SetWarning(eval, os.str());
    }
    return "";  // Just return an empty string
}
//------------------------------------------------------------------------------
// Processess cell array for strvcat
//------------------------------------------------------------------------------
void BuiltInFuncsString::StrvcatHelperCellArray(EvaluatorInterface        eval,
                                                HML_CELLARRAY*            cell,
                                                int                       index,
                                                std::vector<std::string>& out)
{
    int numelem = cell ? cell->Size() : 0;

    std::vector <std::string> outstr;
    if (numelem == 0) return; // Should never get into this situation

    for (int i = 0; i < numelem; ++i)
    {
        Currency cur ((*cell)(i));
        StrvcatHelper(eval, cur, index, out);
    }
}
//------------------------------------------------------------------------------
// Processess matrix for strvcat
//------------------------------------------------------------------------------
void BuiltInFuncsString::StrvcatHelperMatrix(EvaluatorInterface        eval,
                                             const Currency&           in,
                                             int                       index,
                                             std::vector<std::string>& out)
{
    assert(in.IsMatrix());

    const hwMatrix* mtx = in.Matrix();
    if (!mtx) 
    {
        out.push_back("");
        return;
    }
    if (!mtx->IsReal()) throw OML_Error(OML_ERR_NOCOMPLEX, index+1, OML_VAR_VALUE);

    int rows = mtx->M();
    int cols = mtx->N();

    for (int i = 0; i < rows; ++i)
    {
        std::string rowstr("");

        for (int j = 0; j < cols; ++j)
        {
            Currency tmp ((*mtx)(i, j));
            if (tmp.IsScalar())
                rowstr += StrvcatHelperScalar(eval, tmp, index);
        }
        out.push_back (rowstr);
    }
}
//------------------------------------------------------------------------------
// Returns true and tokenizes given string
//------------------------------------------------------------------------------
bool BuiltInFuncsString::hml_strtok(EvaluatorInterface           eval,
                                    const std::vector<Currency>& inputs,
                                    std::vector<Currency>&       outputs)

{
    size_t nargin = (!inputs.empty()) ? inputs.size() : 0;
    if (nargin < 1 || nargin > 2) throw OML_Error(OML_ERR_NUMARGIN);

    // Check if first input is valid
    const Currency& input1 = inputs[0];
    if (!input1.IsString()) throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);

    const hwMatrix* strmtx = input1.Matrix();
    if (!strmtx || !strmtx->IsVector()) 
        throw OML_Error(OML_ERR_STRING_ONEDIMENSION, 1, OML_VAR_TYPE);

    bool horiz = (strmtx->M() == 1); // Determines dimensions of outputs
    
    // Get the delimiter
    std::string delim;
    if (nargin > 1)
    {
        const Currency &input2 = inputs[1];
        if (input2.IsString())
            delim = BuiltInFuncsUtils::GetOrderedStringVal(input2);

        else if (input2.IsScalar())
            delim = BuiltInFuncsUtils::GetValidChar(eval, input2.Scalar());
        
        else if (input2.IsComplex())
        {
            outputs.push_back(input1);
            outputs.push_back(std::string());
        }

        else if (input2.IsMatrix())
        {
            const hwMatrix *m = input2.Matrix();
            assert(m);
            int numelem = m->Size();
            bool isreal = m->IsReal();

            for (int i = 0; i < numelem; ++i)
            {
                if (isreal)
                    delim += BuiltInFuncsUtils::GetValidChar(eval, (*m)(i));
                else
                {
                    hwComplex cplx = m->z(i);
                    if (cplx.IsReal(1.0e-15))
                        delim += BuiltInFuncsUtils::GetValidChar(eval, cplx.Real());
                }
            }
        }
        else throw OML_Error(OML_ERR_STRSCALARCOMPLEXMTX, 2, OML_VAR_TYPE);
    }

    if (delim.empty())
    {
        delim = " ";
    }
    else // \t, \r, \n are interpreted as 2 characters each in oml
    {
        size_t pos = delim.find("\\t"); 
        if (pos != std::string::npos)
        {
            delim.replace(pos, 2, "\t");
        }

        pos = delim.find("\\n"); 
        if (pos != std::string::npos)
        {
            delim.replace(pos, 2, "\n");
        }
        pos = delim.find("\\r"); 
        if (pos != std::string::npos)
        {
            delim.replace(pos, 2, "\r");
        }
    }
    
    std::string str (BuiltInFuncsUtils::GetOrderedStringVal(input1));
 
    // Strip leading delimiters
    while (1)
    {
        if (str.empty()) break;

        size_t idx = str.find_first_of(delim);

        if (idx != 0) break;

        str.erase(str.begin());
    }

    Currency tmp(str);
    const hwMatrix* mtx = tmp.Matrix();
        
    size_t i = str.find_first_of(delim);
    if (i == std::string::npos)
    {
        outputs.push_back(tmp);
        outputs.push_back(std::string());
        return true;
    }

    int len     = static_cast<int>(str.length());
    int rowsTok = 1;
    int rowsRem = 1;
    int colsTok = static_cast<int>(i);
    int colsRem = len - static_cast<int>(i);
    if (!horiz)
    {
        rowsTok = static_cast<int>(i);
        rowsRem = len - static_cast<int>(i);
        colsTok = 1;
        colsRem = 1;
    }

    hwMatrix* tok = EvaluatorInterface::allocateMatrix(rowsTok, colsTok, hwMatrix::REAL);
    hwMatrix* rem = EvaluatorInterface::allocateMatrix(rowsRem, colsRem, hwMatrix::REAL);   

    for (int j = 0; j < i; ++j)
        (*tok)(j) = (double) (*mtx)(j);

    for (int j = static_cast<int>(i); j < mtx->Size(); ++j)
        (*rem)(j - (const int)i) = (double) (*mtx)(j);

    Currency out1(tok);
    out1.SetMask(Currency::MASK_STRING);

    Currency out2(rem);
    out2.SetMask(Currency::MASK_STRING);

    outputs.push_back(out1);
    outputs.push_back(out2);
    return true;
}
//------------------------------------------------------------------------------
// Returns true and returns a string of blanks
//------------------------------------------------------------------------------
bool BuiltInFuncsString::hml_blanks(EvaluatorInterface           eval,
                                    const std::vector<Currency>& inputs,
                                    std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1) throw OML_Error(OML_ERR_NUMARGIN);

    // Check if input is valid
    const Currency& in = inputs[0];

    if (!in.IsInteger()) throw OML_Error(OML_ERR_FINITE_NATURALNUM);

    double val = in.Scalar();
    if (val < 0 || IsInf_T(val) || IsNaN_T(val) || IsNegInf_T (val))
        throw OML_Error(OML_ERR_FINITE_NATURALNUM);

    hwMatrix* mtx = 0;
    int n = static_cast<int>(val);
    std::string str;
    if (n != 0) 
        str.append(static_cast<size_t>(n), ' ');

    outputs.push_back(str);
    return true;
}
//------------------------------------------------------------------------------
// Returns true and reads formatted input from a string
//------------------------------------------------------------------------------
bool BuiltInFuncsString::Sscanf(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
{
    int numargin = inputs.empty() ? 0 : static_cast<int>(inputs.size());
    if (numargin < 2) throw OML_Error(OML_ERR_NUMARGIN);

    // First argument gives the input string to be scanned
    const Currency& cur1 = inputs[0];
    if (!cur1.IsString()) throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);

    std::string in = cur1.StringVal();
    if (in.empty())
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
        return true;
    }

    // Second argument gives format descriptions separated by spaces. 
    const Currency& cur2 = inputs[1];
    if (!cur2.IsString()) throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);

    // Third argument, if available, specifies size for output.
    double rows        = 0;
    double cols        = 0;
    bool   hasSizeMtx  = false;
    bool hasSizeSpec = (numargin > 2);
	if (numargin > 2)
        BuiltInFuncsUtils::GetSizeSpecifications(inputs[2], 3, rows, cols, hasSizeMtx);
        
    // Read formatted input from the string
    bool        validformat = true;
    std::string fmtdesc   (cur2.StringVal());
    std::string validfmts ("%f%g%d%s");
    Currency output  = BuiltInFuncsUtils::GetFormattedInput(in, fmtdesc,
                       validfmts, rows, cols, hasSizeMtx, hasSizeSpec, validformat);
    if (!validformat)
    {
        std::ostringstream os;
        os << "Warning: invalid format specified in argument " << 2
           << "; valid formats are %d, %f, %g and %s";
        BuiltInFuncsUtils::SetWarning(eval, os.str());
    }

	outputs.push_back(output);
    return true;
}
//------------------------------------------------------------------------------
// Returns true if successful in converting string to ascii row matrix
//------------------------------------------------------------------------------
bool BuiltInFuncsString::ToAscii(EvaluatorInterface           eval,
                                 const std::vector<Currency>& inputs,
                                 std::vector<Currency>&       outputs)
{
    if (inputs.empty()) throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input = inputs[0];  // Value to convert
    if (!input.IsString()) throw OML_Error(OML_ERR_STRING, 1, OML_VAR_VALUE);

    std::string strval (input.StringVal());
    if (strval.empty())
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
        return true;
    }

	const hwMatrix* in_mat = input.Matrix();
    assert(in_mat);
    int m = in_mat->M();
    int n = in_mat->N();

    hwMatrix* out = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);

    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < m; ++i)
            (*out)(i, j) = (*in_mat)(i, j);
    }
    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// Returns true and trims leading/trailing spaces from strings/cells
//------------------------------------------------------------------------------
bool BuiltInFuncsString::Strtrim(EvaluatorInterface           eval,
                                 const std::vector<Currency>& inputs,
                                 std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1) throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(StrtrimHelper(eval, inputs[0]));
    return true;
}
//------------------------------------------------------------------------------
// Helper function for strtrim
//------------------------------------------------------------------------------
Currency BuiltInFuncsString::StrtrimHelper(EvaluatorInterface& eval, 
                                           const Currency&     cur)
{
    bool iscell   = cur.IsCellArray();
    bool isstring = cur.IsString();
    if (!iscell && !isstring) throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);

    if (iscell)
    {
        HML_CELLARRAY* cell     = cur.CellArray();
        HML_CELLARRAY* result   = EvaluatorInterface::allocateCellArray(cell->M(), cell->N());
        int            cellsize = cell->Size();
        for (int i = 0; i < cellsize; ++i)
            (*result)(i) = BuiltInFuncsString::StrtrimHelper(eval, (*cell)(i));
        return Currency(result);
    }
    
    // Process strings
    const hwMatrix* str = cur.Matrix();
    int minpos = str->N() - 1;
    int maxpos = 0;

    for (int i = 0; i < str->M(); ++i)
    {
        bool encountered = false;
        for (int j = 0; j < str->N(); ++j)
        {
            if (!isspace((int) (*str)(i,j)))
            {
                maxpos = std::max(maxpos, j + 1);
                if (!encountered)
                {
                    encountered = true;
                    minpos = std::min(minpos, j);
                }
            }
        }
    }

    hwMatrix* result = NULL;
    if (!maxpos)
        result = EvaluatorInterface::allocateMatrix();
    else
    {
        int numrows = str->M();
        result = EvaluatorInterface::allocateMatrix(numrows, maxpos - minpos, hwMatrix::REAL);
        for (int i = 0; i < numrows; ++i)
        {
            for (int j = minpos; j < maxpos; ++j)
                (*result)(i,j - minpos) = (*str)(i,j);
        }
    }

    Currency out(result);
    out.SetMask(Currency::MASK_STRING);
    return out;
}
//------------------------------------------------------------------------------
// Returns true and writes contents of a matrix to a string
//------------------------------------------------------------------------------
bool BuiltInFuncsString::Mat2Str(EvaluatorInterface           eval,
                                 const std::vector<Currency>& inputs,
                                 std::vector<Currency>&       outputs)
{
    int nargin = inputs.empty() ? 0 : static_cast<int>(inputs.size());
    if (nargin < 1) throw OML_Error(OML_ERR_NUMARGIN);
	
	if (nargin > 3) throw OML_Error(OML_ERR_NUMARGIN);

    std::string output ("[");
    const Currency& cur = inputs[0];

    // Check for matrix/scalar/complex as 1x1 matrices are scalars or complex
    bool isscalar  = cur.IsScalar();
    bool iscomplex = cur.IsComplex();
    if (!cur.IsMatrix() && !isscalar && !iscomplex) 
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_TYPE);

    std::string precreal;
    std::string precimag;
    bool        printclass = false;
    if (nargin >= 2)
    {
        const Currency& tmp = inputs[1];
        if (!tmp.IsScalar() && !tmp.IsVector() && !tmp.IsString())
            throw OML_Error(OML_ERR_SCALAR_VECTOR_STRING, 2, OML_VAR_TYPE);

        if (tmp.IsScalar())
        {
            precreal = "%";
            double val = tmp.Scalar();
            if (tmp.IsInteger())
            {
                precreal += "." + std::to_string(static_cast<long long>(val)) + "g";
            }
            else
            {
                std::string strval (tmp.GetValues(eval.GetOutputFormat())); 
                if (strval.find("e") != std::string::npos || 
                    strval.find("E") != std::string::npos)
                {
                    std::ostringstream os;
                    os << std::fixed << val;
                    std::string tmp (os.str());
                    size_t pos = tmp.find(".");
                    if (pos != std::string::npos)
                    {
                        std::string tmp1 = tmp.substr(0, pos);
                        std::string tmp2 = tmp.substr(pos);
                        double decval = atoi(tmp.c_str());
                        decval = ceil(decval);
                        os.str("");
                        os.clear();
                        os << std::fixed << std::setprecision(1) 
                           << tmp1 << decval;
                    }
                    precreal += os.str() + "g";     
                }            
                else
                {
#ifdef OS_WIN
                    precreal += std::to_string(static_cast<long double>(val)) + "g";
#else
                    char buff[1024];
                    sprintf(buff, "%.12g", val);
                    precreal += buff + std::string("g");
#endif
                }
            }
            precimag = precreal;  // Use the same for real/imag
        }
        else if (tmp.IsVector())
        {
            std::vector<double> vec (tmp.Vector());
            if (vec.empty() || vec.size() < 2)
                throw OML_Error(OML_ERR_VECTOR2, 2, OML_VAR_TYPE);
            double rval = vec[0];
            double ival = vec[1];

            precreal = "%";
            if (BuiltInFuncsUtils::IsInt(rval))
                precreal += "." + std::to_string(static_cast<long long>(rval)) + "g";
            else
            {
                Currency curRval (rval);
                std::string strval (curRval.GetValues(eval.GetOutputFormat())); 
                if (strval.find("e") != std::string::npos || 
                    strval.find("E") != std::string::npos)
                {
                    std::ostringstream os;
                    os << std::fixed << rval;
                    std::string tmp (os.str());
                    size_t pos = tmp.find(".");
                    if (pos != std::string::npos)
                    {
                        std::string tmp1 = tmp.substr(0, pos);
                        std::string tmp2 = tmp.substr(pos);
                        double decval = atoi(tmp.c_str());
                        decval = ceil(decval);
                        os.str("");
                        os.clear();
                        os << std::fixed << std::setprecision(1) 
                           << tmp1 << decval;
                    }
                }            
                else
                {
#ifdef OS_WIN
                    precreal += std::to_string(static_cast<long double>(rval)) + "g";
#else
                    char buff[1024];
                    sprintf(buff, "%.12g", rval);
                    precreal += buff + std::string("g");
#endif
                }
            }

            precimag = "%";
            if (BuiltInFuncsUtils::IsInt(ival))
                precimag += "." + std::to_string(static_cast<long long>(ival)) + "g";
            else
            {
                Currency curIval (ival);
                std::string strval (curIval.GetValues(eval.GetOutputFormat())); 
                if (strval.find("e") != std::string::npos || 
                    strval.find("E") != std::string::npos)
                {
                    std::ostringstream os;
                    os << std::fixed << ival;
                    std::string tmp (os.str());
                    size_t pos = tmp.find(".");
                    if (pos != std::string::npos)
                    {
                        std::string tmp1 = tmp.substr(0, pos);
                        std::string tmp2 = tmp.substr(pos);
                        double decval = atoi(tmp.c_str());
                        decval = ceil(decval);
                        os.str("");
                        os.clear();
                        os << std::fixed << std::setprecision(1) 
                           << tmp1 << decval;
                    }
                    precimag += os.str() + "g";
                }            
                else
                {
#ifdef OS_WIN
                    precimag += std::to_string(static_cast<long double>(ival)) + "g";
#else
                    char buff[1024];
                    sprintf(buff, "%.12g", ival);
                    precimag += buff + std::string("g");
#endif
                }
            }
        }
        else if (tmp.IsString())
            printclass = (tmp.StringVal() == "class");
    }

    output += MatrixDisplay::GetOutputValues(cur, 
        eval.GetOutputFormat(), ";", " ", precreal, precimag, -1);

    assert(!output.empty());
    size_t len = output.size();
    if (output[len - 1] == ';')
        output.erase(len - 1);
    output += "]";

    if (!printclass && nargin >= 3)
    {
        const Currency& tmp = inputs[2];
        printclass = (tmp.IsString() && tmp.StringVal() == "class");
    }
    
    if (!printclass)
        outputs.push_back(output);
    else
        outputs.push_back("double(" + output + ")");

    return true;
}
//------------------------------------------------------------------------------
// Returns true after matching/replacing  substrings [regexprep]
//------------------------------------------------------------------------------
bool BuiltInFuncsString::Regexprep(EvaluatorInterface           eval,
                                   const std::vector<Currency>& inputs,
                                   std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.empty() ? 0 : inputs.size();
    if (nargin < 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    BuiltInFuncsString funcs;
    std::vector<std::string> srcvec (funcs.Currency2StringVec(inputs[0], 1));
    std::vector<std::string> patvec (funcs.Currency2StringVec(inputs[1], 2));
    std::vector<std::string> repvec (funcs.Currency2StringVec(inputs[2], 3));
    std::vector<std::string> opts;

    std::regex_constants::match_flag_type replaceflag =
        std::regex_constants::match_default;
    std::regex_constants::syntax_option_type syntaxflag =
        std::regex_constants::ECMAScript;

    // Get other options, if any
    bool emptymatch = false;
    bool ignorecase = false;
    for (int i = 3; i < nargin; ++i)
    {
        if (!inputs[i].IsString())
        {
            throw OML_Error(OML_ERR_STRING, i + 1, OML_VAR_TYPE);
        }
        std::string val (inputs[i].StringVal());
        if (val.empty())
        {
            continue;
        }
            
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);
        if (val == "once")
        {
            replaceflag = std::regex_constants::format_first_only;
        }
        else if (val == "ignorecase")
        {
            syntaxflag = std::regex_constants::ECMAScript | std::regex_constants::icase;
            ignorecase = true;
        }
        else if (val == "emptymatch")
        {
            emptymatch = true;
        }
    }

    if (!emptymatch)
    {
        replaceflag |= std::regex_constants::match_not_null;
    }
    if (srcvec.empty() || patvec.empty() || repvec.empty())
    {
        outputs.push_back(inputs[0]);
        return true;
    }

    std::vector<std::string> out;
    out.reserve(srcvec.size());

    size_t repidx = 0;
    size_t numrep = repvec.size();

    
    {
        for (std::vector<std::string>::const_iterator itr1 = srcvec.begin();
             itr1 != srcvec.end(); ++itr1)
        {
            std::string txt (*itr1);
            if (txt.empty())
            {
                out.push_back("");
                continue;
            }

            for (std::vector<std::string>::const_iterator itr2 = patvec.begin();
                itr2 != patvec.end(); ++itr2)
            {
                std::string pat(*itr2);
                if (pat.empty())
                {
                    continue;
                }

#ifndef OS_WIN
                // * in the pattern causes a zero search. So replace
                // with + which includes the rest of the characters
                size_t len = pat.size();
                std::string tmp;
                for (size_t k = 0; k < len; ++k)
                {
                    char ch = pat[k];
                    if (ch != '*' ||
                        (ch == '*' && ((k > 0 && pat[k - 1] == '\\') || emptymatch)))
                    {
                        continue;
                    }
                    pat[k] = '+';  // Substitute just '*' for '+'
                }

                if (ignorecase)
                {
                    // option 'ignorecase' does not seem to work on Linux
                    len = pat.find_last_of(']');
                    if (len != std::string::npos)
                    {
                        pat.insert(len, "(?-i)");
                    }
                }
#endif
                std::string rep (repvec[repidx++]);
                if (repidx >= numrep)
                {
                    repidx = 0;
                }
                std::string result;
                try
                {
                    std::regex regexpat(pat, syntaxflag);  
                    result = std::regex_replace(txt, regexpat, rep, replaceflag);
                }
                catch (std::regex_error& err)
                {
                    ThrowRegexError(err.code());
                }
                catch (...)
                {
                    throw OML_Error("Error: cannot parse regular expression in argument 2");
                }
                txt = result;
            }
            out.push_back(txt);
        }
    }

    if (out.empty())
    {
        outputs.push_back(inputs[0]);
        return true;
    }

    size_t numout = out.size();
    if (numout == 1 || inputs[0].IsString())
    {
        outputs.push_back(out[0]);
        return true;
    }

    HML_CELLARRAY* cell = EvaluatorInterface::allocateCellArray((int)numout, 1);
    assert(cell);
    for (int i = 0; i < numout; ++i)
    {
        Currency tmp(out[i]);
        tmp.SetMask(Currency::MASK_STRING);
        (*cell)(i, 0) = tmp;
    }
    outputs.push_back(cell);
    return true;
}
//------------------------------------------------------------------------------
// Helper function to get vector of strings from currency
//------------------------------------------------------------------------------
std::vector<std::string> BuiltInFuncsString::Currency2StringVec(const Currency& cur, int idx)
{
    if (!cur.IsString() && !cur.IsCellArray())
        throw OML_Error(OML_ERR_STRING_STRINGCELL, idx, OML_VAR_TYPE);

    std::vector<std::string> out;
    if (cur.IsString())
    {
        out.push_back(cur.StringVal());
        return out;
    }

    HML_CELLARRAY *cell = cur.CellArray();
    assert(cell);
    int numelem = cell->Size();
    out.reserve(numelem);
    int numrows = cell->M();
    int numcols = cell->N();

    for (int j = 0; j < numcols; ++j)
        for (int i = 0; i < numrows; ++i)
    {
        Currency tmp = (*cell)(i, j);
        if (!tmp.IsString())
            throw OML_Error(OML_ERR_STRING_STRINGCELL, idx, OML_VAR_VALUE);
        out.push_back(tmp.StringVal());
    }
    return out;
}
//------------------------------------------------------------------------------
// Throws regex error
//------------------------------------------------------------------------------
void BuiltInFuncsString::ThrowRegexError(std::regex_constants::error_type code)
{
    std::string errmsg;
    switch (code)
    {
        case std::regex_constants::error_collate:
            errmsg = "invalid collating element name";
            break;
        case std::regex_constants::error_ctype:
            errmsg = "invalid character class name";
            break;
        case std::regex_constants::error_escape:
            errmsg = "invalid escape sequence";
            break;
        case std::regex_constants::error_backref:
            errmsg = "invalid back reference";
            break;
        case std::regex_constants::error_brack:
#ifdef LINUX 
            // Square brackets[] is not supported with gcc 4.8
            errmsg = "invalid square brackets - not supported in regular expressions";
#else
            errmsg = "invalid unmatched square brackets";
#endif
            break;
        case std::regex_constants::error_paren:
            errmsg = "invalid unmatched parenthesis";
            break;
        case std::regex_constants::error_brace:
            errmsg = "invalid unmatched curly braces";
            break;
        case std::regex_constants::error_badbrace:
            errmsg = "invalid range in braces";
            break;
        case std::regex_constants::error_range:
            errmsg = "invalid character range";
            break;
        case std::regex_constants::error_space:
        case std::regex_constants::error_stack:
            errmsg = "not enough memory";
            break;
        case std::regex_constants::error_badrepeat:
            errmsg = "invalid section to repeat";
            break;
        case std::regex_constants::error_complexity:
            errmsg = "match was too complex";
            break;
#ifdef OS_WIN
        case std::regex_constants::error_syntax:
            errmsg = "syntax error with regular expression";
            break;
#endif
        default:
            errmsg = "cannot parse regular expression";
            break;
    }

    throw OML_Error("Error: " + errmsg);
}
//------------------------------------------------------------------------------
// Returns true after converting a numerical value to a string
//------------------------------------------------------------------------------
bool BuiltInFuncsString::Num2Str(EvaluatorInterface           eval,
                                 const std::vector<Currency>& inputs,
                                 std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    Currency cur = inputs[0];
    if (!cur.IsMatrixOrString() && !cur.IsScalar() && !cur.IsComplex())
    {
        throw OML_Error(OML_ERR_STRSCALARCOMPLEXMTX, 1, OML_VAR_TYPE);
    }

    if (cur.IsString())
    {
        outputs.push_back(cur);
        return true;
    }

    int precision = -1;
    std::string fmt;

    if (inputs.size() > 1)
    {
        const Currency& tmp = inputs[1];
        if (!tmp.IsInteger() && !tmp.IsString())
        {
            throw OML_Error(OML_ERR_STRING_NATURALNUM, 2, OML_VAR_TYPE);
        }

        if (tmp.IsString())
        {
            fmt = tmp.StringVal();
        }
        else
        {
            precision = static_cast<int>(tmp.Scalar());
            if (precision < 0)
            {
                throw OML_Error(OML_ERR_STRING_NATURALNUM, 2, OML_VAR_VALUE);
            }
        }
    }

    OutputFormat ofmt = (precision >= 0) ? 
                        OutputFormat(1, precision) : *(eval.GetOutputFormat());
         
    BuiltInFuncsString funcs;
    cur.DispOutput();
    if (cur.IsScalar())
    {
        std::string val = funcs.Dbl2Str(cur, fmt, &ofmt, precision);
        outputs.push_back(val);
        return true;
    }
    else if (cur.IsComplex())
    {
        outputs.push_back(funcs.Complex2Str(cur, fmt, &ofmt, precision));
        return true;
    }

    std::string strprec = fmt;
    if (strprec.empty() && precision >= 0)
    {
        strprec = "%";
        if (precision == 0)
        {
            strprec += "0.0lf";
        }
        else
        {
            strprec += "0." + std::to_string(static_cast<long long>(precision)) 
                       + "lf";
        }
    }
    std::string strval = MatrixDisplay::GetOutputValues(cur, &ofmt, "\n", "  ", 
              strprec, strprec, -1);

    if (!strval.empty() && strval[strval.length() - 1] == '\n')
    {
        strval.pop_back();
    }
    outputs.push_back(strval);
    return true;
}
//------------------------------------------------------------------------------
// Converts double to string
//------------------------------------------------------------------------------
std::string BuiltInFuncsString::Dbl2Str(const Currency&     cur,
                                        const std::string&  fmtstr,
                                        const OutputFormat* fmt) const
{
    assert(cur.IsScalar());

    double val = cur.Scalar();
    if (IsNegInf_T (val))       
        return "-Inf";
    else if (IsInf_T(val))
        return "Inf";
    else if (IsNaN_T(val))
        return "NaN";

    if (!fmtstr.empty())
    {
        char   tmp [2048];
        size_t pos = fmtstr.find("d");

        if (pos != std::string::npos) // Print as integer
        {
            sprintf(tmp, fmtstr.c_str(), static_cast<int>(val));
        }
        else
        {
            sprintf(tmp, fmtstr.c_str(), val);
        }
        
        std::string  out = tmp;
        return out;
    }

    // No decimal places
    if (fmt && fmt->GetDecimalPart() == 0 && fmt->GetIntegerPart() != 0) 
    {
        Currency icur(static_cast<int>(val));
        icur.DispOutput();
        return icur.GetOutputString(fmt);
    }

    return cur.GetOutputString(fmt);
}
//------------------------------------------------------------------------------
// Converts double to string
//------------------------------------------------------------------------------
std::string BuiltInFuncsString::Dbl2Str(const Currency&     cur,
                                        const std::string&  fmtstr,
                                        const OutputFormat* fmt,
                                        int                 totaldigits) const
{
    assert(cur.IsScalar());

    double val = cur.Scalar();
    if (IsNegInf_T (val))   
    {
        return "-Inf";
    }
    else if (IsInf_T(val))
    {
        return "Inf";
    }
    else if (IsNaN_T(val))
    {
        return "NaN";
    }

    if (!fmtstr.empty())
    {
        char   tmp [2048];
        size_t pos = fmtstr.find("d");

        if (pos != std::string::npos) // Print as integer
        {
            sprintf(tmp, fmtstr.c_str(), static_cast<int>(val));
        }
        else
        {
            sprintf(tmp, fmtstr.c_str(), val);
        }
        
        std::string  out = tmp;
        return out;
    }

    if (totaldigits > 0)  // Total digits specified
    {
        std::ostringstream os;
        if (BuiltInFuncsUtils::IsInt(val))
        {
            os.setf(static_cast<std::ios_base::fmtflags>(0), std::ios::floatfield);  
        }
        os << std::setprecision(static_cast<std::streamsize>(totaldigits));
        os << val;
        std::string strval (os.str());
        if (!strval.empty())
        {
            size_t pos = strval.find(".");
            if (pos != std::string::npos) // Decimal digits are accurate only to 16th place
            {
                size_t len = strval.size();
                if (len - pos >= 16)
                {
                    strval = strval.substr(0, pos + 16 + 1);
                }
            }
        }
        return strval;
    }

    // No decimal places
    if (fmt && fmt->GetDecimalPart() == 0 && fmt->GetIntegerPart() != 0) 
    {
        Currency icur(static_cast<int>(val));
        icur.DispOutput();
        return icur.GetOutputString(fmt);
    }

    return cur.GetOutputString(fmt);
}
//------------------------------------------------------------------------------
// Converts complex to string
//------------------------------------------------------------------------------
std::string BuiltInFuncsString::Complex2Str(const Currency&     cur,
                                            const std::string&  fmtstr,
                                            const OutputFormat* fmt,
                                            int                 totaldigits) const
{
    assert(cur.IsComplex());

    hwComplex val = cur.Complex();

    // Split the complex number into real, imaginary and sign
    double rval = 0.0;
    double ival = 0.0;
    std::string isign;
    CurrencyDisplay::GetComplexNumberVals(val, rval, ival, isign);

    Currency rcur(rval);
    Currency icur(ival);

    rcur.DispOutput();
    icur.DispOutput();

    if (!fmtstr.empty())
    {
        std::string out = Dbl2Str(rcur, fmtstr, nullptr) + isign + 
                          Dbl2Str(icur, fmtstr, nullptr) + "i";
        return out;
    }
    else if (totaldigits != -1)
    {
        std::string out = Dbl2Str(rcur, fmtstr, fmt, totaldigits) + isign + 
                          Dbl2Str(icur, fmtstr, fmt, totaldigits) + "i";
        return out;
    }

    // No decimal places
    if (fmt && fmt->GetDecimalPart() == 0 && fmt->GetIntegerPart() != 0) 
    {
        std::string out = Dbl2Str(rcur, fmtstr, fmt) + isign + 
                          Dbl2Str(icur, fmtstr, fmt) + "i";
        return out;
    }
    
    std::string out = cur.GetOutputString(fmt);
    return out;
}
//------------------------------------------------------------------------------
// Returns true and creates single matrix from string inputs [str2mat]
//------------------------------------------------------------------------------
bool BuiltInFuncsString::Str2mat(EvaluatorInterface           eval,
                                 const std::vector<Currency>& inputs,
                                 std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    int nargin = static_cast<int>(inputs.size());
    int m      = 0;
    int n      = 0;

    std::vector<Currency> vecCur;
    vecCur.reserve(nargin);

    for (int i = 0; i < nargin; ++i)
    {
        Currency cur = inputs[i];
        if (cur.IsScalar())
        {
            vecCur.push_back(cur);
            m += 1;
            n = std::max(n, 1);
        }
        else if (cur.IsMatrixOrString())
        {
            const hwMatrix* mtx = cur.Matrix();
            assert(mtx);
            vecCur.push_back(cur);
            m += mtx->M();
            n = std::max(n, mtx->N());
        }
        else if (cur.IsCellArray())
        {
            HML_CELLARRAY* cell = cur.CellArray();
            if (!cell)
            {
                continue;
            }
            int cm = cell->M();
            int cn = cell->N();
            for (int ci = 0; ci < cm; ++ci)
            {
                for (int cj = 0; cj < cn; ++cj)
                {
                    Currency tmp = (*cell)(ci, cj);
                    if (!tmp.IsMatrixOrString())
                    {
                        throw OML_Error(OML_ERR_STRING_STRINGCELL, i + 1, OML_VAR_TYPE);
                    }
                    const hwMatrix* mtx = tmp.Matrix();
                    assert(mtx);
                    vecCur.push_back(tmp);
                    m += mtx->M();
                    n = std::max(n, mtx->N());
                }

            }
        }
        else
        {
            throw OML_Error(OML_ERR_STRING_STRINGCELL, i, OML_VAR_TYPE);
        }
    }
	
    std::unique_ptr<hwMatrix> mtx(EvaluatorInterface::allocateMatrix(
                                  m, n, hwMatrix::REAL));
    int row = 0;
    for (std::vector<Currency>::const_iterator itr = vecCur.begin(); 
         itr != vecCur.end(); ++itr)
    {
        Currency cur = *itr;
        if (cur.IsScalar())
        {
            (*mtx)(row, 0) = cur.Scalar();
            for (int j = 1; j < n; ++j)
            {
                (*mtx)(row, j) = ' ';
            }
            ++row;
        }
        else
        {
            const hwMatrix* tmp = cur.Matrix();
            assert(tmp);
            int numrows = tmp->M();
            int numcols = tmp->N();

            for (int i = 0; i < numrows; ++i)
            {
                for (int j = 0; j < numcols; ++j)
                {
                    (*mtx)(row, j) = (*tmp)(i, j);
                }
                for (int j = numcols; j < n; ++j)
                {
                    (*mtx)(row, j) = ' ';
                }
                ++row;
            }
        }
    }

    Currency c (mtx.release());
    c.SetMask(Currency::MASK_STRING);
    outputs.push_back(c);

    return true;
}