/**
* @file BuiltInFuncsString.cpp
* @date November 2015
* Copyright (C) 2015-2022 Altair Engineering, Inc.  
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

    hwMatrix* tok = EvaluatorInterface::allocateMatrix(rowsTok, colsTok, true);
    hwMatrix* rem = EvaluatorInterface::allocateMatrix(rowsRem, colsRem, true);   

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

    hwMatrix* out = EvaluatorInterface::allocateMatrix(m, n, true);

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
    if (!iscell && !isstring)
    {
        throw OML_Error(OML_ERR_STRING_STRINGCELL, 1);
    }

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
        result = EvaluatorInterface::allocateMatrix(numrows, maxpos - minpos, true);
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
                catch (std::regex_error&)
                {
                    
                    // Try with basic flag
                    std::regex_constants::syntax_option_type flag = (ignorecase) ?
                        std::regex_constants::basic :
                        std::regex_constants::basic | std::regex_constants::icase;
                    try
                    {
                        // Negative look behind is like an assertion and is not
                        // supported
                        std::string neglookbehind("(?<!");
                        size_t pos1 = pat.find(neglookbehind);
                        if (pos1 != std::string::npos)
                        {
                            int numopenbracket = 0;
                            size_t startpos = pos1 + neglookbehind.length();
                            std::string newpat(pat.substr(0, pos1));
                            for (size_t k = startpos; k < pat.length(); ++k)
                            {
                                char ch = pat[k];
                                if (ch == '(')
                                {
                                    numopenbracket++;
                                }
                                else if (ch == ')')
                                {
                                    if (numopenbracket > 0)
                                    {
                                        numopenbracket--;
                                        continue;
                                    }
                                    newpat += pat.substr(k + 1);
                                    break;
                                }
                            }
                            if (!newpat.empty())
                            {
                                pat = newpat;
                            }
                        }
                        std::regex regexpat(pat, flag);
                        result = std::regex_replace(txt, regexpat, rep, replaceflag);
                    }
                    catch (std::regex_error& err)
                    {
                        BuiltInFuncsUtils utils;
                        utils.ThrowRegexError(err.code());
                    }
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
// Returns true and creates cell array from an input character vector e.g. string matrix [cellstr]
//------------------------------------------------------------------------------
bool BuiltInFuncsString::CellStr(EvaluatorInterface           eval,
	const std::vector<Currency>& inputs,
	std::vector<Currency>& outputs)
{
	size_t nargin = (!inputs.empty()) ? inputs.size() : 0;
	if (nargin < 1 || nargin > 1)
		throw OML_Error(OML_ERR_NUMARGIN);

	const Currency& cur = inputs[0];
	if (!cur.IsString())
		throw OML_Error(OML_ERR_STRINGVECTOR, 1, OML_VAR_TYPE);

	const hwMatrix* mtx1 = cur.ConvertToMatrix();
	int m = mtx1->M();

	assert(mtx1);
	if (!mtx1 || mtx1->Size() == 0)
	{
		outputs.push_back(EvaluatorInterface::allocateCellArray());
		return true;
	}

#if 0
	HML_CELLARRAY* cell1 = nullptr;
	cell1 = EvaluatorInterface::allocateCellArray(m, 1);
#endif 
	std::unique_ptr<HML_CELLARRAY> cell1(EvaluatorInterface::allocateCellArray(m, 1));

	BuiltInFuncsUtils utils;

	
	for (int i = 0; i < m; ++i)
	{
		std::unique_ptr<hwMatrix> row(EvaluatorInterface::allocateMatrix());
		utils.CheckMathStatus(eval, mtx1->ReadRow(i, *row));
		Currency tmp(row.release());
		tmp.SetMask(Currency::MASK_STRING);
		(*cell1)(i) = tmp;
	}
	outputs.push_back(cell1.release());
	
	return true;
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
                                  m, n, true));
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
//------------------------------------------------------------------------------
// Returns true and converts string to double without eval [str2double]
//------------------------------------------------------------------------------
bool BuiltInFuncsString::Str2Double(EvaluatorInterface           eval,
                                    const std::vector<Currency>& inputs,
                                    std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    bool isString = inputs[0].IsString();
    if (!isString && !inputs[0].IsCellArray())
    {
        outputs.push_back(std::numeric_limits<double>::quiet_NaN());
        return true;  // No error is thrown, instead NaN is returned
    }

    BuiltInFuncsString funcs;

    bool iscomplex = false;
    bool isscalar  = true;
    
    double rval = 0;
    double ival = 0;

    if (isString)
    {
        const hwMatrix* mtx = inputs[0].Matrix();
        if (!mtx || mtx->Size() == 0)
        {
            outputs.push_back(std::numeric_limits<double>::quiet_NaN());
            return true;
        }
        int m = mtx->M();
        if (m == 1)
        {
            funcs.Str2Num(inputs[0].StringVal(), rval, ival, isscalar);
            if (isscalar)
            {
                outputs.push_back(rval);
            }
            else
            {
                outputs.push_back(hwComplex(rval, ival));
            }
            return true;
        }

        // Multiline string
        std::unique_ptr<hwMatrix> out(EvaluatorInterface::allocateMatrix(
            m, 1, std::numeric_limits<double>::quiet_NaN()));

        BuiltInFuncsUtils  utils;

        for (int i = 0; i < m; ++i)
        {
            std::unique_ptr<hwMatrix> row(EvaluatorInterface::allocateMatrix());
            utils.CheckMathStatus(eval, mtx->ReadRow(i, *row));

            Currency rowcur(row.release());
            rowcur.SetMask(Currency::MASK_STRING);

            funcs.Str2Num(rowcur.StringVal(), rval, ival, isscalar);
            if (!iscomplex && !isscalar)
            {
                utils.CheckMathStatus(eval, out->MakeComplex());
                iscomplex = true;
            }
            if (!iscomplex)
            {
                (*out)(i) = rval;
            }
            else
            {
                out->z(i) = hwComplex(rval, ival);
            }
        }
        outputs.push_back(out.release());
    }
    else
    {
        HML_CELLARRAY* cell = inputs[0].CellArray();
        if (!cell || cell->Size() == 0)
        {
            outputs.push_back(std::numeric_limits<double>::quiet_NaN());
            return true;
        }

        int m = cell->M();
        int n = cell->N();

        std::unique_ptr<hwMatrix> mtx(EvaluatorInterface::allocateMatrix(m, n,
            std::numeric_limits<double>::quiet_NaN()));

        BuiltInFuncsUtils  utils;

        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                if (!(*cell)(i, j).IsString())
                {
                    continue;
                }
                funcs.Str2Num((*cell)(i, j).StringVal(), rval, ival, isscalar);
                if (!iscomplex && !isscalar)
                {
                    utils.CheckMathStatus(eval, mtx->MakeComplex());
                    iscomplex = true;
                }

                if (!iscomplex)
                {
                    (*mtx)(i, j) = rval;
                }
                else
                {
                    mtx->z(i, j) = hwComplex(rval, ival);
                }
            }
        }
        outputs.push_back(mtx.release());
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns true successful in converting string to scalar/complex
//------------------------------------------------------------------------------
bool BuiltInFuncsString::Str2Num(const std::string& in,
                                 double&            rval,
                                 double&            ival,
                                 bool&              isscalar)
{
    rval = std::numeric_limits<double>::quiet_NaN();
    ival = std::numeric_limits<double>::quiet_NaN();
    if (in.empty())
    {
        return false;
    }
    else if (in == "i" || in == "j")  // Special case
    {
        rval = 0;
        ival = 1;
        isscalar = false;
        return true;
    }
    std::string tmp(in);

    // Linux and VS2015 will not process numbers with double scientific 
    // notation which use 'D' instead of 'E'
    std::replace(tmp.begin(), tmp.end(), 'D', 'E');
    std::replace(tmp.begin(), tmp.end(), 'd', 'e');
    
    char*  cptr   = nullptr;
    double result = strtod(tmp.c_str(), &cptr);
    std::string imagstr = (cptr) ? cptr : nullptr;
    if (!IsValidStrtodResult(tmp, imagstr, result))
    {
        return false;
    }

    // Remove spaces
    if (!imagstr.empty()) 
    {
        std::string::iterator itr = std::remove(imagstr.begin(), imagstr.end(), ' ');
        imagstr.erase(itr, imagstr.end());
    }

    if (imagstr.empty())
    {
        rval = result;
        return true;   // Conversion succeeded
    }
    char lastch = imagstr[imagstr.length() - 1];
    if (lastch == 'i' || lastch == 'j')
    {
        imagstr.pop_back();
        for (size_t i = 0; i < imagstr.length(); ++i)
        {
            char ch = imagstr[i];
            if (ch == '.' || isdigit(ch) || (i == 0 && (ch == '+' || ch == '-')))
            {
                continue;
            }
            return false; // Conversion failed
        }
        isscalar = false;
        if (imagstr.empty())  // Only a complex number
        {
            ival = result;
            rval = 0;
        }
        else
        {
            char* endptr = nullptr;
            ival = strtod(imagstr.c_str(), &endptr);
            tmp = (endptr) ? endptr : "";
            if (!IsValidStrtodResult(imagstr, tmp, ival))
            {
                ival     = std::numeric_limits<double>::quiet_NaN();
                isscalar = true;
                return false;
            }
            rval = result;
        }
        return true;
    }
    return false;
}
//------------------------------------------------------------------------------
// Returns true if strtod conversion is successful
//------------------------------------------------------------------------------
bool BuiltInFuncsString::IsValidStrtodResult(const std::string& in,
                                             const std::string& end,
                                             double             val)
{
    if (in.empty() || errno == ERANGE)
    {
        return false;
    }
    else if (!IsZero(val) || in != end)
    {
        return true;
    }

    if (in.find_first_of("0") != std::string::npos)
    {
        return true;
    }
    
    return (end.find_first_of("ij") != std::string::npos);
}
//------------------------------------------------------------------------------
// Returns true if given pattern is found in the input [contains]
//------------------------------------------------------------------------------
bool BuiltInFuncsString::Contains(EvaluatorInterface           eval,
                                  const std::vector<Currency>& inputs,
                                  std::vector<Currency>&       outputs)
{
    int nargin = inputs.empty() ? 0 : static_cast<int>(inputs.size());
    if (nargin != 2 && nargin != 4)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    // Validate inputs
    bool isStringTxt = inputs[0].IsString();
    if (!isStringTxt && !inputs[0].IsCellArray())
    {
        throw OML_Error(OML_ERR_STRING_STRINGCELL, 1);
    }

    bool isStringPat = inputs[1].IsString();
    if (!isStringPat && !inputs[1].IsCellArray())
    {
        throw OML_Error(OML_ERR_STRING_STRINGCELL, 2);
    }

    bool ignorecase = false;
    if (nargin == 4)
    {
        if (!inputs[2].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 3);
        }
        std::string val(inputs[2].StringVal());
        if (!val.empty())
        {
            std::transform(val.begin(), val.end(), val.begin(), ::tolower);
        }
        if (val != "ignorecase")
        {
            throw OML_Error(OML_ERR_BAD_STRING, 3);
        }

        if (!inputs[3].IsLogical() && !inputs[3].IsInteger())
        {
            throw OML_Error(OML_ERR_LOGICAL, 4);
        }
        int tmp = static_cast<int>(inputs[3].Scalar());
        if (tmp != 0 && tmp != 1)
        {
            throw OML_Error(OML_ERR_FLAG_01, 4);
        }
        ignorecase = (tmp == 0) ? false : true;
    }

    const hwMatrix* mtx1 = (isStringTxt) ? inputs[0].Matrix() : nullptr;
    const hwMatrix* mtx2 = (isStringPat) ? inputs[1].Matrix() : nullptr;

    if (mtx1 && mtx2 && mtx1->Size() == 1 && mtx2->Size() == 1)
    {
        std::string txt(inputs[0].StringVal());
        std::string pat(inputs[1].StringVal());

        if (ignorecase)
        {
            if (!txt.empty())
            {
                std::transform(txt.begin(), txt.end(), txt.begin(), ::tolower);
            }
            if (!pat.empty())
            {
                std::transform(pat.begin(), pat.end(), pat.begin(), ::tolower);
            }
        }
        bool contains = false;
        if (!txt.empty() && !pat.empty() && txt.find(pat) != std::string::npos)
        {
            contains = true;
        }
        Currency out(contains);
        out.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(out);
        return true;
    }

    HML_CELLARRAY* cell1        = nullptr;
	bool           delete_cell1 = false;

    BuiltInFuncsUtils utils;
    if (inputs[0].IsCellArray())
    {
        cell1 = inputs[0].CellArray();
    }
    else
    {
        if (mtx1 && mtx1->Size())
        {
            int m = mtx1->M();
            cell1 = EvaluatorInterface::allocateCellArray(m, 1);
			delete_cell1 = true;

            for (int i = 0; i < m; ++i)
            {
                std::unique_ptr<hwMatrix> row(EvaluatorInterface::allocateMatrix());
                utils.CheckMathStatus(eval, mtx1->ReadRow(i, *row));
                Currency tmp(row.release());
                tmp.SetMask(Currency::MASK_STRING);
                (*cell1)(i) = tmp;
            }
        }
    }
    if (!cell1)
    {
        Currency cur(false);
        cur.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(cur);
        return true;
    }

    HML_CELLARRAY* cell2        = nullptr;
	bool           delete_cell2 = false;
    if (inputs[1].IsCellArray())
    {
        cell2 = inputs[1].CellArray();
    }
    else
    {
        if (mtx2 && mtx2->Size())
        {
            int m2 = mtx2->M();
            cell2 = EvaluatorInterface::allocateCellArray(m2, 1);
			delete_cell2 = true;

            for (int i = 0; i < m2; ++i)
            {
                std::unique_ptr<hwMatrix> row(EvaluatorInterface::allocateMatrix());
                utils.CheckMathStatus(eval, mtx2->ReadRow(i, *row));
                Currency tmp(row.release());
                tmp.SetMask(Currency::MASK_STRING);
                (*cell2)(i, 0) = tmp;
            }
        }
    }

    int m1 = cell1->M();
    int n1 = cell1->N();

    hwMatrix* out(EvaluatorInterface::allocateMatrix(m1, n1, true));

    out->SetElements(0);
    if (!cell2 || cell2->IsEmpty())
    {
        Currency cur(out);
        cur.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(cur);
        return true;
    }

    int m2 = cell2->M();
    int n2 = cell2->N();

    for (int i1 = 0; i1 < m1; ++i1)
    {
        for (int j1 = 0; j1 < n1; ++j1)
        {
            Currency curTxt = (*cell1)(i1, j1);
            if (!curTxt.IsString())
            {
                throw OML_Error(OML_ERR_STRING_STRINGCELL, 1);
            }
            std::string txt(curTxt.StringVal());
            if (txt.empty())
            {
                continue;
            }
            else if (ignorecase)
            {
                std::transform(txt.begin(), txt.end(), txt.begin(), ::tolower);
            }
            int contains = 0;
            for (int i2 = 0; i2 < m2 && contains == 0; ++i2)
            {
                for (int j2 = 0; j2 < n2 && contains == 0; ++j2)
                {
                    Currency curPat = (*cell2)(i2, j2);
                    if (!curPat.IsString())
                    {
                        throw OML_Error(OML_ERR_STRING_STRINGCELL, 2);
                    }
                    std::string pat(curPat.StringVal());
                    if (pat.empty())
                    {
                        continue;
                    }
                    else if (ignorecase)
                    {
                        std::transform(pat.begin(), pat.end(), pat.begin(), ::tolower);
                    }
                    if (txt.find(pat) != std::string::npos)
                    {
                        contains = 1;
                        break;
                    }
                }
            }
            (*out)(i1, j1) = contains;
        }
    }

	if (delete_cell1)
		delete cell1;

	if (delete_cell2)
		delete cell2;

    Currency cur(out);
    cur.SetMask(Currency::MASK_LOGICAL);
    outputs.push_back(cur);

    return true;
}
//------------------------------------------------------------------------------
// Returns true and strips leading/trailing characters from input [strip]
//------------------------------------------------------------------------------
bool BuiltInFuncsString::Strip(EvaluatorInterface eval,
                               const std::vector<Currency>& inputs,
                               std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    const Currency& cur1 = inputs[0];
    bool isstring1 = cur1.IsString();
    if (!isstring1 && !cur1.IsCellArray())
    {
        throw OML_Error(OML_ERR_STRING_STRINGCELL, 1);
    }

    size_t nargin = inputs.size();
    std::string mode("both");
    if (nargin > 1)
    {
        if (!inputs[1].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 2);
        }
        std::string tmp(inputs[1].StringVal());
        if (tmp.empty())
        {
            throw OML_Error(OML_ERR_BAD_STRING, 2);
        }
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
        if (tmp == "left")
        {
            mode = "left";
        }
        else if (tmp == "right")
        {
            mode = "right";
        }
        else if (tmp != "both")
        {
            throw OML_Error(OML_ERR_BAD_STRING, 2);
        }
    }

    std::vector<std::string> trim;
    bool isstring3 = true;
    if (nargin > 2)
    {
        const Currency& cur3 = inputs[2];
        isstring3 = cur3.IsString();
        if (!isstring3 && !cur3.IsCellArray())
        {
            throw OML_Error(OML_ERR_STRING_STRINGCELL, 3);
        }
        
        if (isstring3)
        {
            std::string tmp (cur3.StringVal());
            if (!tmp.empty())
            {
                trim.push_back(tmp);
            }
        }
        else
        {
            HML_CELLARRAY* cell = cur3.CellArray();
            if (cell)
            {
                int numelements = cell->Size();
                for (int i = 0; i < numelements; ++i)
                {
                    const Currency& elem = (*cell)(i);
                    if (!elem.IsString())
                    {
                        throw OML_Error(OML_ERR_STRING_STRINGCELL, 3);
                    }

                    std::string thistrim (elem.StringVal());
                    if (!thistrim.empty())
                    {
                        trim.push_back(thistrim);
                    }
                }
            }
        }
    }

    if (trim.empty())
    {
        trim.push_back(" ");
    }

    BuiltInFuncsString funcs;

    if (isstring1)
    {
        std::string in(cur1.StringVal());
        
        if (mode == "both" || mode == "left")
        {
            funcs.LeftTrim(in, trim);
        }

        if (mode == "both" || mode == "right")
        {
            funcs.RightTrim(in, trim);
        }
        outputs.push_back(in);
        return true;
    }

    // Cell input
    HML_CELLARRAY* cell1 = cur1.CellArray();
    if (!cell1 || cell1->Size() == 0)
    {
        outputs.push_back(EvaluatorInterface::allocateCellArray());
        return true;
    }

    int m = cell1->M();
    int n = cell1->N();
    std::unique_ptr<HML_CELLARRAY> out(EvaluatorInterface::allocateCellArray(
        m, n));

    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            const Currency& elem1 = (*cell1)(i, j);
            if (!elem1.IsString())
            {
                throw OML_Error(OML_ERR_STRING_STRINGCELL, 1);
            }
            std::string in(elem1.StringVal());


            if (mode == "both" || mode == "left")
            {
                funcs.LeftTrim(in, trim);
            }

            if (mode == "both" || mode == "right")
            {
                funcs.RightTrim(in, trim);
            }
            (*out)(i, j) = in;
        }
    }

    outputs.push_back(out.release());
    return true;
}
//------------------------------------------------------------------------------
// Helper method to left trim input
//------------------------------------------------------------------------------
void BuiltInFuncsString::LeftTrim(std::string& in, 
                                  const std::vector<std::string>& vec)
{
    if (vec.empty())
    {
        return;
    }

    std::vector<std::string> trim(vec);
    size_t insize = in.size();
    while (!trim.empty() && !in.empty())
    {
        for (std::vector<std::string>::iterator itr = trim.begin();
            itr != trim.end() && !in.empty();)
        {
            std::string str(*itr);
            if (str.empty())
            {
                itr = trim.erase(itr);
                continue;
            }

            // Cannot use find_first_not_of as if len > 1 as it will match any 
            // character and not the entire substring
            size_t len = str.length();
            if (len == 1)
            {
                size_t pos = in.find_first_not_of(str);
                if (pos < in.length())
                {
                    in = in.substr(pos);
                    if (in.find(str) == std::string::npos)
                    {
                        itr = trim.erase(itr);
                        continue;
                    }
                    ++itr;
                }
                else
                {
                    itr = trim.erase(itr);
                }
                continue;
            }
            size_t pos = in.find(str);
            if (pos == 0)
            {
                in = in.substr(len);
                itr = trim.erase(itr);
                continue;
            }
            else if (pos > in.length())
            {
                itr = trim.erase(itr);
                continue;
            }
            ++itr;
        }
        if (insize == in.size()) // Nothing new to trim
        {
            break;
        }
        insize = in.size();
    }
}
//------------------------------------------------------------------------------
// Helper method to right trim input
//------------------------------------------------------------------------------
void BuiltInFuncsString::RightTrim(std::string&                    in, 
                                   const std::vector<std::string>& vec)
{
    if (in.empty() || vec.empty())
    {
        return;
    }


    std::vector<std::string> trim(vec);
    size_t insize = in.size();
    while (!trim.empty() && !in.empty())
    {
        for (std::vector<std::string>::iterator itr = trim.begin();
            itr != trim.end() && !in.empty();)
        {
            std::string str(*itr);
            if (str.empty())
            {
                itr = trim.erase(itr);
                continue;
            }

            // Cannot use find_last_not_of as if len > 1 as it will match any 
            // character and not the entire substring
            if (str.length() == 1)
            {
                size_t pos = in.find_last_not_of(str);
                if (pos < in.length())
                {
                    in.erase(pos + 1);
                    if (in.find(str) == std::string::npos)
                    {
                        itr = trim.erase(itr);
                        continue;
                    }
                    ++itr;
                }
                else
                {
                    itr = trim.erase(itr);
                }
                continue;
            }

            size_t pos = in.rfind(str);
            if (pos == std::string::npos)
            {
                itr = trim.erase(itr);
                continue;
            }
            else if (pos + str.length() == in.length())
            {
                in = in.substr(0, pos);
                itr = trim.erase(itr);
                continue;
            }
            ++itr;
        }
        if (insize == in.size()) // Nothing new to trim
        {
            break;
        }
        insize = in.size();
    }
}//------------------------------------------------------------------------------
// Returns true if successful in converting string to scalar/complex
//------------------------------------------------------------------------------
bool BuiltInFuncsString::IsNumber(const std::string& in, Currency& result)
{
    if (in.empty())
    {
        return false;
    }

    double rval = 0;
    double ival = 0;
    bool   isscalar = true;
    BuiltInFuncsString funcs;
    if (funcs.Str2Num(in, rval, ival, isscalar))
    {
        if (isscalar)
        {
            result = Currency(rval);
        }
        else
        {
            result = Currency(hwComplex(rval, ival));
        }
        return true;
    }
    return false;
}
//------------------------------------------------------------------------------
// Does a case-sensitive comparison of inputs [strcmp]
//------------------------------------------------------------------------------
bool BuiltInFuncsString::Strcmp(EvaluatorInterface           eval, 
                                const std::vector<Currency>& inputs, 
                                std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    const Currency& cur1 = inputs[0];
    bool            cellInput1 = cur1.IsCellArray();
    if (!cur1.IsString() && !cellInput1)
    {
        outputs.emplace_back(false);
        return true;
    }
    const Currency& cur2 = inputs[1];
    bool            cellInput2 = cur2.IsCellArray();
    if (!cur2.IsString() && !cellInput2)
    {
        outputs.emplace_back(false);
        return true;
    }

    Currency input1(cur1);
    if (cellInput1 && cur1.CellArray()->Size() == 1)
    {
        HML_CELLARRAY* temp = cur1.CellArray();
        if (!temp || !(*temp)(0).IsString())
        {
            outputs.emplace_back(false);
            return true;
        }
        input1 = (*temp)(0);
        input1.SetMask(Currency::MASK_STRING);
        cellInput1 = false;
    }

    Currency input2(cur2);
    if (cellInput2 && cur2.CellArray()->Size() == 1)
    {
        HML_CELLARRAY* temp = cur2.CellArray();
        if (!temp || !(*temp)(0).IsString())
        {
            outputs.emplace_back(false);
            return 1;
        }
        input2 = (*temp)(0);
        input2.SetMask(Currency::MASK_STRING);
        cellInput2 = false;
    }
    if (cellInput1)  // First input is a cell
    {
        HML_CELLARRAY* cell1 = input1.CellArray();
        int m = cell1->M();
        int n = cell1->N();
        int size = cell1->Size();
        std::unique_ptr<hwMatrix> outmtx(
            EvaluatorInterface::allocateMatrix(m, n, 0.0));

        if (cellInput2)            // Second input is also a cell
        {
            HML_CELLARRAY* cell2 = input2.CellArray();
            if (m != cell2->M() || n != cell2->N())
            {
                std::string msg("Error: invalid input in argument 2; ");
                msg += "dimensions [" + std::to_string(cell2->M()) + "x"
                    + std::to_string(cell2->N()) + "] must match dimensions ["
                    + std::to_string(m) + "x" + std::to_string(n)
                    + "] of argument 1";
                throw OML_Error(msg);
            }
            for (int i = 0; i < size; ++i)
            {
                const Currency& c1 = (*cell1)(i);
                const Currency& c2 = (*cell2)(i);
                if (!c1.IsString() || !c2.IsString())
                {
                    continue;
                }
                const hwMatrix* mtx1 = c1.Matrix();
                const hwMatrix* mtx2 = c2.Matrix();
                if (*mtx1 == *mtx2)
                {
                    (*outmtx)(i) = 1;
                }
            }
        }
        else    // Second input is a string
        {
            for (int i = 0; i < size; ++i)
            {
                const Currency& c1 = (*cell1)(i);
                if (c1.IsString() && (*(c1.Matrix()) == *(input2.Matrix())))
                {
                    (*outmtx)(i) = 1;
                }
            }
        }
        Currency out(outmtx.release());
        out.SetMask(Currency::MASK_LOGICAL);
        outputs.emplace_back(out);
        return true;
    }
    else if (cellInput2) // Second input is a cell
    {
        HML_CELLARRAY* cell2 = input2.CellArray();
        int m = cell2->M();
        int n = cell2->N();
        int size = cell2->Size();
        std::unique_ptr<hwMatrix> outmtx(
            EvaluatorInterface::allocateMatrix(m, n, 0.0));

        for (int i = 0; i < size; ++i)
        {
            const Currency& c1 = (*cell2)(i);
            if (c1.IsString() && (*(c1.Matrix()) == *(input1.Matrix())))
            {
                (*outmtx)(i) = 1;
            }
        }
        Currency out(outmtx.release());
        out.SetMask(Currency::MASK_LOGICAL);
        outputs.emplace_back(out);
        return true;
    }
    else // Both inputs are strings
    {
        if (*(input1.Matrix()) == *(input2.Matrix()))
        {
            outputs.emplace_back(true);
        }
        else
        {
            outputs.emplace_back(false);
        }
    }
    return true;
}
//------------------------------------------------------------------------------
// Does a case-insensitive comparison of inputs [strcmpi]
//------------------------------------------------------------------------------
bool BuiltInFuncsString::Strcmpi(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    std::vector<Currency> lowerin;
    lowerin.reserve(2);

    for (std::vector<Currency>::const_iterator itr = inputs.begin();
        itr != inputs.end(); ++itr)
    {
        const Currency& cur = *itr;
        if (cur.IsString())
        {
            lowerin.emplace_back(ToLowerString(cur.Matrix()));
        }
        else if (cur.IsCellArray())
        {
            std::unique_ptr<HML_CELLARRAY> cell(
                EvaluatorInterface::allocateCellArray(cur.CellArray()));
            int csize = cell->Size();
            for (int i = 0; i < csize; ++i)
            {
                if ((*cell)(i).IsString())
                {
                    (*cell)(i) = ToLowerString((*cell)(i).Matrix());
                }
            }
            lowerin.emplace_back(cell.release());
        }
        else
        {
            outputs.push_back(false);
            return true;
        }
    }

    return Strcmp(eval, lowerin, outputs);
}
//------------------------------------------------------------------------------
// Utility function to convert matrix to a lower case string
//------------------------------------------------------------------------------
Currency BuiltInFuncsString::ToLowerString(const hwMatrix* mtx)
{
    if (mtx)
    {
        std::unique_ptr<hwMatrix> out(EvaluatorInterface::allocateMatrix(mtx));

        int size = out->Size();
        for (int i = 0; i < size; ++i)
        {
            (*out)(i) = tolower(static_cast<int>((*out)(i)));
        }
        Currency cur(out.release());
        cur.SetMask(Currency::MASK_STRING);
        return cur;
    }
    return Currency("");
}