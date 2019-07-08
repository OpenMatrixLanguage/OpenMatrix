/**
* @file BuiltInFuncsConvert.cpp
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

#include "BuiltInFuncsConvert.h"

#include <algorithm>
#include <bitset>
#include <cassert>
#include <iomanip>
#include <memory>  // For std::unique_ptr

#include "BuiltInFuncsUtils.h"
#include "MatrixDisplay.h"
#include "OML_Error.h"

#include "hwMatrix.h"

//------------------------------------------------------------------------------
//! Returns true after converting decimal to hexdecimal string
//------------------------------------------------------------------------------
bool BuiltInFuncsConvert::hml_dec2hex(EvaluatorInterface           eval,
                                      const std::vector<Currency>& inputs,
                                      std::vector<Currency>&       outputs)
{
    if (inputs.empty()) throw OML_Error(OML_ERR_NUMARGIN);

    size_t numInputs = inputs.size();

    const Currency& input = inputs[0];  // Value to convert

    int width = 0; 
    if (numInputs > 1)                  // Optional arg specifying min width
    {
        const Currency& len = inputs[1];
	    bool  isInteger     = len.IsInteger();
        if (!isInteger || len.Scalar() < 0)
            throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_VALUE);
        width = static_cast<int>(len.Scalar());
    }

    // Check if this is a valid scalar
    bool validScalar = (input.IsScalar()) ? 
                       IsValidNaturalNumber(input.Scalar()) : false;

    if (!validScalar && !input.IsMatrix() && !input.IsCellArray())
        throw OML_Error(OML_ERR_NATURALNUM_MATRIX_CELLARRAY, 1, OML_VAR_VALUE);
	            
    if (validScalar)  // Valid Scalar
    {
        std::string hexstr(Dec2HexHelper(input.Scalar(), width));
        outputs.push_back(Currency(hexstr));
        return true;
    }

    // Get hexadecimal strings for matrix or cell array
    std::vector <std::string> hexstrings;

    if (input.IsMatrix())
        hexstrings = Dec2HexHelper(input.Matrix(), width);

    else 
        hexstrings = Dec2HexHelper(input.CellArray(), width);

    // Output will be a string matrix with one row per element, padded with 
    // leading zeros to the width of the largest value
    Currency formattedOutput = BuiltInFuncsUtils::FormatOutput(
        hexstrings, true, true, true, '0');

    outputs.push_back(formattedOutput);
    return true;
}
//------------------------------------------------------------------------------
//! Returns hex string given integer and minimum width to pad element with
//------------------------------------------------------------------------------
std::string BuiltInFuncsConvert::Dec2HexHelper(double val, int width)
{
    if (!IsValidNaturalNumber(val)) return "NaN";
			
    std::ostringstream os;
    os << std::right << std::setfill('0') << std::setw(width);
    os << std::hex << static_cast<unsigned long long>(val);

    std::string hexstr (os.str());
    std::transform(hexstr.begin(), hexstr.end(), hexstr.begin(), ::toupper);

    return hexstr;
}
//------------------------------------------------------------------------------
//! Returns vector of hex strings, given matrix and minimum width
//------------------------------------------------------------------------------
std::vector<std::string> BuiltInFuncsConvert::Dec2HexHelper(const hwMatrix* mtx, int width)
{
    int numrows = mtx ? mtx->M() : 0;
    int numcols = mtx ? mtx->N() : 0;

    std::vector <std::string> hexstrings;
    if (numrows == 0 || numcols == 0 || (mtx && !mtx->IsReal()))
        hexstrings.push_back("Nan");

    else
    {
        hexstrings.reserve(numrows * numcols);
        for (int j = 0; j < numcols; ++j) // Column major for conversion funcs
        {
            for (int i = 0; i < numrows; ++i)
            {
                std::string tmp (Dec2HexHelper((*mtx)(i, j), width));
                hexstrings.push_back(tmp);
            }
        }
    }

    return hexstrings;
}
//------------------------------------------------------------------------------
//! Returns vector of hex strings, given cell array and minimum width
//------------------------------------------------------------------------------
std::vector<std::string> BuiltInFuncsConvert::Dec2HexHelper(HML_CELLARRAY* cell, int width)
{
    int numelem = cell ? cell->Size() : 0;

    std::vector <std::string> hexstrings;
    if (numelem == 0)
    {
        hexstrings.push_back("Nan");
        return hexstrings;
    }
    
    hexstrings.reserve(numelem);
    for (int i = 0; i < numelem; ++i)
    {
        Currency    cur ((*cell)(i));
        
        std::string tmp = cur.IsScalar() ? 
                          Dec2HexHelper(cur.Scalar(), width) : "Nan";
        hexstrings.push_back(tmp);    
    }

    return hexstrings;
}
//------------------------------------------------------------------------------
//! Returns true after converting hexdecimal to decimal
//------------------------------------------------------------------------------
bool BuiltInFuncsConvert::hml_hex2dec(EvaluatorInterface           eval,
                                      const std::vector<Currency>& inputs,
                                      std::vector<Currency>&       outputs)
{
    if (inputs.empty()) throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input = inputs[0];  // Value to convert

    if (!input.IsString() && !input.IsCellArray())
        throw OML_Error(OML_ERR_STRING_MATRIX_CELLARRAY, 1, OML_VAR_VALUE);
	            
    Currency out = Hex2DecHelper(eval, input);
    outputs.push_back(out);
        
    return true;
}
//------------------------------------------------------------------------------
//! Returns currency given hex string
//------------------------------------------------------------------------------
Currency BuiltInFuncsConvert::Hex2DecHelper(const std::string& hexstr)
{
    if (hexstr.empty()) return Currency();

    // Check that each string is a valid 
    size_t len = hexstr.length();
    for (int i = 0; i < len; ++i)
    {
        if (!isxdigit(hexstr[i]))
        {
            Currency result (std::numeric_limits<double>::quiet_NaN());
            result.SetMask(Currency::TYPE_SCALAR);
            return result;
        }
    }

    int val = 0;

    std::stringstream ss;
    ss << std::hex << hexstr;
    ss >> val;

    return Currency(val);
}
//------------------------------------------------------------------------------
//! Returns currency (column matrix) for given input
//------------------------------------------------------------------------------
Currency BuiltInFuncsConvert::Hex2DecHelper(EvaluatorInterface& eval,
                                            const Currency&     input)
{
    std::vector<std::string> hexStr;
    if (input.IsMatrixOrString())
    {
        const hwMatrix* mtx = input.Matrix();
        int rows = mtx ? mtx->M() : 0;
           
        hexStr.reserve(rows);
        for (int i = 0; i < rows; ++i)
        {
            Currency thisrow = BuiltInFuncsUtils::ReadRow(eval, input, i);
            hexStr.push_back (thisrow.StringVal());
        }
    }
    else if (input.IsCellArray())
    {
        HML_CELLARRAY* cell = input.CellArray();
        int            rows = cell ? cell->Size() : 0;
        hexStr.reserve(rows);

        for (int i = 0; i < rows; ++i)
        {
            Currency cur ((*cell)(i));
        
            std::string tmp = cur.IsString() ? cur.StringVal() : "";
            hexStr.push_back(tmp);    
        }
    }

    if (hexStr.empty()) return Currency();

    int numrows = static_cast<int>(hexStr.size());
    hwMatrix* out = EvaluatorInterface::allocateMatrix(numrows, 1, hwMatrix::REAL);
    for (int i = 0; i < numrows; ++i)
    {
        std::string thisrow (hexStr[i]);
        if (thisrow.empty())
            (*out)(i, 0) = std::numeric_limits<double>::quiet_NaN();
        else
        {
            Currency tmp = Hex2DecHelper(thisrow);
            (*out)(i, 0) = tmp.Scalar();
        }
    }
    return Currency(out);
}
//------------------------------------------------------------------------------
// Returns true after converting decimal to binary [dec2bin]
//------------------------------------------------------------------------------
bool BuiltInFuncsConvert::hml_dec2bin(EvaluatorInterface           eval,
                                      const std::vector<Currency>& inputs,
                                      std::vector<Currency>&       outputs)
{
    if (inputs.empty()) throw OML_Error(OML_ERR_NUMARGIN);

    size_t numInputs = inputs.size();

    const Currency& input = inputs[0];  // Value to convert

    int width = 0; 
    if (numInputs > 1)                  // Optional arg specifying min width
    {
        const Currency& len = inputs[1];
	    bool  isInteger     = len.IsInteger();
        if (!isInteger || len.Scalar() < 0)
            throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_VALUE);
        width = static_cast<int>(len.Scalar());
    }

    // Check if this is a valid scalar
    bool validScalar = (input.IsScalar()) ? 
                       IsValidNaturalNumber(input.Scalar()) : false;

    if (!validScalar && !input.IsMatrix() && !input.IsCellArray())
        throw OML_Error(OML_ERR_NATURALNUM_MATRIX_CELLARRAY, 1, OML_VAR_VALUE);
	            
    if (validScalar)  // Valid Scalar
    {
        std::string str(Dec2BinHelper(input.Scalar(), width));
        outputs.push_back(Currency(str));
        return true;
    }

    // Get binary strings for matrix or cell array
    std::vector <std::string> outstr;

    if (input.IsMatrix())
        outstr = Dec2BinHelper(input.Matrix(), width);

    else 
        outstr = Dec2BinHelper(input.CellArray(), width);

    // Output will be a string matrix with one row per element, padded with 
    // leading zeros to the width of the largest value
    Currency out = BuiltInFuncsUtils::FormatOutput(outstr, false, true, true, '0');
    outputs.push_back(out);

    return true;
}
//------------------------------------------------------------------------------
//! Returns binary string given value and minimum width to pad element with
//------------------------------------------------------------------------------
std::string BuiltInFuncsConvert::Dec2BinHelper(double val, int width)
{
    if (!IsValidNaturalNumber(val)) throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_VALUE);

    unsigned long long dec    = static_cast<unsigned long long>(val);
    const size_t numberOfBits = sizeof(dec) * 8;

    std::bitset<numberOfBits> bs(dec);
    std::string tmp (bs.to_string()); // Converts to binary

    // Truncate leading 0's from the result
    size_t npos = tmp.find_first_of("1");
    if (npos == std::string::npos)
        tmp = "0";
    else
        tmp = tmp.substr(npos);

    if (width == 0) return tmp;

    std::ostringstream os;
    os << std::right << std::setfill('0') << std::setw(width) << tmp;
    return os.str();
}
//------------------------------------------------------------------------------
//! Returns vector of binary strings, given matrix and minimum width
//------------------------------------------------------------------------------
std::vector<std::string> BuiltInFuncsConvert::Dec2BinHelper(const hwMatrix* mtx, int width)
{
    int numrows = mtx ? mtx->M() : 0;
    int numcols = mtx ? mtx->N() : 0;

    std::vector <std::string> str;
    if (numrows == 0 || numcols == 0 || (mtx && !mtx->IsReal()))
        throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_VALUE);

    else
    {
        str.reserve(numrows * numcols);
        for (int j = 0; j < numcols; ++j) // Column major for conversion funcs
        {
            for (int i = 0; i < numrows; ++i)
            {
                std::string tmp (Dec2BinHelper((*mtx)(i, j), width));
                str.push_back(tmp);
            }
        }
    }

    return str;
}
//------------------------------------------------------------------------------
//! Returns vector of binary strings, given cell array and minimum width
//------------------------------------------------------------------------------
std::vector<std::string> BuiltInFuncsConvert::Dec2BinHelper(HML_CELLARRAY* cell, int width)
{
    int numelem = cell ? cell->Size() : 0;

    std::vector <std::string> outstr;
    if (numelem == 0) return outstr; // Should never get into this situation
    
    outstr.reserve(numelem);
    for (int i = 0; i < numelem; ++i)
    {
        Currency    cur ((*cell)(i));
        
        std::string tmp = cur.IsScalar() ? 
                          Dec2BinHelper(cur.Scalar(), width) : "";
        outstr.push_back(tmp);    
    }

    return outstr;
}
//------------------------------------------------------------------------------
// Returns true after converting binary string to decimal [bin2dec]
//------------------------------------------------------------------------------
bool BuiltInFuncsConvert::hml_bin2dec(EvaluatorInterface           eval,
                                      const std::vector<Currency>& inputs,
                                      std::vector<Currency>&       outputs)
{
    if (inputs.empty()) throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input = inputs[0];  // Value to convert

    if (!input.IsString() && !input.IsMatrix() && !input.IsCellArray())
        throw OML_Error(OML_ERR_STRING_MATRIX_CELLARRAY, 1, OML_VAR_VALUE);
	            
    Currency out = Bin2DecHelper(eval, input);
    outputs.push_back(out);
        
    return true;
}
//------------------------------------------------------------------------------
//! Returns currency given binary string
//------------------------------------------------------------------------------
Currency BuiltInFuncsConvert::Bin2DecHelper(const std::string& str)
{
    if (str.empty()) 
    {
        return Currency();
    }

    // Check that each string is a valid 
    size_t len = str.length();
    std::string input;
    for (size_t i = 0; i < len; ++i)
    {
        char c = str[i];
        if (c != '0' && c != '1' && c != ' ')
        {
            Currency result (std::numeric_limits<double>::quiet_NaN());
            result.SetMask(Currency::TYPE_SCALAR);
            return result;
        }
        if (c != ' ') 
        {
            input += c;
        }
    }
    long long num = atol(input.c_str());

    double result = 0;
    double i = 0;
    while (num != 0)
    {
        int rem = num % 10;
        num /= 10;
        result += rem * std::pow(2, i);
        ++i;
    }

    return Currency(result);
}
//------------------------------------------------------------------------------
//! Returns currency (column matrix) for given input
//------------------------------------------------------------------------------
Currency BuiltInFuncsConvert::Bin2DecHelper(EvaluatorInterface& eval,
                                            const Currency&     input)
{
    std::vector<std::string> outstr;
    if (input.IsMatrixOrString())
    {
        const hwMatrix* mtx = input.Matrix();
        int rows = mtx ? mtx->M() : 0;
           
        outstr.reserve(rows);
        for (int i = 0; i < rows; ++i)
        {
            Currency thisrow = BuiltInFuncsUtils::ReadRow(eval, input, i);
            outstr.push_back (thisrow.StringVal());
        }
    }
    else if (input.IsCellArray())
    {
        HML_CELLARRAY* cell = input.CellArray();
        int            rows = cell ? cell->Size() : 0;
        outstr.reserve(rows);

        for (int i = 0; i < rows; ++i)
        {
            Currency cur ((*cell)(i));
        
            std::string tmp = cur.IsString() ? cur.StringVal() : "";
            outstr.push_back(tmp);    
        }
    }

    if (outstr.empty()) return Currency();

    int numrows = static_cast<int>(outstr.size());
    hwMatrix* out = EvaluatorInterface::allocateMatrix(numrows, 1, hwMatrix::REAL);
    for (int i = 0; i < numrows; ++i)
    {
        std::string thisrow (outstr[i]);
        if (thisrow.empty())
            (*out)(i, 0) = std::numeric_limits<double>::quiet_NaN();
        else
        {
            Currency tmp = Bin2DecHelper(thisrow);
            (*out)(i, 0) = tmp.Scalar();
        }
    }
    return Currency(out);
}
//------------------------------------------------------------------------------
//! Returns true if given value is a non-negative integer
//------------------------------------------------------------------------------
bool BuiltInFuncsConvert::IsValidNaturalNumber(double val)
{
    if (IsNegInf_T(val) || IsInf_T(val) || IsNaN_T(val) || 
        !BuiltInFuncsUtils::IsInt(val) || val < 0) return false;
    return true;
}
//------------------------------------------------------------------------------
// Returns true and converts bit matrix to vector of integers [bi2de]
//------------------------------------------------------------------------------
bool BuiltInFuncsConvert::Bi2De(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
{
    if (inputs.empty()) 
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    size_t nargin = inputs.size();

    if (!inputs[0].IsMatrix())
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_TYPE);
    }
    const hwMatrix* mtx = inputs[0].Matrix();
    assert(mtx);
    if (!mtx)
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_VALUE);
    }
    if (!mtx->IsReal())
    {
        throw OML_Error(OML_ERR_REALMATRIX, 1, OML_VAR_VALUE);
    }

    double base = 2;
    if (nargin >= 2)
    {
        if (!inputs[1].IsPositiveInteger())
        {
            throw OML_Error(OML_ERR_INVALID_BASE, 2, OML_VAR_TYPE);
        }
        base = inputs[1].Scalar();

        if (base < 2)
        {
            throw OML_Error(OML_ERR_INVALID_BASE, 2);
        }
    }

    bool leftmsb = false;
    if (nargin >= 3)
    {
        if (!inputs[2].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);
        }
        std::string msb (inputs[2].StringVal());
        if (!msb.empty())
        {
            std::transform(msb.begin(), msb.end(), msb.begin(), ::tolower);
        }
        if (msb == "left-msb")
        {
            leftmsb = true;
        }
        else if (msb != "right-msb")
        {
            std::string err ("Error: invalid option in argument 3; must be ");
            err += "either 'right-msb' or 'left-msb'";
            throw OML_Error(err);
        }
    }

    int m = mtx->M();
    int n = mtx->N();
    std::unique_ptr<hwMatrix> out (EvaluatorInterface::allocateMatrix(
                                   m, 1, hwMatrix::REAL));
    for (int i = 0; i < m;  ++i)
    {
        double dec = 0;
        for (int j = 0; j < n; ++j)
        {
            double val = (*mtx)(i, j);
            if (!IsValidNaturalNumber(val))
            {
                throw OML_Error(OML_ERR_FINITE_NATURALNUM, 1, OML_VAR_VALUE);
            }
            if (val > base - 1)
            {
                std::string err ("Error: invalid value in argument 1; must be ");
                err += "between 0 and "; 
                err += std::to_string(static_cast<long long>(base - 1));
                throw OML_Error(err);
            }

            if (!leftmsb)
            {
                dec += val * (std::pow(base, j));
            }
            else
            {
                dec += val * (std::pow(base, n - j - 1));
            }
        }
        (*out)(i, 0) = dec;
    }

    outputs.push_back(out.release());
    return true;
}
//------------------------------------------------------------------------------
// Returns true and converts vector to binary representation [de2bi]
//------------------------------------------------------------------------------
bool BuiltInFuncsConvert::De2Bi(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
{
    if (inputs.empty()) 
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    size_t nargin = inputs.size();

    const Currency& cur1 = inputs[0];
    std::unique_ptr<hwMatrix> mtx = nullptr;
    if (!cur1.IsScalar() && !cur1.IsVector())
    {
        throw OML_Error(OML_ERR_SCALARVECTOR, 1, OML_VAR_TYPE);
    }
    if (cur1.IsScalar())
    {
        if (!cur1.Scalar())
        {
            throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_TYPE);
        }
        double val = static_cast<int>(cur1.Scalar());
        if (!IsValidNaturalNumber(val))
        {
            throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_TYPE);
        }
        mtx.reset(EvaluatorInterface::allocateMatrix(1, 1, val));
    }
    else
    {
        mtx.reset(EvaluatorInterface::allocateMatrix(cur1.Matrix()));
    }

    assert(mtx);
    if (!mtx)
    {
        throw OML_Error(OML_ERR_SCALARVECTOR, 1, OML_VAR_TYPE);
    }

    int numcols = -1;
    if (nargin >= 2)
    {
        const Currency& cur2 = inputs[1];
        if (cur2.IsPositiveInteger())
        {
            numcols = static_cast<int>(cur2.Scalar());
        }
        else if (cur2.IsMatrix())
        {
            const hwMatrix* sizemat = cur2.Matrix();
            if (sizemat->Size() == 1)
            {
                numcols = static_cast<int>((*sizemat)(0));
                if (!IsValidNaturalNumber(numcols))
                {
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_TYPE);
                }
            }
            else if (!sizemat->IsEmpty())
            {
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_TYPE);
            }
        }        
    }

    int base = 2;
    if (nargin >= 3)
    {
        if (!inputs[2].IsPositiveInteger())
        {
            throw OML_Error(OML_ERR_INVALID_BASE, 3, OML_VAR_TYPE);
        }
        base = static_cast<int>(inputs[2].Scalar());

        if (base < 2)
        {
            throw OML_Error(OML_ERR_INVALID_BASE, 3);
        }
    }

    bool leftmsb = false;
    if (nargin >= 4)
    {
        if (!inputs[3].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 4, OML_VAR_TYPE);
        }
        std::string msb (inputs[3].StringVal());
        if (!msb.empty())
        {
            std::transform(msb.begin(), msb.end(), msb.begin(), ::tolower);
        }
        if (msb == "left-msb")
        {
            leftmsb = true;
        }
        else if (msb != "right-msb")
        {
            std::string err ("Error: invalid option in argument 4; must be ");
            err += "either 'right-msb' or 'left-msb'";
            throw OML_Error(err);
        }
    }
    BuiltInFuncsUtils utils;

    int m = mtx->Size();
    int n = std::max(numcols, 1);

    std::unique_ptr<hwMatrix> out (EvaluatorInterface::allocateMatrix(
                                   m, n, hwMatrix::REAL));
    out->SetElements(0);

    for (int i = 0; i < m; ++i)
    {
        int val = static_cast<int>((*mtx)(i));
        if (!IsValidNaturalNumber(val))
        {
            throw OML_Error(OML_ERR_NNINTVECTOR, 1, OML_VAR_VALUE);
        }

        // Convert the input
        std::string strout;
        int itr    = 1;
        while (val != 0)
        {
            int rem = val % base;

            val    /= base;
            itr    *= 10;
            strout += std::to_string(static_cast<long long>(rem));
        }
        if (strout.empty())
        {
            strout = "0";
        }

        if (leftmsb)
        {
            std::reverse(strout.begin(), strout.end());
        }

        // Get the correct number of bits to display if output dimensions are 
        // specified
        int len = static_cast<int>(strout.size());
        if (leftmsb && numcols != -1 && len >= numcols)
        {
            strout = strout.substr(len - numcols);
            len = static_cast<int>(strout.size());

        }

        // Push the output to the matrix
        n = out->N();
        if (numcols == -1 && len > n)
        {
            if (!leftmsb)
            {
                // Zeros will be padded to the right of existing data
                utils.CheckMathStatus(eval, out->Resize(m, len, true));
            }
            if (leftmsb)
            {
                // Zeros will be padded to the left, before existing data
                hwMatrix* tmpmtx = EvaluatorInterface::allocateMatrix(out.get());
                out.reset(EvaluatorInterface::allocateMatrix(m, len, 
                    hwMatrix::REAL));
                out->SetElements(0);
                
                int startcol = len - n;
                for (int ii = 0; ii < i; ++ii)
                {
                    for (int jj = 0; jj < n; ++jj)
                    {
                        (*out)(ii, jj + startcol) = (*tmpmtx)(ii, jj);
                    }
                }
                delete tmpmtx;
                tmpmtx = nullptr;
            }
        }

        n = out->N();
        int startcol = 0;
        if (leftmsb && n > len)
        {
            startcol = n - len;
        }
        for (int j = 0; j < len; ++j)
        {
            if (j + startcol >= n)
            {
                break;
            }
            int ival = int(strout[j]) - 48;  // Value of 0 is 48 in ascii
            (*out)(i, j + startcol) = ival;
        }
    }
    outputs.push_back(out.release());
    return true;
} 
//------------------------------------------------------------------------------
// Returns true and converts vector to matrix [vec2mat]
//------------------------------------------------------------------------------
bool BuiltInFuncsConvert::Vec2Mat(EvaluatorInterface           eval,
                                  const std::vector<Currency>& inputs,
                                  std::vector<Currency>&       outputs)
{
    size_t nargin = (inputs.empty()) ? 0 : inputs.size();

    if (nargin < 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsMatrix())
    {
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_TYPE);
    }
    const hwMatrix* mtx = inputs[0].Matrix();
    assert(mtx);
    if (!mtx || mtx->Size() == 0)
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
        return true;
    }

    if (!inputs[1].IsPositiveInteger())
    {
        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_TYPE);
    }
    int n = static_cast<int>(inputs[1].Scalar());

    double    dval = 0;
    hwComplex cval;
    bool hascomplex = (!mtx->IsReal());
    if (nargin > 2)
    {
        if (mtx->IsReal())
        {
            if (!inputs[2].IsScalar())
            {
                throw OML_Error(OML_ERR_SCALAR, 3, OML_VAR_TYPE);
            }
            dval = inputs[2].Scalar();
        }
        else
        {
            if (inputs[2].IsScalar())
            {
                cval = hwComplex(inputs[2].Scalar(), 0);
            }
            else if (inputs[2].IsComplex())
            {
                cval = inputs[2].Complex();
            }
            else
            {
                throw OML_Error(OML_ERR_SCALAR_COMPLEX, 3, OML_VAR_TYPE);
            }
            
        }
    }

    assert(n > 0);
    int size = mtx->Size();
    int m    = size/n;

    if (size % n > 0)
    {
        ++m;
    }

    std::unique_ptr<hwMatrix> out(EvaluatorInterface::allocateMatrix(m, n, mtx->Type()));

    int paddedvals = 0;
    int k = 0;
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (k < size)
            {
                if (!hascomplex)
                {
                    (*out)(i, j) = (*mtx)(k);
                }
                else
                {
                    out->z(i, j) = mtx->z(k);
                }
            }
            else
            {
                if (!hascomplex)
                {
                    (*out)(i, j) = dval;
                }
                else
                {
                    out->z(i, j) = cval;
                }
                paddedvals ++;
            }
            ++k;
        }
    }
    outputs.push_back(out.release());

    if (eval.GetNargoutValue() > 1)
    {
        outputs.push_back(paddedvals);
    }
    return true;
}