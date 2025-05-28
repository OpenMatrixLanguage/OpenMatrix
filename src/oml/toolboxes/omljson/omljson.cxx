/**
* @file omljson.cxx
* @date September 2023
* Copyright (C) 2023-2024 Altair Engineering, Inc.
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

#include "omljson.h"

#include <algorithm>
#include <cassert>
#include <iomanip>

#include "BuiltInFuncsUtils.h"
#include "OML_Error.h"
#include "StructData.h"

std::string g_indent = "  ";
std::string g_quote  = "\"";

#define TBOXVERSION 2024.1

//------------------------------------------------------------------------------
// Returns library version
//------------------------------------------------------------------------------
double GetToolboxVersion(EvaluatorInterface)
{
    return TBOXVERSION;
}
//------------------------------------------------------------------------------
// Entry point which registers the library with oml
//------------------------------------------------------------------------------
int InitDll(EvaluatorInterface eval)
{
    eval.RegisterBuiltInFunction("jsondecode", &OmlJSON::Decode, FunctionMetaData(-7, 1, "JSON"));
    eval.RegisterBuiltInFunction("jsonencode", &OmlJSON::Encode, FunctionMetaData(-5, 1, "JSON"));

    return 1;
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
OmlJSON::OmlJSON()
    : _checkname (true)
    , _convert   (true)
    , _format    (false)
    , _prefix    ("x")
    , _style     (REPLACESTYLE_UNDERSCORE)
{
}
//------------------------------------------------------------------------------
// Decodes JSON text to oml data and returns true [jsondecode]
//------------------------------------------------------------------------------
bool OmlJSON::Decode(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    else if (!inputs [0].IsString()) 
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

    std::string jstr (inputs [0].StringVal());
    Trim(jstr);
    if (jstr.empty())
    {
        outputs.emplace_back(inputs [0].StringVal());
        return true;
    }

    OmlJSON json;

    int     nargin = static_cast<int>(inputs.size());
    for (int i = 1; i < nargin; ++i)
    {
        if (!inputs [i].IsString()) // option
        {
            throw OML_Error(OML_ERR_STRING, i + 1);
        }

        std::string opt(inputs [i].StringVal());
        std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);

        ++i; // Get the value
        if (i >= nargin)
        {
            throw OML_Error(OML_ERR_MISSING_VALUE, opt, i);
        }

        if (opt == "makevalidname")
        {
            if (!inputs [i].IsInteger())
            {
                throw OML_Error(OML_ERR_LOGICAL, i + 1);
            }
            json._checkname = (inputs [i].Scalar() == 0) ? false : true;
            if (json._checkname && inputs [i].Scalar() != 1)
            {
                throw OML_Error(OML_ERR_LOGICAL, i + 1);
            }
        }
        else if (opt == "prefix")
        {
            if (!inputs [i].IsString())
            {
                throw OML_Error(OML_ERR_STRING, i + 1);
            }
            json._prefix = inputs [i].StringVal();
            if (json._prefix.empty())
            {
                throw OML_Error(OML_ERR_NONEMPTY_STR, i + 1);
            }
            if (!IsValidPrefix(json._prefix))
            {
                throw OML_Error("Error: invalid value in argument " +
                    std::to_string(i + 1) +
                    "; prefix must have only alpha-numeric or '_' character(s)");
            }
        }
        else if (opt == "replacementstyle")
        {
            std::string msg ("Error: invalid option in argument ");
            msg += std::to_string(i + 1) + "; use delete, hex or underscore";
            if (!inputs [i].IsString())
            {
                throw OML_Error(msg);
            }
            std::string val(inputs [i].StringVal());
            std::transform(val.begin(), val.end(), val.begin(), ::tolower);
            if (val == "delete")
            {
                json._style = REPLACESTYLE_DELETE;
            }
            else if (val == "hex")
            {
                json._style = REPLACESTYLE_HEX;
            }
            else if (val == "underscore")
            {
                json._style = REPLACESTYLE_UNDERSCORE;
            }
            else
            {
                if (!opt.empty())
                {
                    msg += "; [" + opt + "]";
                }
                throw OML_Error(msg);
            }
        }
        else
        {
            std::string msg ("Error: invalid option in argument ");
            msg += std::to_string(i) + "; use makevalidname, prefix or replacementstyle";
            if (!opt.empty())
            {
                msg += "; [" + opt + "]";
            }
            throw OML_Error(msg);
        }
    }
    outputs.emplace_back(json.JSON2Currency(jstr));

    return true;
}
//------------------------------------------------------------------------------
// Encodes oml data to JSON text and returns true [jsonencode]
//------------------------------------------------------------------------------
bool OmlJSON::Encode(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    
    OmlJSON json;
    int     nargin = static_cast<int>(inputs.size());
    for (int i = 1; i < nargin; ++i)
    {
        if (!inputs[i].IsString()) // option
        {
            throw OML_Error(OML_ERR_STRING, i + 1);
        }
        std::string opt(inputs [i].StringVal());
        std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);

        if (!(opt == "convertinfandnan" || opt == "prettyprint"))
        {
            std::string msg ("Error: invalid option in argument ");
            msg += std::to_string(i + 1) +
                "; use convertinfandnan or prettyprint";
            if (!opt.empty())
            {
                msg += "; [" + opt + "]";
            }
            throw OML_Error(msg);
        }

        ++i; // Get the value
        if (i >= nargin)
        {
            throw OML_Error(OML_ERR_MISSING_VALUE, opt, i);
        }
        else if (!inputs[i].IsInteger())
        {
            throw OML_Error(OML_ERR_LOGICAL, i + 1);
        }
        bool optval = (inputs [i].Scalar() == 1) ? true : false;
        if (!optval && inputs [i].Scalar() != 0)
        {
            throw OML_Error(OML_ERR_LOGICAL, i + 1);
        }
        if (opt == "convertinfandnan")
        {
            json._convert = optval;
        }
        else if (optval && opt == "prettyprint")
        {
            json._format  = true;
            json._newline = "\n";
        }
    }

    std::string jstr(json.Currency2JSON(inputs [0], 0));
    outputs.emplace_back(jstr);

    return true;
}
//------------------------------------------------------------------------------
// Converts scalar to JSON format string with no indent
//------------------------------------------------------------------------------
std::string OmlJSON::Scalar2JSON(double val, bool islogical) const
{
    if (islogical)
    {
        return (val == 0) ? "false" : "true";
    }
    else if (IsNaN_T(val))
    {
        return _convert ? "null" : "NaN";
    }
    else if (IsInf_T(val))
    {
        return _convert ? "null" : "Infinity";
    }
    else if (IsNegInf_T(val))
    {
        return _convert ? "null" : "-Infinity";
    }

    Currency cur (val);
    if (cur.IsInteger())
    {
        std::ostringstream os;
        os << std::setprecision(std::numeric_limits<double>::max_digits10)
            << val;
        return std::string(os.str());
    }
    std::string curval (cur.GetOutputString(nullptr));
    if (curval.find("e") != std::string::npos || curval.length() <= 5)
    {
        return curval;
    }
    std::ostringstream os;
    os << std::setprecision(std::numeric_limits<double>::max_digits10)
       << val;

    return std::string(os.str());
}
//------------------------------------------------------------------------------
// Converts single dimension string to JSON format string
//------------------------------------------------------------------------------
std::string OmlJSON::String2JSON(const std::string& in) const
{
    size_t len = in.length();
    std::string jstr;
    for (size_t i = 0; i < len; ++i)
    {
        char ch = in [i];
        if (ch == '\n')
        {
            jstr += "\\n";
            //jstr += "n";
        }
        else if (ch == '\t')
        {
            jstr += "\\t";
        }
        else if (ch == '\\')
        {
            jstr += "\\\\";
        }
        else if (ch == '\"')
        {
            jstr += "\\\"";
        }
        else
        {
            jstr += ch;
        }
    }
    return (g_quote + jstr + g_quote);
}
//------------------------------------------------------------------------------
//  Converts 2D matrix to JSON format string
//------------------------------------------------------------------------------
std::string OmlJSON::Matrix2JSON(const Currency& in, int level) const
{
    std::string indent0 (Indent(level - 1));
    const hwMatrix* mtx = in.Matrix();
    if (!mtx || mtx->Size() == 0)
    {
        return indent0 + "[]";
    }

    std::string out     (indent0 + "[");
    std::string indent2 (IncrIndent(indent0));

    if (in.IsVector())
    {
        int msize = mtx->Size();
        if (mtx->IsReal())
        {
            for (int i = 0; i < msize; ++i)
            {
                out += _newline + indent2 + Scalar2JSON((*mtx)(i))
                    + Comma(i, msize);
            }
        }
        else
        {
            for (int i = 0; i < msize; ++i)
            {
                double dval = mtx->z(i).Real();
                out += _newline + indent2 + Scalar2JSON(dval) + Comma(i, msize);
            }
        }
        out += _newline + indent0 + "]";
        return out;
    }
    int m = mtx->M();
    int n = mtx->N();

    std::string indent3 (IncrIndent(indent2));

    if (mtx->IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            out += _newline + indent2 + "[";
            
            for (int j = 0; j < n; ++j)
            {
                out += _newline + indent3 + Scalar2JSON((*mtx)(i, j)) + Comma(j, n);
            }
            out += _newline + indent2 + "]" + Comma(i, m);
        }
    }
    else
    {
        for (int i = 0; i < m; ++i)
        {
            out += _newline + indent2 + "[";

            for (int j = 0; j < n; ++j)
            {
                double dval = mtx->z(i, j).Real();
                out += _newline + indent3 + Scalar2JSON(dval) + Comma(j, n);
            }
            out += _newline + indent2 + "]" + Comma(i, m);
        }

    }
    out += _newline + indent0 + "]";
    return out;
}
//------------------------------------------------------------------------------
//  Converts currency to JSON format string
//------------------------------------------------------------------------------
std::string OmlJSON::Currency2JSON(const Currency& in, int level) const
{
    std::string indent0 (Indent(level - 1));  // Same as parent

    if (in.IsScalar())
    {
        return indent0 + Scalar2JSON(in.Scalar(), in.IsLogical());
    }
    else if (in.IsMultilineString())
    {
        return StringND2JSON(in, level);
    }
    else if (in.IsString())
    {
        return indent0 + String2JSON(in.StringVal());
    }
    else if (in.IsComplex())
    {
        return Complex2JSON(in, level);
    }
    else if (in.IsMatrix())
    {
        return Matrix2JSON(in, level);
    }
    else if (in.IsNothing() || in.IsEmpty())
    {
        return indent0 + "[]";
    }
    else if (in.IsStruct())
    {
        return Struct2JSON(in, level);
    }
    else if (in.IsError())
    {
        return indent0 + String2JSON(in.Message());
    }
    else if (in.IsCellArray())
    {
        return Cell2JSON(in, level);
    }
    throw OML_Error(OML_ERR_DATATYPE, in.GetTypeString(), 1);
}
//------------------------------------------------------------------------------
//  Converts JSON format string to Currency
//------------------------------------------------------------------------------
Currency OmlJSON::JSON2Currency(const std::string& in) const
{
    std::string jstr (in);
    Trim(jstr);

    if (jstr.empty())
    {
        return Currency();
    }
    
    char first = jstr [0];

    Currency output;

    if (first == '\"')        // String
    {
        return JSON2String(jstr);
    }
    else if (JSON2Matrix(jstr, output))
    {
        return output;
    }
    else if (first == '{')    // Struct
    {
        return JSON2Struct(jstr);
    }
    else if (JSON2Scalar(jstr, output))
    {
        return output;
    }
    else if (first == '[') // Cell array or struct array
    {
        if (JSON2StructArray(jstr, output)) // Struct array
        {
            return output;
        }
        return JSON2Cell(jstr);  // Process as a cell
    }
    throw OML_Error(OML_ERR_INVALIDFORMAT, 1);
}
//------------------------------------------------------------------------------
// Converts multiline strings to JSON format string
//------------------------------------------------------------------------------
std::string OmlJSON::StringND2JSON(const Currency& in, int level) const
{
    const hwMatrix* mtx = in.Matrix();
    if (!mtx || mtx->Size() == 0)
    {
        return g_quote + g_quote;
    }

    std::string indent1 (Indent(level - 1));
    std::string indent2 (_format ? indent1 + "  " : "");
    std::string out (indent1 + "[" + _newline);

    int rows = mtx->M();
    int cols = mtx->N();
    for (int i = 0; i < rows; ++i)
    {
        std::string row;
        for (int j = 0; j < cols; ++j)
        {
            if ((*mtx)(i, j) != 0x00)
            {
                row += ((*mtx)(i, j) != 0x00) ?
                    static_cast<unsigned char>((*mtx)(i, j)) : ' ';
            }
        }
        out += indent2 + String2JSON(row) + Comma(i, rows) + _newline;
    }
    out += indent1 + ']';
    return out;
}
//------------------------------------------------------------------------------
// Converts struct to JSON format string
//------------------------------------------------------------------------------
std::string OmlJSON::Struct2JSON(const Currency& in, int level) const
{
    std::string indent1 (Indent(level - 1));

    StructData* sd = in.Struct();
    if (!sd || sd->IsEmpty())
    {
        return indent1 + "{}";
    }
    const std::map<std::string, int> fields (sd->GetFieldNames());
    if (fields.empty() || sd->Size() < 1)
    {
        return indent1 + "{}";
    }

    std::string indent2 (Indent(level + 1));

    int nfields = static_cast<int>(fields.size());
    int size    = sd->Size();
    if (size <= 1)  // Single struct
    {
        std::string out(indent1 + "{");
        int index = 0;
        for (std::map<std::string, int>::const_iterator itr = fields.begin();
             itr != fields.end(); ++itr, ++index)
        {
            std::string fld (itr->first);
            out += _newline + indent2 + String2JSON(fld) + ':'
                + Currency2JSON(sd->GetValue(0, 0, fld), level)
                + Comma(index, nfields);
        }
        out += _newline + indent1 + "}";
        return out;
    }

    // Struct array
    int m = sd->M();
    int n = sd->N();

    std::string indent3 (IncrIndent(indent2));
    std::string out(indent1 + "[");

    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < m; ++i)
        {
            out += _newline + indent2 + "{";
            int index = 0;
            for (std::map<std::string, int>::const_iterator itr = fields.begin();
                 itr != fields.end(); ++itr, ++index)
            {
                std::string fld (itr->first);
                out += _newline + indent3 + String2JSON(fld) + ':'
                    + Currency2JSON(sd->GetValue(i, j, fld), level + 2)
                    + Comma(index, nfields);
            }
            out += _newline + indent2 + "}" + Comma(i, m);
        }
        out += _newline;
        if (j < n - 1)
        {
            out += indent2 + Comma(j, n);
        }
    }
    out += indent1 + ']';
    return out;
}
//------------------------------------------------------------------------------
// Utility to convert to JSON format string to struct
//------------------------------------------------------------------------------
Currency OmlJSON::JSON2Struct(const std::string& jstr) const
{
    std::map<std::string, std::vector<Currency> > data;
    ParseStructData(jstr, data);

    std::unique_ptr<StructData> sd(EvaluatorInterface::allocateStruct());
    std::map<std::string, std::vector<Currency> >::const_iterator itr =  data.begin();
    for (; itr != data.end(); ++itr)
    {
        std::vector<Currency> vals (itr->second);
        sd->SetValue(0, 0, itr->first, 
                    (!vals.empty()) ? vals.back() : Currency());
    }
    return sd.release();
}
//------------------------------------------------------------------------------
// Gets name-value pairs for a struct
//------------------------------------------------------------------------------
void OmlJSON::ParseStructData(const std::string& jstr,
    std::map<std::string, std::vector<Currency> >& data) const
{
    if (jstr.empty())
    {
        return;
    }
    std::string str (jstr);
    Trim(str);
    if (str.empty())
    {
        return;
    }
    else if (str [0] == ',')
    {
        str.erase(str.begin());
    }
    Trim(str);
    if (str.empty() || str [0] != '{')
    {
        return;
    }

    assert(str [0] == '{');

    size_t len = str.size();
    if (len < 2 || str [len - 1] != '}')
    {
        throw OML_Error(OML_ERR_INVALIDFORMAT, "missing }", 1);
    }
    str.erase(str.begin());
    str.pop_back();

    while (!str.empty())
    {
        Trim(str);
        
        // Field name
        size_t pos1 = str.find(':');
        if (pos1 == std::string::npos)
        {
            break;
        }
        std::string name (str.substr(0, pos1));
        
        Currency cur = JSON2Currency(name);
        if (!cur.IsString())
        {
            throw OML_Error(OML_ERR_INVALIDFORMAT, "struct field must be a string", 1);
        }
        name = cur.StringVal();

        str = str.substr(pos1 + 1);
        if (str.empty())
        {
            throw OML_Error(OML_ERR_INVALIDFORMAT, "missing value for field" + name, 1);
        }

        name = FieldName(name);

        std::string val (GetFieldValueString(str));
        cur = JSON2Currency(val);
        std::vector<Currency> vec;
        if (data.empty() || data.find(name) == data.end())
        {
            vec.emplace_back(cur);
            data.emplace(std::make_pair(name, vec));
        }
        else
        {
            vec = data [name];
            vec.emplace_back(cur);
            data [name] = vec;
        }
    }
}
//------------------------------------------------------------------------------
// Utility to convert to JSON format string to matrix
//------------------------------------------------------------------------------
bool OmlJSON::JSON2Matrix(const std::string& in, Currency& cur) const
{
    std::string jstr (in);
    Trim(jstr);
    if (jstr.empty() || jstr == "[]" || IsNull(jstr))
    {
        cur = EvaluatorInterface::allocateMatrix();
        return true;
    }
    else if (jstr [0] != '[')
    {
        return false;
    }

    // Could be a Cell / matrix / struct array
    size_t len = jstr.length();
    if (len < 2 || jstr [len - 1] != ']')
    {
        throw OML_Error(OML_ERR_INVALIDFORMAT, "missing end square bracket", 1);
    }
    jstr = jstr.substr(1, len - 1);
    if (jstr.find(':')  != std::string::npos ||
        jstr.find('{')  != std::string::npos ||
        jstr.find('\"') != std::string::npos ||
        jstr.find("[]") != std::string::npos)
    {
        return false; // Could be a cell / struct
    }

    int rows = static_cast<int>(std::count(jstr.begin(), jstr.end(), '['));
    int cols = 0;
    if (rows == 0)
    {
        // This could be a vector
        rows = static_cast<int>(std::count(jstr.begin(), jstr.end(), ',')) + 1;
        if (rows == 0)
        {
            return false;
        }
        std::vector<double> dvals;
        dvals.reserve(rows);

        Trim(jstr);
        jstr = BuiltInFuncsUtils::LTrim(jstr, "[");
        Trim(jstr);
        jstr = BuiltInFuncsUtils::RTrim(jstr, "]");
        if (jstr.empty() || jstr == "]")
        {
            return EvaluatorInterface::allocateMatrix();
        }
        else if (rows == 1)
        {
            Currency elem;
            if (!JSON2Scalar(jstr, elem) || !elem.IsScalar())
            {
                return false;
            }
            cur = elem.Scalar();
            return true;
        }

        std::vector<std::string> rowvals (SplitRow(jstr));
        if (rows != static_cast<int>(rowvals.size()))
        {
            return false;  // This is a cell array
        }

        try
        {
            for (int j = 0; j < rows; ++j)
            {
                Currency elem;
                if (!JSON2Scalar(rowvals [j], elem) || !elem.IsScalar())
                {
                    return false;
                }
                dvals.emplace_back(elem.Scalar());
            }
            cur = Currency(dvals);
        }
        catch (const OML_Error&)
        {
            return false;  // Process as a cell
        }
        return true;
    }

    // 2D matrix
    std::unique_ptr<hwMatrix> mtx = nullptr;
    int m    = 0;
    while (!jstr.empty())
    {
        Trim(jstr);
        jstr = BuiltInFuncsUtils::LTrim(jstr, "[");
        Trim(jstr);
        if (jstr.empty() || m >= rows && jstr == "]")
        {
            break;
        }

        // Row
        size_t pos = jstr.find(']');
        if (pos == std::string::npos)
        {
            throw OML_Error(OML_ERR_INVALIDFORMAT, "missing closing bracket for matrix row " +
                            std::to_string(m + 1), 1);
        }
        std::string row (jstr.substr(0, pos));
        if (row.empty())
        {
            return false;  // This is a cell array
        }
        std::vector<std::string> rowvals (SplitRow(row));
        if (!mtx)
        {
            if (rowvals.empty())
            {
                return true;
            }
            cols = static_cast<int>(rowvals.size());
            mtx.reset(EvaluatorInterface::allocateMatrix(rows, cols, true));
        }
        else if (cols != static_cast<int>(rowvals.size()))
        {
            return false;  // This is a cell array
        }

        try
        {
            for (int j = 0; j < cols; ++j)
            {
                Currency elem;
                if (!JSON2Scalar(rowvals [j], elem) || !elem.IsScalar())
                {
                    return false;
                }
                (*mtx)(m, j) = elem.Scalar();
            }
        }
        catch (const OML_Error&)
        {
            return false;  // Process as a cell
        }

        ++m;
        jstr = jstr.substr(pos + 1);
        jstr = BuiltInFuncsUtils::LTrim(jstr, ",");
        jstr = BuiltInFuncsUtils::LTrim(jstr, "\n");
    }
   
    cur = mtx.release();
    return true;
}
//------------------------------------------------------------------------------
// Splits a string into a vector of tokens 
//------------------------------------------------------------------------------
std::vector<std::string> OmlJSON::SplitRow(const std::string& input) const
{
    std::string in(input);
    if (in.find(" ") != std::string::npos)
    {
        in = BuiltInFuncsUtils::LTrim(in, " ");
    }

    std::vector<std::string> vec;

    if (in.empty())
    {
        return vec;
    }

    std::string delim(",");
    {
        size_t pos = in.find(delim);
        if (pos == std::string::npos)
        {
            vec.emplace_back(in);
            return vec;
        }
    }

    while (!in.empty())
    {
        size_t pos = in.find(delim);
        if (pos == std::string::npos)
        {
            vec.emplace_back(in);
            break;
        }

        std::string token (in.substr(0, pos));
        if (!token.empty())
        {
            vec.emplace_back(token);
        }
        in = in.substr(pos + 1);
    }
    return vec;
}
//------------------------------------------------------------------------------
// Converts cell to JSON format string and returns it
//------------------------------------------------------------------------------
std::string OmlJSON::Cell2JSON(const Currency& in, int level) const
{
    HML_CELLARRAY* cell = in.CellArray();
    if (!cell || cell->Size() == 0)
    {
        return Indent(level - 1) + "[]";
    }

    int size = cell->Size();

    std::string indent0 (Indent(level - 1));

    std::string out (indent0 + "[" + _newline);
    for (int i = 0; i < size; ++i)
    {
        out += Currency2JSON((*cell)(i), level + 2)
            + Comma(i, size) + _newline;
    }
    out += indent0 + "]";
    return out;
}
//------------------------------------------------------------------------------
// Gets the indent string
//------------------------------------------------------------------------------
std::string OmlJSON::Indent(int level) const
{
    if (!_format || level < 1)
    {
        return std::string();
    }
    std::string indent;
    for (int i = 0; i < level; ++i)
    {
        indent += g_indent;
    }
    return indent;
}
//------------------------------------------------------------------------------
// Gets comma string, if needed
//------------------------------------------------------------------------------
std::string OmlJSON::Comma(int index, int max)
{
    return (index < max - 1) ? "," : "";
}
//------------------------------------------------------------------------------
// Trims space and new line characters from beginning and end
//------------------------------------------------------------------------------
void OmlJSON::Trim(std::string& in)
{
    // Need to trim for spaces both before and after checking for new lines
    in = BuiltInFuncsUtils::LTrim(in, " ");
    in = BuiltInFuncsUtils::LTrim(in, "\t");
    in = BuiltInFuncsUtils::LTrim(in, "\n");
    in = BuiltInFuncsUtils::LTrim(in, " ");
    in = BuiltInFuncsUtils::LTrim(in, "\t");

    in = BuiltInFuncsUtils::RTrim(in, " ");
    in = BuiltInFuncsUtils::RTrim(in, "\t");
    in = BuiltInFuncsUtils::RTrim(in, "\n");
    in = BuiltInFuncsUtils::RTrim(in, " ");
    in = BuiltInFuncsUtils::LTrim(in, "\t");
}
//------------------------------------------------------------------------------
// Returns true if conversion to scalar is valid
//------------------------------------------------------------------------------
bool OmlJSON::JSON2Scalar(const std::string& in, Currency& cur) const
{
    std::string jstr(in);
    Trim(jstr);

    if (jstr.empty())
    {
        cur = EvaluatorInterface::allocateMatrix();
        return true;
    }
    else if (IsNull(jstr))
    {
        cur = Currency(std::numeric_limits<double>::quiet_NaN());
        return true;
    }

    std::transform(jstr.begin(), jstr.end(), jstr.begin(), ::tolower);
    if (jstr == "true")
    {
        cur = Currency(true);
    }
    else if (jstr == "false")
    {
        cur = Currency(false);
    }
    else if (jstr == "nan")
    {
        cur = Currency(std::numeric_limits<double>::quiet_NaN());
    }
    else if (jstr == "infinity")
    {
        cur = Currency(std::numeric_limits<double>::infinity());
    }
    else if (jstr == "-infinity")
    {
        cur = Currency(-std::numeric_limits<double>::infinity());
    }
    else if (jstr.find_first_of("ij") != std::string::npos)
    {
        return false;
    }
    else
    {
        char first = jstr [0];
        if (first == '-' || first == '+' || isdigit(first))
        {
            std::string tmp (in);

            // Linux and VS2015 will not process numbers with double scientific 
            // notation which use 'D' instead of 'E'
            std::replace(tmp.begin(), tmp.end(), 'D', 'E');
            std::replace(tmp.begin(), tmp.end(), 'd', 'e');

            char* dummy = nullptr;
            double value = strtod(jstr.c_str(), &dummy);
            if (errno == ERANGE)
            {
                return false;
            }
            else if (!IsZero(value) || (dummy && std::string(dummy) != jstr) ||
                     (IsZero(value) && jstr.find_first_of("0") != std::string::npos))
            {
                cur = value;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns true if this is a null string
//------------------------------------------------------------------------------
bool OmlJSON::IsNull(const std::string& in)
{
    std::string lower(in);
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    return (lower == "null");
}
//------------------------------------------------------------------------------
// Increments current indent, if formatting
//------------------------------------------------------------------------------
std::string OmlJSON::IncrIndent(const std::string& indent) const
{
    return _format ? indent + g_indent : "";
}
//------------------------------------------------------------------------------
// Converts complex to JSON string
//------------------------------------------------------------------------------
std::string OmlJSON::Complex2JSON(const Currency& in, int level) const
{
    assert(in.IsComplex());

    std::string indent0 (Indent(level - 1));
    std::string indent1 (Indent(level));
    std::string indent2 (IncrIndent(indent1));

    std::string out (indent0 + '[');
    out += _newline + indent2 + Scalar2JSON(in.Complex().Real(), false) +
           _newline + indent1 + "]";
    return out;
}
//------------------------------------------------------------------------------
// Converts 2D struct array to JSON format string 
//------------------------------------------------------------------------------
bool OmlJSON::JSON2StructArray(const std::string& in, Currency& cur) const
{
    if (in.empty() || in.find("{") == std::string::npos)
    {
        return false; // Process as a cell
    }
    
    std::string jstr(in);
    Trim(jstr);
    size_t pos = jstr.find_first_not_of("[");
    if (pos == std::string::npos)
    {
        return false;
    }

    jstr = jstr.substr(pos);
    Trim(jstr);

    if (jstr.empty())
    {
        throw OML_Error(OML_ERR_INVALIDFORMAT, "missing ']'", 1);
    }
    else if (jstr[0] != '{')
    {
        return false;
    }

    if (jstr.back() == ']')
    {
        jstr.pop_back();
    }
    else
    {
        throw OML_Error(OML_ERR_INVALIDFORMAT, "missing ']'", 1);
    }
    Trim(jstr);

    std::vector <std::string> structsdata;

    int npenparen = 0;
    int ncloseparen = 0;

    std::string data;
    size_t len = jstr.length();
    for (size_t i = 0; i < len; ++i)
    {
        char ch = jstr [i];
        if (npenparen == 0 && ncloseparen == 0 && data.empty())
        {
            if (ch == ',' || ch == '\n')
            {
                continue;
            }
        }
        data += ch;
        if (ch == '}')
        {
            ncloseparen ++;
            if (npenparen == ncloseparen)
            {
                Trim(data);
                structsdata.emplace_back(data);
                npenparen = 0;
                ncloseparen = 0;
                data = "";
            }
        }
        else if (ch == '{')
        {
            npenparen ++;
        }
    }

    if (npenparen > ncloseparen)
    {
        throw OML_Error(OML_ERR_INVALIDFORMAT, "missing '}'", 1);
    }
    else if (npenparen < ncloseparen)
    {
        throw OML_Error(OML_ERR_INVALIDFORMAT, "missing '{'", 1);
    }
    else if (structsdata.empty())
    {
        return false;
    }

    int rows = static_cast<int>(structsdata.size());
    std::unique_ptr<StructData> sd (EvaluatorInterface::allocateStruct());
    sd->DimensionNew(rows, 1);

    std::vector<std::string> fields;
    for (int i = 0; i < rows; ++i)
    {
        std::map<std::string, std::vector<Currency> > mapdata;
        ParseStructData(structsdata[i], mapdata);
        if (mapdata.empty())
        {
            break;
        }
        if (i == 0)
        {
            fields.reserve(mapdata.size());
        }
        std::map<std::string, std::vector<Currency> >::const_iterator itr = mapdata.begin();
        for (; itr != mapdata.end(); ++itr)
        {
            std::string name (itr->first);
            if (name.empty())
            {
                return false;
            }
            else if (i == 0)
            {
                fields.emplace_back(name);
            }
            else if (std::find(fields.begin(), fields.end(), name) == fields.end())
            {
                return false;
            }
            std::vector<Currency> vals (itr->second);
            assert(vals.size() == 1);
            sd->SetValue(i, 0, name,
                     (vals.size() == 1) ? vals [0] : Currency());
        }
    }

    cur = sd.release();
    return true;
}
//------------------------------------------------------------------------------
// Gets value of struct fields
//------------------------------------------------------------------------------
std::string OmlJSON::GetFieldValueString(std::string& str) const
{
    size_t pos = str.find(',');
    if (pos == std::string::npos)
    {
        std::string data(str);
        str = "";
        return data;
    }

    std::string data (ContainerValue(str));
    return data;
}
//------------------------------------------------------------------------------
// Gets value between open brackets/parenthesis or quotes
//------------------------------------------------------------------------------
std::string OmlJSON::ContainerValue(std::string& in)
{
    Trim(in);
    if (in.empty() || in == "[]" || in == "{}")
    {
        return in;
    }
    else if (in[0] == ',')
    {
        in.erase(in.begin());
        Trim(in);
    }

    if (in.empty() || in == "[]" || in == "{}")
    {
        return in;
    }
    
    char first    = in [0];
    char last     = first;
    size_t len    = in.length();
    int    nopen  = 0;
    int    nclose = 0;

    std::string value;
    for (size_t i = 0; i < len; ++i)
    {
        char ch = in [i];

        if (ch == '\\' && i + 1 < len - 1)
        {
            int charval = in [i + 1];
            switch (charval)
            {
                case 34:                                    // Double quotes
                case 92:  value += charval;  i++; continue; // backslash
                case 110: value += char(10); i++; continue; // lower n -> \n
                case 116: value += char(9);  i++; continue; // lower t -> \t
                default:  break;
            }
        }
        value += ch;

        if (i == 0)
        {
            if (ch == '\"')
            {
                nopen ++;
                last = ch;
            }
            else if (ch == '[')
            {
                nopen ++;
                last = ']';
            }
            else if (ch == '{')
            {
                nopen ++;
                last = '}';
            }
            else
            {
                size_t pos = in.find(',');
                if (pos != std::string::npos)
                {
                    value = in.substr(0, pos);
                    in = in.substr(pos + 1);
                    return value;
                }
                value = in;
                in = "";
                return value;
            }
        }
        else if ((ch == '[' || ch == '{') && first == ch)
        {
            nopen ++;
        }
        else if (ch == last)
        {
            //if (ch == '\"' && i > 0 && in [i - 1] == '\\')
            //{
            //    if (!val.empty())
            //    {

            //    }
            //    continue; // Could be a quote in the string
            //}
            nclose ++;
            if (nopen == nclose || ch == '\"')
            {
                in = in.substr(i + 1);
                Trim(in);
                if (!in.empty() && in [0] == ',')
                {
                    in.erase(in.begin());
                }
                Trim(in);
                return value;
            }
        }
    }
    return "";
}
//------------------------------------------------------------------------------
// Converts JSON format string to a cell and returns it
//------------------------------------------------------------------------------
Currency OmlJSON::JSON2Cell(const std::string& in) const
{
    std::string jstr(in);
    Trim(jstr);
    if (!jstr.empty())
    {
        if (jstr [0] != '[')
        {
            throw OML_Error(OML_ERR_INVALIDFORMAT, "cannot convert to cell");
        }
        jstr.erase(jstr.begin());
        if (jstr.empty() || jstr[jstr.length() - 1] != ']')
        {
            throw OML_Error(OML_ERR_INVALIDFORMAT, "missing ]");
        }
        jstr.pop_back();
    }
    if (jstr.empty())
    {
        return EvaluatorInterface::allocateCellArray();
    }

    std::vector<Currency> children;
    children.reserve(10);

    while (!jstr.empty())
    {
        std::string data (ContainerValue(jstr));
        Currency    cur = JSON2Currency(data);
        children.emplace_back(cur);
    }

    if (children.empty())
    {
        return EvaluatorInterface::allocateCellArray();
    }

    int ncells = static_cast<int>(children.size());

    std::unique_ptr<HML_CELLARRAY> cell (
        EvaluatorInterface::allocateCellArray(ncells, 1));

    for (int i = 0; i < ncells; ++i)
    {
        (*cell)(i) = children [i];
    }

    return cell.release();
}
//------------------------------------------------------------------------------
// Creates the field name
//------------------------------------------------------------------------------
std::string OmlJSON::FieldName(const std::string& in) const
{
    std::string field(in);
    Trim(field);

    if (!_checkname)
    {
        return field;
    }
    
    BuiltInFuncsUtils::ReplaceAll(field, " ", "");
    if (field.empty())
    {
        return _prefix;
    }

    size_t len = field.length();
    bool addprefix = false;

    if (_style == REPLACESTYLE_UNDERSCORE)
    {
        for (size_t i = 0; i < len; ++i)
        {
            char ch = field [i];
            if (!isalnum(ch) && field [i] != '_')
            {
                field [i] = '_';
                if (i == 0)
                {
                    addprefix = true;
                }
            }
        }
    }
    else if (_style == REPLACESTYLE_DELETE)
    {
        std::string tmp ("");
        for (size_t i = 0; i < len; ++i)
        {
            char ch = field [i];
            if (ch == '_' || isalnum(ch))
            {
                tmp += ch;
            }
            else if (i == 0)
            {
                addprefix = true;
            }
        }
        field = tmp;
    }
    else // _style == REPLACESTYLE_HEX
    {
        std::string tmp("");
        tmp.reserve(len * 4);
        for (size_t i = 0; i < len; ++i)
        {
            char ch = field [i];
            if (ch == '_' || isalnum(ch))
            {
                tmp += ch;
            }
            else
            {
                if (i == 0)
                {
                    addprefix = true;
                }

                std::ostringstream os;
                os << std::hex << "0x" << std::setw(2) << std::setfill('0')
                    << std::uppercase << (0xff & (unsigned int) ch);
                tmp += os.str();
            }
        }
        field = tmp;
    }
    
    if (addprefix)
    {
        field = _prefix + field;
    }
    if (!field.empty() && isdigit(field [0]))
    {
        std::string tmp (field);

        // Linux and VS2015 will not process numbers with double scientific 
        // notation which use 'D' instead of 'E'
        std::replace(tmp.begin(), tmp.end(), 'D', 'E');
        std::replace(tmp.begin(), tmp.end(), 'd', 'e');

        char* cptr = nullptr;
        double value = strtod(tmp.c_str(), &cptr);
        if (errno == ERANGE)
        {
            return field;
        }
        std::string dummy (cptr ? cptr : "");
        if (!IsZero(value) || (cptr && dummy != tmp) ||
           (IsZero(value)  && tmp.find_first_of("0") != std::string::npos))
        {
            return _prefix + field;
        }
    }
    return field;
}
//------------------------------------------------------------------------------
// Returns true if given string is a valid prefix
//------------------------------------------------------------------------------
bool OmlJSON::IsValidPrefix(const std::string& in)
{
    std::string prefix(in);
    Trim(prefix);
    BuiltInFuncsUtils::ReplaceAll(prefix, " ", "");
    if (prefix.empty())
    {
        return false;
    }

    bool hasnonnumeric = false;
    size_t len = prefix.length();
    for (size_t i = 0; i < len; ++i)
    {
        char ch = in [i];
        if (ch == '_')
        {
            hasnonnumeric = true;
        }
        else if (!isalnum(ch))
        {
            return false;
        }
        else if (!isdigit(ch))
        {
            hasnonnumeric = true;
        }
    }
    return hasnonnumeric;
}
//------------------------------------------------------------------------------
// Converts JSON string to string currency
//------------------------------------------------------------------------------
Currency OmlJSON::JSON2String(const std::string& in) const
{
    assert (!in.empty());

    size_t len = in.length();
    if (len < 2 || in [len - 1] != '\"')
    {
        throw OML_Error(OML_ERR_INVALIDFORMAT, "missing '\"'", 1);
    }

    std::string value;
    for (size_t i = 1; i < len - 1; ++i)
    {
        char ch = in [i];
        if (i + 1 < len - 1 && ch == '\\')
        {
            int charval = in [i + 1];
            switch (charval)
            {
                case 34:                                    // Double quotes
                case 92:  value += charval;  i++; continue; // backslash
                case 110: value += char(10); i++; continue; // lower n -> \n
                case 116: value += char(9);  i++; continue; // lower t -> \t
                default:  break;
            }
        }
        value += ch;
    }
    return Currency(value);
}