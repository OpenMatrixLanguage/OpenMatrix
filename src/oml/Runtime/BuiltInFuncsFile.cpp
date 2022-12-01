/**
* @file BuiltInFuncsFile.cpp
* @date March 2016
* Copyright (C) 2016-2022 Altair Engineering, Inc.  
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

#include "BuiltInFuncsFile.h"

#include <algorithm>
#include <cassert>
#include <cctype>        // For std::isdigit
#include <cstdio>        // For std::rename
#include <iomanip>       // For std::setprecision
#include <iostream>
#include <limits>
#include <memory>        // For std::unique_ptr
#include <regex>
#include <sstream>

#ifdef OS_WIN
#    define NOMINMAX
#    include <Windows.h>
#    include <WinBase.h>
#else
#    include <ctype.h>
#    include <sys/stat.h>
#endif

#include "BuiltInFuncsUtils.h"
#include "Evaluator.h"
#include "MatrixDisplay.h"
#include "OML_Error.h"

#include "hwMatrix.h"

//------------------------------------------------------------------------------
// Returns true after reading a text file [textread]
//------------------------------------------------------------------------------
bool BuiltInFuncsFile::Textread(EvaluatorInterface           eval, 
	                            const std::vector<Currency>& inputs, 
			                    std::vector<Currency>&       outputs)
{
    size_t numargin = inputs.empty() ? 0 : inputs.size();
	if (numargin == 0)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    // First argument is file
    const Currency& fileinput = inputs[0];        
    std::string filename;

    if (fileinput.IsString())
    {
        filename = fileinput.StringVal();
    }
    else if (fileinput.IsPositiveInteger())
    {
        int fileId = static_cast<int>(fileinput.Scalar());
        filename = eval.GetFileName(fileId);
    }
    else
    {
        throw OML_Error(OML_ERR_STRING_FILESTREAM, 1);
    }

    filename = BuiltInFuncsUtils::Normpath(filename);
    if (!BuiltInFuncsUtils::FileExists(filename)) 
    {
		throw OML_Error(OML_ERR_FILE_NOTFOUND, 1);
    }
    BuiltInFuncsUtils utils;

    std::string mode("r");
#ifdef OS_WIN
    mode += ", ccs=UTF-8";
#endif
    
    std::FILE* f = utils.FileOpen(filename, mode);
    if (!f)
    {
		throw OML_Error(OML_ERR_FILE_CANNOTOPEN, 1);
    }

    // Get input options - which are the second arguments if specified
    // They could be format options specified as '%d %s', 'headerlines' or 'delimiters'
    std::vector<std::string> formats;
    bool        hasFormatSpec  = false;
    int         numheaderlines = 0;
    std::string delims;
    if (numargin > 1)
    {
        const Currency& in2 = inputs[1];
        if (!in2.IsString()) 
        {
            throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
        }

        std::string in2val (in2.StringVal());
        hasFormatSpec = (in2val != "headerlines" && in2val != "delimiter");

        if (hasFormatSpec)
        {
            bool validformat = true;
            formats = BuiltInFuncsUtils::GetFormats(in2val, "%d%f%g%s", validformat);
            if (!validformat)
            {
                std::ostringstream os;
                os << "Warning: invalid format specified in argument " << 2
                   << "; valid formats are %d, %f, %g and %s";
                BuiltInFuncsUtils::SetWarning(eval, os.str());
            }
        }

        size_t index = hasFormatSpec ? 2 : 1;
        for (; index < numargin; ++index)
        {
            const Currency& c   = inputs[index];
            std::string     val = c.IsString() ? c.StringVal() : "";
            if (val.empty() || (val != "headerlines" && val != "delimiter"))
            {
                std::string msg = "Error: invalid option specified in argument " 
                    + std::to_string(static_cast<long long>(index))
                    + "; valid options are 'headerlines' or 'delimiter'";
                throw OML_Error(msg);
            }

            if (index + 1 >= numargin) 
            {
                throw OML_Error(OML_ERR_NUMARGIN);
            }
            index += 1;  // Get the next value
            const Currency& cval = inputs[index];

            if (val == "headerlines")
            {
                if (!(cval.IsPositiveInteger() || (cval.IsScalar() && cval.Scalar() == 0)))
                {
                    throw OML_Error(OML_ERR_NATURALNUM, static_cast<int>(index + 1));
                }
                numheaderlines = static_cast<int>(cval.Scalar());
            }
            else if (val == "delimiter")
            {
                if (!cval.IsString()) 
                {
                    throw OML_Error(OML_ERR_STRING, static_cast<int>(index + 1));
                }
                delims = cval.StringVal();
            }
        }
    }

    // Set defaults
    if (formats.empty())
    {
        formats.push_back("%lf");
    }

    bool newlinedelim = false;
    if (!delims.empty())
    {
        size_t pos = delims.find("\\t");
        if (pos != std::string::npos)
        {
            delims.replace(pos, 2, "\t");
        }

        pos = delims.find("\\n");
        if (pos != std::string::npos)
        {
            if (delims == "\\n")
            {
                newlinedelim = true;
            }
            else
            {
                std::string tmp = delims.substr(0, pos);
                tmp += delims.substr(pos + 2);
                delims = tmp;
            }
        }
        pos = delims.find("\\r");
        if (pos != std::string::npos)
        {
            if (delims == "\\r")
            {
                newlinedelim = true;
            }
            else
            {
                std::string tmp = delims.substr(0, pos);
                tmp += delims.substr(pos + 2);
                delims = tmp;
            }
        }
    }
    if (delims.empty())
    {
        delims = " ";
    }

    // Pad formats to match outputs
    int numoutputs = eval.GetNargoutValue();
    int numformats = static_cast<int>(formats.size());
    if (numoutputs > numformats)
    {
        formats.reserve(numoutputs);
        while (numformats < numoutputs)
        {
            formats.push_back("%lf");
            numformats++;
        }
    }
    else if (numoutputs < 1)
    {
        numoutputs = 1;
    }

    // Check if there are tabs in the format
    int i = 0;
    for (std::vector<std::string>::iterator itr = formats.begin();
         itr != formats.end(); ++itr, ++i)
    {
        std::string fmt (*itr);
        size_t      pos     = fmt.find("\\t");
        bool        replace = false;
        if (pos != std::string::npos)
        {
            fmt.replace(pos, 2, "\t");
            replace = true;
        }

        pos = fmt.find("%g");
        if (pos == std::string::npos)
        {
            pos = fmt.find("%f");
        }
                
        if (pos != std::string::npos)
        {
            fmt.replace(pos, 2, "%lf");
            replace = true;
        }

        if (replace)
        {
            formats[i] = fmt;
        }
    }
    
    // Read the data
    int linenum = 0;      
    const int tmpbuffsize = 2048;

#ifdef OS_WIN
    wchar_t lineC [tmpbuffsize];
#else
    char lineC[tmpbuffsize];

#endif
    bool hasdata = false;
    while (!eval.IsInterrupt() && !feof(f))
    {
        std::string line;

#ifdef OS_WIN
        if (fgetws(lineC, tmpbuffsize, f) == nullptr)
        {
            break;
        }
        line = utils.WString2StdString(lineC);
#else
        if (fgets(lineC, tmpbuffsize, f) == nullptr)
        {
            break;
        }

        line = lineC;
#endif
		if (linenum < numheaderlines || line.empty())
        {
			linenum++;  // Skip the header/empty lines
			continue; 
		}

#ifdef OS_WIN
        if (!wcsstr(lineC, L"\n") && !wcsstr(lineC, L"\r"))
#else
        if (!strstr(lineC, "\n") && !strstr(lineC, "\r"))
#endif
        {
            while(!feof(f))
            {
#ifdef OS_WIN
                if (fgetws(lineC, tmpbuffsize, f) == nullptr)
#else
                if (fgets(lineC, tmpbuffsize, f) == nullptr)
#endif
                {
                    break;
                }
#ifdef OS_WIN
                line += utils.WString2StdString(lineC);
#else
                line += lineC;
#endif

#ifdef OS_WIN
                if (!wcsstr(lineC, L"\n") && !wcsstr(lineC, L"\r"))
#else
                if (!strstr(lineC, "\n") && !strstr(lineC, "\r"))
#endif
                {
                    break;
                }
            }
        }

        if (line[line.size() - 1] == '\n' || line[line.size() - 1] == '\r')
        {
            line.pop_back();
        }

#ifndef OS_WIN
        utils.StripTrailingNewline(line);
#endif

        if (line.empty())
        {
            linenum++;
            continue;
        }

        if (newlinedelim && formats[0].find("%s") != std::string::npos)
        {
            if (outputs.empty())
            {
                HML_CELLARRAY* cell = EvaluatorInterface::allocateCellArray(1, 1);
                (*cell)(0) = line;
                outputs.push_back(cell);
                hasdata = true;
            }
            else
            {
                assert(outputs[0].CellArray());
                int m = outputs[0].CellArray()->M();
                utils.CheckMathStatus(eval, outputs[0].CellArray()->Resize(m + 1, 1));
                (*(outputs[0].CellArray()))(m) = line;
                hasdata = true;
            }
            linenum++;
            continue;
        }
        int   column = 0;

        std::string tmpdelims = (newlinedelim) ? " " : delims;
        char* tok    = strtok((char *)line.c_str(), tmpdelims.c_str());

        while (tok)
        {
            std::string fmt (formats[column]);
            std::string data (tok);

            bool createcur = (outputs.empty() || outputs.size() < column + 1) ?
                             true : false;

            if (fmt.find("%s") != std::string::npos)
            {
                if (createcur)
                {
                    HML_CELLARRAY* cell = EvaluatorInterface::allocateCellArray(1, 1);
                    (*cell)(0) = data;
                    outputs.push_back(cell);
                }
                else
                {
                    int m = outputs[column].CellArray()->M();
                    utils.CheckMathStatus(eval, 
                        outputs[column].CellArray()->Resize(m + 1, 1));
                    (*outputs[column].CellArray())(m) = data;
                }
                hasdata = true;
            }
            else
            {
                bool hasformat = true;
                std::stringstream tokenizer (data);
                double            val = 0.0;
                if (fmt.find("%d") != std::string::npos)
                {
                    int output = 0;
                    tokenizer >> output;
                    if (tokenizer.fail())
                    {
                        hasformat = false;// This is an error
                    }
                    else
                    {
                        val = static_cast<double>(output);
                    }
                }
                else // Default is float format
                {
                    tokenizer >> val;
                    if (tokenizer.fail())
                    {
                        hasformat = false;// This is an error
                    }
                }

                if (hasformat)
                {
                    hasdata = true;
                }
                if (createcur)
                {
                    outputs.push_back(
                        EvaluatorInterface::allocateMatrix(1, 1, val));
                }
                else
                {
                    hwMatrix* mtx = outputs[column].GetWritableMatrix();
                    int m = mtx->M();
                    utils.CheckMathStatus(eval, mtx->Resize(m + 1, 1));
                    (*mtx)(m, 0) = val;
                }
            }

            tok = strtok(NULL, tmpdelims.c_str());
            column++;
            if (column >= numformats)
            {
                column = 0;
            }
        }
        linenum++;
    }
    
     fclose(f);

    // Create the outputs. 
    if (!hasdata || (!outputs.empty() && outputs.size() != formats.size()))
    {
        outputs.clear();
        if (eval.IsInterrupt())
        {
            return true;
        }        
        throw OML_Error(OML_ERR_FORMAT);
    }
      
    return true;
}
//------------------------------------------------------------------------------
// Returns true after writing a matrix to a file [dlmwrite]
//------------------------------------------------------------------------------
bool BuiltInFuncsFile::Dlmwrite(EvaluatorInterface           eval, 
	                            const std::vector<Currency>& inputs, 
			                    std::vector<Currency>&       outputs)
{
    int nargin = inputs.empty() ? 0 : static_cast<int>(inputs.size());
	if (nargin < 2) throw OML_Error(OML_ERR_NUMARGIN);

    // First argument is file
    const Currency& cur1 = inputs[0];     
    if (!cur1.IsString() && !cur1.IsPositiveInteger())
        throw OML_Error(OML_ERR_STRING_FILESTREAM, 1, OML_VAR_TYPE);

    std::FILE* fp = nullptr;
    std::string filename;
    bool     closefile = true;

    BuiltInFuncsUtils utils;
    if (cur1.IsString())
    {
        filename = cur1.StringVal();
        filename = BuiltInFuncsUtils::Normpath(filename);
    }
    else 
    {
        fp = eval.GetFile(static_cast<int>(cur1.Scalar()));
        if (!fp)
            throw OML_Error(OML_ERR_STRING_FILESTREAM, 1, OML_VAR_TYPE);
        closefile = false;
    }

    // Second argument is matrix
    const Currency& cur2 = inputs[1];
    if (!cur2.IsMatrix() && !cur2.IsScalar() && !cur2.IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_TYPE);

    // Subsequent arguments specify options for writing
    BuiltInFuncsFile fileFuncs;
    std::string rdelim("\n");
    std::string cdelim(" ");
    std::string precision;

    bool append     = false;
    bool hasroffset = false;
    int  roffset    = -1;
    int  coffset    = -1;

    // Process other arguments
    for (int i = 2; i < nargin; ++i)
    {
        int         value = 0;
        const Currency& cur = inputs[i];
        if (cur.IsInteger())
        {
            if (hasroffset)
                coffset = static_cast<int>(cur.Scalar());
            else
            {
                roffset    = static_cast<int>(cur.Scalar());
                hasroffset = true;
            }
            continue;
        }

        if (!cur.IsString())
            throw OML_Error(OML_ERR_STRING_INTEGER, i, OML_VAR_TYPE);

        std::string key = cur.StringVal();
        if (key == "-append")
            append = true;                // No value to process
        else if (key == "append")         // Needs to follow with on/off
        {
            std::string val = fileFuncs.GetStringValue(inputs, nargin, i);
            if (val == "on")
                append = true;
            else if (val == "off")
                append = false;
            else 
                throw OML_Error(OML_ERR_FUNCSWITCH, i + 1, OML_VAR_TYPE);
        }
        else if (key == "delimiter")    // Needs to follow with value
            cdelim = fileFuncs.GetStringValue(inputs, nargin, i);
        else if (key == "newline")
        {
            std::string val = fileFuncs.GetStringValue(inputs, nargin, i);
            if (val == "unix")
                rdelim = "\n";
            else if (val == "pc")
                rdelim = "\r\n";
            else if (val == "mac")
                rdelim = "\r";
            else 
                rdelim = val;
        }
        else if (key == "roffset")
            roffset = fileFuncs.GetIntegerValue(inputs, nargin, i);
        else if (key == "coffset")
            coffset = fileFuncs.GetIntegerValue(inputs, nargin, i);
        else if (key == "precision")
        {
            if (i + 1 >= nargin) throw OML_Error(OML_ERR_NUMARGIN);
            ++i;
                
            const Currency& val = inputs[i];
            if (val.IsString())
                precision = val.StringVal();
            else if (val.IsInteger())
            {
                int digits = static_cast<int>(val.Scalar());
                if (digits < 0) throw OML_Error(OML_ERR_STRING_NATURALNUM, i + 1, OML_VAR_TYPE);
                
                precision = "%." + std::to_string(static_cast<long long>(digits))
                            + "g";
            }
            else
                throw OML_Error(OML_ERR_STRING_NATURALNUM, i + 1, OML_VAR_TYPE);
        }
        else
            cdelim = key;
    }

    // \t, \r, \n are interpreted as 2 characters each in oml
    size_t pos = cdelim.find("\\t"); 
    if (pos != std::string::npos)
    {
        cdelim.replace(pos, 2, "\t");
    }
    pos = cdelim.find("\\n"); 
    if (pos != std::string::npos)
    {
        cdelim.replace(pos, 2, "\n");
    }
    pos = cdelim.find("\\r"); 
    if (pos != std::string::npos)
    {
        cdelim.replace(pos, 2, "\r");
    }

    std::string roffsetstr;
    if (roffset > 0)     // Number of delimiter only rows to add
    {
        const hwMatrix* mtx = cur2.Matrix();
        int   numrows       = mtx ? mtx->M() : 0;
        int   numcols       = mtx ? mtx->N() : 0;

        if (coffset <= 0)
        {
            for (int i = 0; i < roffset; ++i)
            {
                for (int j = 0; j < numcols; ++j)
                    roffsetstr += cdelim;
                roffsetstr += rdelim;
            }
        }
        else
        {
            for (int i = 0; i < roffset; ++i)
            {
                for (int j = 0; j < numcols + coffset - 1; ++j)
                    roffsetstr += cdelim;
                roffsetstr += rdelim;
            }

        }
    }

    if (!fp)
    {
        std::string mode = (append) ? "a" : "w";
        fp = utils.FileOpen(filename, mode);
        if (!fp)
        {
            throw (append) ? OML_Error(OML_ERR_FILE_CANNOTWRITE, 1) :
                OML_Error(OML_ERR_FILE_CANNOTOPEN, 1);
        }
    }

    fflush(fp);
    if (!roffsetstr.empty())
    {
        fwrite(roffsetstr.c_str(), sizeof(char), roffsetstr.size(), fp);
    }
    fflush(fp);


    if (precision.empty())
    {
        MatrixDisplay::WriteNonFormattedOutputValues(cur2, rdelim, cdelim, coffset, fp);
    }
    else
    {
        int skipformat = CurrencyDisplay::GetSkipFormat();
        CurrencyDisplay::SetSkipFormat(-1);
        std::string vals (MatrixDisplay::GetOutputValues(cur2,
            eval.GetOutputFormat(), rdelim, cdelim, precision, precision, coffset));
        CurrencyDisplay::SetSkipFormat(skipformat);

        fflush(fp);
        if (!vals.empty())
        {
            fwrite(vals.c_str(), sizeof(char), vals.size(), fp);
            if (!rdelim.empty())
                fwrite(rdelim.c_str(), sizeof(char), rdelim.size(), fp);
        }
        fflush(fp);
    }
    if (closefile)
    {
        fclose(fp);
    }
    return true;

}
//------------------------------------------------------------------------------
//! Returns string value of a currency, if applicable
//------------------------------------------------------------------------------
std::string BuiltInFuncsFile::GetStringValue(const std::vector<Currency>& inputs,
                                             int                          nargin,
                                             int&                         index) const
{
    assert(nargin == static_cast<int>(inputs.size()));

    if (index + 1 >= nargin) throw OML_Error(OML_ERR_NUMARGIN);
    ++index;
                
    if (!inputs[index].IsString())
    {
        throw OML_Error(OML_ERR_STRING, index + 1, OML_VAR_TYPE);
    }
                
    return inputs[index].StringVal();
}
//------------------------------------------------------------------------------
//! Returns integer value of a currency, if applicable
//------------------------------------------------------------------------------
int BuiltInFuncsFile::GetIntegerValue(const std::vector<Currency>& inputs,
                                      int                          nargin,
                                      int&                         index) const
{
    assert(nargin == static_cast<int>(inputs.size()));

    if (index + 1 >= nargin) throw OML_Error(OML_ERR_NUMARGIN);
    ++index;
                
    const Currency& cur = inputs[index];
    if (!cur.IsInteger()) throw OML_Error(OML_ERR_INTEGER, index, OML_VAR_TYPE);
                
    return static_cast<int>(cur.Scalar());
}
//------------------------------------------------------------------------------
// Returns true after copying files/directories [copyfile]
//------------------------------------------------------------------------------
bool BuiltInFuncsFile::Copyfile(EvaluatorInterface           eval, 
                                const std::vector<Currency>& inputs, 
                                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.empty() ? 0 : inputs.size();
    if (nargin < 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }
    std::string src (inputs[0].StringVal());
    if (src.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 1);
    }

    if (!inputs[1].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
    }
    std::string dst (inputs[1].StringVal());
    if (dst.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 2);
    }

    src = BuiltInFuncsUtils::Normpath(src);
    dst = BuiltInFuncsUtils::Normpath(dst);

    bool forcecopy = false;
    if (nargin > 2)
    {
        if (!inputs[2].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);
        }
        std::string opt(inputs[2].StringVal());
        if (!opt.empty())
        {
            std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);
        }
        if (opt != "f")
        {
            throw OML_Error(OML_ERR_OPTION, 3, OML_VAR_VALUE);
        }
        forcecopy = true;
    }

    bool issrcdir  = BuiltInFuncsUtils::IsDir(src);
    bool isdstcdir = BuiltInFuncsUtils::IsDir(dst);
    bool srxexists = BuiltInFuncsUtils::FileExists(src);
    
    // Additional checks if not force copy otherwise application will wait for
    // user input which should not happen
    if (!forcecopy && srxexists)     
    {        
        std::string dstfile;
        bool        dstexists = BuiltInFuncsUtils::FileExists(dst);

        // Check if destination is a directory, src is a file. In this case,
        // we will need to construct the path
        if (isdstcdir && !issrcdir)
        {
            dstfile   = dst + "/" + BuiltInFuncsUtils::GetBaseName(src);
            dstfile = BuiltInFuncsUtils::Normpath(dstfile);
            dstexists = BuiltInFuncsUtils::FileExists(dstfile);
        }
        
        if (dstexists)
        {
            std::string msg = "Error: cannot overwrite [";
            if (!dstfile.empty())
                msg += BuiltInFuncsUtils::Normpath(dstfile);
            else
                msg += BuiltInFuncsUtils::Normpath(dst);
            msg += "] in argument 2; use 'f' option to force copy";
            throw OML_Error(msg);
        }
    }

    std::string strcmd;
    int         returncode = 0;
#ifdef OS_WIN    
    BuiltInFuncsUtils utils;
    std::wstring str1(utils.StdString2WString(src));
    std::wstring wsrc(utils.GetNormpathW(str1));
    std::wstring str2(utils.StdString2WString(dst));
    std::wstring wdst(utils.GetNormpathW(str2));

    std::wstring wcmd;
    if (srxexists && !issrcdir && !wdst.empty() && !utils.DoesPathExistW(wdst))
    {
        wcmd = L"copy \"" + wsrc + L"\" \"" + wdst + L"\"";
    }
    else
    {
        wcmd = L"xcopy \"" + wsrc + L"\" \"" + wdst + L"\" /I /Q";
        if (issrcdir)
            wcmd += L" /E";
    }

    if (forcecopy)
        wcmd += L" /Y";
    else
        wcmd += L" /-Y";

    returncode = _wsystem(wcmd.c_str());
#else
    // Make sure that the underlying directories exist before copying
    std::string base1(BuiltInFuncsUtils::GetBaseDir(dst));
    if (base1.find(' ') != std::string::npos && base1.find('\"') == std::string::npos)
    {
        base1 = '\"' + base1 + '\"';
    }
    std::string mkdircmd = "mkdir -p " + base1;
    system(mkdircmd.c_str());

    strcmd = "cp -R ";
    if (!forcecopy)
        strcmd += "-n ";

    if (src.find(' ') != std::string::npos && src.find('\"') == std::string::npos)
    {
        src = '\"' + src + '\"';
    }

    if (dst.find(' ') != std::string::npos && dst.find('\"') == std::string::npos)
    {
        dst = '\"' + dst + '\"';
    }


    strcmd +=  src + " " + dst;

    returncode = system(strcmd.c_str());

#endif

    std::string msg;
    int         msgid = 0;
    if (returncode != 0)
    {
        msgid = (errno != 0) ? errno : returncode;

        std::string err (strerror(errno));
        if (!err.empty() && err != "No error")
            msg = err;
        
        if (!msg.empty())
            msg += "\n";

#ifdef OS_WIN
        msg += "Copy failed from [" + utils.WString2StdString(wsrc);
        msg += "] to [" + utils.WString2StdString(wdst) + "]";

#else
        msg += "Copy failed from [" + BuiltInFuncsUtils::Normpath(src);
        msg += "] to [" + BuiltInFuncsUtils::Normpath(dst) + "]";
#endif
        outputs.push_back(0);
        outputs.push_back(msg);
        outputs.push_back(msgid);
    }
    else
    {
        outputs.push_back(1);
        outputs.push_back("");
        outputs.push_back("");
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns true after moving files/directories [movefile]
//------------------------------------------------------------------------------
// #ifdef _DEBUG
#if 1
#    define DEBUG_PRINT(s) { std::cout << s << std::endl; }
#else
#    define DEBUG_PRINT(s)
#endif // 0
bool BuiltInFuncsFile::Movefile(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
{
    std::string debugoutput;    // VSM-6719

    size_t nargin = inputs.empty() ? 0 : inputs.size();
    if (nargin < 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }
    std::string src(inputs[0].StringVal());
    if (src.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 1);
    }
    bool srcwildcard;
    if (src.find_first_of("?*") != std::string::npos)
    {
        srcwildcard = true;
    }
    else
    {
        srcwildcard = false;
    }
    bool srcincludespath;
    std::size_t src_path = src.find_last_of("/\\");
    std::string src_object;  // this will be the src name, stripping any path elements
    if (src_path != std::string::npos)
    {
        // source argument includes path delimiters.  
        // This saves the last name in the path which needs to be checked for non force moves.
        srcincludespath = true;
        src_object = src.substr(src_path + 1, string::npos);
    }
    else
    {
        srcincludespath = false;
        src_object = src;
    }
        




    if (!inputs[1].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
    }
    std::string dst(inputs[1].StringVal());
    if (dst.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 2);
    }

    bool forcecopy = false;
    bool testdisplay = false;
    if (nargin > 2)
    {
        if (!inputs[2].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);
        }
        std::string opt(inputs[2].StringVal());
        if (!opt.empty())
        {
            std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);
        }
        // if (opt != "f")
        // {
        //    throw OML_Error(OML_ERR_OPTION, 3, OML_VAR_VALUE);
        // }
        if (opt.find_first_of('f') != std::string::npos)
        {
            forcecopy = true;
        }
        if (opt.find_first_of('t') != std::string::npos)
        {
            testdisplay = true;
        }
        if (testdisplay)
        {
            if(forcecopy)
                DEBUG_PRINT("FORCECOPY!")
            else
                DEBUG_PRINT("NO forcecopy")
            if(srcwildcard)
                DEBUG_PRINT("WILDCARD!")
            else
                DEBUG_PRINT("NO wildcard")
        }
        // forcecopy = true;
    }

    bool issrcdir = BuiltInFuncsUtils::IsDir(src);
    bool isdstdir = BuiltInFuncsUtils::IsDir(dst);
    bool srxexists = BuiltInFuncsUtils::FileExists(src);
    if (testdisplay)
    {
        debugoutput = "movefile '" + src + "' " + (issrcdir ? " IsDirectory " : " Not Directory");
        DEBUG_PRINT(debugoutput);
        debugoutput = "movefile '" + dst + "' " + (isdstdir ? " IsDirectory " : " Not Directory");
        DEBUG_PRINT(debugoutput);
        debugoutput = "movefile '" + src + "' " + (srxexists ? " Source Exists " : " Source NOT Exists");
        DEBUG_PRINT(debugoutput);
    }


#if 1
    // Additional checks if not force copy otherwise application will wait for
    // user input which should not happen
    if (!forcecopy && srxexists)
    {
        std::string dstfile;
        bool        dstexists = BuiltInFuncsUtils::FileExists(dst);


        // Check if destination is a directory, src is a file. 
        if (isdstdir && !issrcdir)
        {
            dstfile = dst + "/" + BuiltInFuncsUtils::GetBaseName(src);
            dstexists = BuiltInFuncsUtils::FileExists(dstfile);
        }

        if (dstexists && !isdstdir)
        {
            std::string msg = "Error: cannot overwrite [";
            if (!dstfile.empty())
                msg += BuiltInFuncsUtils::Normpath(dstfile);
            else
                msg += BuiltInFuncsUtils::Normpath(dst);
            msg += "] in argument 2; use 'f' option to force move";
            // throw OML_Error(msg);
            outputs.push_back(0);
            outputs.push_back(msg);
            outputs.push_back(0);
        }
    }
#endif 
    std::string strcmd;
    int         returncode = 0;
    bool        runcommand = true;  //  this will be turned off when the system command is run.


#ifdef OS_WIN    
    std::wstring tdisplay1, tdisplay2;
    if (testdisplay)
    {
        tdisplay1 = L"";
        tdisplay2 = L"";
    }
    else
    {
        // suppress output from the wsystem call.
        tdisplay1 = L"@";
        tdisplay2 = L" > nul ";
    }

    BuiltInFuncsUtils utils;
    std::wstring str1(utils.StdString2WString(src));
    std::wstring wsrc(utils.GetNormpathW(str1));
    std::wstring str2(utils.StdString2WString(dst));
    std::wstring wdst(utils.GetNormpathW(str2));
    std::wstring wsrc_object(utils.StdString2WString(src_object));

    std::wstring wcmd, wcmd1, wcmd2;
    // if (srxexists && !issrcdir && !wdst.empty() && !utils.DoesPathExistW(wdst))
    if(!srcwildcard)
    {
        //Windows no wildcard

        if (forcecopy)
        {
            wcmd = tdisplay1 + L"move /Y \"" + wsrc + L"\" \"" + wdst + L"\"" + tdisplay2;
            if(testdisplay)
            {
                std::string command_str(utils.WString2StdString(wcmd));
                debugoutput = "BuiltInFuncsFile::Movefile force " + command_str;
                DEBUG_PRINT(debugoutput);
            }
        }
        else
        {
            if (isdstdir)
            {
                if (srcincludespath)
                {
                    wcmd = tdisplay1 + L"if not exist \"" + wdst + L"\\" + wsrc_object + L"\" move \"" + wsrc + L"\" \"" + wdst + L"\"" + tdisplay2;
                }
                else
                {
                    wcmd = tdisplay1 + L"if not exist \"" + wdst + L"\\" + wsrc + L"\" move \"" + wsrc + L"\" \"" + wdst + L"\"" + tdisplay2;
                }
            }
            else
            {
                wcmd = tdisplay1 + L"if not exist \"" + wdst + L"\" move \"" + wsrc + L"\" \"" + wdst + L"\"" + tdisplay2;
            }
            if (testdisplay)
            {
                std::string command_str(utils.WString2StdString(wcmd));
                debugoutput = "BuiltInFuncsFile::Movefile no force " + command_str;
                DEBUG_PRINT(debugoutput);
            }
        }
 
    }
    else
    {
        // Windows Wildcard

        // FILES AND DIRECTORIES 

        if (forcecopy)
        {
            // wcmd = tdisplay1 + L"for %a in (\"" + wsrc + L"\") do move /Y %a \"" + wdst + L"\\%~nxa\" " + tdisplay2 +
            //     L"& " + tdisplay1 + L"for /D %a in (\"" + wsrc + L"\") do move /Y %a \"" + wdst + L"\\%~nxa\"" + tdisplay2;
            wcmd1 = tdisplay1 + L"for %a in (\"" + wsrc + L"\") do " + tdisplay1 + L"move /Y \"%a\" \"" + wdst + L"\\%~nxa\" " + tdisplay2;
            int returncode1 = _wsystem(wcmd1.c_str());
            wcmd2 = tdisplay1 + L"for /D %a in (\"" + wsrc + L"\") do " + tdisplay1 + L"move /Y \"%a\" \"" + wdst + L"\\%~nxa\"" + tdisplay2;
            int returncode2 = _wsystem(wcmd2.c_str());
            runcommand = false;  // This more complex move requires a pair of sysetm commands.  Do not need to run the standard system call
            returncode = returncode1 || returncode2;
        }
        else
        {
            // When forcecopy is not enabled, we will never overwrite an existing file.
            // For wildcard processing using system call with command for loop.
            // Complexity here is a bit increased because separate "for" loops are requred to process files and directories.
            // wcmd = tdisplay1 + L"for %a in (\"" + wsrc + L"\") do " + tdisplay1 + L"if not exist \"" + wdst + L"\\%~nxa\" move %a \"" + wdst + L"\\%~nxa\" " + tdisplay2 +
            //    L"& " + tdisplay1 + L"for /D %a in (\"" + wsrc + L"\") do " + tdisplay1 + L"if not exist \"" + wdst + L"\\%~nxa\" move %a \"" + wdst + L"\\%~nxa\"" + tdisplay2;
            wcmd1 = tdisplay1 + L"for %a in (\"" + wsrc + L"\") do " + tdisplay1 + L"if not exist \"" + wdst + L"\\%~nxa\" move \"%a\" \"" + wdst + L"\\%~nxa\" " + tdisplay2;
            int returncode1 = _wsystem(wcmd1.c_str());
            wcmd2 = tdisplay1 + L"for /D %a in (\"" + wsrc + L"\") do " + tdisplay1 + L"if not exist \"" + wdst + L"\\%~nxa\" move \"%a\" \"" + wdst + L"\\%~nxa\"" + tdisplay2;
            int returncode2 = _wsystem(wcmd2.c_str());
            runcommand = false;  // This more complex move requires a pair of sysetm commands.  Do not need to run the standard system call
            returncode = returncode1 || returncode2;
        }

        if (testdisplay)
        {
            // std::string command_str3(wcmd.begin(), wcmd.end());
            // debugoutput = "BuiltInFuncsFile::Movefile 2 Windows move command: '" + command_str3 + "'";
            // DEBUG_PRINT(debugoutput);
        }
        
        // End Windows wildcard
    }
    if( runcommand )  // run the standard system call if system calls have not previous been run.
        returncode = _wsystem(wcmd.c_str());
// End of Windows Block
#else
//Linux Block
    strcmd = "mv ";
    if (!forcecopy)
        strcmd += "-n ";  // --no-clobber

    if (src.find(' ') != std::string::npos && src.find('\"') == std::string::npos)
    {
        src = '\"' + src + '\"';
    }

    if (dst.find(' ') != std::string::npos && dst.find('\"') == std::string::npos)
    {
        dst = '\"' + dst + '\"';
    }

    src = BuiltInFuncsUtils::Normpath(src);
    dst = BuiltInFuncsUtils::Normpath(dst);

    strcmd += src + " " + dst;

    if (testdisplay)
    {
        debugoutput = "BuiltInFuncsFile::Movefile  Linux MOVE command: '" + strcmd + "'";
        DEBUG_PRINT(debugoutput);
    }
    returncode = system(strcmd.c_str());
// End of Linux Block
#endif

    std::string msg;
    int         msgid = 0;
    if (returncode != 0)
    {
        msgid = (errno != 0) ? errno : returncode;

        std::string err(strerror(errno));
        if (!err.empty() && err != "No error")
            msg = err;

        if (!msg.empty())
            msg += "\n";

#ifdef OS_WIN
        msg += "Move failed from [" + utils.WString2StdString(wsrc);
        msg += "] to [" + utils.WString2StdString(wdst) + "]";

#else
        msg += "Move failed from [" + BuiltInFuncsUtils::Normpath(src);
        msg += "] to [" + BuiltInFuncsUtils::Normpath(dst) + "]";
#endif
        outputs.push_back(0);
        outputs.push_back(msg);
        outputs.push_back(msgid);
    }
    else
    {
        outputs.push_back(1);
        outputs.push_back("");
        outputs.push_back("");
    }
    return true;
}
//------------------------------------------------------------------------------
//! Returns true after importing data from files [importdata]
//------------------------------------------------------------------------------
bool BuiltInFuncsFile::Importdata(EvaluatorInterface           eval, 
	                              const std::vector<Currency>& inputs, 
			                      std::vector<Currency>&       outputs)
{
    int nargin = (!inputs.empty()) ? static_cast<int>(inputs.size()) : 0;
    if (nargin < 1) 
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }

    std::string fname (inputs[0].StringVal());
    fname = BuiltInFuncsUtils::Normpath(fname);
    if (!BuiltInFuncsUtils::FileExists(fname))
    {
        throw OML_Error(OML_ERR_FILE_NOTFOUND, 1);
    }

    if (BuiltInFuncsUtils::IsDir(fname))
    {
        throw OML_Error(OML_Error(OML_ERR_FILE_NOTFOUND).GetErrorMessage() +
            "; argument 1 is a directory");
    }

    std::string err   ("Error: invalid input in argument 1; ");
    std::string npath ("[" + BuiltInFuncsUtils::Normpath(fname) + "]");


    std::string tboxerr ("unsupported file format " + npath + ".\nTo import data, add toolbox ");
    std::string ext     (BuiltInFuncsUtils::GetFileExtension(fname));
    if (!ext.empty())
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    int nargout = eval.GetNargoutValue();

    if (ext == "mat") // Use mat file binary reader
    {
        FunctionInfo* fi      = nullptr;
        FUNCPTR       funcptr = nullptr;
        eval.FindFunctionByName("load", &fi, &funcptr, NULL);
        if (funcptr)
        {
            std::vector<Currency> funcInputs(inputs);

            outputs = eval.DoMultiReturnFunctionCall(funcptr, funcInputs, 
                        nargin, nargout, true);
            return true;
        }
        throw OML_Error(err + tboxerr + "[omlMatio]");
    }
   
    bool istxtfile     = (ext == "csv" || ext == "txt");
    FUNCPTR xlsfuncptr = nullptr;

#ifdef OS_WIN
    if (IsExtXlsCompatible(ext))
    {
        // Try reading non-textfiles using xlsread
        FunctionInfo* fi = nullptr;
        eval.FindFunctionByName("xlsread", &fi, &xlsfuncptr, NULL);
        if (!istxtfile)
        {
            if (xlsfuncptr)
            {
                std::vector<Currency> funcInputs(inputs);
                outputs = eval.DoMultiReturnFunctionCall(xlsfuncptr, funcInputs, 
                          nargin, nargout, true);
                return true;
            }
            throw OML_Error(err + tboxerr + "[omlxlstoolbox]");
        }
    }
#endif
    
    if (nargin > 1 && inputs[nargin-1].IsString())
    {
        std::string opt (inputs[nargin-1].StringVal());
        if (!opt.empty())
        {
            std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);
            if (opt == "-ascii")
            {
                FunctionInfo* fi      = nullptr;
                FUNCPTR       funcptr = nullptr;
                eval.FindFunctionByName("load", &fi, &funcptr, NULL);
                if (funcptr)
                {
                    std::vector<Currency> funcInputs(inputs);
                    outputs = eval.DoMultiReturnFunctionCall(funcptr, funcInputs, 
                                  nargin, nargout, true);
                    return true;
                }
                throw OML_Error(err + tboxerr + "[omlMatio]");
            }
        }
    }
    
    try
    {
        // Use dlmread if there are no input arguments as it will generate a 
        // matrix with the correct columns. Textread without input arguments will
        // generate a matrix with all values in 1 column
        if (nargin == 1)
        {
            Dlmread(eval, inputs, outputs);
        }
        else
        {
            Textread(eval, inputs, outputs);
        }
        return true;
    }
    catch (OML_Error e)
    {
        if (!istxtfile)
            throw (e);
    }

#ifdef OS_WIN
    // Try excel for text files
    if (istxtfile && xlsfuncptr)
    {
        std::vector<Currency> funcInputs(inputs);
        outputs = eval.DoMultiReturnFunctionCall(xlsfuncptr, funcInputs, nargin,
            nargout, true);
    }
    else
    {
        throw OML_Error(err + "unsupported file format " + npath);
    }
#else
    throw OML_Error(err + "unsupported file format " + npath);
#endif
    
    return true;
}
//------------------------------------------------------------------------------
// Returns true if given file extension is Excel compatible
//------------------------------------------------------------------------------
bool BuiltInFuncsFile::IsExtXlsCompatible(const std::string& ext)
{
    if (ext.empty())
        return false;

    if (ext == "xls" || ext == "xlsx" || ext == "csv" || ext == "xlsm" ||
        ext == "xml" || ext == "xlsb" || ext == "slk" || ext == "xltm" ||
        ext == "xlt" || ext == "txt"  || ext == "dif" || ext == "ods" )
        return true;

    return false;
}
//------------------------------------------------------------------------------
// Returns true if given filename is an Ascii file
//------------------------------------------------------------------------------
bool BuiltInFuncsFile::IsAsciiFile(const std::string& input)
{
    if (input.empty())
        return false;

    std::string name(BuiltInFuncsUtils::Normpath(input));
    FILE *fp = fopen(name.c_str(), "r");
    if (!fp)
        return false;

    // Reads upto newline character to check if there are any non-ascii chars
    bool bAscii = true;
    while (bAscii)
    {
        char c = fgetc(fp);
        if (c == EOF || c == '\n' || c == '\r')
            break;

        bAscii = isascii(c);
    }

    fclose(fp);
    return bAscii;
}
//------------------------------------------------------------------------------
// Returns true after reading a text file/string (textscan command)
//------------------------------------------------------------------------------
bool BuiltInFuncsFile::Textscan(EvaluatorInterface           eval, 
	                            const std::vector<Currency>& inputs, 
			                    std::vector<Currency>&       outputs)
{
	if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    // First argument is file
    bool isString = inputs[0].IsString();
    if (!isString && !inputs[0].IsPositiveInteger())
    {
        throw OML_Error(OML_ERR_STRING_FILESTREAM, 1);
    }
    std::string textstr = (isString) ? inputs[0].StringVal() : "";
    
    int numFormatApplied = 0;
    int fileid = -1;
    std::FILE*        f = nullptr;
    BuiltInFuncsUtils utils;
    BuiltInFuncsFile  funcs;
    if (!isString)
    {
        fileid = utils.GetFileId(eval, inputs[0], 1);
        utils.CheckFileIndex(eval, fileid, 1, false); // Don't read std streams
        f = eval.GetFile(fileid);

        if (!f)
        {
            if (fileid != -1)
            {
                eval.CloseFile(fileid);
            }
            throw OML_Error(OML_ERR_FILE_CANNOTOPEN, 1);
        }
    }

    // Get input options - which are the second arguments if specified
    // They could be format options specified as '%d %s', 'headerlines' or 'delimiters'
    std::vector<std::string> basefmts;
    std::vector<std::string> rawfmts;
    std::vector<bool>        usefmts;
    int         numheaderlines = 0;
    int         repeat         = -1;
    size_t      index          = 1;
    bool        hasRepeat      = false;
    bool        returnOnError  = true;
    std::string delims;
    size_t      numargin = inputs.size();
    std::string whitespace(" ");
    std::string endofline;

    std::string prefix;
    if (numargin > 1)
    {
        const Currency& in2 = inputs[1];
        if (in2.IsInteger())  // Format is not specified, this is repeat param
        {
            repeat = static_cast<int>(in2.Scalar());
            if (repeat < 0 && repeat != -1)
            {
                std::string msg ("Error: invalid option specified in argument 1");
                msg += " for repeat parameter; valid options are -1, 0 or ";
                msg += " positive integers.";
                throw OML_Error(msg);
            }
            index ++;
            hasRepeat = true;
        }
        else if (!in2.IsString()) 
        {
            throw OML_Error(OML_ERR_STRING, 2);
        }
        else
        {
            std::string in2val (in2.StringVal());  // Specifies format
            std::string err;
            if (in2val.find("%") != std::string::npos)
            {
                bool valid = funcs.GetFormats(in2val, basefmts, rawfmts,
                    usefmts, err, prefix);
                if (!valid)
                {
                    if (err.empty())
                    {
                        throw OML_Error(OML_ERR_FORMAT, 2);
                    }
                    throw OML_Error(err + " in argument 2");
                }
                index ++;
            }
        }
    
        for (; index < numargin; ++index)
        {
            const Currency& c = inputs[index];
            if (index == 2 && c.IsInteger() && !hasRepeat)
            {
                repeat = static_cast<int>(c.Scalar());
                if (repeat < 0 && repeat != -1)
                {
                    std::string msg ("Error: invalid option specified in argument 1");
                    msg += " for repeat parameter; valid options are -1, 0 or ";
                    msg += " positive integers.";
                    throw OML_Error(msg);
                }
                continue;
            }
            std::string val = c.IsString() ? c.StringVal() : "";
            if (val.empty())
            {
                throw OML_Error(OML_ERR_OPTION, static_cast<int>(index + 1));
            }
            if (index + 1 >= numargin) 
            {
                throw OML_Error(OML_ERR_NUMARGIN);
            }

            int      cIdx = static_cast<int>(index + 1);  // Get the value
            const Currency& cval = inputs[cIdx];

            std::transform(val.begin(), val.end(), val.begin(), ::tolower);
            if (val == "headerlines")
            {
                if (!(cval.IsPositiveInteger() || (cval.IsScalar() && cval.Scalar() == 0)))
                {
                    throw OML_Error(OML_ERR_NATURALNUM, cIdx + 1);
                }
                numheaderlines = static_cast<int>(cval.Scalar());
            }
            else if (val == "delimiter")
            {
                if (!cval.IsString()) 
                {
                    throw OML_Error(OML_ERR_STRING, cIdx + 1);
                }
                delims = cval.StringVal();
            }
            else if (val == "returnonerror")
            {
                if (!cval.IsLogical() && !cval.IsInteger())  
                {
                    throw OML_Error(OML_ERR_LOGICAL, cIdx + 1);
                }
                int tmp = static_cast<int>(cval.Scalar());
                if (tmp != 0 && tmp != 1)
                {
                    throw OML_Error(OML_ERR_FLAG_01, cIdx + 1);
                }
                returnOnError = (tmp == 0) ? false : true;
            }
            else if (val == "whitespace")
            {
                if (!cval.IsString())
                {
                    throw OML_Error(OML_ERR_STRING, cIdx + 1);
                }
                whitespace = cval.StringVal();
            }
            else if (val == "endofline")
            {
                if (!cval.IsString())
                {
                    throw OML_Error(OML_ERR_STRING, cIdx + 1);
                }
                endofline = cval.StringVal();
                if (endofline == "\r"  || endofline == "\n"  || endofline == "\r\n" ||
                    endofline == "\\r" || endofline == "\\n" || endofline == "\\r\\n")
                {
                    endofline = ""; // This case is handled by default
                }
                else if (!endofline.empty() && endofline.length() > 1)
                {
                    throw OML_Error(
                        "Error: invalid value specified for 'endofline'; must be a single character, '\\r' or '\\r\\n'");
                }
            }
            else
            {
                throw OML_Error(OML_ERR_OPTION, static_cast<int>(index + 1));
            }
            index = cIdx;
        }
    }

    if (basefmts.empty())
    {
        basefmts.push_back("%f");
        rawfmts.push_back("%lf");  // Default option
        usefmts.push_back(true);
    }
    assert(basefmts.size() == rawfmts.size());
    assert(basefmts.size() == usefmts.size());
    int numformats = static_cast<int>(basefmts.size());

    bool newlinedelim = false;
    if (!delims.empty())
    {
        size_t pos = delims.find("\\t");
        if (pos != std::string::npos)
        {
            delims.replace(pos, 2, "\t");
        }

        if (endofline.empty())
        {
            pos = delims.find("\\n");
            if (pos != std::string::npos)
            {
                if (delims == "\\n")
                {
                    newlinedelim = true;
                }
                else
                {
                    std::string tmp = delims.substr(0, pos);
                    tmp += delims.substr(pos + 2);
                    delims = tmp;
                }
            }
            pos = delims.find("\\r");
            if (pos != std::string::npos)
            {
                if (delims == "\\r")
                {
                    newlinedelim = true;
                }
                else
                {
                    std::string tmp = delims.substr(0, pos);
                    tmp += delims.substr(pos + 2);
                    delims = tmp;
                }
            }
        }
        else if (delims.find(endofline) != std::string::npos)
        {
            newlinedelim = true;
        }
    }

    bool excludenewline = false;
    if (numformats == 1)
    {
        std::string tmp(rawfmts[0]);

        if (!endofline.empty())
        {
            std::string tpl("[^" + endofline + "]");
            excludenewline = (tmp == tpl);
        }
        else
        {
            if (tmp == "[^\\n]" || tmp == "[^\\n\\r]" || tmp == "[^\\r]")
            {
                excludenewline = true;
            }
        }
    }

    if (delims.empty())
    {
        delims = " ";
    }
    int numcols =0;
    for (int i = 0; i < numformats; ++i)
    {
        if (usefmts[i])
        {
            numcols++;
        }
    }
    bool hasStringFmt = (std::find(basefmts.begin(), basefmts.end(), "%s") != basefmts.end());
    if (!hasStringFmt && whitespace.empty())
    {
        whitespace = " ";
    }
    else if (hasStringFmt && !whitespace.empty() && whitespace.find(" ") == std::string::npos)
    {
        whitespace += " ";
    }

    // Set the output cell elements with the correct type
    std::unique_ptr<HML_CELLARRAY> cell (EvaluatorInterface::allocateCellArray(
                                            1, numcols));
    for (int i = 0, m = 0; i < numformats && m < numcols; ++i)
    {
        if (!usefmts[i])
        {
            continue;
        }
        std::string fmt (basefmts[i]);
        if (fmt == "%s" || fmt == "%c" || fmt == "custom")
        {
            (*cell)(m) = EvaluatorInterface::allocateCellArray();
        }
        else
        {
            (*cell)(m) = EvaluatorInterface::allocateMatrix();
        }
        ++m;
    }

    if (repeat == 0)   // No data is read
    {
        outputs.push_back(cell.release());  
        return true;
    }
    
    // Read the data
    bool quitLoop = false;
    int  linenum  = 0;   
    char lineC [2048];

    bool resetcol = true;
    int   cellcol       = 0;
    int   column        = 0;
    int   outputrow     = 0;

    while (!quitLoop && !eval.IsInterrupt())
    {
        std::string line;
        if (!f)
        {
            line = textstr;
            quitLoop = true;  // Don't need to loop through if text is processed
        }
        else                  // Reading file
        {
            line = GetLine(f, endofline); 
            if (feof(f))
            {
                quitLoop = true;
            }
        }

		if (linenum < numheaderlines || line.empty())
        {
			linenum++;  // Skip the header/empty lines
            if (quitLoop)
            {
                break;
            }
			continue; 
		}

        if (line.empty())
        {
            linenum++;
            if (quitLoop)
            {
                break;
            }
            continue;
        }

        if (excludenewline || (newlinedelim && basefmts[0] == "%s"))
        {
            if (usefmts[0])
            {
                int m = 0;
                if ((*cell)(0).CellArray()->Size() == 0)
                {
                    (*cell)(0) = EvaluatorInterface::allocateCellArray(1, 1);
                }
                else
                {
                    m = (*cell)(0).CellArray()->M();
                    outputrow = m;
                    utils.CheckMathStatus(eval, 
                        (*cell)(0).CellArray()->Resize(m + 1, 1));
                }
                line = utils.RTrim(line, whitespace); // Trim whitespace
                (*(*cell)(0).CellArray())(m, 0) = line;
            }
            linenum++;
            if (quitLoop)
            {
                break;
            }
            continue;
        }

        // Start of a new line
        if (column >= numformats)
        {
            numFormatApplied++;
            if (repeat > 0 && numFormatApplied >= repeat)
            {
                quitLoop = true;
                break;
            }
        }
        if (resetcol || column >= numformats || cellcol >= numcols)
        {
            cellcol       = 0;
            column        = 0;
        }
        else
        {
            resetcol = true;
        }
        int   columnWithErr = -1;

        // Handle prefix
        std::string line2(line);
        if (!prefix.empty())
        {
            size_t pos2 = line2.find(prefix);
            if (pos2 == 0)
            {
                line2 = line2.substr(prefix.length());
            }
        }

        char* tok = (!line2.empty()) ?
            strtok((char*)line2.c_str(), delims.c_str()) : nullptr;
        std::string strToProcess;

        while (tok || !strToProcess.empty())
        {
            std::string data = (strToProcess.empty()) ? tok : strToProcess;
            data = utils.RTrim(data, whitespace); // Trim whitespace

            if (data.empty())
            {
                quitLoop = true;
                break;
            }

            if (column >= numformats)
            {
                column = 0;
                numFormatApplied++;
                if (repeat > 0 && numFormatApplied >= repeat)
                {
                    quitLoop = true;
                    break;
                }
            }
            if (column == 0)
            {
                data = prefix + data;
            }

            bool use = usefmts[column];
            std::string basefmt (basefmts[column]);
            std::string rawfmt  (rawfmts[column]);
            std::string origfmt (rawfmt);

            assert(!basefmt.empty());
            assert(!rawfmt.empty());

            int result = 0;
            int used   = -1;
            if (basefmt != "custom")
            {
                rawfmt += "%n";  // To give the number of characters still left
            }

            if (basefmt == "%s")
            {
				memset(lineC, 0, sizeof(lineC));
                result = sscanf(data.c_str(), rawfmt.c_str(), lineC, &used);

                std::string out (lineC);
                size_t datalen = data.size();
                size_t numused = (used > 0) ? static_cast<size_t>(used) : 0;
                if (numused > 0 && numused < datalen &&
                    (delims.find(" ")  == std::string::npos ||
                        delims.find("\t") == std::string::npos))
                {
                    char ch = data[numused - 1];
                    if (isspace(ch))
                    {
                        out += ch;
                    }

                    // sscanf will split on whitespaces or tabs, which may not be
                    // in this delimiter list
                    for (size_t i = numused; i < datalen;)
                    {
                        char ch = data[i];
                        if (isspace(ch))
                        {
                            out += ch;
                            ++i;
                            continue;
                        }
                        std::string newstr (data.substr(i));
                        if (newstr.empty())
                        {
                            break;
                        }
                        memset(lineC, 0, sizeof(lineC));
                        int tmp = 0;
                        int newresult = sscanf(newstr.c_str(), rawfmt.c_str(),
                            lineC, &tmp);
                        if (newresult <= 0)
                        {
                            break;
                        }
                        out += lineC;
                        i   += strlen(lineC);
                    }
                }
                if (result > 0)
                {
                    if (use)
                    {
                        int m = 0;
                        if ((*cell)(cellcol).CellArray()->Size() == 0)
                        {
                            (*cell)(cellcol) = EvaluatorInterface::allocateCellArray(1, 1);
                        }
                        else
                        {
                            m = (*cell)(cellcol).CellArray()->M();
                            outputrow = m;
                            utils.CheckMathStatus(eval,
                                (*cell)(cellcol).CellArray()->Resize(m + 1, 1));
                        }
                        (*(*cell)(cellcol).CellArray())(m) = out;
                        cellcol++;
                    }
                    used = static_cast<int>(datalen);
                }
                else
                {
                    columnWithErr = cellcol;
                }
            }
            else if (basefmt == "%c")
            {
                if (strToProcess.empty())
                {
                    strToProcess = data;
                }

                memset(lineC, 0, sizeof(lineC));
                result = sscanf(data.c_str(), rawfmt.c_str(), lineC, &used);
                
                if (result > 0)
                {
                    if (use)
                    {
                        int m = 0;
                        if ((*cell)(cellcol).CellArray()->Size() == 0)
                        {
                            (*cell)(cellcol) = EvaluatorInterface::allocateCellArray(1, 1);
                        }
                        else
                        {
                            m = (*cell)(cellcol).CellArray()->M();
                            outputrow = m;
                            utils.CheckMathStatus(eval,
                                (*cell)(cellcol).CellArray()->Resize(m + 1, 1));
                        }
                        (*(*cell)(cellcol).CellArray())(m) = std::string(lineC);
                        cellcol++;
                    }
                }
                else
                {
                    columnWithErr = cellcol;
                }
            }
            else if (basefmt == "custom")
            {
                std::string matchstr;
                bool foundregexmatch = false;
                bool newlinematch = (rawfmt == "[\\n]" || rawfmt == "[\\r\\n]" || rawfmt == "[\\r]");
                if (newlinematch)
                {
                    foundregexmatch = true;
                }
                else
                {
                    std::regex  re(rawfmt);

                    std::sregex_iterator next(data.begin(), data.end(), re);
                    std::sregex_iterator end = std::sregex_iterator();

                    while (next != end)
                    {
                        foundregexmatch = true;

                        std::smatch match = *next;

                        std::string pre(match.prefix());
                        if (!pre.empty())
                        {
                            strToProcess = match.str();
                            strToProcess += match.suffix();
                            break;
                        }

                        matchstr += match.str();
                        next++;

                        std::string post(match.suffix());
                        if (post.empty())
                        {
                            matchstr = strToProcess;
                            strToProcess = "";
                        }
                    }
                }
                    

                if (foundregexmatch)
                {
                    if (use)
                    {
                        int m = 0;
                        if ((*cell)(cellcol).CellArray()->Size() == 0)
                        {
                            (*cell)(cellcol) = EvaluatorInterface::allocateCellArray(1, 1);
                        }
                        else
                        {
                            m = (*cell)(cellcol).CellArray()->M();
                            outputrow = m;
                            utils.CheckMathStatus(eval,
                                (*cell)(cellcol).CellArray()->Resize(m + 1, 1));
                        }
                        (*(*cell)(cellcol).CellArray())(m) = matchstr;
                        cellcol++;
                    }

                    if (newlinematch)
                    {
                        strToProcess = data;
                        if (numformats == 1)
                        {
                            tok = strtok(nullptr, delims.c_str());
                        }
                            
                        column++;
                        if (column >= numformats)  // End of format. Reset column number
                        {
                            column = 0;
                            numFormatApplied++;
                            if (repeat > 0 && numFormatApplied >= repeat)
                            {
                                quitLoop = true;
                                break;
                            }

                            if (!f)
                            {
                                linenum++;
                                if (repeat > 0 && linenum >= repeat)
                                {
                                    quitLoop = true;
                                    break;
                                }
                            }
                        }
                        if (cellcol >= numcols)
                        {
                            cellcol = 0;
                        }
                        continue;
                    }
                }
                else
                {
                    strToProcess += data;
                    if (numformats == 1)
                    {
                        tok = strtok(nullptr, delims.c_str());
                        if (!tok)
                        {
                            // Last line
                            if (!f)
                            {
                                linenum++;
                                if (repeat > 0 && linenum >= repeat)
                                {
                                    quitLoop = true;
                                    break;
                                }
                            }
                        }
                        continue;
                    }
                }
            }
            else
            {
                double val = std::numeric_limits<double>::quiet_NaN();
                if (basefmt == "%d" || basefmt == "%n")
                {
                    int output = 0;
                    result = sscanf(data.c_str(), rawfmt.c_str(), &output, &used);
                    if (result <= 0)
                    {
                        columnWithErr = cellcol;
                    }
                    else
                    {                        
                        val = static_cast<double>(output);
                        if (origfmt == basefmt && used < static_cast<int>(data.size()))  // Just plain number
                        {
                            used = static_cast<int>(data.size());
                        }
                    }
                }
                else // Default is float format
                {
                    std::string tmpfmt (rawfmt);
                    bool   hasspec = (origfmt.find("%lf") == std::string::npos);
                    size_t pos     = 0;
                    if (hasspec)
                    {
                        // Additional processing needed
                        pos = origfmt.find("lf");
                        if (pos != std::string::npos)
                        {
                            tmpfmt = "%";
                            tmpfmt += origfmt.substr(pos) + "%n";
                        }
                        else
                        {
                            tmpfmt = "%lf%n";
                        }
                    }
                    result = sscanf(data.c_str(), tmpfmt.c_str(), &val, &used);
                    columnWithErr = (result <= 0) ? cellcol : -1;

                    if (columnWithErr == -1 && hasspec)
                    {
                        // Update the value based on the specifiers
                        std::string specstr = (pos != std::string::npos) ?
                                                origfmt.substr(0, pos) : "";
                        specstr += "lf";
                        char tmp[128];
                        sprintf(tmp, specstr.c_str(), val);
                        val = atof(tmp);
                    }
                }
                if (columnWithErr != -1)
                {
                    if (!returnOnError)
                    {
                        std::string msg(OML_Error(OML_ERR_FORMAT).GetErrorMessage());
                        msg += " in field ";
                        msg += std::to_string(static_cast<long long>(columnWithErr + 1));
                        if (basefmt != "custom")
                        {
                            msg += " [" + basefmt + "]";
                        }
                        msg += " to row [";
                        msg += std::to_string(static_cast<long long>(outputrow + 1));
                        msg += "] of the output";
                        throw OML_Error(msg);
                    }
                    quitLoop = true;
                    break;
                }

                if (use)
                {
                    if ((*cell)(cellcol).Matrix()->Size() == 0)
                    {
                        hwMatrix* mtx = EvaluatorInterface::allocateMatrix(1, 1, val);
                        (*cell)(cellcol) = mtx;
                    }
                    else
                    {
                        hwMatrix* mtx = (*cell)(cellcol).GetWritableMatrix();
                        int       m   = mtx->M();
                        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(m + 1, 1));
                        (*mtx)(m, 0) = val;
                    }
                    cellcol++;
                }
            }

            if (basefmt == "custom" && !strToProcess.empty())
            {

            }
            else if (!strToProcess.empty() && used > 0 && 
                used < static_cast<int>(data.size()))
            {
                strToProcess = data.substr(used);
            }
            else
            {
                if (used > 0 && used < static_cast<int>(data.size()) &&
                    strToProcess.empty())
                {
                    data = data.substr(used);
                    strToProcess = data;
                }
                else
                {
                    tok = strtok(nullptr, delims.c_str());
                    strToProcess = "";  // Reset
                }
            }

            column++;
            if (column >= numformats)  // End of format. Reset column number
            {
                column = 0;  
                numFormatApplied++;
                if (repeat > 0 && numFormatApplied >= repeat)
                {
                    quitLoop = true;
                    break;
                }

                if (!f)
                {
                    linenum++;
                    if (repeat > 0 && linenum >= repeat)
                    {
                        quitLoop = true;
                        break;
                    }
                }
            }
            if (cellcol >= numcols)
            {
                cellcol = 0;
            }
        }
        if (f)
        {
            linenum++;
            if (repeat > 0 && linenum >= repeat)
            {
                quitLoop = true;
                break;
            }
        }

        // Check if there is an newline in the next format, which may need
        // to be excluded
        if (!tok)  // End of parsing the line
        {
            if (strstr(rawfmts[column].c_str(), "\\n"))
            {
                if (!usefmts[column])
                {
                    column++;
                    resetcol = false;
                }
                else
                {
                    int m = 0;
                    if ((*cell)(cellcol).CellArray()->Size() == 0)
                    {
                        (*cell)(cellcol) = EvaluatorInterface::allocateCellArray(1, 1);
                    }
                    else
                    {
                        m = (*cell)(cellcol).CellArray()->M();
                        utils.CheckMathStatus(eval,
                            (*cell)(cellcol).CellArray()->Resize(m + 1, 1));
                    }
                    (*(*cell)(cellcol).CellArray())(m) = "";
                    cellcol++;
                }
            }
        }
        else
        {
            resetcol = true;
        }
    } 

    outputs.push_back(cell.release());
    
    return true;
}
//------------------------------------------------------------------------------
// Returns true after reading formatted input from file/string (fscanf)
//------------------------------------------------------------------------------
bool BuiltInFuncsFile::Fscanf(EvaluatorInterface           eval,
	                          const std::vector<Currency>& inputs, 
		                      std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    // First argument specifies the file
    const Currency& cur1 = inputs[0];
    if (!cur1.IsString() && !cur1.IsPositiveInteger()) 
    {
        throw OML_Error(OML_ERR_STRING_FILESTREAM, 1);
    }
    BuiltInFuncsUtils utils;
    int               fileid = utils.GetFileId(eval, cur1, 1);
    utils.CheckFileIndex(eval, fileid, 1, false); // Don't read std streams
    std::FILE* file = eval.GetFile(fileid);

    int numargin = static_cast<int>(inputs.size());
    if (numargin < 2) 
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    // Second argument gives format descriptions separated by spaces. 
    const Currency& cur2 = inputs[1];
    if (!cur2.IsString()) 
    {
        throw OML_Error(OML_ERR_STRING, 2);
    }

    // Third argument, if available, specifies size for output.
    double rows      = 0;
    double cols      = 0;
    bool hasSizeMtx  = false;
    bool hasSizeSpec = (numargin > 2);
	if (numargin > 2)
    {
        utils.GetSizeSpecifications(inputs[2], 3, rows, cols, hasSizeMtx);
    }
    std::string fmtdesc(cur2.StringVal());

    // Get the basefmt and the rawfmts
    std::vector<std::string> basefmts;
    std::vector<std::string> rawfmts;
    std::string err;
    bool validformat = utils.GetFormats(fmtdesc, basefmts, rawfmts, err);
    if (!validformat)
    {
        if (err.empty())
        {
            err = "Error: invalid format specified";
        }
        throw OML_Error(err + " in argument 2");
    }

    // Check if string mask needs to be added and add default format
    bool addStringMask = false;
    if (basefmts.empty())
    {
        basefmts.push_back("%f");
        rawfmts.push_back("%lf");  // Default option
    }
    else if (basefmts.size() == 1)
    {
        addStringMask = (std::find(basefmts.begin(), basefmts.end(), "%s") != 
                            basefmts.end());
    }
    assert(basefmts.size() == rawfmts.size());

    // Read formatted input from the file
    // For each format, read in one value from input and push into matrix
    int sizelimit = static_cast<int>(hasSizeMtx ? rows * cols : rows);
    std::unique_ptr<hwMatrix> mtx(EvaluatorInterface::allocateMatrix(
                                    1, 1, true));
 
    int  row       = 0;
	int  count     = 0;
    bool keepgoing = true;
    char buff [1028];

    // Scan the file
    BuiltInFuncsFile funcs;
    while (keepgoing)
    {
        int i = 0;
        for (std::vector<std::string>::const_iterator itr = basefmts.begin();
                itr != basefmts.end(); ++itr, ++i)
	    {
            std::string basefmt (*itr);
            std::string rawfmt  (rawfmts[i]);
            std::string origfmt (rawfmt);

            assert(!basefmt.empty());
            assert(!rawfmt.empty());

            int result = 0;
            int used   = -1;
            rawfmt += "%n";  // To give the number of characters still left

            bool   iscomplex = false;
            double complex   = 0.0;

            bool   isDbl  = false;
            double dblVal = 0.0;

            bool        isStr = false;
            std::string strVal;

            bool ignoreval = (rawfmt.find("%*") != std::string::npos);

            if (basefmt == "%f" || basefmt == "%g")
	        {
                std::string tmpfmt (rawfmt);
                bool   hasspec = (origfmt.find("%lf") == std::string::npos);
                size_t pos     = 0;
                if (hasspec)
                {
                    // Additional processing needed
                    pos = origfmt.find("lf");
                    if (pos != std::string::npos)
                    {
                        tmpfmt = "%";
                        tmpfmt += origfmt.substr(pos) + "%n";
                    }
                    else
                    {
                        tmpfmt = "%lf%n";
                    }
                }
                result = funcs.ReadDouble(file, tmpfmt, true, iscomplex, dblVal, complex);
                if (ignoreval)
                {
                    continue;
                }
                if (result > 0 && hasspec)
                {
                    // Update the value based on the specifiers
                    std::string specstr = (pos != std::string::npos) ?
                                            origfmt.substr(0, pos) : "";
                    specstr += "lf";
                    char tmp[128];
                    memset(tmp, 0, sizeof(tmp));
                    sprintf(tmp, specstr.c_str(), dblVal);
                    dblVal = atof(tmp);
                }
                isDbl  = true;
                isDbl  = true;
            }
            else if (basefmt == "%d" || basefmt == "%n")
            {
                result = funcs.ReadDouble(file, rawfmt, false, iscomplex, dblVal, complex);
                if (ignoreval)
                {
                    continue;
                }
                isDbl  = true;
            }
            else if (basefmt == "%s")            
            {
                memset(buff, 0, sizeof(buff));
		        result = fscanf(file, rawfmt.c_str(), buff, &used);
                if (ignoreval)
                {
                    continue;
                }
                isStr  = true;
                strVal = buff;
            }
            else if (basefmt == "%c")
            {
                char ch;
                result = fscanf(file, rawfmt.c_str(), &ch, &used);
                if (ignoreval)
                {
                    continue;
                }
                dblVal = static_cast<double>(ch);
                isDbl  = true;
            }

            if (result != 1) 
            {
                keepgoing = false;
                break;
            }

            int newrow = (!isStr) ? 
                            row + 1 : row + static_cast<int>(strVal.size());
            if (iscomplex)
            {
                newrow += 1;
            }

            hwMathStatus stat = mtx->Resize(newrow, 1);
            utils.CheckMathStatus(eval, stat);
            if (isDbl)
            {
                (*mtx)(row) = dblVal;
                row++;
                count++;

                if (iscomplex)  // There is a complex value
                {
                    (*mtx)(row) = complex;
                    row++;     // Don't increase the count as this is part of the number
                }
            }
            else
            {
                size_t len = strVal.size();
                for (size_t i = 0; i < len; ++i)
                {
                    (*mtx)(row) = strVal[i];
                    row++;
                }
                count++;
            }
            if (count == sizelimit)
            {
                keepgoing = false;
            }
        }
	}

    int numvals = mtx->Size();

    // Split the output values vector into the correct size (if any) Right now, 
    // the result is Nx1 (which is the default)
    // Handles cases where inf is passed. Don't cast and use values.size() as 
    // compiler will convert
	int    numrows      = numvals;
	int    numcols      = 1;
	double rawnumcols   = 1.0;
	
    if (hasSizeSpec)                // Has size specifications
    {
        if (hasSizeMtx)             // Has matrix
        {
            numrows    = static_cast<int>(rows);
            rawnumcols = cols;      // This can be infinity
        }
        else if (addStringMask)
        {
            rawnumcols   = numrows;
            numrows      = static_cast<int>((IsInf_T(rows) ? numvals : rows));
        }
        else
        {
			rawnumcols = 1;
            numrows    = static_cast<int>(IsInf_T(rows) ? numvals : rows);
			if (numrows > numvals)  // Don't cast here when comparing
			    numrows = numvals;
        }
    }

	if (IsInf_T(rawnumcols) || IsNegInf_T(rawnumcols) || IsNaN_T(rawnumcols))
    {
        int tmp = (numrows == 0) ? 1 : numrows;
		numcols = (int)numvals / tmp;
    }
	else
		numcols = static_cast<int>(rawnumcols);

	if (numvals < numrows * numcols)
	{
        int tmp = (numrows == 0) ? 1 : numrows;
		if (numvals % tmp == 0)
			numcols = (int)numvals / tmp;
	}

	hwMathStatus stat = mtx->Reshape(numrows, numcols);
    if (stat.IsOk())
    {
        Currency out(mtx.release());
	    if (addStringMask)
        {
		    out.SetMask(Currency::MASK_STRING);
        }   
        outputs.push_back(out);
        return true;
    }

    // There is a mismatch in matrix dimensions, so copy data manually and pad 
    // the matrix
    hwMatrix* tmp = EvaluatorInterface::allocateMatrix(numrows, numcols, 0.0);
    int       tmpSize = tmp->Size();
    for (int i = 0; i < tmpSize && i < numvals; ++i)
    {
        (*tmp)(i) = (*mtx)(i);
    }

    Currency out(tmp);
	if (addStringMask)
    {
		out.SetMask(Currency::MASK_STRING);
    }   

    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// Returns true after renaming a file [rename]
//------------------------------------------------------------------------------
bool BuiltInFuncsFile::Rename(EvaluatorInterface           eval,
	                          const std::vector<Currency>& inputs, 
		                      std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }
    std::string oldname(inputs[0].StringVal());

    if (!inputs[1].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
    }
    std::string newname (inputs[1].StringVal());
    newname = BuiltInFuncsUtils::Normpath(newname);

    int result = 0;
#ifdef OS_WIN
    BuiltInFuncsUtils utils;
    std::wstring wold(utils.StdString2WString(oldname));
    std::wstring wnew(utils.StdString2WString(newname));
    
    result = _wrename(wold.c_str(), wnew.c_str());
#else
    result = std::rename(oldname.c_str(), newname.c_str());
#endif
    outputs.push_back(result);

    std::string msg = (result == 0) ? "" : std::string(strerror(errno));
    outputs.push_back(msg);

    return true;
}
//------------------------------------------------------------------------------
// Returns true if this is a complex number along with the value
//------------------------------------------------------------------------------
bool BuiltInFuncsFile::IsComplex(const std::string& in, double& val)
{
    if (in.empty())
    {
        return false;
    }

    std::string tmp (in);

    // Linux and VS2015 will not process numbers with double scientific 
    // notation which use 'D' instead of 'E'
    std::replace(tmp.begin(), tmp.end(), 'D', 'E');
    std::replace(tmp.begin(), tmp.end(), 'd', 'e');

    char* dummy = nullptr;

	val = strtod(tmp.c_str(), &dummy);
	if (dummy && (*dummy == 'i' || *dummy == 'j' || *dummy == 'I' || 
                  *dummy == 'J'))
    {
        // Check if the character is the last one in the string all the rest
        // are numbers
        size_t len = tmp.length();
        for (size_t i = 0; i < len - 1; ++i)
        {
            char ch = tmp[i];
            if (ch == '.' || isdigit(ch) || (i == 0 && (ch == '+' || ch == '-')))
            {
                continue;
            }
            return false;
        }
        return true;
    }
    return false;
}
//------------------------------------------------------------------------------
// Returns result and scans number from a file
//------------------------------------------------------------------------------
int BuiltInFuncsFile::ReadDouble(FILE*              file, 
                                 const std::string& fmt,
                                 bool               isdouble,
                                 bool&              iscomplex,
                                 double&            realval,
                                 double&            complexval) 
{
    int ival   = 0;
    int used   = 0;
    int result = (isdouble) ? fscanf(file, fmt.c_str(), &realval, &used) :
                              fscanf(file, fmt.c_str(), &ival,    &used);
    if (result != 1)
    {
        return result;
    }

    if (!isdouble)
    {
        realval = static_cast<double>(ival);
    }

    // Check if this is a complex number with/without spaces that is being read
    long pos1 = ftell(file);
    if (pos1 == -1)
    {
        return result;
    }

    // Read the next string to determine if we are processing a complex number
    char buff [1028];
    memset(buff, 0, sizeof(buff));
    int tmpresult = fscanf(file, "%s", buff);
    if (tmpresult != 1)
    {
        fseek(file, pos1, SEEK_SET);   // Rewind as we are not dealing with complex num
        return result;
    }
    long pos2 = ftell(file);
    std::string tmp1(buff);

    if (tmp1.empty())
    {
        fseek(file, pos1, SEEK_SET);   // Rewind as this is not a complex num
        return result;
    }

    bool issignonly = (tmp1 == "+" || tmp1 == "-") ? true : false;
    if (issignonly)  
    {
        // Since this is sign only, check if there is a number after
        memset(buff, 0, sizeof(buff));
        tmpresult = fscanf(file, "%s", buff);
        if (tmpresult != 1)
        {
            fseek(file, pos1, SEEK_SET);   // Rewind as we are not dealing with complex num
            return result;
        }
        tmp1 += buff;
    }

    // Check if there is an i at the end
    if (tmp1.empty())
    {
        fseek(file, pos1, SEEK_SET);   // Rewind as we are not dealing with complex num
        return result;
    }

    // Check if there is an imaginary part
    iscomplex = IsComplex(tmp1, complexval);
    if (iscomplex)
    {
        return result;
    }

    // Check if this has i - which means the first value is a complex number
    std::string lower (tmp1);
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "i" || lower == "j")
    {
        fseek(file, pos2 - 2, SEEK_SET);  // Check if there is a space before i
        char ch = fgetc(file);
        if (!isdigit(ch) || pos2 != -1)
        {
            fseek(file, pos1, SEEK_SET);
        }
        else
        {
            fseek(file, pos2, SEEK_SET);
        }
        return result;
    }

    // Check if there is an imaginary part
    iscomplex = IsComplex(tmp1, complexval);
    if (!iscomplex)
    {
        fseek(file, pos1, SEEK_SET);   // Rewind as this is not a complex num
    }
    return result;
}
//------------------------------------------------------------------------------
// Returns true and reads float/integer data from a file (dlmread)
//------------------------------------------------------------------------------
bool BuiltInFuncsFile::Dlmread(EvaluatorInterface           eval, 
	                           const std::vector<Currency>& inputs, 
			                   std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    // First argument is file
    bool isString = inputs[0].IsString();
    if (!isString && !inputs[0].IsPositiveInteger())
    {
        throw OML_Error(OML_ERR_FILE_FILESTREAM, 1, OML_VAR_TYPE);
    }

    BuiltInFuncsUtils utils;
    std::string       filename;
    std::wstring      wfile;
    int               fileid = -1;
    std::string       ext;

    if (!isString)
    {
        fileid = utils.GetFileId(eval, inputs[0], 1);
        utils.CheckFileIndex(eval, fileid, 1, false); // Don't read std streams

        filename = eval.GetFileName(fileid);
        filename = utils.Normpath(filename);
        wfile = utils.StdString2WString(filename);
    }
    else
    {
        bool exists = false;
#ifdef OS_WIN
        wfile = utils.StdString2WString(inputs[0].StringVal());
        wfile = utils.GetNormpathW(wfile);
        exists = utils.DoesPathExistW(wfile);
        filename = (exists) ? utils.WString2StdString(wfile) : "";
        ext = utils.WString2StdString(utils.GetFileExtensionW(wfile));
#else
        filename = inputs[0].StringVal();
        filename = utils.Normpath(filename);
        exists = utils.DoesPathExist(filename);
#endif
        if (!exists)
        {
            throw OML_Error(OML_ERR_FILE_NOTFOUND, 1);
        }

        ext = utils.GetFileExtension(filename);
    }

    if (IsExtXlsCompatible(ext))
    {
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext != "csv" && ext != "txt")
        {
            std::string err("Error: invalid input in argument 1; cannot read file; ");
            err += "Excel files need to be imported as .csv";
#ifdef OS_WIN
            err += " or read with xlsread command";
#endif
            throw OML_Error(err);
        }
    }

    std::unique_ptr<hwMatrix> mtx = nullptr;

    double remptyval  = 0;
    double iemptyval  = 0;
    bool   hascomplex = false;

    std::string delim;
    std::string strrange;

    bool ignoreMultidelim = false;

    int startrow = -1;
    int startcol = -1;
    int endrow   = -1;
    int endcol   = -1;
    int nargin   = static_cast<int>(inputs.size());
    for (int i = 1; i < nargin; ++i)
    {
        const Currency& cur = inputs[i];
        if (cur.IsString())  // Could be separator/range/emptyval
        {
            std::string val(cur.StringVal());
            if (!val.empty())
            {
                std::string lower(val);
                std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                if (lower == "emptyvalue")
                {
                    ++i;
                    if (i >= nargin)
                    {
                        throw OML_Error(OML_ERR_NUMARGIN);
                    }
                    if (inputs[i].IsScalar())
                    {
                        remptyval = inputs[i].Scalar();
                    }
                    else if (inputs[i].IsComplex())
                    {
                        hwComplex c = inputs[i].Complex();
                        remptyval   = c.Real();
                        iemptyval   = c.Imag();
                        hascomplex  = true;
                    }
                    else
                    {
                        throw OML_Error(OML_ERR_SCALAR_COMPLEX, i + 1, OML_VAR_TYPE);
                    }
                    continue;
                }
                else if (lower == "ignoremultidelim")
                {
                    ++i;
                    if (i >= nargin)
                    {
                        throw OML_Error(OML_ERR_NUMARGIN);
                    }
                    if (inputs[i].IsScalar())
                    {
                        int val = static_cast<int>(inputs[i].Scalar());
                        ignoreMultidelim = (val == 1);
                    }
                    else if (inputs[i].IsString())
                    {
                        std::string tmp(inputs[i].StringVal());
                        if (!tmp.empty())
                        {
                            if (tmp == "on")
                            {
                                ignoreMultidelim = 1;
                            }
                            else if (tmp != "off")
                            {
                                throw OML_Error(OML_ERR_FUNCSWITCH, i + 1);
                            }
                        }
                        else
                        {
                            throw OML_Error(OML_ERR_FUNCSWITCH, i + 1);
                        }

                    }
                    else
                    {
                        throw OML_Error(OML_ERR_FUNCSWITCH, i + 1);
                    }
                    continue;
                }

                else if (lower == "linestoskip" || lower == "headerlines")
                {
                    throw OML_Error(OML_ERR_BAD_STRING, i + 1);
                }
                bool hasalpha = false;
                bool hassep = false;
                if (strstr(val.c_str(), ":") || strstr(val.c_str(), ".."))
                {
                    hassep = true;
                }
                size_t strsize = val.size();
                for (size_t j = 0; j < strsize; ++j)
                {
                    if (isalpha(val[j]))
                    {
                        hasalpha = true;
                        break;
                    }
                }
                if (hasalpha && hassep)
                {
                    BuiltInFuncsFile funcs;
                    std::vector<int> range = funcs.String2VectorRange(val, i);
                    assert(range.size() == 4);
                    startrow = range[0];
                    startcol = range[1];
                    endrow = range[2];
                    endcol = range[3];
                }
                else
                {
                    delim = val;
                }
            }
        }
        else if (cur.IsInteger())
        {
            int val = static_cast<int>(cur.Scalar());
            if (val < 0)
            {
                throw OML_Error(OML_ERR_NATURALNUM, i + 1, OML_VAR_VALUE);
            }
            if (startrow == -1)
            {
                startrow = val;
            }
            else
            {
                startcol = val;
            }
        }
        else if (cur.IsVector())
        {
            std::vector<double> dvec (cur.Vector());
            if (dvec.size() != 4)
            {
                std::string msg ("Error: invalid input; must be a 4-element vector");
                msg += " in argument " + 
                        std::to_string(static_cast<long long>(i + 1));
                throw OML_Error(msg);
            }
            for (int j = 0; j < 4; ++j)
            {
                double tmp = dvec[j];
                if (tmp < 0 || !utils.IsInt(tmp) || IsNegInf_T(tmp) ||
                    IsInf_T(tmp) || IsNaN_T(tmp))
                {
                    std::string msg ("Error: invalid input; vector must contain");
                    msg += "non-negative, finite integers in argument " + 
                            std::to_string(static_cast<long long>(i + 1));
                    throw OML_Error(msg);
                }
            }
            startrow = static_cast<int>(dvec[0]);
            startcol = static_cast<int>(dvec[1]);
            endrow   = static_cast<int>(dvec[2]);
            endcol   = static_cast<int>(dvec[3]);
        }
        else
        {
            throw OML_Error(OML_ERR_SCALAR_VECTOR_STRING, i + 1, OML_VAR_TYPE);
        }
    }

#ifdef OS_WIN
    std::wfstream ifs;
    ifs.open(wfile, std::ios_base::in | std::ios::binary);
#else
    std::ifstream ifs;
    ifs.open(filename, std::ios_base::in | std::ios::binary);
#endif
    if (!ifs)
    {
        throw OML_Error(OML_ERR_FILE_CANNOTOPEN, 1);
    }

    SetDefaultDelimiter(filename, delim);

    // Start reading the data
    try
    {
        bool quitLoop    = false;
        int  currentrow  = 0;
        int  m           = 0;
        int  n           = 0;
        int  maxn        = 0;
        int  rowcapacity = 1000;
        int  colcapacity = 1000;

        // Cannot directly push data into matrix as for large matrices, resizing
        // is very slow
        std::vector< std::vector<double> > data;
        data.reserve(rowcapacity);

        std::vector< std::vector<double> > idata;
        idata.reserve(rowcapacity);

        bool useLinuxLineDelim = false;

#ifdef OS_WIN
        if (ifs.good())
        {
            std::wstring firstline;
            std::getline(ifs, firstline);
            if (firstline.find(L"\r\n") == std::string::npos &&
                firstline.find(L'\r') != std::string::npos)
            {
                useLinuxLineDelim = true;
            }
            ifs.seekg(0);
        }
#else
        if (ifs.good())
        {
            std::string firstline;
            std::getline(ifs, firstline);
            if (firstline.find("\r\n") == std::string::npos &&
                firstline.find('\r')  != std::string::npos)
            {
                useLinuxLineDelim = true;
            }
            ifs.seekg(0);
        }
#endif

        // Read the data
        while (!quitLoop && !eval.IsInterrupt() && !ifs.eof())
        {
            std::string line;
#ifdef OS_WIN
            std::wstring wline;
            if (useLinuxLineDelim)
            {
                std::getline(ifs, wline, L'\r');
                if (ifs.good() && ifs.peek() == L'\n') // Could have mismatched line endings
                {
                    wchar_t nextch;  
                    ifs.get(nextch);
                    wline += nextch;
                }
            }
            else
            {
                std::getline(ifs, wline);
            }
            if (wline.empty() && !ifs.good())
            {
                break;  // End of file not detected on windows for empty file opened in wide mode
            }
            line = utils.WString2StdString(wline);
#else
            if (useLinuxLineDelim)
            {
                std::getline(ifs, line, '\r');
                if (ifs.good() && ifs.peek() == '\n') // Could have mismatched line endings
                {
                    char nextch;
                    ifs.get(nextch);
                    line += nextch;
                }
            }
            else
            {
                std::getline(ifs, line);
            }
            if (line.empty() && !ifs.good())
            {
                break;  
            }
#endif

		    if (startrow != -1 && currentrow < startrow || line.empty())
            {
			    currentrow++;  // Skip the header/empty lines
                if (quitLoop)
                {
                    break;
                }
			    continue; 
		    }
		    if (endrow != -1 && currentrow > endrow)
            {
                break;
		    }

            if (line[line.size() - 1] == '\n')
            {
                line.pop_back();
            }

            if (!line.empty() && line[line.size() - 1] == '\r')
            {
                line.pop_back();
            }
            if (!line.empty())
            {   // Some files could have junk characters
                size_t len = line.length();
                for (size_t i = 0; i < len; ++i)
                {
                    unsigned char ch = line[i];
                    if (isprint(ch))
                    {
                        if (i > 0)
                        {
                            line = line.substr(i);
                        }
                        break;
                    }
                }
            }


            if (line.empty())
            {
                currentrow++;
                if (quitLoop)
                {
                    break;
                }
                continue;
            }
            
            // Start of line parsing
            int currentcol = 0;
            n = 0;
            int used = 0;
            int result = 0;

            if (delim.empty())
            {
                if (std::any_of(line.begin(), line.end(), ::isdigit)) 
                {    
                    double dval = 0;
                    sscanf(line.c_str(), "%lf%n", &dval, &used);
                }
                for (size_t j = static_cast<size_t>(used); j < line.size(); ++j)
                {
                    char ch = line[j];
                    if (isdigit(ch) || ch == '-' || ch == '+' || ch == '.')
                    {
                        break;
                    }
                    else if (isspace(ch) || (ispunct(ch) && ch != '%'))
                    {
                        delim += ch;
                        if (ignoreMultidelim)
                        {
                            break;
                        }
                    }
                    else if (!delim.empty() && isalpha(ch))
                    {
                        break;
                    }
                }
                if (static_cast<size_t>(used) == line.size() && delim.empty())
                {
                    delim = " ";
                }
            }
            
            if (delim.empty())
            {
                throw OML_Error("Error: cannot determine delimter to process file in argument 1");
            }

            bool newlinedelim = false;
            EscapeDelimiters(delim, newlinedelim);
            std::vector<double> thiscol;
            thiscol.reserve(colcapacity);

            std::vector<double> thiscol2;   // Imaginary numbers
            thiscol2.reserve(colcapacity);

            bool iscomment = (line[0] == '%');

            size_t linepos = 0;
            bool quit = false;

            if (ignoreMultidelim)
            {
                line = BuiltInFuncsUtils::LTrim(line, delim);
                line = BuiltInFuncsUtils::RTrim(line, delim);
                if (!line.empty())
                {
                    size_t tmplen = line.length();
                    std::string newstr;
                    newstr += line[0];
                    bool wasdelim = false;
                    for (size_t k = 1; k < tmplen; ++k)
                    {
                        char ch = line[k];
                        if (delim.find(ch) == std::string::npos)
                        {
                            newstr += ch;
                            wasdelim = false;
                        }
                        else if (!wasdelim)
                        {
                            wasdelim = true;
                            newstr += ch;
                        }
                    }
                    line = newstr;
                }
            }
            while (!line.empty() && !quit)
            {
                linepos = line.find(delim);

                if (currentcol < startcol)
                {
                    if (linepos != std::string::npos)
                    {
                        line = line.substr(linepos + delim.length());
                    }
                    currentcol++;
                    continue;
                }
                if (endcol >= 0 && currentcol > endcol)
                {
                    break;
                }

                std::string tok = line.substr(0, linepos);
                double rval = remptyval;
                double ival = iemptyval;
                bool   isscalar = true;
                if (!iscomment)
                {
                    BuiltInFuncsUtils::ReadNumber(tok, rval, ival, isscalar);
                    if (!isscalar)
                    {
                        hascomplex = true;
                    }

                    if (!thiscol.empty() && thiscol.size() >= colcapacity)
                    {
                        colcapacity *= 2;
                        thiscol.reserve(colcapacity);
                        thiscol2.reserve(colcapacity);
                    }
                }
                thiscol.emplace_back(rval);
                thiscol2.emplace_back(ival);
                if (linepos != std::string::npos)
                {
                    line = line.substr(linepos + delim.length());
                }
                else
                {
                    quit = true;
                }
                n++;
                currentcol++;
            }

            if (!data.empty() && data.size() >= rowcapacity)
            {
                rowcapacity *= 2;
                data.reserve(rowcapacity);
            }    

            int colcount = (thiscol.empty()) ? 0 : static_cast<int>(thiscol.size());
            maxn = std::max(maxn, colcount);
            data.push_back(thiscol);
            idata.push_back(thiscol2);  // Imaginary vals
            m++;
            currentrow++;
        }

        ifs.close();

        assert(data.size() == m);
        if (data.size() != m)
        {
            throw OML_Error(OML_ERR_INTERNAL);
        }

        if (!hascomplex)
        {
            mtx.reset(EvaluatorInterface::allocateMatrix(m, maxn, true));
            mtx->SetElements(remptyval);

            for (int i = 0; i < m; ++i)
            {
                std::vector <double> rvals(data[i]);
                if (rvals.empty())
                {
                    continue;
                }
                int n = static_cast<int>(rvals.size());
                for (int j = 0; j < n; ++j)
                {
                    (*mtx)(i, j) = rvals[j];
                }
            }
        }
        else
        {
            mtx.reset(EvaluatorInterface::allocateMatrix(m, maxn, false));
            mtx->SetElements(hwComplex(remptyval, iemptyval));

            assert(data.size() == idata.size());
            for (int i = 0; i < m; ++i)
            {
                std::vector <double> rvals(data[i]);
                std::vector <double> ivals(idata[i]);
                if (rvals.empty())
                {
                    continue;
                }

                int n = static_cast<int>(rvals.size());
                assert(rvals.size() == ivals.size());
                for (int j = 0; j < n; ++j)
                {
                    mtx->z(i, j) = hwComplex(rvals[j], ivals[j]);
                }
            }
        }
    }
    catch (const OML_Error& e)
    {
        if (ifs.is_open())
        {
            ifs.close();
        }
        throw e;
    }
    catch (...)
    {
        if (ifs.is_open())
        {
            ifs.close();
        }
        throw OML_Error(OML_ERR_INTERNAL);
    }
    if (ifs.is_open())
    {
        ifs.close();
    }
    outputs.push_back(mtx.release());
    return true;
}
//------------------------------------------------------------------------------
// Returns a vector range, given a string range in excel format
//------------------------------------------------------------------------------
 std::vector<int> BuiltInFuncsFile::String2VectorRange(const std::string& strRange,
                                                       int                idx)
{
    std::vector<int> range = std::vector<int>(4, -1);

    std::string strbegin;
    std::string strend;
    std::string delim = ":";

    size_t pos = strRange.find(delim);
    if (pos == std::string::npos)
    {
        delim = "..";
        pos = strRange.find(delim);
    }

    if (pos == std::string::npos)
    {
        throw OML_Error(OML_ERR_SPREADSHEET_RANGE, idx + 1);
    }

    strbegin = strRange.substr(0, pos);
    strend   = strRange.substr(pos + delim.length());

    if (strbegin.empty() || strend.empty())
    {
        throw OML_Error(OML_ERR_SPREADSHEET_RANGE, idx + 1);
    }

    int startrow = -1;
    int startcol = -1;
    int endrow   = -1;
    int endcol   = -1;

    OML_Error err(OML_ERR_SPREADSHEET_RANGE, idx + 1);
    GetRowColInfo(strbegin, err.GetErrorMessage(), startrow, startcol);
    GetRowColInfo(strend,   err.GetErrorMessage(), endrow,   endcol);

    range[0] = startrow;
    range[1] = startcol;
    range[2] = endrow;
    range[3] = endcol;
    return range;
}
 //-----------------------------------------------------------------------------
 // Parses given string to get row/column information
 //-----------------------------------------------------------------------------
 void BuiltInFuncsFile::GetRowColInfo(const std::string& str,
                                      const std::string& err,
                                      int&               row,
                                      int&               col)
{
    std::string alpha;
     
    for (size_t i = 0; i < str.length(); ++i)
    {
        char ch = str[i];
        if (isdigit(ch))
        {
            std::string tmp = str.substr(i);
            row = atoi(tmp.c_str()) - 1;
            if (row < 0)
            {
                throw OML_Error(err);
            }
            break;
        }
        else if (isalpha(ch))
        {
            alpha += ch;
        }
        else
        {
            throw OML_Error(err);
        }
    }

    if (alpha.empty())
    {
        throw OML_Error(err);
    }

    col = 0;
    std::transform(alpha.begin(), alpha.end(), alpha.begin(), ::toupper);
    for (size_t i = 0; i < alpha.size(); ++i)
    {
        col *= 26;
        col += alpha[i] - 'A' + 1;
    }

    col -= 1;

    if (col < 0)
    {
        throw OML_Error(err);
    }
 }
 //-----------------------------------------------------------------------------
 // Returns true and opens a file for read/write [fopen]
 //-----------------------------------------------------------------------------
 bool BuiltInFuncsFile::Fopen(EvaluatorInterface           eval,
                              const std::vector<Currency>& inputs,
                              std::vector<Currency>&       outputs)
 {
     if (inputs.empty())
     {
         throw OML_Error(OML_ERR_NUMARGIN);
     }
     size_t nargin = inputs.size();

     const Currency& input = inputs[0];
     if (!input.IsString() && !input.IsPositiveInteger())
     {
         throw OML_Error(OML_ERR_STRING_FILESTREAM, 1);
     }

     BuiltInFuncsUtils utils;
     outputs.reserve(2);
     if (!input.IsString())  // File id
     {
         int fileid = static_cast<int>(input.Scalar());
         utils.CheckFileIndex(eval, fileid, 1, true);
         outputs.push_back(eval.GetFileName(fileid));
         outputs.push_back(eval.GetFileMode(fileid));
         return true;
     }

     // File name
     std::string fname (input.StringVal());
     if (input.IsMultilineString())
     {
         size_t pos = fname.find_first_not_of("\n");
         if (pos != std::string::npos)
         {
             fname = fname.substr(pos);
         }
         if (!fname.empty())
         {
             pos = fname.find("\n");
             if (pos != std::string::npos)
             {
                 fname = fname.substr(0, pos);
             }
         }
     }
     if (fname.empty())
     {
         outputs.push_back(-1);
         outputs.push_back("File name must be a non-empty string");
         return true;
     }
     else if (fname == "all")  // Returns all open file ids
     {
         std::vector<int> indices(eval.GetFileIndices(FIRST_USER_FILE));
         std::unique_ptr<hwMatrix> mtx(EvaluatorInterface::allocateMatrix());
         if (!indices.empty())
         {
             int numvals = static_cast<int>(indices.size());
             mtx.reset(EvaluatorInterface::allocateMatrix(1, numvals, true));
             for (int i = 0; i < numvals; ++i)
             {
                 (*mtx)(i) = indices[i];
             }
         }
         outputs.push_back(mtx.release());
         outputs.push_back(std::string());
         return true;
     }

     std::string mode = "rb"; // Default mode
     if (nargin > 1)
     {
         if (!inputs[1].IsString())
         {
             throw OML_Error(OML_ERR_STRING, 2);
         }
         mode = inputs[1].StringVal();
         if (mode.empty() || mode.length() > 3)
         {
             throw OML_Error(OML_ERR_INVALID_FILE_MODE, 2);
         }
         std::transform(mode.begin(), mode.end(), mode.begin(), ::tolower);
         if (mode[0] != 'r' && mode[0] != 'w' && mode[0] != 'a')
         {
             throw OML_Error(OML_ERR_INVALID_FILE_MODE, 2);
         }
         size_t len = mode.length();
         if (len >= 2)
         {
             char second = mode[1];
             if ((second != '+' && second != 'b'  && second != 't') ||
                 (len == 3      && (second == 'b' || second == 't')))
             {
                 throw OML_Error(OML_ERR_INVALID_FILE_MODE, 2);
             }
         }
         char last = mode[len - 1];
         if (last == 't')
         {
             mode.pop_back();
         }
         else if (mode.find("b") == std::string::npos)
         {
             mode += "b";
         }
     }
     std::FILE* f = utils.FileOpen(utils.Normpath(fname), mode);
     if (f)
     {
         outputs.push_back(eval.AddFile(f, fname, mode));
         outputs.push_back(std::string());
#ifndef OS_WIN
         if (mode[0] == 'a' || mode[0] == 'w' || mode.find("r+") != std::string::npos)
         {
             std::string normpath(utils.Normpath(fname));
             int result = chmod(utils.GetAbsolutePath(normpath).c_str(), S_IRWXU);
             if (result != 0)
             {
                 utils.SetWarning(eval, strerror(errno));
             }
         }
#endif
     }
     else
     {
         outputs.push_back(-1.0);
         outputs.push_back(std::string(strerror(errno)));
     }

     return true;
 }
//------------------------------------------------------------------------------
// Returns true after parsing input and gets formats
// This is very similar to BuiltInFuncsUtils::GetFormats but handles the case of
// '* %d'. Not sure whether this can be implemented across the board for all
// file utilities which use BuiltInFuncsUtils::GetFormats
//------------------------------------------------------------------------------
 bool BuiltInFuncsFile::GetFormats(const std::string& input,
     std::vector<std::string>& basefmts,
     std::vector<std::string>& rawfmts,
     std::vector<bool>& usefmts,
     std::string& err,
     std::string& prefix)
 {
     if (input.empty())
     {
         basefmts.push_back("%f");
         rawfmts.push_back("%lf");  // Default option
         usefmts.push_back(true);
         return true;
     }

     // Handle case where the format string does not start with %
     std::string in(input);
     size_t pos = in.find_first_of("%");
     if (pos > 0 && pos != std::string::npos)
     {
         prefix = input.substr(0, pos);
         in     = input.substr(pos);
     }

     bool firsttoken = true;
     char* tok = strtok((char*)in.c_str(), "%");

     while (tok)
     {
         std::string tmp(tok);

         std::string basefmt;
         std::string rawfmt;
         bool        use = true;
         int  numDecimalPts = 0;
         bool hasOnlyDigits = true;

         size_t len = tmp.size();

         if (tmp[0] == '[' && tmp[len - 1] == ']')
         {
             basefmts.push_back("custom");
             rawfmts.push_back(tok);
             usefmts.push_back(true);
             tok = strtok(NULL, "%");
             continue;
         }
         for (size_t i = 0; i < len; ++i)
         {
             char ch = tmp[i];
             if (i == 0 && ch == '*')
             {
                 use = false;
                 continue;
             }
             else if (!use)
             {
                 if (ch == '[')
                 {
                     basefmt = "custom";
                     rawfmt += ch;
                     continue;
                 }
             }
             if (!basefmt.empty())
             {
                 rawfmt += ch; // Add any characters after the specifier are ok
                 continue;
             }
             else if (ch == 'c' || ch == 'd' || ch == 'e' || ch == 'f' ||
                 ch == 'g' || ch == 'i' || ch == 's' || ch == 'n')
             {
                 basefmt += ch;
                 if (ch == 'f' || ch == 'g')
                 {
                     if (!hasOnlyDigits || numDecimalPts > 1)
                     {
                         rawfmt = "";
                         err = "Error: invalid format %" + tmp;
                         break;
                     }
                     rawfmt += "lf";
                 }
                 else if (ch == 'c' || ch == 'd' || ch == 'n')
                 {
                     if (!hasOnlyDigits || numDecimalPts > 0)
                     {
                         rawfmt = "";
                         err = "Error: invalid format %" + tmp;
                         break;
                     }
                     if (ch == 'n')
                     {
                         ch = 'd';
                     }
                     rawfmt += ch;
                 }
                 else
                 {
                     rawfmt += ch;
                 }
                 continue;
             }
             rawfmt += ch;

             if (ch == '.')
             {
                 numDecimalPts++;
             }
             else if (!std::isdigit(static_cast<unsigned char>(ch)))
             {
                 hasOnlyDigits = false;
             }
         }

         if (basefmt.empty() || rawfmt.empty())
         {
             return false;
         }

         std::string percent = (basefmt != "custom") ? "%" : "";
         basefmts.push_back(percent + basefmt);

         if (firsttoken)
         {
             rawfmts.push_back(prefix + percent + rawfmt);
             firsttoken = false;
         }
         else
         {
             rawfmts.push_back(percent + rawfmt);
         }
         usefmts.push_back(use);

         tok = strtok(NULL, "%");
     }

     // Check if there are tabs in the format
     size_t i = 0;
     for (std::vector<std::string>::iterator itr = rawfmts.begin();
         itr != rawfmts.end(); ++itr, ++i)
     {
         std::string fmt(*itr);
         size_t      pos = fmt.find("\\t");
         bool        replace = false;
         if (pos != std::string::npos)
         {
             fmt.replace(pos, 2, "\t");
             rawfmts[i] = fmt;
         }
     }
     return true;
 }
 //-----------------------------------------------------------------------------
 // Sets default delimiter for certain file types
 //-----------------------------------------------------------------------------
 void BuiltInFuncsFile::SetDefaultDelimiter(const std::string& filename, 
                                            std::string&       delim)
 {
     if (delim.empty())
     {
         BuiltInFuncsUtils utils;
         std::string ext(utils.GetFileExtension(filename));
         if (!ext.empty())
         {
             std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
             if (ext == "csv")
             {
                 delim = ",";
             }
         }
     }
 }
 //-----------------------------------------------------------------------------
 // Returns true and reads a comma separated file [csvread]
 //-----------------------------------------------------------------------------
 bool BuiltInFuncsFile::Csvread(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
 {
     if (inputs.empty())
     {
         throw OML_Error(OML_ERR_NUMARGIN);
     }

     std::vector<Currency> in(inputs);
     in.insert(in.begin() + 1, Currency(","));

     return Dlmread(eval, in, outputs);
 }
 //-----------------------------------------------------------------------------
// Returns true and writes a comma separated file [csvwrite]
//-----------------------------------------------------------------------------
 bool BuiltInFuncsFile::Csvwrite(EvaluatorInterface           eval,
                                 const std::vector<Currency>& inputs,
                                 std::vector<Currency>& outputs)
 {
     if (inputs.empty() || inputs.size() < 2)
     {
         throw OML_Error(OML_ERR_NUMARGIN);
     }

     std::vector<Currency> in(inputs);
     in.insert(in.begin() + 2, Currency(","));

     return Dlmwrite(eval, in, outputs);
 }
 //-----------------------------------------------------------------------------
 // Escapes delimiters
 //-----------------------------------------------------------------------------
 void BuiltInFuncsFile::EscapeDelimiters(std::string& delims, bool& newlinedelim)
 {
     newlinedelim = false;
     if (delims.empty())
     {
         delims = " ";
         return;
     }

     size_t pos = delims.find("\\t");
     if (pos != std::string::npos)
     {
         delims.replace(pos, 2, "\t");
     }

     pos = delims.find("\\n");
     if (pos != std::string::npos)
     {
        if (delims == "\\n")
        {
            newlinedelim = true;
        }
        else
        {
            std::string tmp = delims.substr(0, pos);
            tmp += delims.substr(pos + 2);
            delims = tmp;
        }
    }
    pos = delims.find("\\r");
    if (pos != std::string::npos)
    {
        if (delims == "\\r")
        {
            newlinedelim = true;
        }
        else
        {
            std::string tmp = delims.substr(0, pos);
            tmp += delims.substr(pos + 2);
            delims = tmp;
        }
    }
 }
 //-----------------------------------------------------------------------------
 // Returns a line read from the given file, with the given end of line character
 //-----------------------------------------------------------------------------
 std::string BuiltInFuncsFile::GetLine(std::FILE* fp, const std::string& endofline)
 {
     if (!fp || feof(fp))
     {
         return std::string();
     }

     std::string line;

     if (endofline.empty())
     {
         char lineC[2048];
         while (!feof(fp))  // Loop till we hit \r or \n
         {
             memset(lineC, 0, sizeof(lineC));
             if (fgets(lineC, sizeof(lineC), fp) == nullptr)
             {
                 break;
             }
             line += lineC;
             if (strstr(lineC, "\n") || strstr(lineC, "\r"))
             {
                 break;
             }
         }
         if (!line.empty() && line.back() == '\n')
         {
             line.pop_back();
         }
         if (!line.empty() && line.back() == '\r')
         {
             line.pop_back();
         }
         return line;
     }

     char ref = endofline[0];
     while (!feof(fp))  // Loop till we hit \r or \n
     {
         int ch = fgetc(fp);
         if (ch == EOF || ch == ref)
         {
             break;
         }
         line += ch;
     }
     return line;
 }