/**
* @file BuiltInFuncsUtils.cpp
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

// Begin defines/includes
#include "BuiltInFuncsUtils.h"

#include <cassert>
#include <cctype>  // For std::isdigit
#include <climits>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <memory>  // For std::unique_ptr
#include <sys/stat.h>

#ifdef OS_WIN
#    include <io.h>
#    include <Windows.h>
#else
#    include <dirent.h>
#    include <dlfcn.h>
#    include <glob.h>
#    include <stdlib.h>
#    include <sys/resource.h>
#    include <sys/times.h>
#    include <time.h>
#    include <unistd.h>
#endif

#include "ErrorInfo.h"
#include "Evaluator.h"
#include "MatrixNDisplay.h"
#include "OML_Error.h"

#include <hwMatrix.h>

#ifdef OS_WIN
// Undefine these as Microsoft defines them elsewhere in windows headers
#    undef FormatMessage
#endif
// End defines/includes

//------------------------------------------------------------------------------
// Returns Currency after reading matrix row, assumes input is a matrix or string
//------------------------------------------------------------------------------
Currency BuiltInFuncsUtils::ReadRow(EvaluatorInterface& eval, 
                                    const Currency&     input, 
                                    int                 index)
{
    const hwMatrix* loc = input.Matrix();
    bool alloc = false;

    if (!loc)
    {
        alloc = true;
        loc = EvaluatorInterface::allocateMatrix();
    }

    std::unique_ptr<hwMatrix>row(EvaluatorInterface::allocateMatrix(1, loc->N(), loc->Type()));

    if (index < loc->M())
    {
        hwMathStatus stat = loc->ReadRow(index, *row);
        if (alloc && !stat.IsOk())
        {
            delete loc;
            loc = nullptr;
        }

        CheckMathStatus(eval, stat);
    }

    if (alloc && loc)
    {
        delete loc;
        loc = nullptr;
    }

    Currency out(row.release());
    out.SetMask(input.GetMask());
    return out;
}
//------------------------------------------------------------------------------
// Checks math status and throws error/reports warning if not ok
//------------------------------------------------------------------------------
void BuiltInFuncsUtils::CheckMathStatus(EvaluatorInterface& eval, 
                                        hwMathStatus        mstat)
{
    if (mstat.IsOk()) return;

    if (!mstat.IsWarning())
        throw OML_Error(mstat);

    SetWarning(eval, mstat.GetMessageString());
}
//------------------------------------------------------------------------------
// Sets warning message \todo: Replace hmlwarning method
//------------------------------------------------------------------------------
void BuiltInFuncsUtils::SetWarning(EvaluatorInterface& eval, const std::string& str)
{
	std::string warningstr = eval.FormatMessage(str);
    Currency warn(warningstr);
    warn.DispOutput();
    eval.PrintResult(warn);
    EvaluatorInterface::SetLastWarning(warningstr);
}
//------------------------------------------------------------------------------
// Returns currency which gives formatted output of strings. Strings will be 
//------------------------------------------------------------------------------
Currency BuiltInFuncsUtils::FormatOutput(const std::vector<std::string>& in,
                                         bool                            handleEmptyStr,
                                         bool                            handleNaN,
                                         bool                            alignright,
                                         char                            pad)
{
    int rows = static_cast<int>(in.size());

    // Scan elements to get maxwidth, which will be the number of columns in mtx
    int cols = 0;
    for (int i = 0; i < rows; ++i)
    {
        std::string tmp (in[i]); 

        bool empty = tmp.empty();

        if (!handleEmptyStr && empty) continue;

        int thiswidth = empty ? 3 : static_cast<int>(tmp.size());  // 'Nan' if empty
        if (cols < thiswidth)
            cols = thiswidth;
    }

    // Pad output with leading zeros to width of largest element
    hwMatrix* out = EvaluatorInterface::allocateMatrix(rows, cols, hwMatrix::REAL);
    for (int i = 0; i < rows; ++i)
    {
        std::string tmp (in[i]);
        bool empty = tmp.empty();

        if (empty && handleEmptyStr)
            tmp = "Nan";
                        
        std::ostringstream os;

        (alignright) ? os << std::right : os << std::left;

        if (handleNaN && tmp == "Nan")                     
            os << std::setfill(' ') << std::setw(cols);
        else
            os << std::setfill(pad) << std::setw(cols);
            
        os << tmp;
        std::string val (os.str());
        
        for (int j = 0; j < cols; ++j) 
		    (*out)(i, j) = (unsigned char) val[j];
    }

    Currency result (out);
    result.SetMask(Currency::MASK_STRING);
    return result;
}
//------------------------------------------------------------------------------
// Returns true if the given double is an integer
//------------------------------------------------------------------------------
bool BuiltInFuncsUtils::IsInt(double val)
{
    double intpart      = 0.0;
    double fractionpart = std::modf(val, &intpart);

    if (abs(fractionpart) < std::numeric_limits<double>::epsilon())
        return true;

    return false;
}
//------------------------------------------------------------------------------
// Checks char bounds and returns valid char
//------------------------------------------------------------------------------
unsigned char BuiltInFuncsUtils::GetValidChar(EvaluatorInterface& eval, 
                                              double              val,
                                              bool                showerr)
{
    if (val >= 0 && val <= UCHAR_MAX) return (static_cast<unsigned char>(val));

    if (showerr) throw OML_Error(HW_ERROR_SCALOUTCHARRANGE);

    SetWarning(eval, "Warning: scalar outside of character range");
    return 0;
}
//------------------------------------------------------------------------------
// Gets ordered string. Assumes all elements of mtx are within unsigned char range
//------------------------------------------------------------------------------
std::string BuiltInFuncsUtils::GetOrderedStringVal(const Currency& cur)
{
    const hwMatrix* strmtx = cur.Matrix();
    if (!strmtx) return "";

    int         numelem = static_cast<int>(strmtx->Size());
    std::string returnval;
    for (int i = 0; i < numelem; ++i)
        returnval += static_cast<unsigned char>((*strmtx)(i));

    return returnval;
}
//------------------------------------------------------------------------------
// Parses given string input and gets a vector of format options
//------------------------------------------------------------------------------
std::vector<std::string> BuiltInFuncsUtils::GetFormats(const std::string& formatdesc,
                                                       const std::string& validfmtdesc,
                                                       bool&              validformat)
{
    if (formatdesc.empty()) return std::vector<std::string>();

    std::vector<std::string> formats;
    size_t                   startpos      = 0;
    size_t                   formatdescLen = formatdesc.size();
    std::string              allfmt (formatdesc);
	while (!allfmt.empty())
	{        
        size_t pos = allfmt.find_first_of("%");
        if (pos == std::string::npos) return formats;

        size_t      nextpos = allfmt.find("%", pos + 1);
        std::string rawfmt  = (nextpos == std::string::npos) ?
                               allfmt.substr(pos) :
                               allfmt.substr(pos, nextpos - pos);
        
        if (rawfmt.empty()) return formats;

        // Find the substring before the '%' if any
        std::string prev; 
        if (pos > 0)
            prev = allfmt.substr(0, pos);

        // Find the alpha character after the raw format. We should not allow
        // options like %.5f. This will be treated as %f and we can throw a warning
        std::string fmt;
        size_t      rawfmtSize   = rawfmt.size();
        
        for (size_t i = 0; i < rawfmtSize; ++i)
        {
            char c = rawfmt[i];
            if (i == 0 && c == '%')
            {
                fmt += c;
                continue;
            }
            if (!isalpha(c))
            {
                validformat = false; 
                continue;
            }

            // First alpha character
            std::string thisfmt ("%");
            thisfmt += c;
            if (!validfmtdesc.empty() && validfmtdesc.find(thisfmt) == std::string::npos)
            {
                validformat = false;  // Skip this format like it was not there
                fmt = "";
                break;
            }
            fmt += c;
            if (i < rawfmtSize - 1)
                fmt += rawfmt.substr(i + 1);

            break;
        }

        if (fmt != "%" && !fmt.empty())
            formats.push_back(prev + fmt);
        
        if (nextpos == std::string::npos) break;

        allfmt = allfmt.substr(nextpos);
    }

    return formats;
}
//------------------------------------------------------------------------------
// Utility function which reads currency for sizes, returns true on success
//------------------------------------------------------------------------------
void BuiltInFuncsUtils::GetSizeSpecifications(const Currency& cur, 
                                              int             varindex,
                                              double&         rows, 
                                              double&         cols,
                                              bool&           isMtx)
{
    if (!cur.IsScalar() && !cur.IsMatrix())
        throw OML_Error(OML_ERR_POS_INTEGER_MTX_INF, varindex, OML_VAR_DIM);

    if (cur.IsScalar())
    {
        double val = cur.Scalar();
        if (!cur.IsPositiveInteger() && !IsInf_T(val))
            throw OML_Error(OML_ERR_POS_INTEGER_MTX_INF, varindex, OML_VAR_DIM);

        cols  = 0; 
        rows  = cur.Scalar();
        isMtx = false;
        return;
    }
    
    const hwMatrix* mtx = cur.Matrix();
    if (!mtx || mtx->Size() < 2) 
        throw OML_Error(OML_ERR_POS_INTEGER_MTX_SIZE2, varindex, OML_VAR_DIM);

    rows = (*mtx)(0);
    if (!IsInt(rows) || !(rows > 0) || IsInf_T(rows))
        throw OML_Error(OML_ERR_POSINTEGER, varindex, OML_VAR_DIM);

    cols = (*mtx)(1);
    if (!IsInf_T(cols) && !IsInt(cols) && !(cols > 0))
        throw OML_Error(OML_ERR_POS_INTEGER_MTX_INF, varindex, OML_VAR_DIM);
    isMtx = true;
}
//------------------------------------------------------------------------------
// Returns currency after reading formatted input from a file/string
//------------------------------------------------------------------------------
Currency BuiltInFuncsUtils::GetFormattedInput(const std::string& in,
                                              const std::string& formatdesc,
                                              const std::string& validfmtdesc,
                                              double             rows,
                                              double             cols,
                                              bool               hasSizeMtx,
                                              bool               hasSizeSpec,
                                              bool&              validformat)
{
    assert(!in.empty());

    // Split tokens in format description to get the vector of format options
    std::vector<std::string> formats = GetFormats(formatdesc, validfmtdesc, validformat);

    // For each format, read in one value from input and push into values vector
    BuiltInFuncsUtils utils;
    std::vector<Currency> values;
    int                   sizelimit = (int)(hasSizeMtx ? rows * cols : rows);
    utils.ReadFormattedInput(in, sizelimit, formats, values);

    // Split the output values vector into the correct size (if any) Right now, 
    // the result is Nx1 (which is the default)
    // Handles cases where inf is passed. Don't cast and use values.size() as 
    // compiler will convert
    int    numvals      = values.empty() ? 0 : static_cast<int>(values.size());
	int    numrows      = numvals;
	int    numcols      = 1;
	double rawnumcols   = 1.0;
	bool   stringoutput = false;
	
    if (hasSizeSpec)                // Has size specifications
    {
        if (hasSizeMtx)             // Has matrix
        {
            numrows    = (int)rows;
            rawnumcols = cols;      // This can be infinity
        }
        else if ((formats.size() == 1) && (formats[0].find("%s") != std::string::npos))
        {
            rawnumcols   = numrows;
            numrows      = (int)(IsInf_T(rows) ? numvals : rows);
            stringoutput = true;
        }
        else
        {
			rawnumcols = 1;
            numrows    = (int)(IsInf_T(rows) ? numvals : rows);
			if (numrows > values.size())  // Don't cast here when comparing
			    numrows = numvals;
        }
    }

	if (IsInf_T(rawnumcols) || IsNegInf_T(rawnumcols) || IsNaN_T(rawnumcols))
    {
        int tmp = (numrows == 0) ? 1 : numrows;
		numcols = (int)values.size() / tmp;
    }
	else
		numcols = static_cast<int>(rawnumcols);

	if (values.size() < numrows * numcols)
	{
        int tmp = (numrows == 0) ? 1 : numrows;
		if (values.size() % tmp == 0)
			numcols = (int)values.size() / tmp;
	}


	hwMatrix* mymat = EvaluatorInterface::allocateMatrix(numrows, numcols, 0.0);
    assert(mymat);

	for (int j = 0; j < values.size(); ++j)
	{
		if (j < values.size() && j < mymat->Size())
			(*mymat)(j) = values[j].Scalar();
	}

	Currency result(mymat);
	if (stringoutput)
		result.SetMask(Currency::MASK_STRING);

	return result;
}
//------------------------------------------------------------------------------
// Reads formatted float from string using sscanf, returns true if successul
//------------------------------------------------------------------------------
bool BuiltInFuncsUtils::SscanfHelperFloat(const std::string&     in,
                                          const std::string&     fmt,
                                          std::vector<Currency>& outvals,
                                          std::string&           stringread)
{    
    // Read into double to handle precision issues
    std::string updatedfmt(fmt);
    size_t pos = fmt.find("%g");
    if (pos == std::string::npos)
        pos = fmt.find("%f");                
    if (pos != std::string::npos)
        updatedfmt.replace(pos, 2, "%lf"); 

    double output = 0.0;
    if (sscanf(in.c_str(), updatedfmt.c_str(), &output) != 1) return false;

    outvals.push_back(output);

    std::ostringstream os;
    os << output;
    stringread = os.str();
    if (stringread.find("e") != std::string::npos ||
        stringread.find("E") != std::string::npos)
    {
        os.str("");
        os.clear();

        if (stringread.find("e") != std::string::npos)
        {
            os << std::scientific << output;
        }
        else
        {
            os << std::scientific << std::uppercase << output;
        }
        stringread = os.str();
    }

    // sscanf will ignore decimal point in the case of 5.0000
    if (stringread.find(".") == std::string::npos)
    {   // Read till we get to decimal point
        pos = in.find_first_of(stringread);
        if (pos != std::string::npos)
        {
            size_t len = in.size();
            for (size_t i = pos + stringread.size(); i < len; ++i)
            {
                char c = in[i];
                if (c == '.') 
                {
                    stringread += c;
                    break;
                }
                else if (!isdigit(c)) 
                    break;

                stringread += c;
            }
        }
    }
    // sscanf will read '00' as 0, so we need to take care of this condition
    pos = in.find_first_of(stringread);
    if (pos != std::string::npos)
    {
        size_t len = in.size();
        for (size_t i = pos + stringread.size(); i < len; ++i)
        {
            char c = in[i];
            if (!isdigit(c)) break;
            stringread += c;   
        }
    }
    return true;
}
//------------------------------------------------------------------------------
// Reads formatted string from string using sscanf and returns true if successul
//------------------------------------------------------------------------------
bool BuiltInFuncsUtils::SscanfHelperString(const std::string&     in,
                                           const std::string&     fmt,
                                           std::vector<Currency>& outvals,
                                           std::string&           stringread)
{
    char output [1024];
	if (sscanf(in.c_str(), fmt.c_str(), &output) != 1) return false;

    size_t len = strlen(static_cast<char*>(output));
	for (size_t j = 0; j < len; ++j)
		outvals.push_back(static_cast<int>(output[j]));

    stringread = output;
    return true;
}
//------------------------------------------------------------------------------
// sscanf helper function, reads formatted input from string, returns true if successful
//------------------------------------------------------------------------------
bool BuiltInFuncsUtils::SscanfHelper(std::string&           in,
                                     const std::string&     fmt,
                                     std::vector<Currency>& outvals)
{
    assert(!fmt.empty());

    std::string stringread;

    if (fmt.find("%f") != std::string::npos || fmt.find("%g") != std::string::npos)
    {    // Float format
        bool result = SscanfHelperFloat(in, fmt, outvals, stringread);
        if (!result) return false;
    }
    else if (fmt.find("%d") != std::string::npos) // Integer format
	{
        std::string format  = fmt + "%n";
        int         numread = 0;
        int         output  = 0;

        if (sscanf(in.c_str(), format.c_str(), &output, &numread) != 1) 
        {
            return false;
        }
        outvals.push_back(output);

        // Chop the input string
        if (numread > 0)
        {
            in = in.substr(numread);
        }
        return true;  // Done parsing
    }
    else if (fmt.find("%s") != std::string::npos)  // String format
	{
        bool result = SscanfHelperString(in, fmt, outvals, stringread);
        if (!result) return false;
    }
	
    if (stringread.empty()) return false;

    if (fmt.size() > 2)  // Tack whatever strings are there in format template
        stringread += fmt.substr(2);

    size_t curpos = in.find_first_of(stringread);
    if (curpos !=  std::string::npos)
        curpos += stringread.size();

    if (curpos >= in.size())
    {
        in = "";  // Done parsing
        return false;
    }
        
    in = in.substr(curpos);   // Chop the input string
    return true;
}
//------------------------------------------------------------------------------
// Reads formatted input from string
//------------------------------------------------------------------------------
void BuiltInFuncsUtils::ReadFormattedInput(const std::string&              input,
                                           int                             sizelimit,
                                           const std::vector<std::string>& formats,
                                           std::vector<Currency>&          values)
{
    // For each format, read in one value from input and push into values vector
	int  count = 0;
    bool keepgoing = true;
    std::string in(input);
    BuiltInFuncsUtils utils;
    while (keepgoing && !in.empty())
    {
        size_t oldvaluesSize = values.empty() ? 0 : values.size();
        for (std::vector<std::string>::const_iterator itr = formats.begin();
             itr != formats.end(); ++itr)
	    {
            bool result = utils.SscanfHelper(in, *itr, values);
            if (!result || in.empty())
            {
                keepgoing = false;
                break;
            }
            count++;
            if (count == sizelimit)
                keepgoing = false;
        }

        size_t newvaluesSize = values.empty() ? 0 : values.size();
        if (oldvaluesSize == newvaluesSize) break; // Nothing has been read
	}
}
//------------------------------------------------------------------------------
// Returns true if file exists
//------------------------------------------------------------------------------
bool BuiltInFuncsUtils::FileExists(const std::string& name)
{
    if (name.empty()) return false;

#ifdef OS_WIN
    return (_access(name.c_str(), 0) == 0);
#else
    return (access(name.c_str(), 0) == 0);
#endif
}
//------------------------------------------------------------------------------
// Utility to set slice in an ND matrix
//------------------------------------------------------------------------------
void BuiltInFuncsUtils::SetMatrixNSlice(const hwMatrix* mtx, size_t index, hwMatrixN* lhs)
{
    if (!mtx || !lhs) return;

    hwMatrixN* rhs = ExprTreeEvaluator::Convert2DtoND(mtx);
    if (!rhs) return;

    std::vector<hwSliceArg> slices;
    slices.push_back(hwSliceArg());
    slices.push_back(hwSliceArg());
    slices.push_back((int)index);
    lhs->SliceLHS(slices, *rhs);

    delete rhs;
}
//------------------------------------------------------------------------------
// Returns true if given file path is absolute
//------------------------------------------------------------------------------
bool BuiltInFuncsUtils::IsAbsolutePath(const std::string& path)
{
    if (path.empty()) 
        return false;

    char ch = path[0];
#if OS_WIN
    if (!(ch == '/' || ch == '\\' || (path.length() > 1 && path[1] == ':')))
        return false;
#else
    if (ch != '/')
        return false;
#endif

    return true;
}
//------------------------------------------------------------------------------
// Strips trailing slash
//------------------------------------------------------------------------------
void BuiltInFuncsUtils::StripTrailingSlash(std::string& path)
{
    while (!path.empty())
    {
        char last = path.back();
        if (last == '/' || last == '\\')
            path.erase(path.end() - 1);
        else
            break;  // Reached the end
    }
}
//------------------------------------------------------------------------------
// Returns absolute path
//------------------------------------------------------------------------------
std::string BuiltInFuncsUtils::GetAbsolutePath(const std::string& path)
{
    if (path.empty()) 
        return "";

    char* tmp = NULL;
#ifdef OS_WIN
    if (path.length() == 2 && isalpha(path[0]) && path[1] == ':')
        return path;

    tmp = _fullpath(nullptr, path.c_str(), 0);
#else    
    tmp = realpath(path.c_str(), NULL);
#endif
   
    if (!tmp)
        return "";

    std::string abspath = tmp;
    free(tmp);
    tmp = NULL;

    return abspath;
}
//------------------------------------------------------------------------------
// Returns true if this is the root dir (for windows)
//------------------------------------------------------------------------------
bool BuiltInFuncsUtils::IsRootDir(const std::string& path)
{
#ifdef OS_WIN
    return (path.length() == 2 && isalpha(path[0]) && path[1] == ':');
#else
    return false;
#endif
}
//------------------------------------------------------------------------------
// Gets current working directory
//------------------------------------------------------------------------------
std::string BuiltInFuncsUtils::GetCurrentWorkingDir()
{
    std::string workingdir;

#ifdef OS_WIN
    int   expectedSize = GetCurrentDirectory(0, nullptr);
    char* cwd          = new char[expectedSize+1];
    memset(cwd, 0, sizeof(char)*(expectedSize+1));
    if (!GetCurrentDirectory(expectedSize, cwd))
        throw OML_Error(HW_ERROR_NOTFINDCURWORKDIR);

    workingdir = cwd;
    delete [] cwd;
    cwd = NULL;
#else
    char *cwd = getcwd(nullptr, 0);
    if (!cwd)
        throw OML_Error(HW_ERROR_NOTFINDCURWORKDIR);
    workingdir = cwd;
    free(cwd);
#endif

    return workingdir;
}
//------------------------------------------------------------------------------
// Gets file names matching pattern in current directory
//------------------------------------------------------------------------------
std::vector<std::string> BuiltInFuncsUtils::GetMatchingFiles(const std::string& pattern)
{
    std::vector<std::string> files;
#ifdef OS_WIN
    WIN32_FIND_DATA data;
    HANDLE          handle = FindFirstFile(pattern.c_str(), &data);

    while (handle != INVALID_HANDLE_VALUE)
    {
        files.push_back(std::string(data.cFileName));

        if (!FindNextFile(handle, &data))
            break;
    }
#else
    glob_t data;
    int retcode = glob(pattern.c_str(), GLOB_TILDE, nullptr, &data);

    if (retcode)
    {
        globfree(&data);
        if (retcode == GLOB_NOSPACE)
            throw OML_Error(HW_ERROR_OUTMEM);
        else if (retcode == GLOB_ABORTED)
            throw OML_Error(HW_ERROR_READ);
    }
    else
    {
        for (size_t i = 0; i < data.gl_pathc; ++i)
        {
            files.push_back(basename(data.gl_pathv[i]));
        }
        globfree(&data);
    }
#endif

    return files;
}
//------------------------------------------------------------------------------
// Returns normalized path for the operating system
//------------------------------------------------------------------------------
std::string BuiltInFuncsUtils::Normpath(const std::string& path)
{
    if (path.empty()) return path;

    std::string normpath(path);

#ifdef OS_WIN
    std::replace(normpath.begin(), normpath.end(), '/', '\\');
#else
    std::replace(normpath.begin(), normpath.end(), '\\', '/');
#endif

    return normpath;
}
//------------------------------------------------------------------------------
// Returns environment variable
//------------------------------------------------------------------------------
std::string BuiltInFuncsUtils::GetEnv(const std::string& name)
{
    if (name.empty()) return "";

    char* cenv = getenv(name.c_str());
    
    return (cenv ? std::string(cenv) : "");
}
//------------------------------------------------------------------------------
// True for double values like Nan/Inf/-Inf, which ignore standard prinf format
//------------------------------------------------------------------------------
bool BuiltInFuncsUtils::NeedsSpecialPrintfFormat(double val)
{
    return (IsNaN_T(val) || IsInf_T(val) || IsNegInf_T (val));
}
//------------------------------------------------------------------------------
// Gets string for doubles like Nan/Inf/-Inf, which ignore standard prinf format
//------------------------------------------------------------------------------
std::string BuiltInFuncsUtils::GetSpecialPrintfFormat(const std::string& format,
                                                      double             val)
{
    std::string strtoprint;
    
    if (IsNaN_T(val))     
        strtoprint = "NaN";
    else if (IsInf_T(val))     
        strtoprint = "Inf";
    else if (IsNegInf_T (val)) 
        strtoprint = "-Inf";
    else
        return "";  // Does not need special formatting

    size_t pos = (!format.empty()) ? format.find('%') : std::string::npos;
    if (pos == std::string::npos)
        return strtoprint;

    // Get updated format which will replace %f,%g etc with %s
    std::string strformat;
    
    // Copy whatever string is before '%'
    if (pos != 0)
        strformat = format.substr(0, pos);

    strformat += "%s";

    // Copy whatever characters come after the alpha char in the format
    size_t len = format.length();
    size_t index = pos + 1;
    while (index < len)
    {
        if (isalpha(format[index])) break;
        index++;
    }

    if (index + 1 < len)
        strformat += format.substr(index+1);

    char buf[1028];
    sprintf(buf, strformat.c_str(), strtoprint.c_str());
    strtoprint = buf;
    return strtoprint;
}
//------------------------------------------------------------------------------
// Returns true if the given absolute path is a directory
//------------------------------------------------------------------------------
bool BuiltInFuncsUtils::IsDir(const std::string& path)
{
    if (path.empty()) 
    {
        return false;
    }
    std::string dir (Normpath(path));
    StripTrailingSlash(dir);
    struct stat file_stat;
    if (stat(dir.c_str(), &file_stat) == -1)
    {
        return false;
    }
    int returncode = 0;

#ifdef OS_WIN
    returncode = (file_stat.st_mode & _S_IFDIR);
#else
    returncode = S_ISDIR(file_stat.st_mode);
#endif

    if (returncode == 0)
    {
        return false;
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns base directory of the given path/file
//------------------------------------------------------------------------------
std::string BuiltInFuncsUtils::GetBaseDir(const std::string& path)
{
    std::string basedir;
    if (!path.empty())
    {
        size_t pos = path.find_last_of("/\\");
        if (pos != std::string::npos)
            basedir = path.substr(0, pos);
    }

    if (basedir.empty())
        basedir = GetCurrentWorkingDir();

    std::string tmp = GetAbsolutePath(basedir);
    if (!tmp.empty())
        basedir = tmp;

    return basedir;
}
//------------------------------------------------------------------------------
// Returns the base name if a file or the directory name for the given path
//------------------------------------------------------------------------------
std::string BuiltInFuncsUtils::GetBaseName(const std::string& path)
{
    if (path.empty()) return "";

    size_t pos = path.find_last_of("/\\");
    if (pos == std::string::npos)
        return path;

    return path.substr(pos+1);
}
//------------------------------------------------------------------------------
// Returns extension for a given file name
//------------------------------------------------------------------------------
std::string BuiltInFuncsUtils::GetFileExtension(const std::string& filename)
{
    if (filename.empty()) return "";

    size_t pos = filename.find_last_of(".");
    if (pos == std::string::npos)
        return "";

    return filename.substr(pos+1);
}
//------------------------------------------------------------------------------
// Returns true if there is a toolbox function of the given name
//------------------------------------------------------------------------------
bool BuiltInFuncsUtils::IsToolboxFunction(EvaluatorInterface eval,
                                          const std::string& name)
{
    if (name.empty()) 
        return false;

    FunctionInfo *fi   = nullptr;
    FUNCPTR       fptr = nullptr;
    eval.FindFunctionByName(name, &fi, &fptr);

    if (fptr)
        return true;

    return false;
}
//------------------------------------------------------------------------------
// Throws an error if the given file index is invalid
//------------------------------------------------------------------------------
void BuiltInFuncsUtils::CheckFileIndex(EvaluatorInterface  eval, 
                                       int                 fileid, 
                                       int                 argindex,
                                       bool                checkStdStreams)
{

    if (fileid < (checkStdStreams ? 0 : FIRST_USER_FILE) || 
        fileid >= eval.GetNumFiles()                     || 
        !eval.GetFile(fileid))
    {
        std::string msg ("Error: invalid file stream id [");
        msg += std::to_string(static_cast<long long>(fileid)) + "]";
        if (argindex >= 0)
        {
            msg += " specified in argument " + 
                   std::to_string(static_cast<long long>(argindex));
        }
        throw OML_Error(msg);
    }
}
//------------------------------------------------------------------------------
// Gets file id from input currency
//------------------------------------------------------------------------------
int BuiltInFuncsUtils::GetFileId(EvaluatorInterface eval,
                                 const Currency&    cur,
                                 int                idx)
{
    if (cur.IsString())
    {
        std::string nameInput = cur.StringVal();
        int         numfiles  = eval.GetNumFiles();
        for (int i = 0; i < numfiles; i++)
        {
            std::string fname = eval.GetFileName(i);
            if (!fname.empty() && fname == nameInput)
            {
                return i;
            }
        }
        return -1;
    }

    if (!cur.IsPositiveInteger())
    {
        throw OML_Error(OML_ERR_INTEGER, idx, OML_VAR_FILEID);
    }
    return static_cast<int>(cur.Scalar());
}
//------------------------------------------------------------------------------
// Returns true after parsing input and gets formats
//------------------------------------------------------------------------------
bool BuiltInFuncsUtils::GetFormats(EvaluatorInterface        eval,
                                   const std::string&        in,
                                   std::vector<std::string>& basefmts,
                                   std::vector<std::string>& rawfmts,
                                   bool                      checkSpec)
{
    if (in.empty())
    {
        basefmts.push_back("%f");
        rawfmts.push_back("%lf");  // Default option
        return true;
    }

    char* tok = strtok((char *)in.c_str(), "%");
    while (tok)
    {
        std::string tmp (tok);
        if (tmp.empty())
        {
            continue;
        }
        std::string basefmt;
        std::string rawfmt;

        int  numDecimalPts = 0;
        bool hasOnlyDigits = true;

        size_t len = tmp.size();
        for (size_t i = 0; i < len; ++i)
        {
            char ch = tmp[i];
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
                    if (checkSpec &&  (!hasOnlyDigits || numDecimalPts > 1))
                    {
                        rawfmt  = "";
                        SetWarning(eval, "Warning: ignoring incorrect format specifier");
                    }
                    rawfmt += "lf";
                }
                else if (ch == 'c' || ch == 'd' || ch == 'n')
                {
                    if (checkSpec &&  (!hasOnlyDigits || numDecimalPts > 0))
                    {
                        rawfmt  = "";
                        SetWarning(eval, "Warning: ignoring incorrect format specifier");
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

            if (checkSpec)
            {
                if (ch == '.')
                {
                    numDecimalPts++;
                }
                else if (!std::isdigit(static_cast<unsigned char>(ch)))
                {
                    hasOnlyDigits = false;
                }
            }
        }

        if (basefmt.empty() || rawfmt.empty())
        {
            return false;
        }

        basefmts.push_back("%" + basefmt);
        rawfmts.push_back("%"  + rawfmt);

        
        tok = strtok(NULL, "%");
    }

    // Check if there are tabs in the format
    size_t i = 0;
    for (std::vector<std::string>::iterator itr = rawfmts.begin();
         itr != rawfmts.end(); ++itr, ++i)
    {
        std::string fmt (*itr);
        size_t      pos     = fmt.find("\\t");
        bool        replace = false;
        if (pos != std::string::npos)
        {
            fmt.replace(pos, 2, "\t");
            rawfmts[i] = fmt;
        }
    }
    return true;
}
//------------------------------------------------------------------------------
// Adds trailing slash
//------------------------------------------------------------------------------
void BuiltInFuncsUtils::AddTrailingSlash(std::string& path)
{
    if (!path.empty())
    {
        if (path[path.size() - 1] != '\\' && path[path.size() - 1] != '/')
        {
            path += "/";
        }
        path = Normpath(path);
    }
}
//------------------------------------------------------------------------------
// Gets relative path
//------------------------------------------------------------------------------
std::string BuiltInFuncsUtils::GetRelativePath(const std::string& path,
                                               const std::string& basedir)
{
    size_t pos = path.find(basedir);
    if (pos != std::string::npos)
    {
        std::string relpath = path.substr(pos + basedir.length());
        return relpath;
    }
    return path;
}
//------------------------------------------------------------------------------
// Strips trailing newline
//------------------------------------------------------------------------------
void BuiltInFuncsUtils::StripTrailingNewline(std::string& str)
{
    while (!str.empty())
    {
        char last = str.back();
        if (last == '\n' || last == '\\r')
        {
            str.erase(str.end() - 1);
        }
        else
        {
            break;  // Reached the end
        }
    }
}
//------------------------------------------------------------------------------
// True if pagination environment is enabled
//------------------------------------------------------------------------------
bool BuiltInFuncsUtils::IsPaginationEnvEnabled()
{
    const char* paginationEnv = getenv("OML_PAGINATE");
    if (!paginationEnv) 
    {
        return true;  // By default it is enabled
    }

    if (paginationEnv == "0") 
    {
        return false;
    }

    std::string strEnv (paginationEnv);
    std::transform(strEnv.begin(), strEnv.end(), strEnv.begin(), ::tolower);

    if (strEnv == "false" || strEnv == "0") 
    {
        return false;
    }
    return true;
}
