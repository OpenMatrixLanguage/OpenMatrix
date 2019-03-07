/**
* @file BuiltInFuncsCore.cpp
* @date February 2016
* Copyright (C) 2016-2018 Altair Engineering, Inc.  
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

#include "BuiltInFuncsCore.h"

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <time.h>

#ifdef OS_WIN
#   include <Windows.h>
#else
#   include <unistd.h>
#   include <glob.h>
#   include <dirent.h>
#   include <dlfcn.h>
#   include <sys/resource.h>
#   include <sys/times.h>
#endif

#include "BuiltInFuncsUtils.h"
#include "CurrencyDisplay.h"
#include "ErrorInfo.h"
#include "Evaluator.h"
#include "FunctionInfo.h"
#include "Interpreter.h"
#include "OutputFormat.h"
#include "OML_Error.h"
#include "SignalHandlerBase.h"
#include "StructData.h"

#include "hwMatrixN.h"

#include "OMLInterface.h"

typedef int (*initModuleFP)(EvaluatorInterface);
typedef int (*initAltFP)(OMLInterface*);

//#define BUILTINFUNCSCORE_DBG 1  // Uncomment to print debug info
#if defined (BUILTINFUNCSCORE_DBG) && (OS_WIN) 
#   define BUILTINFUNCSCORE_DBG_ERR(s) {                                \
        std::cerr << "Error: Source [" << s << "]";                     \
        std::cerr << ", Error [" << GetLastError() << "]" << std::endl; }
#  else
#    define BUILTINFUNCSCORE_DBG_ERR(s) 0
#endif   
// End defines/includes

//------------------------------------------------------------------------------
// Returns true after getting user input. [input command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Input(EvaluatorInterface           eval, 
	                         const std::vector<Currency>& inputs, 
			                 std::vector<Currency>&       outputs)
{
    if (inputs.empty()) throw OML_Error(OML_ERR_NUMARGIN);

    size_t numInputs = inputs.size();

    const Currency &inputCur = inputs[0];
	if (!inputCur.IsString()) throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);

    std::string userInput;
	std::string prompt (inputCur.StringVal());
    std::string type;
    if (numInputs > 1)
    {
        const Currency &typeCur = inputs[1];
        if (typeCur.IsString())
            type = typeCur.StringVal();
    }
	
    Interpreter interp(eval);  // Create a new interp as this is in middle of an eval
    eval.OnUserInput(prompt, type, userInput);

    Currency result = interp.DoString(userInput);
    // Return empty matrix if there is an error or user did no input anything
    if (result.IsError() || result.IsNothing() || result.IsEmpty())
        result = eval.allocateMatrix(); 

	outputs.push_back(result);
	return true;
}
//------------------------------------------------------------------------------
// Returns true after pausing a script [pause command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Pause(EvaluatorInterface           eval, 
                             const std::vector<Currency>& inputs, 
                             std::vector<Currency>&       outputs)
{
    size_t numargs = inputs.empty() ? 0 : inputs.size();
    if (numargs > 1) throw OML_Error(OML_ERR_NUMARGIN);

    std::string msg;
    bool        wait     = true;
    double      waittime = 0.0;

    if (numargs == 0) // Pause till user presses a key
        msg += "Execution is paused. Press any key in the OML Command Window to continue...";

    else              // Pause for number of seconds specified by user
    {
        const Currency& in = inputs[0];
        if (!in.IsScalar() || !(in.Scalar() > 0.0)) 
            throw OML_Error(OML_ERR_POSITIVE_SCALAR, 1, OML_VAR_VALUE);

        waittime = in.Scalar();

        std::ostringstream os;
        os << "Execution is paused for " << waittime << " second(s)...";
        msg += os.str();
        wait = false;
    }

    eval.OnPauseStart(msg, wait);
    BuiltInFuncsCore funcs;
    if (!wait)
    {
        funcs.SleepFunc(waittime); 
        eval.OnPauseEnd();
    }

    return true;
}
//------------------------------------------------------------------------------
// Sleeps for the given amount of seconds
//------------------------------------------------------------------------------
void BuiltInFuncsCore::SleepFunc(double seconds) const
{
#ifdef OS_WIN
        Sleep((DWORD)(seconds * 1000)); // sleep specified amount of time
#else
        usleep(seconds * 1000000);
#endif
}
//------------------------------------------------------------------------------
// Returns true after specifying format options [format command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Format(EvaluatorInterface           eval,
                              const std::vector<Currency>& inputs,
                              std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.empty() ? 0 : inputs.size();
    
    // If there are no inputs, format is reset
    OutputFormat *format = new OutputFormat(*eval.GetOutputFormat());
    if (!format)
        format = new OutputFormat();

    assert(format);
    format->Reset();     // Reset format to get default values

    if (nargin == 0)     // Just return the default format
    {
        outputs.push_back(format);
        return true;
    }

    const Currency& in1 = inputs[0];
    bool  isString      = in1.IsString();
    if (!isString && !in1.IsPositiveInteger() && in1.Scalar() !=  static_cast<int>(0)) 
        throw OML_Error(OML_ERR_STRING_NATURALNUM, 1, OML_VAR_TYPE);

    if (isString) // Set precision
    {
        std::string precisionStr (in1.StringVal());
        if (precisionStr == "long")  
            format->SetLongFormat();

        else if (precisionStr == "short")
            format->SetShortFormat();
        else
            throw OML_Error(OML_ERR_FORMAT, 1, OML_VAR_TYPE);

        if (nargin >= 2) // Set format flag
        {
            const Currency& in2 = inputs[1];
            if (!in2.IsString()) throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);

            std::string type (in2.StringVal());
            if (type == "e")  
                format->SetFlags(std::ios::scientific);
            else if (type == "E")
                format->SetFlags(std::ios::scientific | std::ios_base::uppercase);
            else
                throw OML_Error(OML_ERR_FORMAT, 2, OML_VAR_TYPE);
        }

        outputs.push_back(format);
        return true;
    }
    if (nargin < 2) throw OML_Error(OML_ERR_NUMARGIN);
        
    const Currency& in2 = inputs[1];
    if (!in2.IsPositiveInteger() && in2.Scalar() < 0) 
    {
        throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_TYPE);
    }
    format->SetCustomFormat(static_cast<int>(in1.Scalar()), static_cast<int>(in2.Scalar()));

    outputs.push_back(format);
    return true;
}
//------------------------------------------------------------------------------
// Gets number of bytes used by Currency
//------------------------------------------------------------------------------
int BuiltInFuncsCore::GetsBytesUsed(const Currency& cur) const
{
    BuiltInFuncsCore funcs;
    if (cur.IsScalar())
        return (sizeof(double));
    
    else if (cur.IsComplex())
        return (2 * sizeof(double));

    else if (cur.IsString())
    {
        const hwMatrix* mtx = cur.Matrix();
        return (mtx ? mtx->Size() * sizeof(unsigned char) : 0);
    }
    else if (cur.IsMatrix())
    {
        const hwMatrix* mtx = cur.Matrix();
        if (!mtx) return 0;

        int realsize = mtx->Size() * sizeof(double);
        return (mtx->IsRealData() ? realsize : realsize * 2);
    }
    else if (cur.IsCellArray())
    {
        HML_CELLARRAY *cell = cur.CellArray();
        if (!cell) return 0;

        int sz    = cell->Size();
        int bytes = 0;
        for (int i = 0; i < sz; ++i)
            bytes += funcs.GetsBytesUsed((*cell)(i));
        
        return bytes;
    }
    else if (cur.IsStruct())
    {
        int         bytes = 0;
        StructData* sd    = cur.Struct();
        if (!sd) return 0;

        const std::map<std::string, int> fieldnames = sd->GetFieldNames();
        for (std::map<std::string, int>::const_iterator itr = fieldnames.begin(); 
             itr != fieldnames.end(); ++itr)
        {
            for (int i = 0; i < sd->M(); ++i)
            {
                for (int j = 0; j < sd->N(); ++j)
                    bytes += funcs.GetsBytesUsed(sd->GetValue(i, j, itr->first));
            }
        }
        return bytes;
    }
    else if (cur.IsBoundObject())
        return sizeof(cur.BoundObject());
    
    else if (cur.IsNDMatrix())
    {
        const hwMatrixN* mtxn  = cur.MatrixN();
        assert(mtxn);
        int bytes = mtxn->Size() * sizeof(double);
        return (mtxn->IsReal() ? bytes : bytes * 2);
    }

    return 0;
}
//------------------------------------------------------------------------------
// Returns true, gets information of variables defined in current scope [whos command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Whos(EvaluatorInterface           eval, 
	                        const std::vector<Currency>& inputs, 
			                std::vector<Currency>&       outputs)
{
    std::vector<std::string> names;  // vector of names
    if (inputs.empty())
        names = eval.GetVariableNames();
    else
    {
        int i = 1;
        for (std::vector<Currency>::const_iterator itr = inputs.begin();
             itr != inputs.end(); ++itr, ++i)
        {
            const Currency &input = *itr;
            if (!input.IsString()) 
                throw OML_Error(OML_ERR_STRING, i, OML_VAR_TYPE);

            names.push_back(input.StringVal());
        }
    }

    std::vector<std::string> attrs;   // vector of attributes            
    std::vector<std::string> dims;    // vector of sizes
    std::vector<std::string> bytes;   // vector of bytes
    std::vector<std::string> types;   // vector of class names
    
    size_t numvars = !names.empty() ? names.size() : 0;
    attrs.reserve(numvars);
    dims.reserve(numvars);
    bytes.reserve(numvars);
    types.reserve(numvars);

    // Check if the output needs to be printed or if it is returned as a struct
    StructData* sd = (eval.GetNargoutValue() != 0) ? new StructData() : 0;

    // Max widths for formatting table columns
    size_t wattr  = 4;
    size_t wname  = 4;
    size_t wsize  = 4;
    size_t wbytes = 5; 
    size_t wclass = 5;

    BuiltInFuncsCore funcs;
    std::vector<std::string>::const_iterator itr = names.begin();
    for (int i = 0; itr != names.end(); ++itr, ++i)
    {
        std::string name (*itr);
        Currency    var      = eval.GetValue(name);
        bool        isreal   = true;
        int         numbytes = funcs.GetsBytesUsed(var);
        std::string sz ("1x1");
        std::string cl (funcs.GetCurrencyClass(var));

        if (var.IsComplex())
            isreal = (var.Complex().Imag() == 0.0);

        else if (var.IsString())
        {
            const hwMatrix*  mtx = var.Matrix();
            assert(mtx);
            if (mtx)
            {                
                std::vector<int> dims;
                dims.push_back(mtx->M());
                dims.push_back(mtx->N());
                sz       = funcs.DimensionToString(dims);
                numbytes = mtx->Size() * sizeof(unsigned char);
            }
        }       
        else if (var.IsMatrix())
        {
            const hwMatrix*  mtx = var.Matrix();
            assert(mtx);
            if (mtx)
            {                
                std::vector<int> dims;
                dims.push_back(mtx->M());
                dims.push_back(mtx->N());
                sz       = funcs.DimensionToString(dims);
                isreal   = mtx->IsRealData();
                numbytes = isreal ? (mtx->Size() * sizeof(double)) :
                                    (mtx->Size() * sizeof(double)) * 2;
            }
        }
        else if (var.IsCellArray())
        {
            assert(var.CellArray());

            std::vector<int> dims;
            dims.push_back(var.CellArray()->M());
            dims.push_back(var.CellArray()->N());
            sz = funcs.DimensionToString(dims);
        }
        else if (var.IsStruct())
        {
            assert(var.Struct());

            std::vector<int> dims;
            dims.push_back(var.Struct()->M());
            dims.push_back(var.Struct()->N());
            sz = funcs.DimensionToString(dims);
        }
        else if (var.IsNDMatrix())
        {
            const hwMatrixN* mtxn = var.MatrixN();
            sz     = funcs.DimensionToString(mtxn->Dimensions());
            isreal = mtxn->IsReal();
        }

        if (sd)
        {
            // Set the value in the struct and don't print result
            sd->SetValue(i, 0, "complex", isreal ? 0.0 : 1.0);
            sd->SetValue(i, 0, "name",    name);
            sd->SetValue(i, 0, "bytes",   numbytes);
            sd->SetValue(i, 0, "class",   cl);
            sd->SetValue(i, 0, "size",    sz);
            continue;
        }

        // Add to results that need to be printed
        std::string attr = isreal ? "" : "c"; //\todo: Implement other attrs
        std::string byte = std::to_string(static_cast<long long>(numbytes));

        attrs.push_back(attr);
        dims.push_back(sz);
        bytes.push_back(byte);
        types.push_back(cl);
        
        if (!name.empty() && wname < name.size())
            wname = name.size();

        if (!sz.empty() && wsize < sz.size())
            wsize = sz.size();

        if (!byte.empty() && wbytes < byte.size())
            wbytes = byte.size();

        if (!cl.empty() && wclass < cl.size())
            wclass = cl.size();
    }

    if (sd)
    {
        outputs.push_back(sd);
        return true;
    }

    std::string b1 ("====");
    std::string b2 ("=====");

    std::stringstream ss;
    // Labels
    ss << " " << std::left  << std::setfill(' ') << std::setw(wattr)  << "Attr";
    ss << " " << std::left  << std::setfill(' ') << std::setw(wname)  << "Name";
    ss << " " << std::left  << std::setfill(' ') << std::setw(wsize)  << "Size";
    ss << " " << std::right << std::setfill(' ') << std::setw(wbytes) << "Bytes";
    ss << " " << std::left  << std::setfill(' ') << std::setw(wclass) << "Class";
    ss << std::endl;

    // Borders
    ss << " " << std::left  << std::setfill(' ') << std::setw(wattr)  << b1; 
    ss << " " << std::left  << std::setfill(' ') << std::setw(wname)  << b1; 
    ss << " " << std::left  << std::setfill(' ') << std::setw(wsize)  << b1; 
    ss << " " << std::right << std::setfill(' ') << std::setw(wbytes) << b2;
    ss << " " << std::left  << std::setfill(' ') << std::setw(wclass) << b2;
    ss << std::endl;

    // Data
    for (size_t i = 0; i < numvars; ++i)
    {
        ss << " " << std::left  << std::setfill(' ') << std::setw(wattr)  << attrs[i];
        ss << " " << std::left  << std::setfill(' ') << std::setw(wname)  << names[i];
        ss << " " << std::left  << std::setfill(' ') << std::setw(wsize)  << dims[i];
        ss << " " << std::right << std::setfill(' ') << std::setw(wbytes) << bytes[i];
        ss << " " << std::left  << std::setfill(' ') << std::setw(wclass) << types[i];
        ss << std::endl;
    }
    ss << std::endl;
    Currency out (ss.str());
    out.DispOutput();
    eval.PrintResult(out);

    return true;
}
//------------------------------------------------------------------------------
// Converts dimensions to string - whos helper method
//------------------------------------------------------------------------------
std::string BuiltInFuncsCore::DimensionToString( const std::vector<int>& dims) const
{
    std::string sz;
    for (std::vector<int>::const_iterator itr = dims.begin(); 
        itr != dims.end(); ++itr)
    {
        if (!sz.empty())
            sz += "x";

        sz += std::to_string(static_cast<long long>(*itr));
    }
    return sz;
}
//------------------------------------------------------------------------------
// Returns true after getting class information [class command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Class(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (eval.GetNargoutValue() > 1)
        throw OML_Error(OML_ERR_NUMARGOUT);

    BuiltInFuncsCore funcs;
    std::string strclass (funcs.GetCurrencyClass(inputs[0]));
    if (strclass.empty())
        throw OML_Error(HW_ERROR_UNKOWNTYPE);

    outputs.push_back(strclass);
    return true;
}
//------------------------------------------------------------------------------
// Gets class for given currency
//------------------------------------------------------------------------------
std::string BuiltInFuncsCore::GetCurrencyClass(const Currency& cur) const
{
    if (cur.IsScalar())
        return (cur.IsLogical() ? "logical" : "double");
    
    else if (cur.IsString())
        return "char";
    
    else if (cur.IsComplex())
        return "double";

    else if (cur.IsMatrix())
        return (cur.IsLogical() ? "logical" : "double");

    else if  (cur.IsNDMatrix())
        return (cur.IsLogical() ? "logical" : "double");

    else if (cur.IsFunctionHandle())
        return "function_handle";
    
    else if (cur.IsCellArray())
        return "cell";
    
    else if (cur.IsStruct())
        return "struct";
    
    else if (cur.IsObject() || cur.IsBoundObject())
		return cur.GetClassname();

    return "";
}
//------------------------------------------------------------------------------
// Returns true after getting class information [class command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::IsA(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (eval.GetNargoutValue() > 1)
        throw OML_Error(OML_ERR_NUMARGOUT);

	Currency    obj          = inputs[0];
	Currency    class_type   = inputs[1];
	std::string target_class;

	if (class_type.IsString())
		target_class = class_type.StringVal();
	else
		throw OML_Error("");

	Currency val(0.0);
	val.SetMask(Currency::MASK_LOGICAL);

	if (obj.IsObject())
	{
		if (obj.GetClassname() == target_class)
			val.ReplaceScalar(1.0);
	}
	else if (obj.IsMatrix() || obj.IsScalar() || obj.IsNDMatrix() || obj.IsComplex())
	{
		if ((target_class == "numeric") || (target_class == "float"))
			val.ReplaceScalar(1.0);
	}
		
	outputs.push_back(val);

    return true;
}
//------------------------------------------------------------------------------
// Dynamically loads a library and returns handle to it
//------------------------------------------------------------------------------
void* BuiltInFuncsCore::DyLoadLibrary(const std::string& name)
{
#ifdef OS_WIN
    return (void*)LoadLibrary(name.c_str());
#else  // OS_UNIX
    return dlopen(name.c_str(), RTLD_GLOBAL|RTLD_LAZY);
#endif  // OS_UNIX
}
//------------------------------------------------------------------------------
// Returns function pointer from a dynamically loaded library
//------------------------------------------------------------------------------
void* BuiltInFuncsCore::DyGetFunction(void* handle, const std::string& name)
{
    if (!handle || name.empty())
        return nullptr;

#ifdef OS_WIN
    return GetProcAddress((HINSTANCE)handle, name.c_str());
#else  
    return dlsym(handle, name.c_str());
#endif
}
//------------------------------------------------------------------------------
// Release dynamically loaded library
//------------------------------------------------------------------------------
void BuiltInFuncsCore::DyFreeLibrary(void* handle)
{
#ifdef OS_WIN
    FreeLibrary((HINSTANCE)handle);
#else  // OS_UNIX
    dlclose(handle);
#endif  // OS_UNIX
}


//------------------------------------------------------------------------------
// Returns true after adding a toolbox [addtoolbox command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::AddToolbox(EvaluatorInterface           eval, 
	                              const std::vector<Currency>& inputs, 
			                      std::vector<Currency>&       outputs)
{
    if (inputs.empty()) throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& in = inputs[0];
    if (!in.IsString())
        OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);

    std::string importDll (in.StringVal());

#ifndef OS_WIN
    importDll = "lib" + importDll + ".so";
#endif

     void *vResult = BuiltInFuncsCore::DyLoadLibrary(importDll.c_str());

    if (!vResult)
    {
        BUILTINFUNCSCORE_DBG_ERR(importDll);
	    throw OML_Error(OML_ERR_INVALID_DLL);
    }

    void* symbol = BuiltInFuncsCore::DyGetFunction(vResult, "InitDll");
    if (!symbol)  // Check if this is being loaded using oml wrappers
    {
        symbol = BuiltInFuncsCore::DyGetFunction(vResult, "InitDllOmlWrap");
    }
    else // Additional version check for toolboxes
    {
        // Adding additional check only for omlziptoolbox now
        if (!DyGetFunction(vResult, "GetToolboxVersion"))
        {
            throw OML_Error(OML_MSG_INVALID_VERSION + std::string(" in ") +
                            importDll);  
        }
    }

	if (!symbol)	
	{
		symbol = BuiltInFuncsCore::DyGetFunction(vResult, "InitToolbox");

		if (symbol)
		{
			initAltFP fp = (initAltFP)symbol;
			OMLInterfaceImpl impl(&eval);

			eval.SuspendFuncListUpdate();   

			fp(&impl);

			eval.UnsuspendFuncListUpdate();
			eval.OnUpdateFuncList();
			return true;
		}
	}

    initModuleFP fp = (initModuleFP)symbol;
    if (!fp)
    {
        BUILTINFUNCSCORE_DBG_ERR(importDll);   
        throw OML_Error(OML_ERR_INVALID_INITDLL);
    }
    // Suspend updating funcs till all functions in the toolbox are registered
    eval.SuspendFuncListUpdate();   

    fp(eval);

    // Unsuspend and add new functions registered in toolbox
    eval.UnsuspendFuncListUpdate();
    eval.OnUpdateFuncList();
    return true;
}
//------------------------------------------------------------------------------
// Returns true after displaying a count of functions [funccount command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Funccount(EvaluatorInterface           eval, 
                                 const std::vector<Currency>& inputs, 
                                 std::vector<Currency>&       outputs)
{
    outputs.push_back(eval.GetFunctionCount());
    return true;
}
//------------------------------------------------------------------------------
// Returns true after displaying a list of functions along with total count [funclist command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Funclist(EvaluatorInterface           eval, 
                                const std::vector<Currency>& inputs, 
                                std::vector<Currency>&       outputs)
{
    std::vector<std::string> funclist (eval.GetFunctionNames());
    std::string cmpstr;

    if (!inputs.empty())
    {
        const Currency& in = inputs[0];
        if (!in.IsString()) throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);

        cmpstr = in.StringVal();
    }

    int         numfuncs = 0;
    bool        addspace = false;
    bool        cmpname  = (!cmpstr.empty());
    std::string out;
    for (std::vector<std::string>::const_iterator itr = funclist.begin();
         itr != funclist.end(); ++itr)
    {
        std::string thisfunc (*itr);
        if (cmpname && thisfunc.find(cmpstr) == std::string::npos) continue;
        
        if (addspace)
            out += " ";

        out += thisfunc;
        numfuncs ++;
        addspace = true;
    }

    out += "\nFunctions displayed: " + std::to_string(static_cast<long long>(numfuncs));

    Currency currOut (out);
    currOut.DispOutput();
    outputs.push_back(currOut);
    return true;
}
//------------------------------------------------------------------------------
// Returns true after displaying a list of variables in the interpreter [varlist command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Varlist(EvaluatorInterface           eval, 
                               const std::vector<Currency>& inputs, 
                               std::vector<Currency>&       outputs)
{
    std::vector<std::string> varlist (eval.GetVariableNames());

    std::string msg;
    for (std::vector<std::string>::const_iterator itr = varlist.begin();
         itr != varlist.end(); ++itr)
    {
        std::string var (*itr);
        if (!var.empty())
            msg += var + " ";
    }
    Currency cur (msg);
    cur.DispOutput();
    outputs.push_back(cur);
    return true;
}
//------------------------------------------------------------------------------
// Gets build number from the given file
//------------------------------------------------------------------------------
std::string BuiltInFuncsCore::GetBuildNumber(const std::string& versionfile)
{
    if (versionfile.empty()) return "";

    std::ifstream ifs;
    ifs.open(versionfile);
    if (!ifs) 
        return "";

    std::string buildnumber;
    std::string search ("BuildNumber value = \"");
    while (1)
    {
        std::string thisline;
        if (std::getline(ifs, thisline).eof()) break;

        if (thisline.empty()) continue;

        size_t pos = thisline.find(search);
        if (pos == std::string::npos) continue;

        buildnumber = thisline.substr(pos + search.length());

        // Strip trailing quote
        if (!buildnumber.empty())
        {
            pos = buildnumber.find("\"");
            if (pos != std::string::npos)
                buildnumber = buildnumber.substr(0, pos);
        }
        break;
    }
    ifs.close();

    if (!buildnumber.empty())
       return (" Build " + buildnumber);

    return "";
}
//------------------------------------------------------------------------------
// Returns true after getting current file name being processed [omlfilename command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Omlfilename(EvaluatorInterface           eval, 
                                   const std::vector<Currency>& inputs, 
                                   std::vector<Currency>&       outputs)
{
    std::string option;

    size_t nargin = inputs.size();
    if (nargin > 0)
    {
        Currency cur = inputs[0];
        if (!cur.IsString())
            throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
        option = cur.StringVal();
        if (option != "fullpath" && option != "fullpathext")
            throw OML_Error(OML_ERR_OPTION, 1, OML_VAR_VALUE);
    }

    std::string currfile (eval.GetCurrentFilename());
    bool        isbeingedited = false;
    if (currfile.empty() || currfile == "dummy")
    {
        SignalHandlerBase* handler = eval.GetSignalHandler();
        assert(handler);
        currfile = handler ? handler->GetEditorFileName() : "";
        isbeingedited = true;
    }

    if (currfile.empty() || (!isbeingedited && currfile == "dummy"))
    {
        outputs.push_back("");
        return true;
    }

    std::string filename;      // Only the file name
    std::string fullpath;      // Path with file name and no extension
    std::string fullpathext;   // Path with file name and extension
    
    // Check if there is an absolute path for the file name
    currfile   = BuiltInFuncsUtils::Normpath(currfile);
    size_t pos = currfile.find_last_of("/\\");
    if (pos == std::string::npos)
    {
        fullpathext = BuiltInFuncsUtils::GetCurrentWorkingDir() + "/" + currfile;
        filename    = currfile;
    }
    else
    {
        fullpathext = currfile;
        filename = currfile.substr(pos+1);
    }

    fullpathext = BuiltInFuncsUtils::Normpath(fullpathext);
    if (option.empty())
    {
        pos = filename.find_last_of(".");
        if (pos != std::string::npos)
            filename = filename.substr(0, pos);
        outputs.push_back(filename);
    }
    else if (option == "fullpath")
    {
        pos = fullpathext.find_last_of(".");
        if (pos != std::string::npos)
            fullpathext = fullpathext.substr(0, pos);
        outputs.push_back(fullpathext);
    }
    else
        outputs.push_back(fullpathext);

    return true;
}
//------------------------------------------------------------------------------
// Returns true after checking if given variable/file/path/function exists [exist command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Exist(EvaluatorInterface           eval, 
                             const std::vector<Currency>& inputs, 
                             std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.empty() ? 0 : inputs.size();
    if (nargin < 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    Currency c1 = inputs[0];
    if (!c1.IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }

    std::string name (c1.StringVal());
    if (name.empty())
    {
        outputs.push_back(0);  // Name does not exist
        return true;
    }

    std::string type;
    bool        notype     = true;
    int         returncode = 0; // Does not exist
    if (nargin > 1)
    {
        Currency c2 = inputs[1];
        if (!c2.IsString())
        {
            throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
        }
        type = c2.StringVal();
        if (type == "var" || type == "builtin" || type == "file" || type == "dir")
        {
            notype = false; // Valid options
        }
        else if (type == "class")
        {
            type = "";      // \todo: Ignored for now, implement later
        }
        else
        {
            BuiltInFuncsUtils::SetWarning(eval, 
                "Warning: invalid option [" + type + "] in argument 2");
            type = "";
        }
    }

    // Check if this is a variable
    if (notype || type == "var")
    {
        returncode = (eval.Contains(name)) ? 1 : 0;   // Variable
        if (returncode != 0 || !notype)
        { 
            outputs.push_back(returncode);  
            return true;
        }
    }

    FunctionInfo *fi = nullptr;
    if (notype || type == "builtin")
    {
        FUNCPTR fptr = nullptr;
        eval.FindFunctionByName(name, &fi, &fptr);
        returncode = fptr ? 5 : 0;    // Built-in function
        if (returncode != 0 || !notype)
        {
            outputs.push_back(returncode);  
            return true;
        }
    }

    BuiltInFuncsUtils::StripTrailingSlash(name);
#ifdef OS_WIN
    if (name.size() == 2 && name[1] == ':')
        name += "/";            // Add trailing slash if this is a drive letter
#endif

	std::string absIfPath = BuiltInFuncsUtils::GetAbsolutePath(name);
    if (notype || type == "dir" || type == "file")
    {
        returncode = (BuiltInFuncsUtils::IsDir(absIfPath)) ? 7 : 0; // Directory
        if (returncode != 0 || (!notype && type != "file"))
        {
            outputs.push_back(returncode);  
            return true;
        }
    }

    if (notype || type == "file")
    {
        returncode = BuiltInFuncsUtils::FileExists(absIfPath) ? 2 : 0; // File
        if (returncode != 0 || !notype)
        {
            outputs.push_back(returncode);  
            return true;
        }
    }

    if (fi)
    {
        outputs.push_back(103); // Any other function
        return true;
    }

    outputs.push_back(returncode);
    return true;
}
//------------------------------------------------------------------------------
// Returns true and returns a description of type of the given input [type command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Type(EvaluatorInterface           eval, 
                            const std::vector<Currency>& inputs, 
                            std::vector<Currency>&       outputs)
{
    if (inputs.empty())
        throw OML_Error(OML_ERR_NUMARGIN);

    size_t nargin = inputs.size();

    // Check for quiet option, input validity
    std::map<int, std::string> in;
    int  i     = 1;
    bool quiet = false;
    for (std::vector<Currency>::const_iterator itr = inputs.begin(); 
         itr != inputs.end(); ++itr, ++i)
    {
        Currency tmp = *itr;
        if (!tmp.IsString())
            throw OML_Error(OML_ERR_STRING, i, OML_VAR_TYPE);

        std::string val (tmp.StringVal());
        if (val.empty())
            throw OML_Error(OML_ERR_NONEMPTY_STR, i, OML_VAR_VALUE);

        const hwMatrix* mat = tmp.Matrix();
        if (mat && mat->M() != 1)
            throw OML_Error(OML_ERR_ONEROW_STRING, i, OML_VAR_VALUE);

        if (val == "-q" || val == "-Q")
            quiet = true;

        else
            in[i] = val;
    }
    int nargout = eval.GetNargoutValue();
    if (in.empty())
    {
        if (nargout > 0)
            outputs.push_back(EvaluatorInterface::allocateCellArray());
        return true;
    }

    HML_CELLARRAY* cell = (nargout > 0) ? 
        EvaluatorInterface::allocateCellArray((int)nargin, 1) :
        nullptr;
    
    i = 0;
    std::string reservedtokens ("[]{}();,%\"");
    for (std::map<int, std::string>::const_iterator itr = in.begin();
         itr != in.end(); ++itr, ++i)
    {
        std::string str = itr->second;
        std::string type;
        if (eval.IsOperator(str))
            type = "'" + str + "' is an operator";

        else if (eval.IsKeyword(str))
            type = "'" + str + "' is a keyword";

        else if (eval.IsStdFunction(str))
            type = "'" + str + "' is a built-in function";
       
        else if (eval.IsUserFunction(str))
            type = "'" + str + "' is a user-defined function";
        
        else if (reservedtokens.find(str) != std::string::npos)
            type = "'" + str + "' is a reserved token";

        else
        {
            const Currency& val = eval.GetValue(str);
            if (!val.IsNothing())
               type = "'" + str + "' is a variable\n" + val.GetOutputString(eval.GetOutputFormat());
            else 
            {
                std::string fname (str);
                if (BuiltInFuncsUtils::IsDir(fname))
                {
                    throw OML_Error("Error: invalid input in argument " +
                        std::to_string(static_cast<long long>(itr->first)) + 
                        "; [" + BuiltInFuncsUtils::Normpath(fname) + 
                        "] is a directory and is not defined");
                }

                if (!BuiltInFuncsUtils::FileExists(fname) && 
                     BuiltInFuncsUtils::GetFileExtension(fname).empty())
                    fname += ".oml";
                bool foundfunc = eval.FindFileInPath(fname, fname);
                if (!quiet && foundfunc)
                    type = "'" + str + "' is a function defined in [" + 
                           BuiltInFuncsUtils::Normpath(fname) + "]";

                if (!foundfunc && !BuiltInFuncsUtils::FileExists(fname))
                {
                    throw OML_Error("Error: invalid input in argument "    +
                        std::to_string(static_cast<long long>(itr->first)) + 
                        "; '" + str + "' is not defined");
                }
                std::ifstream ifs(fname.c_str());
                if (ifs.fail()) 
                {
                    throw OML_Error("Error: invalid input in argument "    +
                        std::to_string(static_cast<long long>(itr->first)) + 
                        "; unable to open file for reading [" + 
                        BuiltInFuncsUtils::Normpath(fname) + "]");
                }
                // Read the entire file
                while (!ifs.eof())
                {
                    std::string thisline;
                    std::getline(ifs, thisline);

                    if (!type.empty())
                        type += "\n";
                    type += thisline;
                }
                ifs.close();
            }
        }

        Currency tmp(type);

        if (cell)
            (*cell)(i, 0) = tmp;
        else
        {
            tmp.DispOutput();
            eval.PrintResult(tmp);
        }
    }

    if (cell)
        outputs.push_back(cell);
    return true;
}
//------------------------------------------------------------------------------
// Enables/disables pagination (omlpaginate)
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::OmlPaginate(EvaluatorInterface           eval, 
                                   const std::vector<Currency>& inputs, 
                                   std::vector<Currency>&       outputs)
{
    size_t numargin = inputs.empty() ? 0 : inputs.size();
	if (numargin != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    Currency cur = inputs[0];    
    if (!cur.IsLogical() && !cur.IsInteger())  
    {
        throw OML_Error(OML_ERR_LOGICAL, 1, OML_VAR_TYPE);
    }
    int tmp = static_cast<int>(cur.Scalar());
    if (tmp != 0 && tmp != 1)
    {
        throw OML_Error(OML_ERR_FLAG_01, 1, OML_VAR_VALUE);
    }

    bool val = (tmp == 0) ? false : true;
    CurrencyDisplay::SetPaginate(val);
    return true;
}
//------------------------------------------------------------------------------
// Sleeps for the given time period (sleep)
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::OmlSleep(EvaluatorInterface           eval, 
                             const std::vector<Currency>& inputs, 
                             std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    Currency cur = inputs[0];
    if (!cur.IsScalar())
    {
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);
    }

    BuiltInFuncsCore funcs;
    funcs.SleepFunc(cur.Scalar());
    return true;
}
