/**
* @file BuiltInFuncsCore.cpp
* @date February 2016
* Copyright (C) 2016-2023 Altair Engineering, Inc.  
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
#include <memory>
#include <sstream>
#include <time.h>

#ifdef OS_WIN
#   include <Windows.h>
#   include <psapi.h>    // MS SDK header, PPROCESS_MEMORY_COUNTERS
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

#include "hwMatrixN_NMKL.h"
#include "hwMatrixS_NMKL.h"

#include "OMLInterface.h"

typedef int (*initModuleFP)(EvaluatorInterface);
typedef void (*exitModuleFP)();
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

std::map<std::string, void*> BuiltInFuncsCore::libs;

//------------------------------------------------------------------------------
// Returns true after getting user input. [input command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Input(EvaluatorInterface           eval, 
	                         const std::vector<Currency>& inputs, 
			                 std::vector<Currency>&       outputs)
{
    size_t nargin = (!inputs.empty()) ? inputs.size() : 0;
    if (nargin < 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    else if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }
    
    std::string cmd;
    SignalHandlerBase* handler = eval.GetSignalHandler();
    if (handler && handler->IsWaiting(cmd))
    {
        BuiltInFuncsUtils::SetWarning(eval, 
            "Warning: cannot run [input] as [" + cmd + "] is running");
        outputs.emplace_back(EvaluatorInterface::allocateMatrix());
        return true;
    }

    std::string userInput;
	std::string prompt (inputs[0].StringVal());

    bool       stringInput = (nargin > 1 && inputs[1].IsMatrixOrString());
    std::string type       = (stringInput) ? "s" : "";
	
    // Create a new interp as this is in middle of an eval
    Interpreter interp(eval);  
    eval.OnUserInput(prompt, type, userInput);

    Currency result = interp.DoString(userInput);

    if (stringInput)
    {
        if (result.IsNothing() || result.IsEmpty() || result.IsError())
        {
            result = Currency(userInput); // Just show what the user typed
        }
    }
    else if (result.IsNothing() || result.IsEmpty() || result.IsError())
    {
        // Return empty matrix if there is an error or user did no input anything
        result = EvaluatorInterface::allocateMatrix();
    }

	outputs.emplace_back(result);
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

        if (eval.GetVerbose() > 0)
        {
            msg += "Execution is paused for ";
            if (in.IsInteger())
            {
                msg += std::to_string(static_cast<long long>(waittime));
            }
            else
            {
                msg += std::to_string(static_cast<long double>(waittime));
            }
            msg += " second(s)...";
        }
        wait = false;
    }

    if (wait)
    {
        std::string cmd;
        SignalHandlerBase* handler = eval.GetSignalHandler();
        if (handler && handler->IsWaiting(cmd))
        {
            BuiltInFuncsUtils::SetWarning(eval,
                "Warning: cannot run [pause] as [" + cmd + "] is running");
            return true;
        }
    }
    if (!msg.empty())
    {
        Currency tmp(msg);
        tmp.DispOutput();
        eval.PrintResult(tmp);
        msg = "";
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
        throw OML_Error(OML_ERR_STRING_NATURALNUM, 1);

    if (isString) // Set precision
    {
        std::string precisionStr (in1.StringVal());
        if (precisionStr == "long")  
            format->SetLongFormat();

        else if (precisionStr == "short")
            format->SetShortFormat();
        else
            throw OML_Error(OML_ERR_FORMAT, 1);

        if (nargin >= 2) // Set format flag
        {
            const Currency& in2 = inputs[1];
            if (!in2.IsString()) throw OML_Error(OML_ERR_STRING, 2);

            std::string type (in2.StringVal());
            if (type == "e")  
                format->SetFlags(std::ios::scientific);
            else if (type == "E")
                format->SetFlags(std::ios::scientific | std::ios_base::uppercase);
            else
                throw OML_Error(OML_ERR_FORMAT, 2);
        }

        outputs.push_back(format);
        return true;
    }
    if (nargin < 2) throw OML_Error(OML_ERR_NUMARGIN);
        
    const Currency& in2 = inputs[1];
    if (!in2.IsPositiveInteger() && in2.Scalar() < 0) 
    {
        throw OML_Error(OML_ERR_NATURALNUM, 2);
    }
    format->SetCustomFormat(static_cast<int>(in1.Scalar()), static_cast<int>(in2.Scalar()));

    outputs.push_back(format);
    return true;
}
//------------------------------------------------------------------------------
// Gets number of bytes used by Currency
//------------------------------------------------------------------------------
unsigned long long BuiltInFuncsCore::GetsBytesUsed(const Currency& cur) const
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
        unsigned long long bytes = 0;
        for (int i = 0; i < sz; ++i)
            bytes += funcs.GetsBytesUsed((*cell)(i));
        
        return bytes;
    }
	else if (cur.IsNDCellArray())
	{
		HML_ND_CELLARRAY* cells = cur.CellArrayND();
		if (!cells) return 0;

		int sz = cells->Size();
		unsigned long long bytes = 0;
		for (int i = 0; i < sz; ++i)
			bytes += funcs.GetsBytesUsed((*cells)(i));

		return bytes;
	}
    else if (cur.IsStruct())
    {
        unsigned long long  bytes = 0;
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
// Returns true, gets information of variables defined in current scope [whos]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Whos(EvaluatorInterface           eval, 
	                        const std::vector<Currency>& inputs, 
			                std::vector<Currency>&       outputs)
{
    std::vector<std::string> names;  // vector of names
    if (inputs.empty())
    {
        names = eval.GetVariableNames();
    }
    else
    {
        int i = 1;
        names.reserve(inputs.size());
        for (std::vector<Currency>::const_iterator itr = inputs.begin();
             itr != inputs.end(); ++itr, ++i)
        {
            const Currency &input = *itr;
            if (!input.IsString()) 
                throw OML_Error(OML_ERR_STRING, i, OML_VAR_TYPE);

            names.emplace_back(input.StringVal());
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
        Currency var = eval.GetValue(name);
        if (var.IsNothing())
        {
            // Create a second interpreter so outputs don't get mixed up
            Interpreter interp2(eval); 
            const Currency& var2 = interp2.DoString(name); // Try and evaluate
            if (!(var2.IsNothing() || var2.IsEmpty() || var2.IsError()))
            {
                var = var2;
            }
        }
        bool        isreal   = true;
        size_t      numbytes = funcs.GetsBytesUsed(var);
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
		else if (var.IsNDCellArray())
		{
			HML_ND_CELLARRAY* cells = var.CellArrayND();
			sz = funcs.DimensionToString(cells->Dimensions());
		}
        else if (var.IsSparse())
        {
            const hwMatrixS* mtx = var.MatrixS();
            if (mtx)
            {
                sz = std::to_string(static_cast<long long>(mtx->M())) + "x" +
                    std::to_string(static_cast<long long>(mtx->N()));
                isreal = mtx->IsReal();
                numbytes = isreal ? (mtx->NNZ() * sizeof(double)) :
                    (mtx->NNZ() * sizeof(double)) * 2;

            }
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

        sz += std::to_string(*itr);
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
        throw OML_Error(HW_ERROR_UNKNOWNTYPE);

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

	else if (cur.IsNDMatrix())
		return (cur.IsLogical() ? "logical" : "double");

	else if (cur.IsFunctionHandle())
		return "function_handle";

	else if (cur.IsCellArray())
		return "cell";

    else if  (cur.IsNDMatrix() || cur.IsSparse())
        return (cur.IsLogical() ? "logical" : "double");

	else if (cur.IsStruct())
		return "struct";

	else if (cur.IsObject() || cur.IsBoundObject())
		return cur.GetClassname();

	else if (cur.IsPointer())
		return cur.Pointer()->GetClassname();

    return "";
}
//------------------------------------------------------------------------------
// Returns true after getting class information [isa command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::IsA(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs, 
                           std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[1].IsString() && !inputs[1].IsCellArray())
    {
        throw OML_Error(OML_ERR_STRING_STRINGCELL, 2);
    }
    std::unique_ptr<HML_CELLARRAY> type = nullptr;
    if (inputs[1].IsCellArray())
    {
        type.reset(EvaluatorInterface::allocateCellArray(inputs[1].CellArray()));
    }
    else
    {
        type.reset(EvaluatorInterface::allocateCellArray(1,1));
        (*type)(0) = inputs[1].StringVal();
    }
    std::unique_ptr<HML_CELLARRAY> out(EvaluatorInterface::allocateCellArray(
        type->M(), type->N()));

    int numtypes = type->Size();
    
    const Currency& cur1 = inputs[0];

	std::string classname = cur1.GetClassname();
	bool  isobj = false;

	if (classname.length())
		isobj = true;

    for (int i = 0; i < numtypes; ++i)
    {
        Currency element = (*type)(i);
        if (!element.IsString())
        {
            throw OML_Error(OML_ERR_STRING_STRINGCELL, 2);
        }
        std::string cmp(element.StringVal());
        bool result = false;
        
        if (isobj)
        {
            result = (classname == cmp);
            if (!result) // Check if this is a derived class
            {
                result = eval.IsA(cur1, cmp);
            }
        }
        else if (!cmp.empty())
        {
            std::transform(cmp.begin(), cmp.end(), cmp.begin(), ::tolower);

            if (cmp == "integer")
            {
                result = (cur1.IsInteger() || cur1.IsLogical());
            }
            else if (cmp == "float" || cmp == "double" || cmp == "numeric")
            {
                result = (cur1.IsScalar()   || cur1.IsMatrix() || 
                          cur1.IsNDMatrix() || cur1.IsSparse() || 
                          cur1.IsComplex()  || cur1.IsLogical());
            }
			else if (cmp == "struct")
			{
				result = cur1.IsStruct();
			}
        }
        Currency tmp(result);
        tmp.SetMask(Currency::MASK_LOGICAL);
        (*out)(i) = tmp;
    }

    if (out->Size() == 1)
    {
        outputs.push_back(Currency((*out)(0)));
    }
    else
    {
        outputs.push_back(out.release());
    }
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
    return dlopen(name.c_str(), RTLD_LAZY);
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
// Get filename of a dynamically loaded library
//------------------------------------------------------------------------------
std::string BuiltInFuncsCore::DyGetLibraryPath(void* handle)
{
	char buff[2048];

	if (handle)
	{
#ifdef OS_WIN
		GetModuleFileName((HMODULE)handle, buff, sizeof(buff));
#else  // OS_UNIX
#include <dlfcn.h>
		Dl_info dli;
	    	dladdr(handle, &dli);
		strcpy(buff, dli.dli_fname);
#endif  // OS_UNIX
	}

	return buff;
}
//------------------------------------------------------------------------------
// Release dynamically loaded library
//------------------------------------------------------------------------------
void BuiltInFuncsCore::DyFreeLibrary(void* handle)
{
    if (!handle)
        return;

#ifdef OS_WIN
    FreeLibrary((HINSTANCE)handle);
#else  // OS_UNIX
    dlclose(handle);
#endif  // OS_UNIX
}
//------------------------------------------------------------------------------
// Loads the "Finalize" function and executes it
//------------------------------------------------------------------------------
void BuiltInFuncsCore::Finalize(void* handle)
{
	if (!handle)
		return;

	void* symbol = DyGetFunction(handle, "Finalize");

	if (symbol)
	{
		exitModuleFP fp = (exitModuleFP)symbol;
		fp();
	}
}

//------------------------------------------------------------------------------
// Returns true after adding a toolbox [addtoolbox]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::AddToolbox(EvaluatorInterface           eval,
	const std::vector<Currency>& inputs,
	std::vector<Currency>&       outputs)
{
	if (inputs.empty()) throw OML_Error(OML_ERR_NUMARGIN);

	const Currency& in = inputs[0];
	if (!in.IsString())
		throw OML_Error(OML_ERR_STRING, 1);

	std::string importDll(in.StringVal());

    eval.SuspendFuncListUpdate();
    // Any code updates should be added to AddLibraryImpl as this is used both
    // by the addtoolbox command and the LoadDefaultLibraries method which is
    // used for loading several libraries at once by the GUI
    try
    {
        AddLibraryImpl(eval, in.StringVal());
    }
    catch (const OML_Error& err)
    {
        // Unsuspend
        eval.UnsuspendFuncListUpdate();
        throw err;
    }

    // Unsuspend and add new functions registered in toolbox
    eval.UnsuspendFuncListUpdate();
    eval.OnUpdateFuncList();

    return true;
}

bool BuiltInFuncsCore::RegisterLibrary(EvaluatorInterface           eval,
	const std::vector<Currency>& inputs,
	std::vector<Currency>& outputs)
{
	if (inputs.empty()) throw OML_Error(OML_ERR_NUMARGIN);

	const Currency& in = inputs[0];
	if (!in.IsString())
		throw OML_Error(OML_ERR_STRING, 1);

	std::vector<Currency> ins;
	std::vector<Currency> outs;
	ins.push_back(in);

	eval.SetDLLHierarchy(inputs[1].StringVal());
	return BuiltInFuncsCore::AddToolbox(eval, ins, outs);
}
//------------------------------------------------------------------------------
// Removes a library/toolbox [removelibrary
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::RemoveToolbox(EvaluatorInterface           eval,
	                                 const std::vector<Currency>& inputs,
	                                 std::vector<Currency>&       outputs)
{
	if (inputs.empty()) throw OML_Error(OML_ERR_NUMARGIN);

	const Currency& in = inputs[0];
	if (!in.IsString())
		throw OML_Error(OML_ERR_STRING, 1);

	std::string removeDll = in.StringVal();

	bool success = eval.RemoveLibrary(removeDll);

#ifndef OS_WIN
    if (!success)
    {
        removeDll = removeDll + ".so";
        success = eval.RemoveLibrary(removeDll);
    }

    if (!success)
    {
        removeDll = "lib" + removeDll;
        success = eval.RemoveLibrary(removeDll);
    }
#endif

	if (success)
	{
		void* handle = libs[removeDll];

		if (handle)
        {
            Finalize(handle);
			DyFreeLibrary(handle);
        }
		libs.erase(removeDll);
	}
	else
	{
		// do we want a warning here?
	}

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
// Returns true after displaying a list of keywords in the interpreter [varlist command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Keywordlist(EvaluatorInterface           eval,
	const std::vector<Currency>& inputs,
	std::vector<Currency>&       outputs)
{
	std::vector<std::string> varlist(eval.GetKeywords());

	std::string msg;
	for (std::vector<std::string>::const_iterator itr = varlist.begin();
		itr != varlist.end(); ++itr)
	{
		std::string var(*itr);
		if (!var.empty())
			msg += var + " ";
	}
	Currency cur(msg);
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
// Returns true after getting current file name being processed [omlfilename]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Omlfilename(EvaluatorInterface           eval, 
                                   const std::vector<Currency>& inputs, 
                                   std::vector<Currency>&       outputs)
{
    std::string option;
    bool        isfullpath = false;

    if (!inputs.empty())
    {
        if (!inputs[0].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
        }
        option = inputs[0].StringVal();
        std::transform(option.begin(), option.end(), option.begin(), ::tolower);
        isfullpath = (option == "fullpath");
        if (!isfullpath && option != "fullpathext")
        {
            throw OML_Error(
                "Error: invalid option in argument 1; must be 'fullpath' or 'fullpathext'");
        }
    }

    std::string currfile(eval.GetCurrentFilename());
    bool isdummyfile = (currfile == "dummy");

    bool        isbeingedited = false;
    if (currfile.empty() || isdummyfile)
    {
        SignalHandlerBase* handler = eval.GetSignalHandler();
        assert(handler);
        currfile = handler ? handler->GetEditorFileName() : "";
        isbeingedited = true;
    }

    if (currfile.empty() || (!isbeingedited && isdummyfile))
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
// Returns true after getting current file name being processed [omlfilename command]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Omllinenumber(EvaluatorInterface           eval,
	const std::vector<Currency>& inputs,
	std::vector<Currency>& outputs)
{
	outputs.push_back(eval.GetCurrentLinenumber());
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

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }

    std::string name (inputs[0].StringVal());
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
        if (!inputs[1].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
        }
        type = inputs[1].StringVal(); 
        if (!type.empty())
        {
            std::transform(type.begin(), type.end(), type.begin(), ::tolower);
        }
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
        // VSM-8200 try/catch to reduce unwanted error for exist check.
        std::string error, formaterr;
        try 
        {
            eval.FindFunctionByName(name, &fi, &fptr, NULL);
        }
        catch (OML_Error& e)
        {
            error = e.GetErrorMessage();
            std::size_t suppress_error_check = error.find("cannot have executable statements");
            // suppress thrown error when this error is NOT located
            if (suppress_error_check == std::string::npos)
                throw e;
            // If a rejected file is not a function file, check if name matches a standard function to adjust the return.
            if (eval.IsStdFunction(name))
                returncode = 5;
        }
        if(returncode == 0)
            returncode = fptr ? 5 : 0;    // Built-in function
        if (returncode != 0 || !notype)
        {
            outputs.push_back(returncode);  
            return true;
        }

        if (fi)
        {
            outputs.push_back(103); // Any other function
            return true;
        }
    }

    BuiltInFuncsUtils::StripTrailingSlash(name);
#ifdef OS_WIN
    if (name.size() == 2 && name[1] == ':')
        name += "/";            // Add trailing slash if this is a drive letter
#endif

    std::string absIfPath(BuiltInFuncsUtils::Normpath(name));
	absIfPath = BuiltInFuncsUtils::GetAbsolutePath(absIfPath);

    bool filefound = false;
    if (notype || type == "file" || type == "dir")
    {
        returncode = BuiltInFuncsUtils::IsDir(absIfPath) ? 7 : 0; // Dir in current directory, return immediately
        if (returncode != 0)
        {
            outputs.push_back(returncode);
            return true;
        }
        if (BuiltInFuncsUtils::FileExists(absIfPath))
        {
            filefound = true; // File in current directory, but check for directory 
        }
        // Search the OML search path for a matching directory or file name.
        // Return 7 if a directory is found in the current directory or anywhere in the path.
        // Return 2 if a file is found in the currect directory or anywhere in the path.
        std::vector<std::string> paths = eval.GetPaths();
        std::vector<std::string>::const_iterator iter = paths.cbegin();
        std::string FileSearchPath;
        while (iter != paths.cend())
        {
            FileSearchPath = *iter++ + "/" + name;
            if (BuiltInFuncsUtils::IsDir(FileSearchPath))
            {
                returncode = 7;  
                if (notype || type == "dir")
                {
                    // immediately push back the results when a directory is found.
                    outputs.push_back(returncode);
                    return true;
                }
            }
            // If a file is located save the result, which will be returned unless a directory is located.
            if(BuiltInFuncsUtils::FileExists(FileSearchPath))
            {
                filefound = true; // File in current directory, but check for directory 
            }
        }
        if (filefound)  // if a file was located in the search, 2 will be returned here.
        {
            returncode = 2;
            outputs.push_back(returncode);
            return true;
        }
    }


    outputs.push_back(returncode);
    return true;
}
//------------------------------------------------------------------------------
// Returns true and returns a description of type of the given input [type]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Type(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    // Go through inputs and check for quiet option, validity
    bool quiet   = false;
    int nargin   = static_cast<int>(inputs.size());
    int  numcols = nargin;
    for (int i = 0; i < nargin; ++i)
    {
        if (!inputs[i].IsString())
        {
            throw OML_Error(OML_ERR_STRING, i + 1);
        }

        std::string val(inputs[i].StringVal());
        if (val.empty())
        {
            throw OML_Error(OML_ERR_NONEMPTY_STR, i + 1);
        }
        if (val == "-q" || val == "-Q")
        {
            quiet = true;
            numcols--;
        }
    }

    int nargout = eval.GetNargoutValue();

    if (numcols <= 0)  // Only quiet options specified
    {
        if (nargout > 0)
        {
            outputs.push_back(EvaluatorInterface::allocateCellArray());
        }
        return true;
    }
     
    std::unique_ptr<HML_CELLARRAY> cell = nullptr;
    if (nargout > 0)
    {
        cell.reset(EvaluatorInterface::allocateCellArray(1, numcols));
    }

    std::string reservedtokens("[]{}();,%\"");
    const OutputFormat* fmt = eval.GetOutputFormat();
    BuiltInFuncsUtils utils;

    std::string cwd = utils.GetCurrentWorkingDir();
    for (int i = 0, cellindex = 0; i < nargin; ++i)
    {
        std::string str (inputs[i].StringVal());
        if (quiet && str == "-q" || str == "-Q")
        {
            continue;
        }

        std::string type;
        if (eval.IsOperator(str))
        {
            type = str + " is an operator";
        }
        else if (eval.IsKeyword(str))
        {
            type = str + " is a keyword";
        }
        else if (eval.IsStdFunction(str))
        {
            type = str + " is a built-in function";
        }
        else if (eval.IsUserFunction(str))
        {
            type = str + " is a user-defined function";
        }
        else if (reservedtokens.find(str) != std::string::npos)
        {
            type = str + " is a reserved token";
        }
        else
        {
            const Currency& val = eval.GetValue(str);
            if (!val.IsNothing())
            {
                if (!quiet)
                {
                    type = str + " is a variable\n";
                }
                Currency elementval = eval.GetValue(str);
                elementval.DispOutput();
                type += elementval.GetOutputString(fmt);
            }
            else
            {
                std::string fname(str);
                if (utils.IsDir(str))
                {
                    throw OML_Error("Error: invalid input in argument " +
                        std::to_string(static_cast<long long>(i + 1)) + "; [" + 
                        utils.Normpath(str) + "] is a directory and is not defined");
                }

                bool exists = utils.FileExists(fname);
                std::string ext = utils.GetFileExtension(fname);
                if (!exists && ext.empty())
                {
                    fname += ".oml";
                }
                // Make sure we send only the relative path
                std::string rpath = fname;
                if (rpath.find(cwd) != std::string::npos)
                {
                    rpath = utils.GetRelativePath(fname, cwd);

                }
                bool foundfunc = eval.FindFileInPath(rpath, fname);
                if (foundfunc)
                {
                    if (!quiet && ext.empty())
                    {
                        type = str + " is a function defined in [" +
                            utils.Normpath(fname) + "]";
                    }
                }
                else if (!foundfunc && !exists)
                {
                    throw OML_Error("Error: invalid input in argument " +
                        std::to_string(static_cast<long long>(i + 1)) +
                        "; '" + str + "' is not defined");
                }

                fname = BuiltInFuncsUtils::Normpath(fname);

#ifdef OS_WIN
                BuiltInFuncsUtils utils;
                std::wstring wtype;
                std::wstring wname = utils.StdString2WString(fname);
                FILE* f = _wfopen(wname.c_str(), L"r, ccs=UTF-8");
                if (!f)
                {
                    throw OML_Error(OML_ERR_FILE_CANNOTREAD, i + 1);
                }
                struct _stat fstat;
                _wstat(wname.c_str(), &fstat);
                long fsize = fstat.st_size;

                // Read entire file contents in to memory
                if (fsize > 0)
                {
                    wtype.resize(fsize);
                    size_t numread = fread(&(wtype.front()), sizeof(wchar_t), fsize, f);
                    wtype.resize(numread);
                    wtype.shrink_to_fit();
                }
                fclose(f);
                type = utils.WString2StdString(wtype);
#else
                std::ifstream ifs(fname.c_str());
                if (ifs.fail())
                {
                    throw OML_Error(OML_ERR_FILE_CANNOTREAD, i + 1);
                }
                // Read the entire file
                while (!ifs.eof())
                {
                    std::string thisline;
                    std::getline(ifs, thisline);

                    if (!type.empty())
                    {
                        type += "\n";
                    }
                    type += thisline;
                }
                ifs.close();

#endif
            }
        }

        if (cell)
        {
            (*cell)(0, cellindex++) = type;
        }
        else
        {
            Currency display(type);
            display.DispOutput();
            eval.PrintResult(display);
        }
    }

    if (cell)
    {
        outputs.push_back(cell.release());
    }

    return true;
}
//------------------------------------------------------------------------------
// Sets pagination mode (omlpaginate)
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::OmlPaginate(EvaluatorInterface           eval,
    const std::vector<Currency>& inputs,
    std::vector<Currency>& outputs)
{
    size_t numargin = inputs.empty() ? 0 : inputs.size();
    if (numargin != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    Currency cur = inputs[0];
    if (cur.IsLogical() || cur.IsInteger())
    {
        int tmp = static_cast<int>(cur.Scalar());
        if (tmp >= CurrencyDisplay::PAGINATE_OFF ||
            tmp <= CurrencyDisplay::PAGINATE_INTERACTIVE)
        {
            CurrencyDisplay::SetPaginate((CurrencyDisplay::PAGINATE)tmp);
            return true;
        }
        throw OML_Error(OML_ERR_BAD_STRING, 1);
    }
    else if (cur.IsString())
    {
        std::string tmp(cur.StringVal());
        if (!tmp.empty())
        {
            BuiltInFuncsUtils utils;
            utils.SetPaginate(tmp);
            return true;
        }
        throw OML_Error(OML_ERR_BAD_STRING, 1);
    }
    throw OML_Error(OML_ERR_STRING_NATURALNUM, 1);
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
//------------------------------------------------------------------------------
// Returns true and prints memory usage (memoryuse)
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::MemoryUse(EvaluatorInterface           eval,
                                 const std::vector<Currency>& inputs,
                                 std::vector<Currency>&       outputs)
{
    std::string out;

#ifdef OS_WIN
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
    {
        out += "\tProcessPageFaultCount: " +
            std::to_string(static_cast<long long>(pmc.PageFaultCount));
        out += "\n\tProcessPeakWorkingSetSize: " +
            std::to_string(static_cast<long long>(pmc.PeakWorkingSetSize));
        out += "\n\tProcessWorkingSetSize: " +
            std::to_string(static_cast<long long>(pmc.WorkingSetSize));
        out += "\n\tProcessQuotaPeakPagedPoolUsage: " +
            std::to_string(static_cast<long long>(pmc.QuotaPeakPagedPoolUsage));
        out += "\n\tProcessQuotaPagedPoolUsage: " +
            std::to_string(static_cast<long long>(pmc.QuotaPagedPoolUsage));
        out += "\n\tProcessQuotaPeakNonPagedPoolUsage: " +
            std::to_string(static_cast<long long>(pmc.QuotaPeakNonPagedPoolUsage));
        out += "\n\tProcessQuotaNonPagedPoolUsage: " +
            std::to_string(static_cast<long long>(pmc.QuotaNonPagedPoolUsage));
        out += "\n\tProcessPagefileUsage: " +
            std::to_string(static_cast<long long>(pmc.PagefileUsage));
        out += "\n\tProcessPeakPagefileUsage: " +
            std::to_string(static_cast<long long>(pmc.PeakPagefileUsage));
    }
#endif

    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// Clears screen (clc command)
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Clc(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    CurrencyDisplay::ClearLineCount();

#if OS_WIN
    system("CLS");
#else
    system("clear");
#endif

    return true;
}
//------------------------------------------------------------------------------
// Sets limit which causes printing without formatting [skipformat]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::SkipFormat(EvaluatorInterface           eval,
                                  const std::vector<Currency>& inputs,
                                  std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsInteger())
    {
        throw OML_Error(OML_ERR_INTEGER, 1);
    }
    CurrencyDisplay::SetSkipFormat(static_cast<int>(inputs[0].Scalar()));
    return true;
}
//------------------------------------------------------------------------------
// Displays the currency output [display]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Display(EvaluatorInterface           eval,
                               const std::vector<Currency>& inputs,
                               std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    eval.PrintResult(inputs[0]);
    return true;
}
//------------------------------------------------------------------------------
// Exits the application [exit/quit]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Exit(EvaluatorInterface           eval,
    const std::vector<Currency>& inputs,
    std::vector<Currency>& outputs)
{
    int returnCode = EXIT_SUCCESS;
    if (!inputs.empty())
    {
        if (!inputs[0].IsInteger())
        {
            throw OML_Error(OML_ERR_INTEGER, 1);
        }

        returnCode = static_cast<int>(inputs[0].Scalar());
    }
    eval.SetQuit(true);
    eval.OnSaveOnExit(returnCode);
    return true;
}
//------------------------------------------------------------------------------
// Helper method to load default libraries. This can be done at startup to save
// time.
//------------------------------------------------------------------------------
void BuiltInFuncsCore::LoadDefaultLibraries(EvaluatorInterface              eval,
                                            const std::vector<std::string>& tboxes,
                                            const std::vector<std::string>& scripts)
{
    eval.SuspendFuncListUpdate(); // Suspend function list update till done
    
    for (std::vector<std::string>::const_iterator itr = tboxes.begin();
        itr != tboxes.end(); ++itr)
    {
        std::string name(*itr);
        try
        {
            AddLibraryImpl(eval, name, false);
        }
        catch (const OML_Error&) // Load silently
        {
        }
    }
    std::vector<Currency> scriptname = { Currency() };
    for (std::vector<std::string>::const_iterator itr = scripts.begin();
        itr != scripts.end(); ++itr)
    {
        scriptname[0] = *itr;
        try
        {
            eval.CallFunction("run", scriptname);
        }
        catch (const OML_Error&) // Load silently
        {
        }
    }

    eval.UnsuspendFuncListUpdate(); // Update function list
    eval.OnUpdateFuncList();
}
//------------------------------------------------------------------------------
// Executes keyboard command to allow the user to evaluate strings [keyboard]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::Keyboard(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
{
    SignalHandlerBase* parentHandler = eval.GetSignalHandler();
    if (!parentHandler || parentHandler->IsInBatchMode())
    {
        return true;
    }
    std::string cmd;
    if (parentHandler->IsWaiting(cmd))
    {
        BuiltInFuncsUtils::SetWarning(eval,
            "Warning: cannot run [keyboard] as [" + cmd + "] is running");
        return true;
    }

    std::string prompt("debug> ");
    if (!inputs.empty())
    {
        if (!inputs[0].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 1);
        }
        prompt = inputs[0].StringVal();
    }

    std::string type;
    

    // Don't to interactive pagination
    CurrencyDisplay::PAGINATE cachedPaginate = CurrencyDisplay::GetPaginate();
    bool resetPagination = CurrencyDisplay::IsPaginateInteractive();
    if (resetPagination)
    {
        CurrencyDisplay::SetPaginate(CurrencyDisplay::PAGINATE_ON);
    }

    std::vector<Currency> kinput = { Currency() };
    while (!eval.IsInterrupt())
    {
        std::string userInput;
        eval.OnUserInput(prompt, type, userInput);

        std::string lower(userInput);
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower == "return" || lower == "dbcont")
        {
            break;  // End of keyboard
        }
        else if (lower == "keyboard" || lower.find("keyboard") == 0)
        {
            continue;
        }

        kinput[0] = userInput;
        try
        {
            Currency result = eval.CallFunction("eval", kinput);
            eval.PushResult(result);
        }
        catch (const OML_Error& e)
        {
            Currency cur(Currency::TYPE_ERROR, e.GetErrorMessage());
            eval.PushResult(cur);
        }
    }
    if (resetPagination)
    {
        CurrencyDisplay::SetPaginate(cachedPaginate);
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns true if given data is a scalar [isfloat]
//------------------------------------------------------------------------------
bool BuiltInFuncsCore::IsFloat(EvaluatorInterface           eval,
                               const std::vector<Currency>& inputs,
                               std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    const Currency& cur = inputs[0];
    if (cur.IsLogical() || cur.IsCellArray() || cur.IsNDCellArray() ||
        cur.IsCellList()  || cur.IsString()  || cur.IsStruct())
    {
        outputs.emplace_back(false);
    }
    else if (cur.IsScalar() || cur.IsComplex()  || cur.IsEmpty() ||
             cur.Matrix()   || cur.IsNDMatrix() || cur.IsSparse())
    {
        outputs.emplace_back(true);
    }
    else
    {
        outputs.emplace_back(false);
    }
    return true;
}
//------------------------------------------------------------------------------
// Helper method for adding a library - used addtoolbox command and adding libraries in GUI
//------------------------------------------------------------------------------
void BuiltInFuncsCore::AddLibraryImpl(EvaluatorInterface eval, 
                                      const std::string& dllpath,
                                      bool warn)
{
    std::string importDll(dllpath);
    void* vResult = BuiltInFuncsCore::DyLoadLibrary(importDll.c_str());

#ifndef OS_WIN
    if (!vResult)
    {
        importDll = importDll + ".so";
        vResult = BuiltInFuncsCore::DyLoadLibrary(importDll.c_str());
    }

    if (!vResult)
    {
        importDll = "lib" + importDll;
        vResult = BuiltInFuncsCore::DyLoadLibrary(importDll.c_str());
    }
#endif

    if (!vResult)
    {
        throw OML_Error(OML_ERR_INVALID_DLL);
    }

    libs[importDll] = vResult;

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
            // Setting a warning here so we know which dll has an issue
            if (warn)
            {
                BuiltInFuncsUtils::SetWarning(eval,
                    "Warning: Invalid version in " + importDll);
            }
            throw OML_Error(OML_ERR_INVALID_VERSION, 1);
        }
    }

    if (!symbol)
    {
        symbol = BuiltInFuncsCore::DyGetFunction(vResult, "InitToolbox");

        if (symbol)
        {
            initAltFP fp = (initAltFP)symbol;
            OMLInterfaceImpl* impl = new OMLInterfaceImpl(&eval);

            eval.RegisterDLL(vResult);
#ifdef OS_WIN
            eval.SetDLLContext(DyGetLibraryPath(vResult));
#else
            eval.SetDLLContext(DyGetLibraryPath(symbol));
#endif
            fp(impl);
            eval.SetDLLContext("");

            eval.BCIGarbageCollect();
            return;
        }
    }

    initModuleFP fp = (initModuleFP)symbol;
    if (!fp)
    {
        throw OML_Error(OML_ERR_INVALID_INITDLL);
    }
    eval.RegisterDLL(vResult);

#ifdef OS_WIN
    eval.SetDLLContext(DyGetLibraryPath(vResult));
#else
    eval.SetDLLContext(DyGetLibraryPath(symbol));
#endif
    fp(eval);
    eval.SetDLLContext("");
}