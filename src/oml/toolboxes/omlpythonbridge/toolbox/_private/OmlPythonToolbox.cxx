/**
* @file OmlPythonToolbox.cxx
* @date December, 2017
* Copyright (C) 2015-2020 Altair Engineering, Inc.  
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

#include "OmlPythonToolbox.h"
#include "OmlPythonBridge.h"
#include "EvaluatorInt.h" 
#include "OML_Error.h"
#define TBOXVERSION 2020
// End defines/includes

bool Oml_EvalPythonFile(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    bool results = false;
    size_t nargin = inputs.size();

    if (nargin != 1)
        throw OML_Error("Valid syntax is \"evalpythonfile filename\".");
	
    if (!(inputs[0].IsString()))
        throw OML_Error("File name must be string.");

    if(OmlPythonBridge::GetInstance(eval))
    {
        OmlPythonBridge::GetInstance(eval)->SetLastErrorMessage("");
        results = OmlPythonBridge::GetInstance(eval)->EvalPythonFile(inputs[0].StringVal());
        if (results)
		{
			outputs.push_back(1);
		}
		else
		{
			outputs.push_back(0);
		}
		outputs.push_back(OmlPythonBridge::GetInstance(eval)->GetLastErrorMessage());
    }
    else
    {
		outputs.push_back(0);
		outputs.push_back("Python Interpreter is not avialable.");
    }
    return results;
}

bool Oml_EvalPythonScript(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    bool results = false;	
    size_t nargin = inputs.size();

    if (nargin != 1)
        throw OML_Error("Valid syntax is \"evalpythonfile script\".");

    if (!(inputs[0].IsString()))
        throw OML_Error("Script must be string.");

    if(OmlPythonBridge::GetInstance(eval))
    {
        OmlPythonBridge::GetInstance(eval)->SetLastErrorMessage("");
        results = OmlPythonBridge::GetInstance(eval)->EvalPythonScript(inputs[0].StringVal());
        if (results)
		{
			outputs.push_back(1);
		}
		else
		{
			outputs.push_back(0);
		}
		outputs.push_back(OmlPythonBridge::GetInstance(eval)->GetLastErrorMessage());
    }
    else
    {
		outputs.push_back(0);
		outputs.push_back("Python Interpreter is not avialable.");
    }
    return results;
}

bool Oml_GetPythonVariable(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    bool results = false;
    size_t nargin = inputs.size();

    if (nargin != 1)
        throw OML_Error("Valid syntax is \"getpythonvar python_var\".");
	
    if (inputs[0].IsEmpty())
        throw OML_Error("Python variable name should not be empty.");

    if (!(inputs[0].IsString()))
        throw OML_Error("Python variable name must be string.");
    
    if(OmlPythonBridge::GetInstance(eval))
    {
        OmlPythonBridge::GetInstance(eval)->SetLastErrorMessage("");
        results = OmlPythonBridge::GetInstance(eval)->GetPythonVariable(inputs[0].StringVal(),outputs);
        if (results)
		{
			outputs.push_back(1);
		}
		else
		{
			outputs.push_back("");
			outputs.push_back(0);
		}
		outputs.push_back(OmlPythonBridge::GetInstance(eval)->GetLastErrorMessage());
    }
    else
    {
        outputs.push_back("");
		outputs.push_back(0);
		outputs.push_back("Python Interpreter is not avialable.");
    }
    return results;
}

bool Oml_SetPythonVariable(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    bool results = false;
    size_t nargin = inputs.size();

    if (nargin != 2)
        throw OML_Error("Valid syntax is \"exporttopython oml_var python_var\".");


    if (!(inputs[0].IsString() || inputs[0].IsScalar() || inputs[0].IsComplex() || inputs[0].IsMatrix() 
        || inputs[0].IsNDMatrix() || inputs[0].IsStruct() || inputs[0].IsCellArray() || inputs[0].IsSparse()))
        throw OML_Error("Supported OML types for export are numbers, strings, complex, struct, Cell, matrix, ndmatrix and sparse matrix.");


    if (inputs[1].IsEmpty())
        throw OML_Error("Python variable name should not be empty.");
	
    if (!(inputs[1].IsString()))
        throw OML_Error("Python variable name must be string.");
    
    if(OmlPythonBridge::GetInstance(eval))
    {
        OmlPythonBridge::GetInstance(eval)->SetLastErrorMessage("");
        results = OmlPythonBridge::GetInstance(eval)->SetPythonVariable(inputs[1].StringVal(),inputs[0]);
        if (results)
		{
			outputs.push_back(1);
		}
		else
		{
			outputs.push_back(0);
		}
		outputs.push_back(OmlPythonBridge::GetInstance(eval)->GetLastErrorMessage());
    }
    else
    {
		outputs.push_back(0);
		outputs.push_back("Python Interpreter is not avialable.");
    }

    return results;
}


//! Register Python interface methods
extern "C" OMLPYTHONBRIDGE_DECLS int InitDll(EvaluatorInterface evl)
{
    evl.RegisterBuiltInFunction("evalpythonfile", Oml_EvalPythonFile, FunctionMetaData(1, 1, "PythonBridgeOmlCommands"));
    evl.RegisterBuiltInFunction("evalpythonscript", Oml_EvalPythonScript, FunctionMetaData(1, 1, "PythonBridgeOmlCommands"));
    evl.RegisterBuiltInFunction("getpythonvar", Oml_GetPythonVariable, FunctionMetaData(1, 1, "PythonBridgeOmlCommands"));
	evl.RegisterBuiltInFunction("exporttopython", Oml_SetPythonVariable, FunctionMetaData(2, 2, "PythonBridgeOmlCommands"));

    return 0;
}

//------------------------------------------------------------------------------
// Returns toolbox version
//------------------------------------------------------------------------------
extern "C" OMLPYTHONBRIDGE_DECLS double GetToolboxVersion(EvaluatorInterface eval)
{
    return TBOXVERSION;
}
