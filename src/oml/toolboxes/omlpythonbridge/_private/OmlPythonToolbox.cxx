/////////////////////////////////////////////////////////////////////////
 // File : OmlPythonToolbox.cxx
 // Copyright (c) 2017 solidThinking Inc.  All Rights Reserved
 // Contains trade secrets of solidThinking, Inc.  Copyright notice
 // does not imply publication.  Decompilation or disassembly of this
 // software is strictly prohibited.
 ////////////////////////////////////////////////////////////////////////

// Begin defines/includes

#include "OmlPythonToolbox.h"
#include "OmlPythonBridge.h"
#include "EvaluatorInt.h" 
#include "OML_Error.h"

// End defines/includes

bool Oml_EvalPythonFile(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    bool results = false;
    size_t nargin = inputs.size();

    if (nargin != 1)
        throw OML_Error("Valid syntax is \"evalpythonfile filename\".");
	
    if (!(inputs[0].IsString()))
        throw OML_Error("File name must be string.");

    if(OmlPythonBridge::GetInstance())
    {
        OmlPythonBridge::GetInstance()->ClearLastErrorMessage();
        results = OmlPythonBridge::GetInstance()->EvalPythonFile(inputs[0].StringVal());
        if (results)
		{
			outputs.push_back(1);
		}
		else
		{
			outputs.push_back(0);
		}
		outputs.push_back(OmlPythonBridge::GetInstance()->GetLastErrorMessage());
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

    if(OmlPythonBridge::GetInstance())
    {
        OmlPythonBridge::GetInstance()->ClearLastErrorMessage();
        results = OmlPythonBridge::GetInstance()->EvalPythonScript(inputs[0].StringVal());
        if (results)
		{
			outputs.push_back(1);
		}
		else
		{
			outputs.push_back(0);
		}
		outputs.push_back(OmlPythonBridge::GetInstance()->GetLastErrorMessage());
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
    
    if(OmlPythonBridge::GetInstance())
    {
        OmlPythonBridge::GetInstance()->ClearLastErrorMessage();
        results = OmlPythonBridge::GetInstance()->GetPythonVariable(inputs[0].StringVal(),outputs);
        if (results)
		{
			outputs.push_back(1);
		}
		else
		{
			outputs.push_back("");
			outputs.push_back(0);
		}
		outputs.push_back(OmlPythonBridge::GetInstance()->GetLastErrorMessage());
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
        || inputs[0].IsNDMatrix() || inputs[0].IsStruct() || inputs[0].IsCellArray()))
        throw OML_Error("Supported OML types for export are numbers, strings, complex, struct, Cell and matrix, ndmatrix.");


    if (inputs[1].IsEmpty())
        throw OML_Error("Python variable name should not be empty.");
	
    if (!(inputs[1].IsString()))
        throw OML_Error("Python variable name must be string.");
    
    if(OmlPythonBridge::GetInstance())
    {
        OmlPythonBridge::GetInstance()->ClearLastErrorMessage();
        results = OmlPythonBridge::GetInstance()->SetPythonVariable(inputs[1].StringVal(),inputs[0]);
        if (results)
		{
			outputs.push_back(1);
		}
		else
		{
			outputs.push_back(0);
		}
		outputs.push_back(OmlPythonBridge::GetInstance()->GetLastErrorMessage());
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
