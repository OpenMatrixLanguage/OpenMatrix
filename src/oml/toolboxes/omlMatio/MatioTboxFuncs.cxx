/**
* @file MatioTboxFuncs.cxx
* @date November 2015
* Copyright (C) 2015-2020 Altair Engineering, Inc.  
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

#include "MatioTboxFuncs.h"

#include <cassert>
#include <fstream>
#include <memory>

#include "omlmatio.h"

#include "BuiltInFuncsUtils.h"
#include "OML_Error.h"
#include "StructData.h" 

#include "matio.h"

#define TBOXVERSION 2021

// Ascii files
// Returns true after loading file in ascii format
bool LoadTxtFile(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
// Returns true after saving file in ascii format
bool SaveAsciiFile(EvaluatorInterface              eval, 
                   const std::string&              filename,
                   const std::vector<std::string>& vars);

// Debugging
void PrintMatioFileVersion(mat_t*);    // Prints version to stdout for debugging
void SetMatioMessage(int, char*);      // Sets error from the matio library

//# define OMLMATIO_DBG 1  // Uncomment to print debug info
#ifdef OMLMATIO_DBG


#    define OMLMATIO_PRINT(str, m) { std::cout << str << m << std::endl; }
#    define OMLMATIO_SHOWERRORS { Mat_LogInitFunc("OML", SetMatioMessage); }// Crash in release with wide strings
#else
#    define OMLMATIO_PRINT(str, m) 0
#    define OMLMATIO_SHOWERRORS  0
#endif
//------------------------------------------------------------------------------
// Helper class for loading variables
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Entry point which registers load/save functions with oml
//------------------------------------------------------------------------------
int InitDll(EvaluatorInterface eval)
{
    eval.RegisterBuiltInFunction("load", &OmlLoad, FunctionMetaData(-1 ,1, 
                                 "FileIO"));
    eval.RegisterBuiltInFunction("save", &OmlSave, FunctionMetaData(-1, 0, 
                                 "FileIO"));
    return 1;
}
//------------------------------------------------------------------------------
// Returns true after loading the given file using MATIO library [load command]
//------------------------------------------------------------------------------
bool OmlLoad(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    int nargin = static_cast<int>(inputs.size());

	std::vector<std::string> target_variables;
    target_variables.reserve(nargin);

	if (nargin > 1)
	{
        if (!inputs[nargin - 1].IsString())
        {
            throw OML_Error(OML_ERR_STRING, nargin);
        }

    	std::string lastinput (inputs[nargin-1].StringVal());
        if (!lastinput.empty())
        {
            std::string lower(lastinput);
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            if (lower == "-ascii")
            {
                return LoadTxtFile(eval, inputs, outputs);
            }
        }

		for (int j = 1; j < nargin; ++j)
		{
            if (!inputs[j].IsString())
            {
                throw OML_Error(OML_ERR_STRING, j + 1);
            }
    		target_variables.push_back(inputs[j].StringVal());
		}
	}

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

	int nargout = eval.GetNargoutValue();
    std::unique_ptr<StructData> out_sd = nullptr;

    if (nargout)
    {
        out_sd.reset(EvaluatorInterface::allocateStruct());
    }
    std::string filename(inputs[0].StringVal());
    BuiltInFuncsUtils utils;
#ifdef OS_WIN
    if (utils.HasWideChars(filename))
    {
        throw OML_Error(OML_ERR_UNICODE_FILENAME, 1);
    }
#endif

    OMLMATIO_SHOWERRORS; // Causes a crash in release if there are wide strings

	mat_t* m = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
    if (!m)
    {
        try  // Try reading as an txt before quitting
        {
            return LoadTxtFile(eval, inputs, outputs);
        }
        catch (const OML_Error& e)
        {
            throw e;
        }
        catch (...)
        {
            throw OML_Error(OML_ERR_FILE_CANNOTREAD, 1);
        }
	}

    PrintMatioFileVersion(m);  // Prints file version to stdout for debugging

    size_t varIdx  = 0;
    size_t numvars = (target_variables.empty()) ? 0 : target_variables.size();
    bool   hasvars = (!target_variables.empty());

    OmlMatio omlMatio(filename, eval.GetVerbose(), OmlMatio::MATFILEVERSION_73);

    bool loadedVars = false;
    while (1)
    {
        matvar_t* var = nullptr;
        if (!hasvars)
        {
            var = Mat_VarReadNext(m);
        }
        else if (varIdx < numvars)
        {
            var = Mat_VarRead(m, target_variables[varIdx].c_str());
            varIdx++;
        }
           
        if (!var)
        {
            if (!hasvars || varIdx >= numvars)
            {
                break;   // End of file and no variables specified
            }

            continue;
        }

        try
        {
            Currency cur = omlMatio.MatVarToCurrency(var, eval);
            if (!cur.IsNothing())
            {
                std::string name(omlMatio.GetName(var));
                if (!out_sd)
                {
                    eval.SetValue(name, cur);
                }
                else
                {
                    out_sd->addField(name);
                    out_sd->SetValue(0, 0, name, cur);
                }
                loadedVars = true;
            }
        }
        catch (const OML_Error& e)
        {
            Mat_VarFree(var);
            var = nullptr;

            Mat_Close(m);
            m = nullptr;

            std::string warn(omlMatio.GetWarning());
            if (!warn.empty())
            {
                BuiltInFuncsUtils::SetWarning(eval, warn);
            }

            throw e;
        }
        Mat_VarFree(var);
        var = nullptr;
    }
    
	Mat_Close(m);
    m = nullptr;

    if (out_sd)
    {
        outputs.push_back(out_sd.release());
    }

    if (!loadedVars)
    {
        omlMatio.AddWarning("No variables loaded from [" +
            BuiltInFuncsUtils::Normpath(filename) + "]");
    }

    BuiltInFuncsUtils::SetWarning(eval, omlMatio.GetWarning());
	return true;
}
//------------------------------------------------------------------------------
// Returns true after saving the given file using MATIO library [save]
//------------------------------------------------------------------------------
bool OmlSave(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    int nargin = static_cast<int>(inputs.size());

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }

    std::string filename (inputs[0].StringVal());
    BuiltInFuncsUtils utils;

#ifdef OS_WIN
    if (utils.HasWideChars(filename))
    {
        throw OML_Error(OML_ERR_UNICODE_FILENAME, 1);
    }
#endif
    bool isascii = false;  // Saves as binary files by default

    mat_ft            version     = MAT_FT_MAT73;         // Default version
    matio_compression compression = MAT_COMPRESSION_ZLIB; // Default compression
    int               cmpidx      = -1;

    std::vector<std::string> target_variables;
    target_variables.reserve(nargin);
	
    for (int i = 1; i < nargin; ++i)
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

        if (val[0] != '-')  // This is a variable
        {
            target_variables.push_back(val);
            continue;
        }

        // Processing matio options
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);
        if (val == "-v7" || val == "-7" || val == "-v7.3" || val == "-7.3")
        {
            version = MAT_FT_MAT73;
        }
        else if (val == "-v5" || val == "-5")
        {
            version = MAT_FT_MAT5;
        }
        else if (val == "-v4" || val == "-4")
        {
            std::string msg = "Error: unsupported format. Valid formats are ";
            msg += "-v5 and -v7 in argument " + 
                    std::to_string(static_cast<long long>(i + 1));
            throw OML_Error(msg);

        }
        else if (val == "-nozip")
        {
            compression = MAT_COMPRESSION_NONE;
            cmpidx      = i + 1;
        }
        else if (val == "-ascii")
        {
            isascii = true;
        }
        else
        {
            throw OML_Error(OML_ERR_OPTIONVAL, i + 1, OML_VAR_VALUE);
        }
    }

    if (isascii)
    {
        return SaveAsciiFile(eval, filename, target_variables);
    }
		
    OmlMatio::MATFILEVERSION omlMatVer = (version == MAT_FT_MAT5) ?
        OmlMatio::MATFILEVERSION_5 : OmlMatio::MATFILEVERSION_73;
    OmlMatio omlMatio(filename, eval.GetVerbose(), omlMatVer);

    if (version != MAT_FT_MAT5 && compression != MAT_COMPRESSION_ZLIB)
    {
        std::string msg = "Invalid compression specified for version;";
        msg += " ignoring option in argument " + 
               std::to_string(static_cast<long long>(cmpidx));
        omlMatio.AddWarning(msg);
        compression = MAT_COMPRESSION_ZLIB;
    }

    OMLMATIO_SHOWERRORS; // Causes a crash in release if there are wide strings

    // Creates the file to save
	mat_t* m = Mat_CreateVer(filename.c_str(), nullptr, version);
	if (!m)
    {
        BuiltInFuncsUtils::SetWarning(eval, omlMatio.GetWarning());
        throw OML_Error("Error: invalid value in argument 1; cannot create file ["
            + BuiltInFuncsUtils::Normpath(filename) + "]");
    }


    bool hastargets   = (!target_variables.empty());
    bool handle_check = eval.HasBuiltin("ishandle");
    std::vector<std::string> varnames (eval.GetVariableNames());

    int savedVars = false;
	for (std::vector<std::string>::const_iterator iter = varnames.begin(); 
         iter != varnames.end(); ++iter)
	{
        std::string name (*iter);
		if (hastargets)
		{
			std::vector<std::string>::iterator iter2 = std::find(
                target_variables.begin(), target_variables.end(), name);
			if (iter2 == target_variables.end())
            {
				continue;
            }
		}

		Currency cur = eval.GetValue(name);
        if (handle_check)
        {
            std::vector<Currency> temp;
            temp.push_back(cur);

            Currency is_handle = eval.CallFunction("ishandle", temp);

            if (cur.IsScalar())
            {
                if (cur.IsPositiveInteger() || cur.Scalar() == 0.0)
                {
                    is_handle = Currency(false);
                }
            }
            if (is_handle.IsLogical() && is_handle.Scalar() == 1.0)
            {
                continue;
            }
        }

        matvar_t* var = nullptr;
        try
        {
            var = omlMatio.CurrencyToMatVar(cur.GetOutputName().c_str(), cur);
            if (var)
            {
                savedVars = true;
                Mat_VarWrite(m, var, compression);
                Mat_VarFree(var);
                var = nullptr;
            }
        }
        catch (const OML_Error& e)
        {
            if (var)
            {
                Mat_VarFree(var);
                var = nullptr;

            }
            Mat_Close(m);
            m = nullptr;

            BuiltInFuncsUtils::SetWarning(eval, omlMatio.GetWarning());
            throw(e);
        }
    }

	Mat_Close(m);
    m = nullptr;

    if (hastargets)
    {
        bool hasvars = (!(varnames.empty()));
        for (std::vector<std::string>::const_iterator itr = target_variables.begin();
            itr != target_variables.end(); ++itr)
        {
            if (!hasvars || 
                std::find(varnames.begin(), varnames.end(), *itr) == varnames.end())
            {
                omlMatio.AddWarning("Ignoring [" + *itr + "]; variable not found");
            }
        }
    }
    if (!savedVars)
    {
        omlMatio.AddWarning("No variables saved in [" + 
            BuiltInFuncsUtils::Normpath(filename) + "]");
    }
    BuiltInFuncsUtils::SetWarning(eval, omlMatio.GetWarning());

	return true;
}
//-----------------------------------------------------------------------------
// Returns true after reading an ascii file
//-----------------------------------------------------------------------------
bool LoadTxtFile(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs)
{
    assert(!inputs.empty());
	
    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }
		
    std::string file (inputs[0].StringVal());	

	std::ifstream ifs;
	ifs.open(file, std::ifstream::in);

	if (!ifs.good())
    {
        throw OML_Error(OML_ERR_FILE_CANNOTREAD, 1);
    }
	std::vector<std::vector<double>> elements;

	while (!ifs.eof())
	{
        std::string s;
        std::getline(ifs, s);

        if (s.empty() || s[0] == '%')
        {
            continue;
        }

#ifndef OS_WIN
        if (s[s.size() - 1] == '\r')
        {
            s.erase(s.size() - 1);
            if (s.empty())
                continue;
        }
#endif
        std::vector<double> row_elements;
        double val = 0.0;

		std::size_t prev = 0;
		std::size_t pos  = 0;

		while ((pos = s.find_first_of(" ,%\t", prev)) != std::string::npos)
		{
			if (s[pos] == '%')
			{
				prev = s.length();
				break;
			}

			if (pos > prev)
			{
				std::string substr = s.substr(prev, pos-prev);

				try
				{
					val = std::stod(substr);
				}
				catch (...)
				{
                    throw OML_Error(OML_ERR_FILE_CANNOTREAD, 1);
				}

				row_elements.push_back(val);
			}
			prev = pos+1;
		}

		if (prev < s.length())
		{
			std::string substr = s.substr(prev, std::string::npos);

			try
			{
				val = std::stod(substr);
			}
			catch (...)
			{
                throw OML_Error(OML_ERR_FILE_CANNOTREAD, 1);
			}

			row_elements.push_back(val);
		}

		if (row_elements.size())
			elements.push_back(row_elements);
	}

	int num_rows = static_cast<int>(elements.size());
	int num_cols = (num_rows) ? static_cast<int>(elements[0].size()) : 0;

    std::unique_ptr<hwMatrix> mtx(EvaluatorInterface::allocateMatrix(
        num_rows, num_cols, hwMatrix::REAL));

	for (int j=0; j<num_rows; j++)
	{
		if (elements[j].size() != num_cols)
			throw OML_Error(OML_ERR_PLOT_DIM_NOT_MATCH);

		for (int k=0; k<num_cols; k++)
			(*mtx)(j,k) = (elements[j])[k];
	}

    BuiltInFuncsUtils utils;

    Currency result(mtx.release());
    std::string name(utils.GetBaseName(file));
    if (!name.empty())
    {
        size_t pos = name.find_last_of(".");
        if (pos != std::string::npos)
        {
            name = name.substr(0, pos);
        }
    }

    if (!name.empty())
    {
        std::replace(name.begin(), name.end(), ' ', '_');
    }

    eval.SetValue(name, result);
	outputs.push_back(result);
	return true;
}
//------------------------------------------------------------------------------
// Returns true after saving ascii file
//------------------------------------------------------------------------------
bool SaveAsciiFile(EvaluatorInterface              eval, 
                   const std::string&              filename,
                   const std::vector<std::string>& vars)
{
    bool hastargets = true;
    std::vector<std::string> varnames  (vars);
    std::vector<std::string> scopevars (eval.GetVariableNames());
    if (varnames.empty())
    {
        varnames   = scopevars;
        hastargets = false;
    }

    bool hashandle  = eval.HasBuiltin("ishandle");

    std::FILE* fp = fopen(filename.c_str(), "w");
    if (!fp)
    {
        throw OML_Error("Error: Cannot open file: " + 
            BuiltInFuncsUtils::Normpath(filename));
    }

    for (std::vector<std::string>::const_iterator itr = varnames.begin(); 
         itr != varnames.end(); ++itr)
    {
        std::string name (*itr);
        if (name.empty())
        {
            continue;
        }

        Currency cur = eval.GetValue(name);

        if (hashandle)
        {
			std::vector<Currency> temp;
			temp.push_back(cur);

			Currency ishandle = eval.CallFunction("ishandle", temp);
			if (cur.IsScalar())
			{
				if (cur.IsPositiveInteger() || cur.Scalar() == 0.0)
                {
					ishandle = Currency(false);
                }
			}
			if (ishandle.IsLogical() && ishandle.Scalar() == 1.0)
            {
				continue;
            }
        }

        if (cur.IsScalar())
        {
            fprintf(fp, "%g\n", cur.Scalar());
        }
        else if (cur.IsMatrixOrString())
        {
            const hwMatrix* mtx = cur.Matrix();
            if (!mtx) 
            {
                continue;
            }
            if (!mtx->IsRealData())
            {
                std::string msg = " cannot save variable [" + name + 
                                  "] as matrix value is not real";
                if (hastargets)
                {
                    fclose(fp);
                    throw OML_Error("Error:" + msg);
                }
                else
                {
                    BuiltInFuncsUtils::SetWarning(eval,"Warning:" + msg);
                    continue;
                }
            }
            int m = mtx->M();
            int n = mtx->N();
			for (int j = 0; j < m; ++j)
			{
				for (int k = 0; k < n; ++k)
                {
					fprintf(fp, "%g ", (*mtx)(j,k));
                }
				fprintf(fp, "\n");
			}
        }
        else
        {
            std::string msg = " cannot save variable [" + name + 
                                "] as it is not a scalar, real matrix or string";
            if (hastargets)
            {
                fclose(fp);
                if (hastargets)
                {
                    if (scopevars.empty() || 
                        std::find(scopevars.begin(), scopevars.end(), name) == scopevars.end())
                    {
                        throw OML_Error("Error: cannot save variable [" + name +
                            "] as it does not exist in the current session");
                    }
                }
                throw OML_Error("Error:" + msg);
            }
            else
            {
                BuiltInFuncsUtils::SetWarning(eval, "Warning:" + msg);
                continue;
            }
        }
    }

    fclose(fp);
	return true;
}
//------------------------------------------------------------------------------
// Returns toolbox version
//------------------------------------------------------------------------------
double GetToolboxVersion(EvaluatorInterface eval)
{
    return TBOXVERSION;
}
//------------------------------------------------------------------------------
// Prints mat file version to stdout for debugging
//------------------------------------------------------------------------------
void PrintMatioFileVersion(mat_t* m)
{
#ifdef OMLMATIO_DBG
    assert(m);
    if (!m)
    {
        return;
    }
    mat_ft ver = Mat_GetVersion(m);
    switch (ver)
    {
    case MAT_FT_MAT73: OMLMATIO_PRINT("Matio file v", 7.3); break;
    case MAT_FT_MAT5:  OMLMATIO_PRINT("Matio file v", 5);   break;
    case MAT_FT_MAT4:  OMLMATIO_PRINT("Matio file v", 4);   break;
    default:           OMLMATIO_PRINT("Matio file version undefined ", ver); break;
    }
#endif
}
//------------------------------------------------------------------------------
// Prints messages from matio library
//------------------------------------------------------------------------------
void SetMatioMessage(int level, char* msg)
{
#ifdef OMLMATIO_DBG
    if (!msg)
    {
        return;
    }

    std::string strmsg(msg);
    if (level & MATIO_LOG_LEVEL_CRITICAL)
    {
        throw OML_Error("Matio critical error; " + strmsg);
    }
    else if (level & MATIO_LOG_LEVEL_ERROR)
    {
        throw OML_Error("Matio error; " + strmsg);
    }
    else if (level & MATIO_LOG_LEVEL_WARNING)
    {
        OMLMATIO_PRINT("Matio warning; " + strmsg, 0);
    }
    else if (level & MATIO_LOG_LEVEL_DEBUG)
    {
        OMLMATIO_PRINT("Matio debug; " + strmsg, 0);
    }
    else if (level & MATIO_LOG_LEVEL_MESSAGE)
    {
        OMLMATIO_PRINT("Matio debug; " + strmsg, 0);
    }
#endif
}
