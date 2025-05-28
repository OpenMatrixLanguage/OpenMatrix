/**
* @file MatioTboxFuncs.cxx
* @date November 2015
* Copyright (C) 2015-2024 Altair Engineering, Inc.  
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
#include "MatrixDisplay.h"
#include "OML_Error.h"
#include "StructData.h" 

#include "matio.h"

#define TBOXVERSION 2024.1

// Ascii files
// Returns true after loading file in ascii format
bool LoadTxtFile(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
// Returns true after saving file in ascii format
bool SaveAsciiFile(EvaluatorInterface, const std::string&, const std::vector<std::string>&, bool);
// Gets the string version of mat file format
std::string MatioFileVersion(mat_t*);
// Imports variables in ascii file
Currency ImportVar(const std::string&, const std::string&, const std::vector <std::vector<double> >&);
//# define OMLMATIO_DBG 1  // Uncomment to print debug info
#ifdef OMLMATIO_DBG
#    define OMLMATIO_PRINT(str, m) { std::cout << str << m << std::endl; }
#    define OMLMATIO_SHOWERRORS { Mat_LogInitFunc("OML", SetMatioMessage); }// Crash in release with wide strings
#else
#    define OMLMATIO_PRINT(str, m) 0
#    define OMLMATIO_SHOWERRORS  0
#endif
void PrintMatioFileVersion(mat_t*);    // Prints version to stdout for debugging
void SetMatioMessage(int, char*);      // Sets error from the matio library
bool SaveStructFields(EvaluatorInterface, OmlMatio*, mat_t*, const Currency&, bool, matio_compression); // Saves struct fields as variables

//------------------------------------------------------------------------------
// Entry point which registers load/save functions with oml
//------------------------------------------------------------------------------
int InitDll(EvaluatorInterface eval)
{
    eval.RegisterBuiltInFunction("load", &OmlLoad, FunctionMetaData(-1 ,1, "FileIO"));
    eval.RegisterBuiltInFunction("save", &OmlSave, FunctionMetaData(-1, 0, "FileIO"));
    return 1;
}
//------------------------------------------------------------------------------
// Returns true after loading the given file using MATIO library [load]
//------------------------------------------------------------------------------
bool OmlLoad(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    else if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

    std::string filename(inputs[0].StringVal());
    if (filename.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 1);
    }

    BuiltInFuncsUtils utils;
    filename = utils.Normpath(filename);
    std::string ext(utils.GetFileExtension(filename));
    if (ext.empty())
    {
        filename += ".mat";
    }
    if (!utils.FileExists(filename))
    {
        if (!eval.FindFileInPath(filename, filename))
        {
            throw OML_Error(OML_ERR_FILE_NOTFOUND, filename, 1);
        }
    }

    int nargin = static_cast<int>(inputs.size());

	std::vector<std::string> target_variables;
    target_variables.reserve(nargin);

	for (int j = 1; j < nargin; ++j)
	{
        if (!inputs[j].IsString())
        {
            throw OML_Error(OML_ERR_STRING, j + 1);
        }
        std::string val (inputs[j].StringVal());
        if (val.empty())
        {
            continue;
        }

        if (val [0] != '-')  // This is a variable
        {
            target_variables.emplace_back(val);
            continue;
        }
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);
        if (val == "-ascii")
        {
            std::vector<Currency> intxt(inputs);
            if (!intxt.empty())
            {
                intxt [0] = filename; // Use the updated filename
            }
            return LoadTxtFile(eval, intxt, outputs);
        }
        else
        {
            throw OML_Error(OML_ERR_OPTIONVAL, val, j + 1);
        }
	}

	int nargout = eval.GetNargoutValue();
    std::unique_ptr<StructData> out_sd = nullptr;

    if (nargout)
    {
        out_sd.reset(EvaluatorInterface::allocateStruct());
    }

    OMLMATIO_SHOWERRORS; // Causes a crash in release if there are wide strings

	mat_t* m = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
    if (!m)
    {
        try  // Try reading as an txt before quitting
        {
            std::vector<Currency> intxt(inputs);
            if (!intxt.empty())
            {
                intxt[0] = filename; // Use the updated filename
            }
            return LoadTxtFile(eval, intxt, outputs);
        }
        catch (const OML_Error& e)
        {
            throw OML_Error(e.GetErrorMessage());
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

    OmlMatio omlMatio(filename, eval.GetVerbose(), OmlMatio::MATFILEVERSION_5);

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

            throw OML_Error(e.GetErrorMessage());
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
bool OmlSave(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    else if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }

    std::string filename (inputs[0].StringVal());
    BuiltInFuncsUtils::CheckSpecialCharsFileName(filename, 1);
    filename = BuiltInFuncsUtils::Normpath(filename);

    bool isascii = false;  // Saves as binary files by default
    bool append  = false;
    std::string strversion;

    mat_ft            version     = MAT_FT_MAT5;         // Default version
    matio_compression compression = MAT_COMPRESSION_ZLIB; // Default compression
    int               cmpidx      = -1;
    int               nargin      = static_cast<int>(inputs.size());

    std::vector<std::string> target_variables;
    target_variables.reserve(nargin);

    bool savingsStructVars = false;
    std::vector<std::string> structvars;
    structvars.reserve(nargin);
	
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
            if (!savingsStructVars)
            {
                target_variables.emplace_back(val);
            }
            else
            {
                structvars.emplace_back(val);
            }
            continue;
        }
        // Processing matio options
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);

        if (val == "-struct")
        {
            savingsStructVars = true;
            continue;
        }
        val.erase(val.begin());
        if (val.empty())
        {
            throw OML_Error(OML_ERR_OPTIONVAL, "-", i + 1);
        }


        if (val == "v7" || val == "7")
        {
            BuiltInFuncsUtils::SetWarning(eval,
                "Warning: unsupported format in argument " +
                std::to_string(i + 1) + "; saving file in v5 format");
            version    = MAT_FT_MAT5;
            strversion = val;
        }
        else if (val == "v7.3" || val == "7.3")
        {
            version    = MAT_FT_MAT73;
            strversion = val;
        }
        else if (val == "v5" || val == "5")
        {
            version    = MAT_FT_MAT5;
            strversion = val;
        }
        else if (val == "v4" || val == "4")
        {
            throw OML_Error("Error: unsupported format in argument " +
                std::to_string(i + 1) + ". Valid formats are -v5(default) and -v7.3");
        }
        else if (val == "nozip")
        {
            compression = MAT_COMPRESSION_NONE;
            cmpidx      = i + 1;
        }
        else if (val == "ascii")
        {
            isascii = true;
        }
        else if (val == "append")
        {
            append = true;
        }
        else
        {
            throw OML_Error(OML_ERR_OPTIONVAL, "-" + val, i + 1);
        }
    }

    if (isascii)
    {
        return SaveAsciiFile(eval, filename, target_variables, append);
    }
		
    OmlMatio::MATFILEVERSION omlMatVer = (version == MAT_FT_MAT73) ?
        OmlMatio::MATFILEVERSION_73 : OmlMatio::MATFILEVERSION_5;

    OmlMatio omlMatio(filename, eval.GetVerbose(), omlMatVer);

    if (version != MAT_FT_MAT5 && compression != MAT_COMPRESSION_ZLIB)
    {
        std::string msg = "Invalid compression specified for version;";
        msg += " ignoring option in argument " + std::to_string(cmpidx);
        omlMatio.AddWarning(msg);
        compression = MAT_COMPRESSION_ZLIB;
    }

    OMLMATIO_SHOWERRORS; // Causes a crash in release if there are wide strings

    // Creates the file to save
	mat_t* m = (!append) ? 
        Mat_CreateVer(filename.c_str(), nullptr, version) :
        Mat_Open(filename.c_str(), MAT_ACC_RDWR);
	if (!m)
    {
        BuiltInFuncsUtils::SetWarning(eval, omlMatio.GetWarning());
        throw OML_Error(OML_ERR_FILE_CANNOTWRITE, filename, 1);
    }

    if (append) // Check version and file existence
    {
        if (BuiltInFuncsUtils::FileExists(filename))
        {
            mat_ft fileversion = Mat_GetVersion(m);
            if (!(fileversion == MAT_FT_MAT73 || fileversion == MAT_FT_MAT5))
            {
                std::string err ("Error: unsupported file format");
                err += (fileversion == MAT_FT_MAT4) ? " [4]" : "";
                err += " in argument 1; valid formats are [5] and [7.3]";
                Mat_Close(m);
                m = nullptr;
                throw OML_Error(err);
            }
            else if (fileversion != version && !strversion.empty())
            {
                omlMatio.SetVersion((fileversion == MAT_FT_MAT73) ?
                    OmlMatio::MATFILEVERSION_73 : OmlMatio::MATFILEVERSION_5);
                omlMatio.AddWarning("Incompatible format [" + strversion +
                    "]; saving format [" + MatioFileVersion(m) + "]");
            }
        }
    }

    bool hasstructvars = false;
    if (savingsStructVars)
    {
        if (structvars.empty())
        {
            Mat_Close(m);
            m = nullptr;

            throw OML_Error(OML_ERR_OPTION,
                "must have at least one struct variable name when using '-struct'");
        }
        hasstructvars = true;
    }

    bool hastargets = (!target_variables.empty());
    bool saved  = false;

    std::vector<std::string> varnames (eval.GetVariableNames());
    for (std::vector<std::string>::const_iterator iter = varnames.begin(); 
         iter != varnames.end(); ++iter)
	{
        std::string name (*iter);
		if (hastargets && 
            std::find(target_variables.begin(), target_variables.end(), name) == target_variables.end())
		{
    		continue;
		}
        else if (hasstructvars &&
            std::find(structvars.begin(), structvars.end(), name) == structvars.end())
        {
            continue;
        }

		const Currency& cur = eval.GetValue(name);
        if (hasstructvars && !cur.IsStruct())
        {
            omlMatio.AddWarning("Ignoring [" + name + "]; must be a struct");
            continue;
        }
        else if (hasstructvars)
        {
            saved = SaveStructFields(eval, &omlMatio, m, cur, append, compression);
            continue;
        }
        int result = 0;
        matvar_t* var = (append) ? 
            Mat_VarRead(m, name.c_str()): nullptr;
        if (var)
        {
            try
            {
                Currency matcur = omlMatio.MatVarToCurrency(var, eval);

                std::vector<Currency> inputs2;
                inputs2.emplace_back(matcur);
                inputs2.emplace_back(cur);
                Currency isequal = eval.CallFunction("isequal", inputs2);
                if (isequal.IsLogical() && isequal.Scalar() == 1)
                {
                    // Nothing to update
                    saved = true;
                    continue;
                }

                result = Mat_VarDelete(m, name.c_str());
                Mat_VarFree(var);
                var = nullptr;
                if (result != 0)
                {
                    omlMatio.AddWarning("Ignoring [" + name + "]; " +
                        omlMatio.MatioError(result));
                    continue;
                }
            }
            catch (const OML_Error&)
            {
                omlMatio.AddWarning("Ignoring [" + name + "]; " +
                                    omlMatio.MatioError(result));
                Mat_VarFree(var);
                var = nullptr;
                Mat_Close(m);
                m = nullptr;

                BuiltInFuncsUtils::SetWarning(eval, omlMatio.GetWarning());
                throw;
            }
        }

        try
        {
            var = omlMatio.CurrencyToMatVar(cur.GetOutputName().c_str(), cur);
            if (var)
            {
                int result = Mat_VarWrite(m, var, compression);
                if (result != 0)
                {
                    omlMatio.AddWarning("Ignoring [" + name + "]; " +
                                        omlMatio.MatioError(result));
                }
                else
                {
                    saved = true;
                }
                Mat_VarFree(var);
            }
            var = nullptr;
        }
        catch (const OML_Error&)
        {
            if (var)
            {
                Mat_VarFree(var);
                var = nullptr;
            }
            
            Mat_Close(m);
            m = nullptr;

            BuiltInFuncsUtils::SetWarning(eval, omlMatio.GetWarning());
            throw;
        }
    }

	Mat_Close(m);
    m = nullptr;

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
    for (std::vector<std::string>::const_iterator itr = structvars.begin();
            itr != structvars.end(); ++itr)
    {
        if (!hasvars ||
            std::find(varnames.begin(), varnames.end(), *itr) == varnames.end())
        {
            omlMatio.AddWarning("Ignoring [" + *itr + "]; variable not found");
        }
    }

    if (!saved)
    {
        omlMatio.AddWarning("No variables saved in [" + filename + "]");
    }
    eval.RefreshPathCache();
    
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
		
    std::string file(inputs[0].StringVal());
    file = BuiltInFuncsUtils::Normpath(file);

    if (!BuiltInFuncsUtils::FileExists(file))
    {
        if (!eval.FindFileInPath(file, file))
        {
            throw OML_Error(OML_ERR_FILE_NOTFOUND, 1);
        }
    }

	std::ifstream ifs;
	ifs.open(file, std::ifstream::in);

	if (!ifs.good())
    {
        throw OML_Error(OML_ERR_FILE_CANNOTREAD, 1);
    }
	std::vector< std::vector<double> > elements;
    std::vector<double>                row_elements;

    std::string basename(BuiltInFuncsUtils::GetBaseName(file));
    if (!basename.empty())
    {
        size_t pos = basename.find_last_of(".");
        if (pos != std::string::npos)
        {
            std::string tmp1(basename);
            basename = tmp1.substr(0, pos);
        }
    }
    if (!basename.empty())
    {
        std::replace(basename.begin(), basename.end(), ' ', '_');
    }
    
    int idx = 0;
    std::string matname;
    std::string mattype;
	while (!ifs.eof())
	{
        std::string s;
        std::getline(ifs, s);

#ifndef OS_WIN
        if (!s.empty() && s [s.size() - 1] == '\r')
        {
            s.erase(s.size() - 1);
        }
#endif
        s = BuiltInFuncsUtils::LTrim(s, " ");
        s = BuiltInFuncsUtils::RTrim(s, " ");
        if (s.empty())
        {
            if (!elements.empty())
            {
                std::string name (matname);
                if (name.empty())
                {
                    name = basename;
                    if (idx > 0)
                    {
                        name += "-" + std::to_string(idx);
                    }
                    idx++;
                }
                Currency var = ImportVar(name, mattype, elements);
                eval.SetValue(name, var);
                outputs.emplace_back(var);
            }
            elements.clear();
            row_elements.clear();
            matname = "";
            mattype = "";
            continue;
        }
        else if (s[0] == '%')
        {
            continue;
        }
        else if (s.find("# name: ") == 0)
        {
            matname = s.substr(std::string("# name: ").length());
            continue;
        }
        else if (s.find("# type: ") == 0)
        {
            mattype = s.substr(std::string("# type: ").length());
            continue;
        }
        else if (s[0] == '#')
        {
            continue;
        }

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
                    throw OML_Error(OML_ERR_INVALIDFORMAT, 1);
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
                throw OML_Error(OML_ERR_INVALIDFORMAT, 1);
			}

			row_elements.push_back(val);
		}

        if (!row_elements.empty())
        {
            elements.push_back(row_elements);
            row_elements.clear();
        }
	}

    if (!elements.empty())
    {
        std::string name (matname);
        if (name.empty())
        {
            name = basename;
            if (idx > 0)
            {
                name += "-" + std::to_string(idx);
            }
            idx++;
        }
        Currency var = ImportVar(name, mattype, elements);        
        eval.SetValue(name, var);
        outputs.emplace_back(var);
    }
	return true;
}
//------------------------------------------------------------------------------
// Returns true after saving ascii file
//------------------------------------------------------------------------------
bool SaveAsciiFile(EvaluatorInterface              eval, 
                   const std::string&              filename,
                   const std::vector<std::string>& uservars,
                   bool                            append)
{
    std::string mode = (append) ? "a" : "w";
    std::FILE* fp = fopen(filename.c_str(), mode.c_str());
    if (!fp)
    {
        throw OML_Error(OML_ERR_FILE_CANNOTOPEN, filename, 1);
    }

    bool hastargets = !(uservars.empty());
    std::vector<std::string> omlvars (eval.GetVariableNames());
    if (hastargets && omlvars.empty())
    {
        fclose(fp);
        throw OML_Error("Error: cannot save [" + filename +
            "]; variables do not exist in the current session");
    }

    std::vector<std::string> varnames ((hastargets) ? uservars : omlvars);

    for (std::vector<std::string>::const_iterator itr = varnames.begin(); 
         itr != varnames.end(); ++itr)
    {
        std::string name (*itr);
        if (name.empty())
        {
            continue;
        }

        const Currency& cur = eval.GetValue(name);

        std::string header;
        if (append)
        {
            header += "\n# name: " + name + "\n# type: " + cur.GetTypeString() + "\n";
        }

        if (cur.IsScalar())
        {
            fprintf(fp, "%s%s\n", header.c_str(),
                CurrencyDisplay::NonFormattedDoubleToString(cur.Scalar()).c_str());
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
                if (hastargets)
                {
                    fclose(fp);
                    throw OML_Error(OML_ERR_REALMATRIX, name);
                }
                else
                {
                    BuiltInFuncsUtils::SetWarning(eval,"Warning: cannot save [" 
                        + name + "]; must be a scalar, real matrix or string");
                    continue;
                }
            }
            if (!header.empty())
            {
                fprintf(fp, "%s", header.c_str());
            }
            MatrixDisplay::WriteNonFormattedOutputValues(cur, "\n", " ", -1, fp);
        }
        else
        {
            std::string msg = " cannot save [" + name + 
                                "]; must be a scalar, real matrix or string";
            if (hastargets)
            {
                fclose(fp);
                if (std::find(omlvars.begin(), omlvars.end(), name) == omlvars.end())
                {
                    throw OML_Error("Error: cannot save [" + name +
                        "]; does not exist in the current session");
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
//------------------------------------------------------------------------------
// Gets mat file version
//------------------------------------------------------------------------------
std::string MatioFileVersion(mat_t* m)
{
    if (m)
    {
        mat_ft ver = Mat_GetVersion(m);
        switch (ver)
        {
            case MAT_FT_MAT73: return "v7.3";
            case MAT_FT_MAT5:  return "v5";
            case MAT_FT_MAT4:  return "v4";
            default:           return "";
        }
    }
    return "";
}
//------------------------------------------------------------------------------
// Imports a variable in an ascii file
//------------------------------------------------------------------------------
Currency ImportVar(const std::string& name,
                   const std::string&                  mattype,
                   const std::vector <std::vector<double> >& vals)
{
    int rows = static_cast<int>(vals.size());
    int cols = (rows > 0) ? static_cast<int>(vals [0].size()) : 0;

    Currency var;
    if (rows == 1 && cols == 1)
    {
        var = (vals [0]) [0];
    }
    else
    {
        std::unique_ptr<hwMatrix> mtx(EvaluatorInterface::allocateMatrix(
            rows, cols, true));

        for (int j = 0; j < rows; j++)
        {
            if (vals [j].size() != cols)
                throw OML_Error(OML_ERR_PLOT_DIM_NOT_MATCH);

            for (int k = 0; k < cols; k++)
                (*mtx)(j, k) = (vals [j]) [k];
        }

        var = mtx.release();
    }

    if (mattype == "char" || mattype == "string")
    {
        var.SetMask(Currency::MASK_STRING);
    }

    var.SetOutputName(name);
    return var;
}
//------------------------------------------------------------------------------
// Saves fields in the structs as individual variables
//------------------------------------------------------------------------------
bool SaveStructFields(EvaluatorInterface eval, 
                      OmlMatio*          omlMatio,
                      mat_t*             m, 
                      const Currency&    parent, 
                      bool               append,
                      matio_compression  compression)
{
    const StructData* sd = parent.Struct();
    assert(sd);
    if (!sd)
    {
        return false;
    }

    bool saved = false;
    const std::map<std::string, int> fields (sd->GetFieldNames());
    for (std::map<std::string, int>::const_iterator itr = fields.begin();
         itr != fields.end(); ++itr)
    {
        std::string name (itr->first);
        const Currency& value = sd->GetValue(0, 0, name);
        matvar_t* var = (append) ? Mat_VarRead(m, name.c_str()) : nullptr;
        if (var)
        {
            const Currency& matcur = omlMatio->MatVarToCurrency(var, eval);
            std::vector<Currency> inputs2;
            inputs2.emplace_back(matcur);
            inputs2.emplace_back(value);
            Currency isequal = eval.CallFunction("isequal", inputs2);
            if (isequal.IsLogical() && isequal.Scalar() == 1)
            {
                // Nothing to update
                saved = true;
                continue;
            }

            int result = Mat_VarDelete(m, name.c_str());
            Mat_VarFree(var);
            var = nullptr;
            if (result != 0)
            {
                omlMatio->AddWarning("Ignoring [" + name + "]; " +
                                    omlMatio->MatioError(result));
            }
        }
        var = omlMatio->CurrencyToMatVar(name.c_str(), value);
        if (var)
        {
            int result = Mat_VarWrite(m, var, compression);
            if (result != 0)
            {
                omlMatio->AddWarning("Ignoring [" + name + "]; " +
                                    omlMatio->MatioError(result));
            }
            else
            {
                saved = true;
            }
            Mat_VarFree(var);
        }
    }
    return saved;
}
