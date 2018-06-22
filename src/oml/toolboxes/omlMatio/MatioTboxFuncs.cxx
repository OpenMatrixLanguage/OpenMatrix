/**
* @file MatioTboxFuncs.cxx
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

#include "MatioTboxFuncs.h"

#include <cassert>
#include <fstream>

#include "BuiltInFuncsUtils.h"
#include "EvaluatorInt.h" 
#include "OML_Error.h"
#include "StructData.h" 
#include "matio.h"

// Returns true after loading file in ascii format
bool ReadASCIIFile(EvaluatorInterface           eval, 
                   const std::vector<Currency>& inputs, 
                   std::vector<Currency>&       outputs);
// Returns true after saving file in ascii format
bool SaveAsciiFile(EvaluatorInterface              eval, 
                   const std::string&              filename,
                   const std::vector<std::string>& vars);

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

void log_helper(int level, char* message)
{
	throw OML_Error(message);
}

Currency MatVarToCurrency(matvar_t* var)
{
	if (var->class_type == MAT_C_DOUBLE)
	{
		if (var->rank == 1) // not sure about this
		{
		}
		else if (var->rank == 2) // matrix or vector or scalar
		{
			if ((var->dims[0] == 1) & (var->dims[1] == 1))
			{
				if (var->isComplex) // complex
				{
					mat_complex_split_t *complex_data = (mat_complex_split_t*)var->data;
					double *rp = (double*)complex_data->Re;
					double *ip = (double*)complex_data->Im;
					return hwComplex(*rp, *ip);
				}
				else // scalar
				{
					double* my_val = (double*)var->data;
					return *my_val;
				}
			}
			else // matrix
			{ 
				if (var->isComplex) // complex matrix
				{
					int num_rows = static_cast<int>(var->dims[0]);
					int num_cols = static_cast<int>(var->dims[1]);

					hwMatrix* mat = EvaluatorInterface::allocateMatrix(num_rows, num_cols, hwMatrix::COMPLEX);
							                    
					mat_complex_split_t *complex_data = (mat_complex_split_t*)var->data;
					double *rp = (double*)complex_data->Re;
					double *ip = (double*)complex_data->Im;

					for (int j=0; j<num_cols; j++)
					{
						for (int k=0; k<num_rows; k++)
						{
							int     idx  = num_rows*j + k;
							mat->z(k, j) = hwComplex(*(rp+idx), *(ip+idx));
						}
					}

					return mat;
				}
				else
				{
					int num_rows = static_cast<int>(var->dims[0]);
					int num_cols = static_cast<int>(var->dims[1]);

					hwMatrix* mat = EvaluatorInterface::allocateMatrix(num_rows, num_cols, hwMatrix::REAL);

					for (int j=0; j<num_cols; j++)
					{
						for (int k=0; k<num_rows; k++)
						{
							int     idx  = num_rows*j + k;
							double* next = (double*)var->data + idx;
							(*mat)(k, j) = *next;
						}
					}
										
					return mat;
				}
			}
		}
		else // ND matrix
		{
			int size = 1;
            int rank = var->rank;
           
            std::vector<int> dims;
            dims.reserve(rank);

			for (int j = 0; j < rank; ++j)
			{
                int val = static_cast<int>(var->dims[j]);
				dims.push_back(val);
				size *= val;
			}

			hwMatrixN* mat_n = NULL;

			if (var->isComplex)
			{
				mat_n = new hwMatrixN(dims, hwMatrixN::COMPLEX);

				mat_complex_split_t *complex_data = (mat_complex_split_t*)var->data;
				double *rp = (double*)complex_data->Re;
				double *ip = (double*)complex_data->Im;

				for (int k=0; k<size; k++)
					mat_n->z(k) = hwComplex(rp[k], ip[k]);
			}
			else
			{
				mat_n = new hwMatrixN(dims, hwMatrixN::REAL);

				double* data = (double*)var->data;

				for (int k=0; k<size; k++)
					(*mat_n)(k) = data[k];
			}

			return mat_n;
		}
	}
	else if (var->class_type == MAT_C_UINT8)
	{
		if (var->rank == 1) // not sure about this
		{
		}
		else if (var->rank == 2) // matrix or vector or scalar
		{
			if ((var->dims[0] == 1) & (var->dims[1] == 1))
			{
				unsigned char* my_val = (unsigned char*)var->data;
				Currency ret((double)*my_val);

				if (var->isLogical)
					ret.SetMask(Currency::MASK_LOGICAL);

				return ret;				
			}
			else // matrix
			{ 
				int num_rows = static_cast<int>(var->dims[0]);
				int num_cols = static_cast<int>(var->dims[1]);

				hwMatrix* mat = EvaluatorInterface::allocateMatrix(num_rows, num_cols, hwMatrix::REAL);

				for (int j=0; j<num_cols; j++)
				{
					for (int k=0; k<num_rows; k++)
					{
						int     idx  = num_rows*j + k;
						char* next = (char*)var->data + idx;
						(*mat)(k, j) = (double)*next;
					}
				}
								
				Currency ret(mat);

				if (var->isLogical)
					ret.SetMask(Currency::MASK_LOGICAL);

				return ret;
			}
		}
		else
		{
            Mat_VarFree(var);
			throw OML_Error(OML_ERR_UNSUPPORTDIM);
		}
	}
	else if (var->class_type == MAT_C_CHAR)
	{
		if (var->rank == 2)
		{
			int num_rows = static_cast<int>(var->dims[0]);
			int num_cols = static_cast<int>(var->dims[1]);

			if ((num_rows == 1) && (num_cols != 0))
			{
				char* text = (char*)var->data;
				return text;
			}
			else
			{
				return "";
			}
		}
	}
	else if (var->class_type == MAT_C_CELL)
	{
		if (var->rank == 2)
		{
			HML_CELLARRAY* cells = EvaluatorInterface::allocateCellArray(
                static_cast<int>(var->dims[0]), static_cast<int>(var->dims[1]));
			matvar_t**     vals  = (matvar_t**)var->data;

			for (int j=0; j<cells->Size(); j++)
				(*cells)(j) = MatVarToCurrency(vals[j]);

			return cells;
		}
	}
	else if (var->class_type == MAT_C_STRUCT)
	{
		if (var->rank == 2)
		{
			StructData* sd = new StructData;
			sd->Dimension(static_cast<int>(var->dims[0]) - 1, 
                          static_cast<int>(var->dims[1]) - 1);

			char* const* fields = Mat_VarGetStructFieldnames(var);
			int num_fields = Mat_VarGetNumberOfFields(var);

			for (int j=0; j<num_fields; j++)
			{
				std::string field = fields[j];
				sd->addField(field);

				for (int k=0; k<var->dims[0]; k++)
				{
					for (int m=0; m<var->dims[1]; m++)
					{
						int index = k* static_cast<int>(var->dims[1]) + m;
						matvar_t* temp = Mat_VarGetStructFieldByName(var, field.c_str(), index);
						Currency cur = MatVarToCurrency(temp);
						sd->SetValue(k, m, field, cur);
					}
				}
			}

			return sd;
		}
	}

	return Currency(-1.0, Currency::TYPE_NOTHING);
}
//------------------------------------------------------------------------------
// Returns true after loading the given file using MATIO library [load command]
//------------------------------------------------------------------------------
bool OmlLoad(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    int nargin = (!inputs.empty()) ? static_cast<int>(inputs.size()) : 0;
    if (nargin < 1) 
        throw OML_Error(OML_ERR_NUMARGIN);

	std::vector<std::string> target_variables;
    target_variables.reserve(nargin);

	if (nargin > 1)
	{
		Currency lastcur = inputs[nargin-1];
        if (!lastcur.IsString())
            throw OML_Error(OML_ERR_STRING, nargin, OML_VAR_TYPE);

    	std::string last_input = inputs[nargin-1].StringVal();
    	if ((last_input == "-ascii") || (last_input == "-ASCII"))
			return ReadASCIIFile(eval, inputs, outputs);

		for (int j=1; j<nargin; ++j)
		{
            Currency cur = inputs[j];
            if (!cur.IsString())
                throw OML_Error(OML_ERR_STRING, j + 1, OML_VAR_TYPE);

    		target_variables.push_back(cur.StringVal());
		}
	}

	if (!inputs[0].IsString())
		throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);

	int nargout        = eval.GetNargoutValue();
	StructData* out_sd = NULL;

	if (nargout)
		out_sd = new StructData;

	Mat_LogInitFunc("OML", log_helper);

	std::string filename = inputs[0].StringVal();

	mat_t* m = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
    if (!m)
    {
        throw OML_Error("Error: invalid value in argument 1; cannot read file ["
            + BuiltInFuncsUtils::Normpath(filename) + "]");
	}

	if (m)
	{
		while (1)
		{
			matvar_t* var = Mat_VarReadNext(m);

			if (!var)
				break;

			if (target_variables.size())
			{
				std::vector<std::string>::iterator iter = std::find(target_variables.begin(), target_variables.end(), std::string(var->name));

				if (iter == target_variables.end())
                {
                    Mat_VarFree(var);
					continue;
                }
			}

			if (nargout == 0)
			{
				eval.SetValue(var->name, MatVarToCurrency(var));
			}
			else
			{
				std::string field = var->name;
				out_sd->addField(field);
				out_sd->SetValue(0, 0, var->name, MatVarToCurrency(var));
			}

			Mat_VarFree(var);
		}
	}

	Mat_Close(m);

	if (nargout == 1)
		outputs.push_back(out_sd);

	return true;
}
//------------------------------------------------------------------------------
// Converts currency to matvar
//------------------------------------------------------------------------------
matvar_t* CurrencyToMatVar(const char* name, const Currency& cur)
{
	matvar_t* var = nullptr;
	size_t    dims[2];

	int flags = 0;
	if (cur.IsLogical())
    {
		flags = MAT_F_LOGICAL;
    }

	if (cur.IsScalar())
	{
		dims[0]    = 1;
		dims[1]    = 1;
        double val = cur.Scalar();
		matvar_t* mvar = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, 
                                       dims, &val, flags);
        return mvar;
	}
	else if (cur.IsMatrix())
	{
		const hwMatrix* mtx = cur.Matrix();
		if (!mtx)
        {
			mtx = EvaluatorInterface::allocateMatrix();
        }
		dims[0] = mtx->M();
		dims[1] = mtx->N();

		if (mtx->IsReal())
		{
			var = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, (void*) mtx->GetRealData(), flags);
		}
		else
		{
			mat_complex_split_t t;
			double* temp_real = new double [mtx->Size()];
			double*	temp_imag = new double [mtx->Size()];

			for (int j=0; j<mtx->Size(); j++)
			{
				hwComplex cplx = mtx->z(j);
				temp_real[j] = cplx.Real();
				temp_imag[j] = cplx.Imag();
			}

			t.Re = temp_real;
			t.Im = temp_imag;

			var = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &t, MAT_F_COMPLEX);

			delete [] temp_real;
			delete [] temp_imag;
		}
	}
	else if (cur.IsNDMatrix())
	{
		const hwMatrixN* mtxn = cur.MatrixN();

		if (!mtxn)
			mtxn = EvaluatorInterface::allocateMatrixN();

		std::vector<int> my_dims = mtxn->Dimensions();

		size_t* dims = new size_t [my_dims.size()];

		for (int j=0; j<my_dims.size(); j++)
			dims[j] = my_dims[j];

		if (mtxn->IsReal())
		{
            int     matsize   = mtxn->Size();
			double* temp_real = new double [matsize];

			for (int j = 0; j < matsize; ++j)
            {
				temp_real[j] = (*mtxn)(j);
            }

            int nsize = (my_dims.empty()) ? 0 : static_cast<int>(my_dims.size());

			var = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, nsize, dims, 
                                temp_real, flags);

			delete [] temp_real;
            temp_real = nullptr;
		}
		else
		{
			mat_complex_split_t t;
			double* temp_real = new double [mtxn->Size()];
			double*	temp_imag = new double [mtxn->Size()];

			for (int j=0; j<mtxn->Size(); j++)
			{
				hwComplex cplx = mtxn->z(j);
				temp_real[j] = cplx.Real();
				temp_imag[j] = cplx.Imag();
			}

			t.Re = temp_real;
			t.Im = temp_imag;

			var = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 
                  static_cast<int>(my_dims.size()), dims, &t, MAT_F_COMPLEX);

			delete [] temp_real;
			delete [] temp_imag;
		}
	}
	else if (cur.IsString())
	{
		std::string str = cur.StringVal();
		const char* ptr = str.c_str();
		dims[0] = 1;
		dims[1] = str.length()+1; // for the NULL terminator

		var = Mat_VarCreate(name, MAT_C_CHAR, MAT_T_UINT8, 2, dims, (void*)ptr, flags);
	}
	else if (cur.IsComplex())
	{
		hwComplex cplx = cur.Complex();
		double real = cplx.Real();
		double imag = cplx.Imag();

		dims[0] = 1;
		dims[1] = 1;

		mat_complex_split_t t;
		t.Re = &real;
		t.Im = &imag;

		var = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, (void*)&t, MAT_F_COMPLEX);
	}
	else if (cur.IsCellArray())
	{
		HML_CELLARRAY* cells = cur.CellArray();
        assert(cells);
		matvar_t** matvar = (matvar_t**)malloc((cells->Size()+1)*sizeof(matvar_t*));

		for (int j=0; j<cells->Size(); j++)
			matvar[j] = CurrencyToMatVar("cell", (*cells)(j));

		matvar[cells->Size()] = NULL;

		dims[0] = cells->M();
		dims[1] = cells->N();

		var = Mat_VarCreate(name, MAT_C_CELL, MAT_T_CELL, 2, dims, (void*)matvar, flags);

		free(matvar);
	}
	else if (cur.IsStruct())
	{
		StructData* sd = cur.Struct();

		std::map<std::string, int>           field_names = sd->GetFieldNames();
		std::map<std::string, int>::iterator iter;

		int    num_fields = static_cast<int>(field_names.size());
		const char** temp_fields = new const char* [num_fields];

		int count = 0;

		for (iter = field_names.begin(); iter != field_names.end(); iter++)
		{
			std::string field = iter->first;
			temp_fields[count] = strdup(field.c_str());
			count++;
		}

		size_t     dims[2];
		dims[0] = sd->M();
		dims[1] = sd->N();

		var = Mat_VarCreateStruct(name, 2, dims, temp_fields, num_fields);

		for (int j=0; j<sd->M(); j++)
		{
			for (int k=0; k<sd->N(); k++)
			{
				for (iter = field_names.begin(); iter != field_names.end(); iter++)
				{
					std::string field = iter->first;
					matvar_t* temp_var = CurrencyToMatVar(field.c_str(), sd->GetValue(j, k, field));
					int index = k+j*sd->N();
					Mat_VarSetStructFieldByName(var, field.c_str(), index, temp_var);
				}
			}
		}

		for (int j=0; j<num_fields; j++)
			free((void*)temp_fields[j]);
		delete [] temp_fields;
	}

	return var;
}
//------------------------------------------------------------------------------
// Returns true after saving the given file using MATIO library [save command]
//------------------------------------------------------------------------------
bool OmlSave(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    int nargin = (!inputs.empty()) ? static_cast<int>(inputs.size()) : 0;
    if (nargin < 1) 
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    Currency cur1 = inputs[0];
    if (!cur1.IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }
    std::string filename (cur1.StringVal());
    bool        isascii = false;  // Saves as binary files by default

    mat_ft            version     = MAT_FT_MAT73;         // Default version
    matio_compression compression = MAT_COMPRESSION_ZLIB; // Default compression
    int               cmpidx      = -1;

    std::vector<std::string> target_variables;
    target_variables.reserve(nargin);
	
    for (int i = 1; i < nargin; ++i)
    {
        Currency cur = inputs[i];
        if (!cur1.IsString())
        {
            throw OML_Error(OML_ERR_STRING, i + 1, OML_VAR_TYPE);
        }
        std::string val (cur.StringVal());
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
		
    if (version != MAT_FT_MAT5 && compression != MAT_COMPRESSION_ZLIB)
    {
        std::string msg = "Warning: invalid compression specified for version;";
        msg += " ignoring option in argument " + 
               std::to_string(static_cast<long long>(cmpidx));
        BuiltInFuncsUtils::SetWarning(eval, msg);
        compression = MAT_COMPRESSION_ZLIB;
    }

    // Initializes logging functions
	Mat_LogInitFunc("OML", nullptr);

    // Creates the file to save
	mat_t* m = Mat_CreateVer(filename.c_str(), nullptr, version);
	if (!m)
    {
        throw OML_Error("Error: invalid value in argument 1; cannot create file ["
            + BuiltInFuncsUtils::Normpath(filename) + "]");
    }

    bool hastargets   = (!target_variables.empty());
	bool handle_check = eval.HasBuiltin("ishandle");
    std::vector<std::string> varnames (eval.GetVariableNames());

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

		matvar_t* var = CurrencyToMatVar(cur.GetOutputName().c_str(), cur);
		if (var)
		{
			Mat_VarWrite(m, var, compression);
			Mat_VarFree(var);
		}
	}

	Mat_Close(m);
	return true;
}

bool ReadASCIIFile(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	std::string file;
	
	if (inputs[0].IsString())
		file = inputs[0].StringVal();
	else
		throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);

	std::ifstream ifs;
	ifs.open(file, std::ifstream::in);

	if (!ifs.good())
    {
        throw OML_Error("Error: invalid value in argument 1; cannot read file ["
            + BuiltInFuncsUtils::Normpath(file) + "]");
    }
	char buffer[1024];

	std::vector<std::vector<double>> elements;

	while (1)
	{
        memset(buffer, 0, sizeof(buffer));
		ifs.getline(buffer, sizeof(buffer));

		int len = static_cast<int>(strlen(buffer));
		if (!len)
        {
			break;
        }

		std::string s = buffer;

		std::vector<double> row_elements;

		double val;

		if (s.length() == 0)
			continue;

		if (s[0] == '%')
			continue;

#ifndef OS_WIN
        if (s[s.size() - 1] == '\r')
        {
            s.erase(s.size() - 1);
            if (s.empty())
                continue;
        }
#endif

		std::size_t prev = 0;
		std::size_t pos;

		while ((pos = s.find_first_of(" ,%", prev)) != std::string::npos)
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
                    throw OML_Error("Error: invalid file in argument 1; cannot read ["
                        + BuiltInFuncsUtils::Normpath(file) + "]");
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
                throw OML_Error("Error: invalid file in argument 1; cannot read ["
                    + BuiltInFuncsUtils::Normpath(file) + "]");
			}

			row_elements.push_back(val);
		}

		if (row_elements.size())
			elements.push_back(row_elements);
	}

	int num_rows = static_cast<int>(elements.size());
	int num_cols = (num_rows) ? static_cast<int>(elements[0].size()) : 0;

	hwMatrix* result = new hwMatrix(num_rows, num_cols, hwMatrix::REAL);

	for (int j=0; j<num_rows; j++)
	{
		if (elements[j].size() != num_cols)
			throw OML_Error(HW_ERROR_INCOMPDIM);

		for (int k=0; k<num_cols; k++)
			(*result)(j,k) = (elements[j])[k];
	}

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
