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

#include "BuiltInFuncsUtils.h"
#include "EvaluatorInt.h" 
#include "FunctionInfo.h"
#include "OML_Error.h"
#include "StructData.h" 

#include "hwMatrixN.h"
#include "hwMatrixS.h"

#include "matio.h"

#define TBOXVERSION 2020

// Returns true after loading file in ascii format
bool ReadASCIIFile(EvaluatorInterface           eval, 
                   const std::vector<Currency>& inputs, 
                   std::vector<Currency>&       outputs);
// Returns true after saving file in ascii format
bool SaveAsciiFile(EvaluatorInterface              eval, 
                   const std::string&              filename,
                   const std::vector<std::string>& vars);
// Converts currency to matvar
matvar_t* CurrencyToMatVar(const char*     name, 
                           const Currency& cur,
                           mat_ft          version,
                           std::string&    warn);

Currency MatVarToCurrency(matvar_t*, std::string&);   // Var to currency
Currency MatCell2Currency(matvar_t*, std::string&);   // Cell to currency
Currency MatChar2Currency(matvar_t*, std::string&);   // Char array to currency
Currency MatSparse2Currency(matvar_t*, std::string&);  // Sparse mtx to currency

matvar_t* Currency2MatSparse(const char*, const Currency&); // Sparse mtx to matvar_t

std::string GetTypeString(matvar_t*);                 // Gets type description

Currency HandleInvalidDims(matvar_t*, std::string&);  // Sets warning and returns empty currency)
Currency HandleInvalidRank(matvar_t*, std::string&);  // Sets warning and returns empty currency


// Prints matio version for debugging
void PrintMatioVersion();
// Prints mat file version to stdout for debugging
void PrintMatioFileVersion(mat_t* m);
// Sets error from the matio library
void SetMatioMessage(int   level,
                     char* msg);
//# define OMLMATIO_DBG 1  // Uncomment to print debug info
#ifdef OMLMATIO_DBG
#    define OMLMATIO_PRINT(str, m) { std::cout << str << m << std::endl; }
#else
#    define OMLMATIO_PRINT(str, m) 0
#endif

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
// Converts matvar_t to currency
//------------------------------------------------------------------------------
Currency MatVarToCurrency(matvar_t* var, std::string& warn)
{
    assert(var);
    if (!var)
    {
        return Currency(-1.0, Currency::TYPE_NOTHING);
    }

    // \todo Convert all to a switch command
    matio_classes type = var->class_type;
    int           rank = var->rank;
    std::string msg("Cannot read variable [");
    if (var->name)
    {
        msg += var->name;
    }
    msg += "]";
    std::string newline = (!warn.empty()) ? "\n" : "";

    if (type == MAT_C_EMPTY) // 0  => Empty array
    {
        return EvaluatorInterface::allocateMatrix();
    }
    else if (type == MAT_C_CELL)   // 1  => Cell array
    {
        return MatCell2Currency(var, warn);
    }
    else if (type == MAT_C_CHAR)   // 4  => Character array
    {
        return MatChar2Currency(var, warn);
    }
    else if (type == MAT_C_OPAQUE || type == MAT_C_FUNCTION)
    {
        warn += newline + msg + "; unsupported type [" + GetTypeString(var) + "]";
        return Currency(-1.0, Currency::TYPE_NOTHING);
    }

    if (rank == 1 && (type == MAT_C_DOUBLE || type == MAT_C_SINGLE  ||
        type == MAT_C_UINT8 || type == MAT_C_UINT32 || MAT_C_UINT64 ||
        type == MAT_C_SPARSE))
    {
        return HandleInvalidRank(var, warn);
    }

	if (type == MAT_C_DOUBLE)
	{
        int num_rows = static_cast<int>(var->dims[0]);
        int num_cols = static_cast<int>(var->dims[1]);

        if (rank == 2) // matrix or vector or scalar
		{
            if (num_rows == 1 && num_cols == 1) // scalar or complex
			{
                if (!var->isComplex)
                {
                    double* my_val = (double*)var->data;
                    return *my_val; // scalar

                }
                mat_complex_split_t* cdata = (mat_complex_split_t*)var->data;
                if (cdata)
                {
                    double *rp = (double*)cdata->Re;
                    double *ip = (double*)cdata->Im;
                    return hwComplex(*rp, *ip);  // Complex
                }
                warn += newline + msg + "; empty complex data";
                return Currency(-1.0, Currency::TYPE_NOTHING);
			}

            hwMatrix::DataType mattype = (var->isComplex) ?
                                         hwMatrix::COMPLEX : hwMatrix::REAL;
            hwMatrix* mat = EvaluatorInterface::allocateMatrix(
                            num_rows, num_cols, mattype);

            if (!var->isComplex)
            {
                for (int j = 0; j < num_cols; ++j)
                {
                    for (int k = 0; k < num_rows; ++k)
                    {
                        int idx = num_rows * j + k;
                        double* next = (double*)var->data + idx;
                        (*mat)(k, j) = *next;
                    }
                }
                return mat;  // Real matrix
            }						                    
            mat_complex_split_t* cdata = (mat_complex_split_t*)var->data;
            if (cdata)
            {
                double *rp = (double*)cdata->Re;
                double *ip = (double*)cdata->Im;

                for (int j = 0; j < num_cols; ++j)
                {
                    for (int k = 0; k < num_rows; ++k)
                    {
                        int     idx = num_rows * j + k;
                        mat->z(k, j) = hwComplex(*(rp + idx), *(ip + idx));
                    }
                }
            }
            return mat;  // Complex matrix
		}
		// ND matrix
		int size = 1;
        std::vector<int> dims;
        dims.reserve(rank);
		for (int j = 0; j < rank; ++j)
		{
            int val = static_cast<int>(var->dims[j]);
			dims.push_back(val);
			size *= val;
		}

        hwMatrixN::DataType mattype = (var->isComplex) ?
            hwMatrixN::COMPLEX : hwMatrixN::REAL;

	    hwMatrixN* mat_n = new hwMatrixN(dims, mattype);
		if (var->isComplex)
		{
            mat_complex_split_t* cdata = (mat_complex_split_t*)var->data;
            if (cdata)
            {
                double *rp = (double*)cdata->Re;
                double *ip = (double*)cdata->Im;

                for (int k = 0; k < size; ++k)
                {
                    mat_n->z(k) = hwComplex(rp[k], ip[k]);
                }
            }
            return mat_n;
        }
		double* data = (double*)var->data;
        for (int k = 0; k < size; ++k)
        {
            (*mat_n)(k) = data[k];
        }
		return mat_n;
	}
    else if (type == MAT_C_UINT8 || type == MAT_C_INT64 ||
        type == MAT_C_UINT64 || type == MAT_C_INT32 || type == MAT_C_INT16 ||
        type == MAT_C_UINT16 || type == MAT_C_INT8)
    {
        if (rank == 2) // matrix or vector or scalar
        {
            if (var->dims[0] == 1 && var->dims[1] == 1)
            {
                unsigned char* my_val = (unsigned char*)var->data;
                Currency ret((double)*my_val);

                if (var->isLogical)
                {
                    ret.SetMask(Currency::MASK_LOGICAL);
                }
                return ret;
            }
            else // matrix
            {
                int num_rows = static_cast<int>(var->dims[0]);
                int num_cols = static_cast<int>(var->dims[1]);

                hwMatrix* mat = EvaluatorInterface::allocateMatrix(
                    num_rows, num_cols, hwMatrix::REAL);

                for (int j = 0; j < num_cols; ++j)
                {
                    for (int k = 0; k < num_rows; ++k)
                    {
                        int     idx = num_rows * j + k;
                        char* next = (char*)var->data + idx;
                        (*mat)(k, j) = (double)*next;
                    }
                }

                Currency ret(mat);
                if (var->isLogical)
                {
                    ret.SetMask(Currency::MASK_LOGICAL);
                }
                return ret;
            }
        }
        else
        {
            throw OML_Error(OML_ERR_UNSUPPORTDIM);
        }
    }
	else if (type == MAT_C_UINT32)
	{
        if (rank == 2) // matrix or vector or scalar
		{
            char *data1 = (char*)var->data;
            if (var->dims[0] == 1 && var->dims[1] == 1)
			{
                unsigned long val = *(mat_uint32_t*)data1;
				Currency ret(static_cast<int>(val));

                if (var->isLogical)
                {
                    ret.SetMask(Currency::MASK_LOGICAL);
                }
				return ret;				
			}
			else // matrix
			{ 
				int num_rows = static_cast<int>(var->dims[0]);
				int num_cols = static_cast<int>(var->dims[1]);

				hwMatrix* mat = EvaluatorInterface::allocateMatrix(
                                num_rows, num_cols, hwMatrix::REAL);

                for (int j = 0; j < num_cols; ++j)
				{
                    for (int k = 0; k < num_rows; ++k)
					{
                        int     idx = num_rows * j + k;
                        unsigned long next = *((mat_uint32_t*)data1 + idx);
                        (*mat)(k, j) = static_cast<int>(next);
					}
				}
								
				Currency ret(mat);
				if (var->isLogical)
                {
                    ret.SetMask(Currency::MASK_LOGICAL);
                }
				return ret;
			}
		}
		else
		{
			throw OML_Error(OML_ERR_UNSUPPORTDIM);
		}
	}
	else if (type == MAT_C_STRUCT)
	{
		if (rank == 2)
		{
			StructData* sd = new StructData;
			sd->Dimension(static_cast<int>(var->dims[0]) - 1, 
                          static_cast<int>(var->dims[1]) - 1);

			char* const* fields = Mat_VarGetStructFieldnames(var);
			int num_fields = Mat_VarGetNumberOfFields(var);

            for (int j = 0; j < num_fields; ++j)
            {
                std::string field(fields[j]);
                sd->addField(field);

                for (int k = 0; k < var->dims[0]; ++k)
                {
                    for (int m = 0; m <var->dims[1]; ++m)
                    {
                        int index = k* static_cast<int>(var->dims[1]) + m;
                        matvar_t* temp = Mat_VarGetStructFieldByName(
                            var, field.c_str(), index);
                        if (!temp)
                        {
                            continue;
                        }
                        Currency cur = MatVarToCurrency(temp, warn);
                        sd->SetValue(k, m, field, cur);
                    }
                }
            }

            return sd;
        }
	}
    else if (type == MAT_C_SINGLE)
    {
        if (rank == 2) // matrix or vector or scalar
        {
            if (var->dims[0] == 1 && var->dims[1] == 1)
            {
                if (var->isComplex) // complex
                {
                    mat_complex_split_t* cdata = (mat_complex_split_t*)var->data;
                    if (cdata)
                    {
                        float *rp = (float*)cdata->Re;
                        float *ip = (float*)cdata->Im;
                        return hwComplex(*rp, *ip);
                    }
                    warn += newline + msg + "; empty complex data";
                    return Currency(-1.0, Currency::TYPE_NOTHING);
                }
                else // scalar
                {
                    float* my_val = (float*)var->data;
                    return Currency(static_cast<double>(*my_val));
                }
            }
            else // matrix
            {
                int num_rows = static_cast<int>(var->dims[0]);
                int num_cols = static_cast<int>(var->dims[1]);

                if (var->isComplex) // complex matrix
                {
                    hwMatrix* mat = EvaluatorInterface::allocateMatrix(
                        num_rows, num_cols, hwMatrix::COMPLEX);

                    mat_complex_split_t* cdata = (mat_complex_split_t*)var->data;
                    if (cdata)
                    {
                        float *rp = (float*)cdata->Re;
                        float *ip = (float*)cdata->Im;

                        for (int j = 0; j < num_cols; ++j)
                        {
                            for (int k = 0; k < num_rows; ++k)
                            {
                                int     idx = num_rows * j + k;
                                mat->z(k, j) = hwComplex((double)*(rp + idx),
                                    (double)*(ip + idx));
                            }
                        }
                    }
                    return mat;
                }
                else
                {
                    hwMatrix* mat = EvaluatorInterface::allocateMatrix(
                        num_rows, num_cols, hwMatrix::REAL);

                    for (int j = 0; j<num_cols; j++)
                    {
                        for (int k = 0; k<num_rows; k++)
                        {
                            int     idx = num_rows * j + k;
                            float* next = (float*)var->data + idx;
                            (*mat)(k, j) = static_cast<double>(*next);
                        }
                    }

                    return mat;
                }
            }
        }
        else // ND matrix
        {
            int size = 1;

            std::vector<int> dims;
            dims.reserve(rank);

            for (int j = 0; j < rank; ++j)
            {
                int val = static_cast<int>(var->dims[j]);
                dims.push_back(val);
                size *= val;
            }

            hwMatrixN* mat_n = nullptr;

            if (var->isComplex)
            {
                mat_n = new hwMatrixN(dims, hwMatrixN::COMPLEX);

                mat_complex_split_t* cdata = (mat_complex_split_t*)var->data;
                if (cdata)
                {
                    float *rp = (float*)cdata->Re;
                    float *ip = (float*)cdata->Im;
                    for (int k = 0; k < size; ++k)
                    {
                        mat_n->z(k) = hwComplex(rp[k], ip[k]);
                    }
                }
            }
            else
            {
                mat_n = new hwMatrixN(dims, hwMatrixN::REAL);
                float* data = (float*)var->data;
                for (int k = 0; k < size; ++k)
                {
                    (*mat_n)(k) = static_cast<double>(data[k]);
                }
            }
            return mat_n;
        }
    }
    else if (type == MAT_C_SPARSE)  // Sparse matrix
    {
        return MatSparse2Currency(var, warn);
    }

    warn += newline + msg + "; type [" + GetTypeString(var) + "]";
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
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

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
                return ReadASCIIFile(eval, inputs, outputs);
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

    Mat_LogInitFunc("OML", SetMatioMessage);

	std::string filename (inputs[0].StringVal());
    std::string warn;

	mat_t* m = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
    if (!m)
    {
        throw OML_Error(OML_ERR_FILE_CANNOTREAD, 1);
	}

    PrintMatioFileVersion(m);  // Prints file version to stdout for debugging

    bool checktargets = (target_variables.empty()) ? false : true;
    while (1)
    {
        matvar_t* var = Mat_VarReadNext(m);
        if (!var)
        {
            break;
        }
        if (!var->name)
        {
            Mat_VarFree(var);
            continue;
        }

        std::string name(var->name);

        if (checktargets)
        {
            std::vector<std::string>::iterator iter = std::find(
                target_variables.begin(), target_variables.end(), name);
            if (iter == target_variables.end())
            {
                Mat_VarFree(var);
                continue;
            }
        }

        try
        {
            if (!out_sd)
            {
                eval.SetValue(name, MatVarToCurrency(var, warn));
            }
            else
            {
                out_sd->addField(name);
                out_sd->SetValue(0, 0, name, MatVarToCurrency(var, warn));
            }
        }
        catch (const OML_Error& e)
        {
            Mat_VarFree(var);
            Mat_Close(m);
            throw e;
        }
        Mat_VarFree(var);
    }

	Mat_Close(m);

    if (out_sd)
    {
        outputs.push_back(out_sd.release());
    }

    if (!warn.empty())
    {
        std::string newline;
        if (warn.find("\n") != std::string::npos)
        {
            newline = "\n";
        }
        BuiltInFuncsUtils::SetWarning(eval, "Warning: " + newline + warn);
    }
	return true;
}
//------------------------------------------------------------------------------
// Converts currency to matvar
//------------------------------------------------------------------------------
matvar_t* CurrencyToMatVar(const char*     name,
                           const Currency& cur,
                           mat_ft          version,
                           std::string&    warn)
{
    if (cur.IsFunctionHandle())  // No code to read/write function handle in matio
    {
        if (!warn.empty())
        {
            warn += "\n";
        }
        warn += "Cannot save function handle";
        if (name)
        {
            warn += " [" + std::string(name) + "]";
        }
        return nullptr;
    }

    size_t    dims[2] = { 1, 1 };
    int       flags = (cur.IsLogical()) ? MAT_F_LOGICAL : 0;
    matvar_t* var = nullptr;

    if (cur.IsScalar())
    {
        double val = cur.Scalar();
        var = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2,
            dims, &val, flags);
        return var;
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
            var = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims,
                (void*)mtx->GetRealData(), flags);
        }
        else
        {
            mat_complex_split_t t;
            double* temp_real = new double[mtx->Size()];
            double*	temp_imag = new double[mtx->Size()];

            for (int j = 0; j < mtx->Size(); ++j)
            {
                hwComplex cplx = mtx->z(j);
                temp_real[j] = cplx.Real();
                temp_imag[j] = cplx.Imag();
            }

            t.Re = temp_real;
            t.Im = temp_imag;

            var = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &t,
                MAT_F_COMPLEX);

            delete[] temp_real;
            delete[] temp_imag;

            temp_real = nullptr;
            temp_imag = nullptr;
        }
        return var;
    }
    else if (cur.IsNDMatrix())
    {
        const hwMatrixN* mtxn = cur.MatrixN();
        if (!mtxn)
        {
            mtxn = EvaluatorInterface::allocateMatrixN();
        }
        std::vector<int> my_dims = mtxn->Dimensions();

        size_t* dims = new size_t[my_dims.size()];

        for (int j = 0; j<my_dims.size(); j++)
            dims[j] = my_dims[j];

        if (mtxn->IsReal())
        {
            int     matsize = mtxn->Size();
            double* temp_real = new double[matsize];

            for (int j = 0; j < matsize; ++j)
            {
                temp_real[j] = (*mtxn)(j);
            }

            int nsize = (my_dims.empty()) ? 0 : static_cast<int>(my_dims.size());

            var = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, nsize, dims,
                temp_real, flags);

            delete[] temp_real;
            temp_real = nullptr;
        }
        else
        {
            mat_complex_split_t t;
            double* temp_real = new double[mtxn->Size()];
            double*	temp_imag = new double[mtxn->Size()];

            for (int j = 0; j<mtxn->Size(); j++)
            {
                hwComplex cplx = mtxn->z(j);
                temp_real[j] = cplx.Real();
                temp_imag[j] = cplx.Imag();
            }

            t.Re = temp_real;
            t.Im = temp_imag;

            var = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE,
                static_cast<int>(my_dims.size()), dims, &t, MAT_F_COMPLEX);

            delete[] temp_real;
            delete[] temp_imag;
        }
    }
    else if (cur.IsString())
    {
        std::string data(cur.StringVal());
        dims[0] = 1;
        dims[1] = data.length();
        if (version == MAT_FT_MAT5)
        {
            var = Mat_VarCreate(name, MAT_C_CHAR, MAT_T_UTF8, 2, dims,
                (void*)(data.c_str()), flags);
        }
        else
        {
            var = Mat_VarCreate(name, MAT_C_CHAR, MAT_T_UINT8, 2, dims,
                (void*)(data.c_str()), flags);
        }
        return var;
    }
    else if (cur.IsComplex())
    {
        hwComplex cplx = cur.Complex();
        double real = cplx.Real();
        double imag = cplx.Imag();

        mat_complex_split_t t;
        t.Re = &real;
        t.Im = &imag;

        var = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, (void*)&t, MAT_F_COMPLEX);
    }
    else if (cur.IsCellArray())
    {
        HML_CELLARRAY* cells = cur.CellArray();
        assert(cells);
        matvar_t** matvar = (matvar_t**)malloc((cells->Size() + 1) * sizeof(matvar_t*));

        for (int j = 0; j < cells->Size(); ++j)
        {
            Currency element = (*cells)(j);
            std::string elemname(element.GetOutputName());
            if (elemname.empty())
            {
                elemname = "cell";
            }
            bool isFunctionHandle = element.IsFunctionHandle();
            if (!isFunctionHandle)
            {
                matvar[j] = CurrencyToMatVar(elemname.c_str(), element, version, warn);
            }
            if (!matvar[j])
            {
                if (!warn.empty())
                {
                    warn += "\n";
                }
                if (!isFunctionHandle)
                {
                    warn += "Unrecognized data type in ";
                }
                else
                {
                    warn += "Cannot save function handle in ";
                }
                warn += std::string(name) + "("
                    + std::to_string(static_cast<long long>(j + 1));

                size_t doubledims[2] = { 1, 1 };
                double dummyval = std::numeric_limits<double>::quiet_NaN();
                matvar[j] = Mat_VarCreate(elemname.c_str(), MAT_C_DOUBLE,
                    MAT_T_DOUBLE, 2, doubledims, &dummyval, 0);
            }
        }

        matvar[cells->Size()] = nullptr;

        dims[0] = cells->M();
        dims[1] = cells->N();

        var = Mat_VarCreate(name, MAT_C_CELL, MAT_T_CELL, 2, dims, (void*)matvar, flags);

        free(matvar);
        matvar = nullptr;
    }
    else if (cur.IsStruct())
    {
        StructData* sd = cur.Struct();
        assert(sd);
        if (!sd)
        {
            return nullptr;
        }
        std::map<std::string, int> field_names(sd->GetFieldNames());
        int    num_fields = static_cast<int>(field_names.size());
        if (num_fields > 4091 && version == MAT_FT_MAT73)  // Limit in matio 1.5.12
        {
            std::string msg("Invalid number of fields; number of fields must");
            msg += " be lesser than 4092 in version 7.3;";
            msg += " save file with version 'v5'";
            throw OML_Error(msg);
        }

        const char** temp_fields = new const char*[num_fields];

        std::map<std::string, int>::const_iterator iter = field_names.begin();
        for (int count = 0; iter != field_names.end(); ++iter, ++count)
        {
            temp_fields[count] = iter->first.c_str();
        }

        dims[0] = sd->M();
        dims[1] = sd->N();

        var = Mat_VarCreateStruct(name, 2, dims, temp_fields, num_fields);

        for (int j = 0; j<sd->M(); j++)
        {
            for (int k = 0; k<sd->N(); k++)
            {
                for (iter = field_names.begin(); iter != field_names.end(); iter++)
                {
                    matvar_t* temp_var = nullptr;
                    std::string field (iter->first);
                    Currency element = sd->GetValue(j, k, field);
                    bool isFunctionHandle = element.IsFunctionHandle();
                    if (!isFunctionHandle)
                    {
                        temp_var = CurrencyToMatVar(field.c_str(),
                                   element, version, warn);
                    }
                    if (!temp_var)
                    {
                        if (!warn.empty())
                        {
                            warn += "\n";
                        }
                        if (!isFunctionHandle)
                        {
                            warn += "Unrecognized data type in ";
                        }
                        else
                        {
                            warn += "Cannot save function handle in ";
                        } 
                        warn += std::string(name) + "("
                            + std::to_string(static_cast<long long>(j + 1)) + ","
                            + std::to_string(static_cast<long long>(k + 1)) + ")"  
                            + "." + field;

                        size_t doubledims[2] = { 1, 1 };
                        double dummyval = std::numeric_limits<double>::quiet_NaN();
                        temp_var = Mat_VarCreate(field.c_str(), MAT_C_DOUBLE,
                            MAT_T_DOUBLE, 2, doubledims, &dummyval, 0);
                    }
                    if (temp_var)
                    {
                        int index = k + j * sd->N();
                        Mat_VarSetStructFieldByName(var, field.c_str(), index, temp_var);
                    }
                }
            }
        }

        delete[] temp_fields;
        temp_fields = nullptr;
    }
    else if (cur.IsSparse())
    {
        return Currency2MatSparse(name, cur);
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
		
    if (version != MAT_FT_MAT5 && compression != MAT_COMPRESSION_ZLIB)
    {
        std::string msg = "Warning: invalid compression specified for version;";
        msg += " ignoring option in argument " + 
               std::to_string(static_cast<long long>(cmpidx));
        BuiltInFuncsUtils::SetWarning(eval, msg);
        compression = MAT_COMPRESSION_ZLIB;
    }

    PrintMatioVersion();                     // Prints version - for debugging
    Mat_LogInitFunc("OML", SetMatioMessage); // Initializes logging functions

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

    std::string warn;
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
        if (cur.IsFunctionHandle())
        {
            BuiltInFuncsUtils::SetWarning(eval, "Warning: ignoring [" + name + "]; saving function handles is not supported");
            continue;
        }
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
            var = CurrencyToMatVar(cur.GetOutputName().c_str(), cur, version, warn);
            if (var)
            {
                Mat_VarWrite(m, var, compression);
                Mat_VarFree(var);
            }
        }
        catch (const OML_Error& e)
        {
            if (var)
            {
                Mat_VarFree(var);
            }
            Mat_Close(m);
            throw(e);
        }
    }

	Mat_Close(m);

    if (!warn.empty())
    {
        BuiltInFuncsUtils::SetWarning(eval, "Warning: " + warn);
    }
	return true;
}
//-----------------------------------------------------------------------------
// Returns true after reading an ascii file
//-----------------------------------------------------------------------------
bool ReadASCIIFile(EvaluatorInterface           eval, 
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
        throw OML_Error("Error: invalid value in argument 1; cannot read file ["
            + BuiltInFuncsUtils::Normpath(file) + "]");
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
//------------------------------------------------------------------------------
// Returns toolbox version
//------------------------------------------------------------------------------
double GetToolboxVersion(EvaluatorInterface eval)
{
    return TBOXVERSION;
}
//------------------------------------------------------------------------------
// Prints matio version for debugging
//------------------------------------------------------------------------------
void PrintMatioVersion()
{
#ifdef OMLMATIO_DBG
    int matversion[3];
    Mat_GetLibraryVersion(matversion, matversion + 1, matversion + 2);
    OMLMATIO_PRINT("Matio major version: ", matversion[0]);
    OMLMATIO_PRINT("Matio minor version: ", matversion[1]);
    OMLMATIO_PRINT("Matio release level: ", matversion[2]);
#endif
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
}
//------------------------------------------------------------------------------
// Gets type description of matvar
//------------------------------------------------------------------------------
std::string GetTypeString(matvar_t* var)
{
    if (!var)
    {
        return "NULL";
    }

    switch (var->class_type)
    {
        case MAT_C_EMPTY:    return "Empty matrix";
        case MAT_C_CELL:     return "Cell array";
        case MAT_C_STRUCT:   return "Struct";
        case MAT_C_OBJECT:   return "Object";
        case MAT_C_CHAR:     return "Character array";
        case MAT_C_SPARSE:   return "Sparse matrix";
        case MAT_C_DOUBLE:   return "Double-precision";
        case MAT_C_SINGLE:   return "Single-precision";
        case MAT_C_INT8:     return "Signed 8-bit integer";
        case MAT_C_UINT8:    return "Unsigned 8-bit integer";
        case MAT_C_INT16:    return "Signed 16-bit integer";
        case MAT_C_UINT16:   return "Unsigned 16-bit integer";
        case MAT_C_INT32:    return "Signed 32-bit integer";
        case MAT_C_UINT32:   return "Unsigned 32-bit integer";
        case MAT_C_INT64:    return "Signed 64-bit integer";
        case MAT_C_UINT64:   return "Unsigned 64-bit integer";
        case MAT_C_FUNCTION: return "Function handle";
        case MAT_C_OPAQUE:   return "Opaque";
        default: break;
    }

    return "Unknown: " + std::to_string(static_cast<long long>(var->class_type));
}
//------------------------------------------------------------------------------
// Converts mat cell array to currency
//------------------------------------------------------------------------------
Currency MatCell2Currency(matvar_t* var, std::string& warn)
{
    assert(var);
    assert(var->class_type == MAT_C_CELL);

    if (var->rank != 2)
    {
        return HandleInvalidRank(var, warn);
    }

    if (!var->dims)
    {
        return HandleInvalidDims(var, warn);
    }

    int m = static_cast<int>(var->dims[0]);
    int n = static_cast<int>(var->dims[1]);

    HML_CELLARRAY* cells = EvaluatorInterface::allocateCellArray(m, n);
    matvar_t**     vals = (matvar_t**)var->data;

    for (int j = 0; j < cells->Size(); ++j)
    {
        (*cells)(j) = MatVarToCurrency(vals[j], warn);
    }
    return cells;
}
//------------------------------------------------------------------------------
// Converts matvar of character array to currency
//------------------------------------------------------------------------------
Currency MatChar2Currency(matvar_t* var, std::string& warn)
{
    assert(var);
    assert(var->class_type == MAT_C_CHAR);

    if (var->rank != 2)
    {
        return HandleInvalidRank(var, warn);
    }

    if (!var->dims)
    {
        return HandleInvalidDims(var, warn);
    }

    int m = static_cast<int>(var->dims[0]);
    int n = static_cast<int>(var->dims[1]);

    if (var->data_type == MAT_T_UINT16 || var->data_type == MAT_T_UTF16)
    {
        std::vector<std::string> rowvals;
        rowvals.reserve(n);

        BuiltInFuncsUtils utils;
        char tmp[256];
        const mat_uint16_t *data = (const mat_uint16_t*)var->data;
        for (int i = 0; i < m; ++i)
        {
            std::string row;
            for (int j = 0; j < n; ++j)
            {
                const mat_uint16_t c = data[j * m + i];
                memset(tmp, 0, sizeof(tmp));
                if (c <= 0x7f)
                {
                    sprintf(tmp, "%c", c);
                }
                else if (c <= 0x7FF)
                {
                    sprintf(tmp, "%c%c", 0xC0 | (c >> 6), 0x80 |
                        (c & 0x3F));
                }
                else
                {
                    sprintf(tmp, "%c%c%c", 0xE0 | (c >> 12), 0x80 |
                        ((c >> 6) & 0x3F), 0x80 | (c & 0x3F));
                }
                row += tmp;
            }
            rowvals.push_back(row);
        }
        return utils.FormatOutput(rowvals, false, false, false, ' ');
    }

    std::unique_ptr<hwMatrix> mtx(EvaluatorInterface::allocateMatrix(
                                  m, n, hwMatrix::REAL));
    const char* tmp = (const char*)var->data;
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            (*mtx)(i, j) = (char)tmp[j * m + i];
        }
    }
    Currency result(mtx.release());
    result.SetMask(Currency::MASK_STRING);
    return result;
}
//------------------------------------------------------------------------------
// Sets warning and returns empty currency
//------------------------------------------------------------------------------
Currency HandleInvalidRank(matvar_t* var, std::string& warn)
{
    if (!var)
    {
        return Currency(-1.0, Currency::TYPE_NOTHING);
    }

    if (!warn.empty())
    {
        warn += "\n";
    }
    warn += "Cannot read variable [";
    if (var->name)
    {
        warn += var->name;
    }
    warn += "]; cannot process rank [" + std::to_string(
        static_cast<long long>(var->rank)) + "] for " + GetTypeString(var);
    return Currency(-1.0, Currency::TYPE_NOTHING);
}
//------------------------------------------------------------------------------
// Sets warning and returns empty currency
//------------------------------------------------------------------------------
Currency HandleInvalidDims(matvar_t* var, std::string& warn)
{
    if (!var)
    {
        return Currency(-1.0, Currency::TYPE_NOTHING);
    }
    if (!warn.empty())
    {
        warn += "\n";
    }
    warn += "Cannot read variable [";
    if (var->name)
    {
        warn += var->name;
    }
    warn += "]; dimensions are empty for " + GetTypeString(var);
    return Currency(-1.0, Currency::TYPE_NOTHING);
}
//------------------------------------------------------------------------------
// Converts matvar of sparse array to currency
//------------------------------------------------------------------------------
Currency MatSparse2Currency(matvar_t* var, std::string& warn)
{
    assert(var);
    assert(var->class_type == MAT_C_SPARSE);


    mat_sparse_t* sparse = (mat_sparse_t*)var->data;
    if (!sparse || !sparse->ir || !sparse->jc || !sparse->data || sparse->njc <= 1 )
    {
        return new hwMatrixS;
    }

    int m = static_cast<int>(var->dims[0]);
    int n = static_cast<int>(var->dims[1]);

    assert(sparse->njc == n + 1);

    // First n elements of jc
    std::vector<int> begincount(sparse->jc, sparse->jc + n);

    // Last n elements of jc
    std::vector<int> endcount(sparse->jc + 1, sparse->jc + sparse->njc);

    if (!var->isComplex)  // Real
    {
        return new hwMatrixS(m, n, begincount.data(), endcount.data(),
            sparse->ir, (double*)sparse->data);
    }

    std::vector<hwComplex> vec;

    mat_complex_split_t* cdata = (mat_complex_split_t*)sparse->data;
    if (cdata)
    {
        double* rp = (double*)cdata->Re;
        double* ip = (double*)cdata->Im;

        vec.reserve(sparse->ndata);
        for (int j = 0; j < sparse->ndata; ++j)
        {
            vec.push_back(hwComplex(*(rp + j), *(ip + j)));
        }
    }

    hwMatrixS* mtx =  new hwMatrixS(m, n, begincount.data(), endcount.data(),
            sparse->ir, vec.data());

    return mtx;
}
//------------------------------------------------------------------------------
// Converts currency to matvar_t for sparse matrices
//------------------------------------------------------------------------------
matvar_t* Currency2MatSparse(const char* name, const Currency& cur)
{
    size_t dims[2] = { 1, 1 };
    int    flags = (cur.IsLogical()) ? MAT_F_LOGICAL : 0;

    const hwMatrixS* spm = cur.MatrixS();
    if (!spm || spm->IsEmpty())
    {
        std::unique_ptr<hwMatrix> mtx(EvaluatorInterface::allocateMatrix());
        return Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims,
                (void*)mtx.get()->GetRealData(), flags);
    }

    int m   = spm->M();
    int n   = spm->N();
    int nnz = spm->NNZ();

    dims[0] = m;  // Number of total rows
    dims[1] = n;  // Number of total columns

    // Populate the sparse matio structure
    mat_sparse_t sparse;
    sparse.nzmax = nnz;  // Max non-zero elements

    // Array of size nzmax where ir[k] is the row of data[k].
    sparse.ir = (mat_int32_t*)spm->rows();
    sparse.nir = nnz; // Number of elements in ir

    // Array size n + 1 with jc[k] being the index into ir / data of the
    // first non - zero element for row k. Need to combine first n elements
    // of pointerB with last element of pointerE
    std::vector<int> jc(spm->pointerB(), spm->pointerB() + n);
    jc.push_back(*(spm->pointerE() + n - 1));
    sparse.jc = (mat_int32_t*)jc.data();
    sparse.njc = n + 1;

    sparse.ndata = nnz;                       // Number of values

    int rank = 2;

    matvar_t* tmp = nullptr;
    if (spm->IsReal())   // Real values
    {
        sparse.data  = (void*)spm->GetRealData(); // Array of values

        tmp = Mat_VarCreate(name, MAT_C_SPARSE, MAT_T_DOUBLE, 2, dims, &sparse, flags);
        if (tmp) // Recommended way of writing in matio
        {
            matvar_t* matvar2 = Mat_VarCreate(name, MAT_C_SPARSE,
                MAT_T_DOUBLE, 2, dims, tmp->data, flags);
            Mat_VarFree(tmp);
            return matvar2;
        }
        return Mat_VarCreate(name, MAT_C_SPARSE, MAT_T_DOUBLE, rank, dims, 
                &sparse, flags);
    }

    double* realvec = new double[nnz];
    double* imagvec = new double[nnz];

    hwComplex* cdata = const_cast<hwComplex*>(spm->GetComplexData());
    for (int i = 0; i < nnz; ++i)
    {
        hwComplex cplx = cdata[i];
        realvec[i] = cplx.Real();
        imagvec[i] = cplx.Imag();
    }

    mat_complex_split_t t{ nullptr, nullptr };
    t.Re = (void*)realvec;
    t.Im = (void*)imagvec;

    sparse.data = &t;

    tmp = Mat_VarCreate(name, MAT_C_SPARSE, MAT_T_DOUBLE, 2, dims, &sparse, MAT_F_COMPLEX);
    if (tmp) // Recommended way of writing in matio
    {
        matvar_t* matvar2 = Mat_VarCreate(name, MAT_C_SPARSE, MAT_T_DOUBLE, 2, dims, tmp->data, MAT_F_COMPLEX);
        Mat_VarFree(tmp);

        delete[] realvec;
        delete[] imagvec;

        realvec = nullptr;
        imagvec = nullptr;

        return matvar2;
    }

    matvar_t* matvar2 =  Mat_VarCreate(name, MAT_C_SPARSE, MAT_T_DOUBLE, rank, dims,
        &sparse, MAT_F_COMPLEX);

    delete[] realvec;
    delete[] imagvec;

    realvec = nullptr;
    imagvec = nullptr;

    return matvar2;
}
