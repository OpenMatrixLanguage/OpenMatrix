/**
* @file omlmatio.cxx
* @date August 2020
* Copyright (C) 2020-2022 Altair Engineering, Inc.
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

#include "omlmatio.h"

#include <cassert>

#include "BuiltInFuncsUtils.h"
#include "EvaluatorInt.h" 
#include "FunctionInfo.h"
#include "OML_Error.h"
#include "StructData.h" 

#include "CellND.cc"
#include "hwMatrixN_NMKL.h"
#include "hwMatrixS_NMKL.h"

#include "matio.h"

static std::string s_omlInternalClass = "_omlInternalClassName";

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
OmlMatio::OmlMatio(const std::string& filename, int verbose, MATFILEVERSION ver)
    : _idx(1)
    , _verbose(verbose)
    , _name   ("empty")
    , _version(ver)
{
    if (!filename.empty())
    {
        // This is in case the varible does not have a name defined
        _name = BuiltInFuncsUtils::GetBaseName(filename);
        size_t pos = _name.find(".");
        if (pos != std::string::npos)
        {
            _name = _name.substr(0, pos);
        }
    }
    _name += '_';
}
void OmlMatio::Reset()
{
    _idx = 1;
    _name = "empty_";
}
void OmlMatio::SetNamePrefix(const std::string& name)
{
    if (!name.empty())
    {
        // This is in case the varible does not have a name defined
        _name = BuiltInFuncsUtils::GetBaseName(name);
        size_t pos = _name.find(".");
        if (pos != std::string::npos)
        {
            _name = _name.substr(0, pos);
        }
        _name += '_';
    }
}
//------------------------------------------------------------------------------
// Returns empty currency and adds invalid dimension warning
//------------------------------------------------------------------------------
Currency OmlMatio::HandleInvalidDims(matvar_t* var, const std::string& name)
{
    AddWarning("Cannot read variable [" + name + "]; dimensions are [] for "
        + GetTypeString(var));

    return Currency(-1.0, Currency::TYPE_NOTHING);
}
//------------------------------------------------------------------------------
// Returns empty currency and adds warning about invalid rank
//------------------------------------------------------------------------------
Currency OmlMatio::HandleInvalidRank(matvar_t* var, const std::string& name)
{
    std::string msg("Cannot read variable [" + name + "]; cannot process rank [");
    if (var)
    {
        msg += std::to_string(static_cast<long long>(var->rank));
    }
    
    AddWarning(msg + "] for " + GetTypeString(var));

    return Currency(-1.0, Currency::TYPE_NOTHING);
}
//------------------------------------------------------------------------------
// Gets type description of matvar
//------------------------------------------------------------------------------
std::string OmlMatio::GetTypeString(matvar_t* var)
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
// Converts 2D/ND cell array to currency
//------------------------------------------------------------------------------
Currency OmlMatio::CellArray2Currency(matvar_t* var, EvaluatorInterface eval)
{
    assert(var);
    assert(var->class_type == MAT_C_CELL);

    
    std::string name = (var->name) ? var->name : "";
    if (!var->dims)
    {
        return HandleInvalidDims(var, name);
    }

    int rank = var->rank;
    if (var->rank < 2)
    {
        return HandleInvalidRank(var, name);
    }

    matvar_t** vals = (matvar_t**)var->data;
    Currency result;
    if (var->rank == 2)   // 2D cell array
    {
        int m = static_cast<int>(var->dims[0]);
        int n = static_cast<int>(var->dims[1]);

        std::unique_ptr<HML_CELLARRAY> cells(
            EvaluatorInterface::allocateCellArray(m, n));

        int cellsize = cells->Size();
        if (vals)
        {
            for (int j = 0; j < cellsize; ++j)
            {
                (*cells)(j) = MatVarToCurrency(vals[j], eval);
            }
        }
        result = (cells.release());
    }
    else
    {
        // ND cell array
        std::vector<int> dims;
        dims.reserve(rank);
        int numcells = 1;
        for (int i = 0; i < rank; ++i)
        {
            int val   = static_cast<int>(var->dims[i]);
            numcells *= val;
            dims.push_back(val);
        }

        std::unique_ptr<HML_ND_CELLARRAY> ndcell(new HML_ND_CELLARRAY(dims, 
            HML_ND_CELLARRAY::REAL));
        if (vals)
        {
            for (int i = 0; i < numcells; ++i)
            {
                (*(ndcell.get()))(i) = MatVarToCurrency(vals[i], eval);
            }
        }
        result = ndcell.release();
    }
    return result;
}
//------------------------------------------------------------------------------
// Converts char array to currency
//------------------------------------------------------------------------------
Currency OmlMatio::Char2Currency(matvar_t* var)
{
    assert(var);
    assert(var->class_type == MAT_C_CHAR);

    std::string name = (var->name) ? var->name : "";
    if (var->rank != 2)
    {
        return HandleInvalidRank(var, name);
    }
    if (!var->dims)
    {
        return HandleInvalidDims(var, name);
    }

    int m = static_cast<int>(var->dims[0]);
    int n = static_cast<int>(var->dims[1]);

    Currency result;
    if (var->data_type == MAT_T_UINT16 || var->data_type == MAT_T_UTF16)
    {
        std::vector<std::string> rowvals;
        rowvals.reserve(n);

        BuiltInFuncsUtils utils;
        const mat_uint16_t *data = (const mat_uint16_t*)var->data;
        for (int i = 0; i < m; ++i)
        {
            char tmp[256];

            std::wstring wrow;
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
                wrow += utils.StdString2WString(tmp);
            }
            rowvals.push_back(BuiltInFuncsUtils::WString2StdString(wrow));
        }
        result = utils.FormatOutput(rowvals, false, false, false, ' ');
    }
    else
    {
        std::unique_ptr<hwMatrix> mtx(EvaluatorInterface::allocateMatrix(
            m, n, true));
        const char* tmp = (const char*)var->data;
        if (tmp)
        {
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    (*mtx)(i, j) = (char)tmp[j * m + i];
                }
            }
        }
        result = mtx.release();
    }
    result.SetMask(Currency::MASK_STRING);
    return result;
}
//------------------------------------------------------------------------------
// Converts matvar of sparse array to currency
//------------------------------------------------------------------------------
Currency OmlMatio::Sparse2Currency(matvar_t* var)
{
    assert(var);
    assert(var->class_type == MAT_C_SPARSE);


    mat_sparse_t* sparse = (mat_sparse_t*)var->data;
    if (!sparse || !sparse->ir || !sparse->jc || !sparse->data || sparse->njc <= 1)
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
            (mat_int32_t *)sparse->ir, (double*)sparse->data);
    }

    std::vector<hwComplex> vec;

    mat_complex_split_t* cdata = (mat_complex_split_t*)sparse->data;
    if (cdata)
    {
        double* rp = (double*)cdata->Re;
        double* ip = (double*)cdata->Im;

        vec.reserve(sparse->ndata);
        int numdata = static_cast<int>(sparse->ndata);
        for (int j = 0; j < numdata; ++j)
        {
            vec.push_back(hwComplex(*(rp + j), *(ip + j)));
        }
    }

    hwMatrixS* mtx = new hwMatrixS(m, n, begincount.data(), endcount.data(),
        (mat_int32_t *)sparse->ir, vec.data());

    return mtx;
}
//------------------------------------------------------------------------------
// Gets name, creating a default one if name is not defined in matvar_t
//------------------------------------------------------------------------------
std::string OmlMatio::GetName(matvar_t* var,bool warn)
{
    assert(var);
    std::string name;
    if (var && var->name)
    {
        name = var->name;
    }
    
    if (name.empty())
    {
        name = _name + std::to_string(static_cast<long long>(_idx++));
        if (warn)
        {
        AddWarning("Variable name is [], creating name [" + name + "]");
    }
    }
    return name;
}
//------------------------------------------------------------------------------
// Utility to add to warning(s)
//------------------------------------------------------------------------------
void OmlMatio::AddWarning(const std::string& txt)
{
    assert(!txt.empty());
    _warn += "         " + txt + '\n';
}
//------------------------------------------------------------------------------
// Converts matvar_t to currency
//------------------------------------------------------------------------------
Currency OmlMatio::MatVarToCurrency(matvar_t* var, EvaluatorInterface eval)
{
    if (!var)
    {
        return Currency(-1.0, Currency::TYPE_NOTHING);
    }
    matio_classes classType = var->class_type;

    switch (classType)
    {
        case MAT_C_EMPTY:     // [0] Empty array
            return EvaluatorInterface::allocateMatrix();

        case MAT_C_CELL:      // [1] 2D/ND cell arrays
            return CellArray2Currency(var, eval);

        case MAT_C_STRUCT:    // [2] Structure     
            return Struct2Currency(var, eval);

        case MAT_C_CHAR:      // [4] Character array
            return Char2Currency(var);

        case MAT_C_SPARSE:    // [5] Sparse array
            return Sparse2Currency(var);

        case MAT_C_DOUBLE:    // [6] Double
        case MAT_C_SINGLE:    // [7] Single-precision
            return Double2Currency(var);

        case MAT_C_INT8:      // [8] Signed    8-bit integer
        case MAT_C_UINT8:     // [9] Unsigned  8-bit integer
        case MAT_C_INT16:     // [9] Signed   16-bit integer
        case MAT_C_UINT16:    // [9] Unsigned 16-bit integer
        case MAT_C_INT32:     // [9] Signed   32-bit integer
        case MAT_C_UINT32:    // [9] Unsigned 32-bit integer
        case MAT_C_INT64:     // [9] Signed   64-bit integer
        case MAT_C_UINT64:    // [9] Unsigned 64-bit integer
        {
            if (!var->isComplex)
            {
                return IntReal2Currency(var);
            }
            return IntComplex2Currency(var);
        }
        case MAT_C_OBJECT:    // [3]  Object 
        case MAT_C_FUNCTION:  // [16] Function 
        case MAT_C_OPAQUE:    // [17] Opaque
        default: 
        {
            std::string name = (var->name) ? var->name : "";
            return HandleInvalidType(var, name);
        }
    }

    return Currency(-1.0, Currency::TYPE_NOTHING);
}
//------------------------------------------------------------------------------
// Returns empty currency and adds warning about unsupported type
//------------------------------------------------------------------------------
Currency OmlMatio::HandleInvalidType(matvar_t* var, const std::string& name)
{
    AddWarning("Cannot read variable [" + name + "]; unsupported type [" +
             GetTypeString(var) + "]");

    return Currency(-1.0, Currency::TYPE_NOTHING);
}
//------------------------------------------------------------------------------
// Converts structure to currency
//------------------------------------------------------------------------------
Currency OmlMatio::Struct2Currency(matvar_t* var, EvaluatorInterface eval)
{
    assert(var);
    assert(var->class_type == MAT_C_STRUCT);

    std::string name = (var->name) ? var->name : "";
    if (!var->dims)
    {
        HandleInvalidDims(var, name);
    }

    int rank = var->rank;
    if (rank != 2) //\todo: Support more than 2 ranks?
    {
        return HandleInvalidRank(var, name);
    }

    int nfields         = Mat_VarGetNumberOfFields(var);
    char* const* fields = Mat_VarGetStructFieldnames(var);

    std::unique_ptr<StructData> sd(EvaluatorInterface::allocateStruct());

    int m = static_cast<int>(var->dims[0]);
    int n = static_cast<int>(var->dims[1]);
    sd->Dimension(m - 1, n - 1);

    std::string classname;

    for (int j = 0; j < nfields; ++j)
    {
        std::string field = (fields[j]) ? fields[j] : "";

        // No string sensitive comparison
        bool isclass(field == s_omlInternalClass);
        if (!isclass)
        {
            sd->addField(field); // Add as a field only if this is not a class
        }

        for (int row = 0; row < m; ++row)
        {
            for (int col = 0; col < n; ++col)
            {
                int index = row * n + col;
                matvar_t* temp = Mat_VarGetStructFieldByName(
                    var, field.c_str(), index);
                if (!temp)
                {
                    continue;
                }
                if (!isclass)
                {
                    sd->SetValue(row, col, field, MatVarToCurrency(temp, eval));
                }
                else
                {
                    Currency curclassname = MatVarToCurrency(temp, eval);
                    if (!curclassname.IsString())
                    {
                        AddWarning("Cannot load variable [" + name +
                            "]; reserved field [" + s_omlInternalClass +
                            "] must be a string");
                        continue;
                    }
                    classname = curclassname.StringVal();
                }
            }
        }
    }

    if (classname.empty())
    {
        return sd.release();
    }

    // Class object
    FunctionInfo* fi      = nullptr;
    FUNCPTR       funcptr = nullptr;
    eval.FindFunctionByName(classname, &fi, &funcptr, nullptr);
    if (!fi && !funcptr)
    {
        AddWarning("Cannot load variable [" + name +
            "]; cannot find class definition [" + classname + "]");
        return Currency(-1.0, Currency::TYPE_NOTHING);
    }

    Currency result(sd.release());
    result.SetClass(classname);
    
    return result;
}
//------------------------------------------------------------------------------
// Converts double to currency
//------------------------------------------------------------------------------
Currency OmlMatio::Double2Currency(matvar_t* var)
{
    assert(var);
    assert(var->class_type == MAT_C_DOUBLE || var->class_type == MAT_C_SINGLE);

    std::string name = (var->name) ? var->name : "";
    if (!var->dims)
    {
        HandleInvalidDims(var, name);
    }

    int rank = var->rank;
    if (rank < 2) 
    {
        return HandleInvalidRank(var, name);
    }

    int rows = static_cast<int>(var->dims[0]);
    int cols = static_cast<int>(var->dims[1]);

    if (rank == 2) // 2D matrix or vector or scalar
    {
        if (rows == 1 && cols == 1) // scalar or complex
        {
            return Number2Currency(var, name);
        }

        return Matrix2Currency(var, name);
    }
    // ND matrix
    long msize = 1;
    std::vector<int> dims;
    dims.reserve(rank);
    for (int j = 0; j < rank; ++j)
    {
        int val = static_cast<int>(var->dims[j]);
        dims.push_back(val);
        msize *= static_cast<long>(val);
    }

    hwMatrixN::DataType mattype = (var->isComplex) ?
        hwMatrixN::COMPLEX : hwMatrixN::REAL;

    std::unique_ptr<hwMatrixN> mat_n(new hwMatrixN(dims, mattype));
    if (var->isComplex)
    {
        mat_complex_split_t* cdata = (mat_complex_split_t*)var->data;
        if (cdata)
        {
            double* rp = nullptr;
            double* ip = nullptr;
            if (var->class_type == MAT_C_DOUBLE)
            {
                rp = (double*)cdata->Re;
                ip = (double*)cdata->Im;
                for (long k = 0; k < msize; ++k)
                {
                    mat_n->z(k) = hwComplex(rp[k], ip[k]);
                }

            }
            else 
            {
                float* rp = (float*)cdata->Re;
                float* ip = (float*)cdata->Im;
                for (long k = 0; k < msize; ++k)
                {
                    mat_n->z(k) = hwComplex(rp[k], ip[k]);
                }
            }
        }
    }
    else
    {
        if (var->class_type == MAT_C_DOUBLE)
        {
            double* data = (double*)var->data;
            for (long k = 0; k < msize; ++k)
            {
                (*mat_n)(k) = data[k];
            }
        }
        else
        {
            float* data = (float*)var->data;
            for (long k = 0; k < msize; ++k)
            {
                (*mat_n)(k) = (double)data[k];
            }
        }
    }
    Currency result(mat_n.release());
    return result;
}
//------------------------------------------------------------------------------
// Converts double number to currency
//------------------------------------------------------------------------------
Currency OmlMatio::Number2Currency(matvar_t* var, const std::string& name)
{
    assert(var);
    assert(var->class_type == MAT_C_DOUBLE || var->class_type == MAT_C_SINGLE);
    assert(var->rank == 2);
    assert(var->dims && var->dims[0] == 1 && var->dims[1] == 1);

    Currency result = Currency(-1.0, Currency::TYPE_NOTHING);
    bool     hasval = false;

    if (!var->isComplex)    // scalar
    {
        if (var->class_type == MAT_C_DOUBLE)
        {
            double* val = (double*)var->data;
            if (val)
            {
                result = Currency(*val);
                if (var->isLogical)
                {
                    result.SetMask(Currency::MASK_LOGICAL);
                }
                hasval = true;
            }
        }
        else
        {
            float* val = (float*)var->data;
            if (val)
            {
                result = Currency((double)*val);
                hasval = true;
            }
        }
    }
    else  // Complex
    {
        mat_complex_split_t* cdata = (mat_complex_split_t*)var->data;
        if (cdata)
        {
            hasval = true;
            if (var->class_type == MAT_C_DOUBLE)
            {
                double *rp = (double*)cdata->Re;
                double *ip = (double*)cdata->Im;
                result = hwComplex(*rp, *ip);  // Complex
            }
            else
            {
                float *rp = (float*)cdata->Re;
                float *ip = (float*)cdata->Im;
                result = hwComplex((double)*rp, (double)*ip);  // Complex

            }
        }
    }
    if (!hasval)
    {
        AddWarning("Cannot read variable [" + name + "]; no data");
    }
    return result;
}
//------------------------------------------------------------------------------
// Converts matrix to currency
//------------------------------------------------------------------------------
Currency OmlMatio::Matrix2Currency(matvar_t* var, const std::string& name)
{
    assert(var);
    assert(var->class_type == MAT_C_DOUBLE || var->class_type == MAT_C_SINGLE);
    assert(var->rank == 2);
    assert(var->dims);

    int rows = static_cast<int>(var->dims[0]);
    int cols = static_cast<int>(var->dims[1]);

    std::unique_ptr<hwMatrix> mat(EvaluatorInterface::allocateMatrix(
        rows, cols, !static_cast<bool>(var->isComplex)));

    if (!var->isComplex)
    {
        if (var->class_type == MAT_C_DOUBLE)
        {
            for (int j = 0; j < cols; ++j)
            {
                for (int k = 0; k < rows; ++k)
                {
                    int idx = rows * j + k;
                    double* next = (double*)var->data + idx;
                    (*mat)(k, j) = *next;
                }
            }
        }
        else
        {
            for (int j = 0; j < cols; ++j)
            {
                for (int k = 0; k < rows; ++k)
                {
                    int idx = rows * j + k;
                    float* next = (float*)var->data + idx;
                    (*mat)(k, j) = *next;
                }
            }
        }
        Currency result(mat.release());
        return result;  // Real matrix
    }
    mat_complex_split_t* cdata = (mat_complex_split_t*)var->data;
    assert(cdata);
    if (cdata)
    {
        if (var->class_type == MAT_C_DOUBLE)
        {
            double* rp = (double*)cdata->Re;
            double* ip = (double*)cdata->Im;

            for (int j = 0; j < cols; ++j)
            {
                for (int k = 0; k < rows; ++k)
                {
                    int     idx = rows * j + k;
                    mat->z(k, j) = hwComplex(*(rp + idx), *(ip + idx));
                }
            }
        }
        else
        {
            float* rp = (float*)cdata->Re;
            float* ip = (float*)cdata->Im;

            for (int j = 0; j < cols; ++j)
            {
                for (int k = 0; k < rows; ++k)
                {
                    int     idx = rows * j + k;
                    mat->z(k, j) = hwComplex(*(rp + idx), *(ip + idx));
                }
            }
        }
    }
    Currency result(mat.release());
    return result;  // Complex matrix
}
//------------------------------------------------------------------------------
// Converts numeric data to currency
//------------------------------------------------------------------------------
Currency OmlMatio::IntReal2Currency(matvar_t* var)
{
    assert(var);

    std::string name = (var->name) ? var->name : "";
    if (!var->dims)
    {
        HandleInvalidDims(var, name);
    }

    int rank = var->rank;
    if (rank < 2)
    {
        return HandleInvalidRank(var, name);
    }

    Currency result;
    if (rank == 2) // matrix or vector or scalar, different from matdump
    {
        if (var->dims[0] == 1 && var->dims[1] == 1)
        {
            int val = 0;
            switch (var->data_type)
            {
                case MAT_T_INT8:
                {
                    val = *(mat_int8_t*)var->data;
                    break;
                }
                case MAT_T_UINT8:
                {
                    val = *(mat_uint8_t*)var->data;
                    break;
                }
                case MAT_T_INT16: 
                {
                    if (var->isLogical)
                    {
                        unsigned char* myval = (unsigned char*)var->data;
                        val = (int)*myval;
                    }
                    else
                    {
                        val = *(mat_int16_t*)var->data;
                    }
                    break;
                }
                case MAT_T_UINT16: 
                {
                    val = *(mat_uint16_t*)var->data;
                    break;
                }
                case MAT_T_INT32: 
                {
                    val = *(mat_int32_t*)var->data;
                    break;
                }
                case MAT_T_UINT32: 
                {
                    char *data1 = (char*)var->data;
                    unsigned long tmp = *(mat_uint32_t*)data1;
                    val = static_cast<int>(tmp);
                    break;
                }
                default:
                {
                    unsigned char* my_val = (unsigned char*)var->data;
                    val = (int)*my_val;
                    break;
                }
            }
            result = val;
        }
        else
        {
            int rows = static_cast<int>(var->dims[0]);
            int cols = static_cast<int>(var->dims[1]);

            std::unique_ptr<hwMatrix> mat(EvaluatorInterface::allocateMatrix(
                rows, cols, true));

            for (int j = 0; j < cols; ++j)
            {
                for (int k = 0; k < rows; ++k)
                {
                    int     idx = rows * j + k;
                    char* next = (char*)var->data + idx;
                    (*mat)(k, j) = (double)*next;
                }
            }

            result = mat.release();
        }

    }
    if (var->isLogical)
    {
        result.SetMask(Currency::MASK_LOGICAL);
    }
    return result;
}
//------------------------------------------------------------------------------
// Converts numeric data to currency
//------------------------------------------------------------------------------
Currency OmlMatio::IntComplex2Currency(matvar_t* var)
{
    assert(var);

    std::string name = (var->name) ? var->name : "";
    if (!var->dims)
    {
        HandleInvalidDims(var, name);
    }

    int rank = var->rank;
    if (rank < 2)
    {
        return HandleInvalidRank(var, name);
    }

    Currency result;
    if (rank == 2) // matrix or vector or scalar, different from matdump
    {
        mat_complex_split_t* cdata = (mat_complex_split_t*)var->data;

        if (var->dims[0] == 1 && var->dims[1] == 1)
        {
            if (cdata)
            {
                double *rp = (double*)cdata->Re;
                double *ip = (double*)cdata->Im;
                result = hwComplex(*rp, *ip);  // Complex
            }
        }
        else
        {
            int rows = static_cast<int>(var->dims[0]);
            int cols = static_cast<int>(var->dims[1]);

            std::unique_ptr<hwMatrix> mat(EvaluatorInterface::allocateMatrix(
                rows, cols, false));

            double* rp = (double*)cdata->Re;
            double* ip = (double*)cdata->Im;

            for (int j = 0; j < cols; ++j)
            {
                for (int k = 0; k < rows; ++k)
                {
                    int     idx = rows * j + k;
                    mat->z(k, j) = hwComplex(*(rp + idx), *(ip + idx));
                }
            }

            result = mat.release();
        }

    }
    if (var->isLogical)
    {
        result.SetMask(Currency::MASK_LOGICAL);
    }
    return result;
}
//------------------------------------------------------------------------------
// Converts currency to matvar
//------------------------------------------------------------------------------
matvar_t* OmlMatio::CurrencyToMatVar(const char* name, const Currency& cur)
{

    size_t    dims[2] = { 1, 1 };
    int       flags = (cur.IsLogical()) ? MAT_F_LOGICAL : 0;
    matvar_t* var = nullptr;

    if (cur.IsScalar())
    {
        double val = cur.Scalar();
        if (!cur.IsSingle())
        {
            return Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims,
                &val, flags);
        }
        float f = static_cast<float>(val);
        return Mat_VarCreate(name, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims, &f, flags);
    }
    else if (cur.IsMatrix())
    {
        return Currency2MatMatrix(name, cur);
    }
    else if (cur.IsNDMatrix())
    {
        return Currency2MatMatrixND(name, cur);
    }
    else if (cur.IsString())
    {
        std::string data(cur.StringVal());
        dims[0] = 1;

        BuiltInFuncsUtils utils;
        if (!utils.HasWideChars(data))
        {
            dims[1] = data.length();
            matio_types type = (_version == MATFILEVERSION_5) ?
                MAT_T_UTF8 : MAT_T_UINT8;
            return Mat_VarCreate(name, MAT_C_CHAR, type, 2, dims,
                (void*)(data.c_str()), 0);
        }
        std::wstring wstr(utils.StdString2WString(data));
        dims[1] = wstr.length();
        matio_types type = (_version == MATFILEVERSION_5) ?
            MAT_T_UTF16 : MAT_T_UINT8;

        return Mat_VarCreate(name, MAT_C_CHAR, type, 2, dims,
               (void*)(wstr.c_str()), 0);
    }
    else if (cur.IsComplex())
    {
        hwComplex cplx = cur.Complex();
        if (!cur.IsSingle())
        {
            double real = cplx.Real();
            double imag = cplx.Imag();
            mat_complex_split_t t = { &real, &imag };
            return Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims,
                (void*)&t, MAT_F_COMPLEX);
        }
        float real = static_cast<float>(cplx.Real());
        float imag = static_cast<float>(cplx.Imag());
        mat_complex_split_t t = { &real, &imag };
        return Mat_VarCreate(name, MAT_C_SINGLE, MAT_T_SINGLE, 2, dims,
            (void*)&t, MAT_F_COMPLEX);
    }
    else if (cur.IsCellArray())
    {
        return Currency2MatCell(name, cur);
    }
    else if (cur.IsStruct())
    {
        return Currency2MatStruct(name, cur);
    }
    else if (cur.IsSparse())
    {
        return Currency2MatSparse(name, cur);
    }
    else if (cur.IsNDCellArray())
    {
        return Currency2MatCellND(name, cur);
    }
    else if (cur.IsObject())
    {
        return Currency2MatObj(name, cur);
    }
    else
    {
        AddWarning("Ignoring [" + std::string(name) + "]; type [" +
            cur.GetTypeString() + "] not supported");
    }

    return var;
}
//------------------------------------------------------------------------------
// Converts currency to mat var of type nd cell
//------------------------------------------------------------------------------
matvar_t* OmlMatio::Currency2MatCellND(const char*     name,
    const Currency& cur)
{
    HML_ND_CELLARRAY* cell = cur.CellArrayND();
    std::vector<int> dims((cell) ? cell->Dimensions() : std::vector<int>());
    if (dims.empty())
    {
        AddWarning("Cannot save variable [" + std::string(name) +
            "]; dimensions are []");
        return nullptr;
    }

    int cellsize = cell->Size();
    if (cellsize < 3)
    {
        std::unique_ptr<HML_CELLARRAY> tmp(EvaluatorInterface::allocateCellArray());
        cell->ConvertNDto2D(*tmp.get());
        return CurrencyToMatVar(name, tmp.release());
    }

    int  rank = static_cast<int>(dims.size());

    std::vector<size_t> dims_data;
    dims_data.reserve(rank);
    for (int i = 0; i < rank; ++i)
    {
        dims_data.push_back(static_cast<size_t>(dims[i]));
    }

    matvar_t* master = Mat_VarCreate(name, MAT_C_CELL, MAT_T_CELL, rank,
        dims_data.data(), nullptr, 0);


    for (int i = 0; i < cellsize; ++i)
    {
        Currency element = (*cell)(i);
        matvar_t* elementvar = CurrencyToMatVar(element.GetOutputName().c_str(),
            element);
        if (elementvar)
        {
            Mat_VarSetCell(master, i, elementvar);
        }
        else
        {
            Mat_VarSetCell(master, i, GetDummyMatVar(element.GetOutputName()));
            if (_verbose)
            {
                AddWarning("Error creating cell element [" + element.GetOutputName()
                    + "] in [" + std::string(name) + "]");
            }
        }
    }

    return master;
}
//------------------------------------------------------------------------------
// Converts currency to matvar_t for sparse matrices
//------------------------------------------------------------------------------
matvar_t* OmlMatio::Currency2MatSparse(const char* name, const Currency& cur)
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

    int m = spm->M();
    int n = spm->N();
    int nnz = spm->NNZ();

    dims[0] = m;  // Number of total rows
    dims[1] = n;  // Number of total columns

                  // Populate the sparse matio structure
    mat_sparse_t sparse;
    sparse.nzmax = nnz;  // Max non-zero elements

                         // Array of size nzmax where ir[k] is the row of data[k].
    sparse.ir = (mat_uint32_t*)spm->rows();
    sparse.nir = nnz; // Number of elements in ir

                      // Array size n + 1 with jc[k] being the index into ir / data of the
                      // first non - zero element for row k. Need to combine first n elements
                      // of pointerB with last element of pointerE
    std::vector<int> jc(spm->pointerB(), spm->pointerB() + n);
    jc.push_back(*(spm->pointerE() + n - 1));
    sparse.jc = (mat_uint32_t*)jc.data();
    sparse.njc = n + 1;

    sparse.ndata = nnz;                       // Number of values

    int rank = 2;

    matvar_t* tmp = nullptr;
    if (spm->IsReal())   // Real values
    {
        sparse.data = (void*)spm->GetRealData(); // Array of values

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

    matvar_t* matvar2 = Mat_VarCreate(name, MAT_C_SPARSE, MAT_T_DOUBLE, rank, dims,
        &sparse, MAT_F_COMPLEX);

    delete[] realvec;
    delete[] imagvec;

    realvec = nullptr;
    imagvec = nullptr;

    return matvar2;
}
//------------------------------------------------------------------------------
// Converts currency to matvar of type matrix
//------------------------------------------------------------------------------
matvar_t* OmlMatio::Currency2MatMatrix(const char* name, const Currency& cur)
{
    assert(cur.IsMatrix());

    const hwMatrix* mtx = cur.Matrix();

    matio_classes ctype = MAT_C_DOUBLE;
    matio_types   dtype = MAT_T_DOUBLE;
    if (cur.IsSingle())
    {
        ctype = MAT_C_SINGLE;
        dtype = MAT_T_SINGLE;
    }
    if (!mtx || mtx->Size() == 0)
    {
        // If matrix is not created, there is a crash with Mat_VarCreate
        std::unique_ptr<hwMatrix> empty(EvaluatorInterface::allocateMatrix());

        size_t dims[2] = { 0, 0 };
        
        return Mat_VarCreate(name, MAT_C_EMPTY, dtype, 2, dims,
                (void*)empty.get()->GetRealData(), 0);
    }

    size_t dims[2] = { static_cast<size_t>(mtx->M()), static_cast<size_t>(mtx->N()) };
    if (mtx->IsReal())
    {
        int flags = (cur.IsLogical()) ? MAT_F_LOGICAL : 0;
        return Mat_VarCreate(name, ctype, dtype, 2, dims, 
            (void*)mtx->GetRealData(), flags);
    }
    int matsize = mtx->Size();
    

    if (!cur.IsSingle())
    {
        std::vector<double> realdata;
        std::vector<double> imagdata;
        realdata.reserve(matsize);
        imagdata.reserve(matsize);

        for (int j = 0; j < matsize; ++j)
        {
            hwComplex cplx = mtx->z(j);
            realdata.push_back(cplx.Real());
            imagdata.push_back(cplx.Imag());
        }

        mat_complex_split_t t;
        t.Re = realdata.data();
        t.Im = imagdata.data();
        return Mat_VarCreate(name, ctype, dtype, 2, dims, &t, MAT_F_COMPLEX);
    }
    std::vector<float> realdata;
    std::vector<float> imagdata;
    realdata.reserve(matsize);
    imagdata.reserve(matsize);

    for (int j = 0; j < matsize; ++j)
    {
        hwComplex cplx = mtx->z(j);
        realdata.push_back(static_cast<float>(cplx.Real()));
        imagdata.push_back(static_cast<float>(cplx.Imag()));
    }
    mat_complex_split_t t;
    t.Re = realdata.data();
    t.Im = imagdata.data();

    return Mat_VarCreate(name, ctype, dtype, 2, dims, &t, MAT_F_COMPLEX);
}
//------------------------------------------------------------------------------
// Converts currency to matvar of type matrix ND
//------------------------------------------------------------------------------
matvar_t* OmlMatio::Currency2MatMatrixND(const char* name, const Currency& cur)
{
    assert(cur.IsNDMatrix());

    const hwMatrixN* mtxn = cur.MatrixN();
    if (!mtxn || mtxn->Size() == 0)
    {
        // If matrix is not created, there is a crash with Mat_VarCreate
        std::unique_ptr<hwMatrix> empty(EvaluatorInterface::allocateMatrix());

        size_t dims[2] = { 0, 0 };
        return Mat_VarCreate(name, MAT_C_EMPTY, MAT_T_DOUBLE, 2, dims,
            (void*)empty.get()->GetRealData(), 0);
    }

    std::vector<int> my_dims = mtxn->Dimensions();
    int rank = static_cast<int>(my_dims.size());

    std::vector<size_t> dims;
    dims.reserve(rank);
    for (int j = 0; j < rank; ++j)
    {
        dims.push_back(static_cast<size_t>(my_dims[j]));
    }

    if (mtxn->IsReal())
    {
        int flags = (cur.IsLogical()) ? MAT_F_LOGICAL : 0;
        return Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, rank, 
            dims.data(), (void*)mtxn->GetRealData(), flags);
    }
    int matsize = mtxn->Size();

    std::vector<double> realdata;
    std::vector<double> imagdata;
    realdata.reserve(matsize);
    imagdata.reserve(matsize);

    for (int j = 0; j < matsize; ++j)
    {
        hwComplex cplx = mtxn->z(j);
        realdata.push_back(cplx.Real());
        imagdata.push_back(cplx.Imag());
    }
    mat_complex_split_t t;
    t.Re = realdata.data();
    t.Im = imagdata.data();

    return Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, rank, dims.data(),
        &t, MAT_F_COMPLEX);
}
//------------------------------------------------------------------------------
// Creates a dummy variable of value NaN
//------------------------------------------------------------------------------
matvar_t* OmlMatio::GetDummyMatVar(const std::string& name)
{
    double val = std::numeric_limits<double>::quiet_NaN();
    size_t dims[2] = { 1, 1 };
    return Mat_VarCreate(name.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, 
        &val, 0);
}
//------------------------------------------------------------------------------
// Converts currency to matvar of type cell
//------------------------------------------------------------------------------
matvar_t* OmlMatio::Currency2MatCell(const char* name, const Currency& cur)
{
    assert(cur.IsCellArray());

    HML_CELLARRAY* cells = cur.CellArray();
    assert(cells);
    if (!cells)
    {
        size_t dims[2] = { 0, 0 };
        return Mat_VarCreate(name, MAT_C_CELL, MAT_T_CELL, 2, dims, nullptr, 0);
    }

    size_t dims[2] = { static_cast<size_t>(cells->M()),
                       static_cast<size_t>(cells->N()) };

    matvar_t* master = Mat_VarCreate(name, MAT_C_CELL, MAT_T_CELL, 2, dims,
        nullptr, 0);

    int cellsize = cells->Size();

    for (int j = 0; j < cellsize; ++j)
    {
        Currency element = (*cells)(j);
        std::string elemname(element.GetOutputName());
        if (elemname.empty())
        {
            elemname = "cell";
        }
        matvar_t* elementvar = CurrencyToMatVar(elemname.c_str(), element);

        if (!elementvar)
        {
            elementvar = GetDummyMatVar(elemname);
            if (_verbose)
            {
                std::string msg("Ignoring [" + std::string(name) + "(" +
                    std::to_string(static_cast<long long>(j + 1)) + ")]; ");
                AddWarning(msg);
            }
        }
        Mat_VarSetCell(master, j, elementvar);
    }

    return master;
}
//------------------------------------------------------------------------------
// Converts currency to matvar of type struct
//------------------------------------------------------------------------------
matvar_t* OmlMatio::Currency2MatStruct(const char* name, const Currency& cur)
{
    assert(cur.IsStruct());
    StructData* sd = cur.Struct();
    if (!sd)
    {
        matvar_t var1;
        size_t dims[2] = { 0, 1 };
        return Mat_VarCreate(name, MAT_C_STRUCT, MAT_T_STRUCT, 2, dims, &var1, 0);
    }

    std::map<std::string, int> names(sd->GetFieldNames());
    int    nfields = static_cast<int>(names.size());
    if (nfields > 4091 && _version == MATFILEVERSION_73)  // Limit in matio v7.3
    {
        // This happens with matio v1.5.17 as well. 
        // \todo: Check if this is still valid when matio version is updated
        std::string msg("Ignoring [" + std::string(name) + "]; ");
        msg += "number of struct fields cannot exceed 4091 in mat version 7.3;";
        msg += " save [" + std::string(name) + "] with 'v5' option";
        AddWarning(msg);
        return nullptr;
    }

    int m = sd->M();
    int n = sd->N();
    size_t dims[2] = { static_cast<size_t>(m),  static_cast<size_t>(n) };

    std::vector <const char*> fields;
    fields.reserve(nfields);
    std::map<std::string, int>::const_iterator iter = names.begin();
    for (int count = 0; iter != names.end(); ++iter, ++count)
    {
        fields.push_back(iter->first.c_str());
    }

    matvar_t* var = Mat_VarCreateStruct(name, 2, dims, fields.data(), nfields);

    for (int j = 0; j < m; ++j)
    {
        for (int k = 0; k < n; ++k)
        {
            for (iter = names.begin(); iter != names.end(); ++iter)
            {
                std::string field(iter->first);
                Currency element = sd->GetValue(j, k, field);
                matvar_t* temp_var = CurrencyToMatVar(field.c_str(), element);
                if (!temp_var)
                {
                    temp_var = GetDummyMatVar(field);
                    if (_verbose)
                    {
                        std::string msg("Ignoring [" + std::string(name) + "(" +
                            std::to_string(static_cast<long long>(j + 1)) + "," +
                            std::to_string(static_cast<long long>(k + 1)) + ")" +
                            "." + field + "]; ");
                        AddWarning(msg);
                    }
                }
                int index = k + j * sd->N();
                Mat_VarSetStructFieldByName(var, field.c_str(), index, temp_var);
            }
        }
    }

    return var;
}
//------------------------------------------------------------------------------
// Gets warning(s)
//------------------------------------------------------------------------------
std::string OmlMatio::GetWarning()
{
    BuiltInFuncsUtils utils;
    utils.StripTrailingNewline(_warn);
    if (_warn.empty())
    {
        return _warn;
    }

    if (_warn.find('\n') != std::string::npos) // Multiple errors
    {
        return "Warnings:\n" + _warn;
    }

    _warn = utils.LTrim(_warn);
    return "Warning: " + _warn;
}
//------------------------------------------------------------------------------
// Converts currency to matvar of type object
//------------------------------------------------------------------------------
matvar_t* OmlMatio::Currency2MatObj(const char* name, const Currency& cur)
{
    assert(cur.IsObject());
    const StructData* sd = cur.Struct();
    if (!sd)
    {
        matvar_t var1;
        size_t dims[2] = { 0, 1 };
        return Mat_VarCreate(name, MAT_C_STRUCT, MAT_T_STRUCT, 2, dims, &var1, 0);
    }

    std::map<std::string, int> names(sd->GetFieldNames());
    if (!names.empty() && names.find(s_omlInternalClass) != names.end())
    {
        std::string msg("Ignoring [" + std::string(name) + 
            "]; cannot save objects with a reserved field [" +
            s_omlInternalClass + "]; ");
        AddWarning(msg);
        return nullptr;
    }

    std::unique_ptr<StructData> sd2(new StructData(*sd));

    sd2->SetValue(0, 0, s_omlInternalClass, cur.GetClassname());

    return Currency2MatStruct(name, Currency(sd2.release()));
}
//------------------------------------------------------------------------------
// Adds warning about invalid dimensions
//------------------------------------------------------------------------------
void OmlMatio::AddInvalidDimsWarning(matvar_t* var, const std::string& name)
{
    AddWarning("Cannot read variable [" + name + "]; dimensions are [] for "
        + GetTypeString(var));
}
//------------------------------------------------------------------------------
// Adds warning about invalid rank
//------------------------------------------------------------------------------
void OmlMatio::AddInvalidRankWarning(matvar_t* var, const std::string& name)
{
    std::string msg("Cannot read variable [" + name + "]; cannot process rank [");
    if (var)
    {
        msg += std::to_string(static_cast<long long>(var->rank));
    }

    AddWarning(msg + "] for " + GetTypeString(var));
}
//------------------------------------------------------------------------------
// Helper method to get children
//------------------------------------------------------------------------------
std::vector<matvar_t*> OmlMatio::GetChildren(EvaluatorInterface eval, matvar_t* var)
{
    std::vector<matvar_t*> children;
    if (!var || !(var->class_type == MAT_C_CELL || var->class_type == MAT_C_STRUCT))
    {
        return children;
    }
    std::string name = (var->name) ? var->name : "";
    if (!var->dims)
    {
        AddInvalidDimsWarning(var, name);
        return children;
    }

    int rank = var->rank;
    if (rank < 2)
    {
        AddInvalidRankWarning(var, name);
        return children;
    }

    if (var->class_type == MAT_C_CELL)
    {
        matvar_t** vals = (matvar_t**)var->data;
        if (vals)
        {

            int numcells = 1;
            for (int i = 0; i < var->rank; ++i)
            {
                int val = static_cast<int>(var->dims[i]);
                numcells *= val;
            }

            children.reserve(numcells);
            for (int i = 0; i < numcells && !eval.IsInterrupt(); ++i)
            {
                children.emplace_back(vals[i]);
            }
        }
    }
    else
    {   
        // Struct
        if (rank != 2)
        {
            AddInvalidRankWarning(var, name);
            return children;
        }

        int          m       = static_cast<int>(var->dims[0]);
        int          n       = static_cast<int>(var->dims[1]);
        if (m == 1 && n == 1)
        {
        int          nfields = Mat_VarGetNumberOfFields(var);
            children.reserve(nfields);
            for (int i = 0; i < nfields && !eval.IsInterrupt(); ++i)
        {
                matvar_t* child = Mat_VarGetStructFieldByIndex(var, i, 0);
                if (child)
                {
                    children.emplace_back(child);
                }
            }
        }
        else
        {
            int numchildren = m * n;
            children.reserve(numchildren);

            for (int i = 0; i < numchildren && !eval.IsInterrupt(); ++i)
            {
                matvar_t* child = Mat_VarGetStructsLinear(var, i, 1, 1, 1);
                if (child)
                {
                    children.emplace_back(child);
                }
            }
        }
    }
    return children;
}