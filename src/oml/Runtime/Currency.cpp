/**
* @file Currency.h
* @date August 2013
* Copyright (C) 2013-2021 Altair Engineering, Inc.  
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

#include "Currency.h"

#include "FunctionInfo.h"
#include "OMLTree.h"
#include "CellDisplay.h"
#include "CellNDisplay.h"
#include "CurrencyDisplay.h"
#include "MatrixDisplay.h"
#include "MatrixNDisplay.h"
#include "OutputFormat.h"
#include "SparseDisplay.h"
#include "StringDisplay.h"
#include "StructData.h"
#include "StructDisplay.h"
#include "GeneralFuncs.h"

#include "hwMatrixN.h"
#include "hwMatrixS.h"

#include "ExprCppTreeLexer.h"

#include <cassert>
#include <sstream>
#include <iomanip>

#ifndef OS_WIN
#define sprintf_s sprintf
#endif

StringManager Currency::vm;
StringManager Currency::pm;
std::set<void*> Currency::ignore_cow_pointers;

bool Currency::_experimental = false;

Currency::Currency(double val): type(TYPE_SCALAR),  mask(MASK_DOUBLE), out_name(NULL), 
    _display(0), _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.value = val;
}

Currency::Currency(int val): type(TYPE_SCALAR), mask(MASK_DOUBLE), out_name(NULL), _display(0), 
    _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.value = val;
}

Currency::Currency(size_t val): type(TYPE_SCALAR), mask(MASK_DOUBLE), out_name(NULL), _display(0), 
    _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.value = (double)val;
}

Currency::Currency(double val, int in_type): type(CurrencyType(in_type)), mask(MASK_DOUBLE), out_name(NULL),
    _display(0), _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.value = val;
}

Currency::Currency(double val, std::string info): type(TYPE_ERROR), mask(MASK_NONE), out_name(NULL), _display(0), 
    _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	message = new std::string;
	*message = info;
	data.value = val;
}

Currency::Currency(bool logical_val) : type(TYPE_SCALAR), mask(MASK_LOGICAL), out_name(NULL), _display(0), 
    _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.value = logical_val ? 1 : 0;
}

Currency::Currency(const char* in_str): type(TYPE_MATRIX), mask(MASK_STRING), out_name(NULL),
    _display(0), _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	size_t length = strlen(in_str);

    data.mtx = ExprTreeEvaluator::allocateMatrix(length ? 1 : 0, (int)length, hwMatrix::REAL);
	for (size_t i = 0; i < length; i++)
	{
		unsigned char next_char = (unsigned char)in_str[i];
		(*data.mtx)((int)i) = next_char;

		if (next_char >= 192)
			_is_utf8 = true;
	}
}

Currency::Currency(const std::string& str): type(TYPE_MATRIX), mask(MASK_STRING), out_name(NULL),
    _display(0), _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
    data.mtx = ExprTreeEvaluator::allocateMatrix(str.length() ? 1 : 0, (int)str.length(), hwMatrix::REAL);
	for (unsigned int i = 0; i < str.length(); i++)
	{
		unsigned char next_char = (unsigned char)str[i];
		(*data.mtx)(i) = next_char;

		if (next_char >= 192)
			_is_utf8 = true;
	}
}

Currency::Currency(const std::vector<double>& in_data): type(TYPE_MATRIX), mask(MASK_DOUBLE), out_name(NULL),
    _display(0), _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{ 
	if (!in_data.empty())
	{
		double* data_ptr = new double[in_data.size()];
		memcpy(data_ptr, &in_data.front(), sizeof(double) * in_data.size());
		data.mtx = ExprTreeEvaluator::allocateMatrix(1, (int)in_data.size(), (void*)data_ptr, hwMatrix::REAL);
		data.mtx->OwnData(true);
	}
	else
	{
		data.mtx = ExprTreeEvaluator::allocateMatrix(1, 0, hwMatrix::REAL);
	}
}

Currency::Currency(const std::vector<std::string>& in_data) : type(TYPE_CELLARRAY), mask(MASK_DOUBLE), out_name(NULL),
_display(0), _outputType(OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	if (!in_data.empty())
	{
		data.cells = ExprTreeEvaluator::allocateCellArray((int)in_data.size(), 1);

		for (int j = 0; j < in_data.size(); ++j)
			(*data.cells)(j) = in_data[j];
	}
	else
	{
		// not sure if this should be 0x1 or 0x0
		data.cells = ExprTreeEvaluator::allocateCellArray();
	}
}

Currency::Currency(hwMatrix* in_data): type (TYPE_MATRIX), mask(MASK_DOUBLE), out_name(NULL),
    _display(0), _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.mtx = in_data;
}

Currency::Currency(hwMatrixN* in_data) : type(TYPE_ND_MATRIX), mask(MASK_DOUBLE), out_name(NULL),
_display(0), _outputType(OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.mtxn = in_data;

	if (data.mtxn)
	{
		if (data.mtxn->Dimensions().size() <= 2)
		{
			hwMatrix* temp = EvaluatorInterface::allocateMatrix();
			data.mtxn->ConvertNDto2D(*temp, false);

			DeleteMatrixN(data.mtxn);
			data.mtx = temp;
			type = TYPE_MATRIX;
		}
	}
}

Currency::Currency(hwMatrixS* in_data) : type(TYPE_SPARSE), mask(MASK_DOUBLE), out_name(NULL),
_display(0), _outputType(OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.mtxs = in_data;
}

Currency::Currency(const hwComplex& cplx): type (TYPE_COMPLEX), mask(MASK_DOUBLE), out_name(NULL),
    _display(0), _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.complex = new hwComplex(cplx);
}

Currency::Currency(): type(TYPE_MATRIX), mask(MASK_DOUBLE), out_name(NULL),
    _display(0), _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.mtx = NULL; 
}

Currency::Currency(HML_CELLARRAY* cell_array): type(TYPE_CELLARRAY), mask(MASK_NONE), out_name(NULL),
    _display(0), _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.cells = cell_array;
}

Currency::Currency(HML_ND_CELLARRAY* cell_array) : type(TYPE_ND_CELLARRAY), mask(MASK_NONE), out_name(NULL),
_display(0), _outputType(OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.cells_nd = cell_array;

	if (data.cells_nd)
	{
		if (data.cells_nd->Dimensions().size() <= 2)
		{
			HML_CELLARRAY* temp = EvaluatorInterface::allocateCellArray();
			data.cells_nd->ConvertNDto2D(*temp, false);

			DeleteCellsN(data.cells_nd);
			data.cells = temp;
			type = TYPE_CELLARRAY;
		}
	}
}

Currency::Currency(FunctionInfo* fi): type(TYPE_FUNCHANDLE), mask(MASK_NONE), out_name(NULL),
    _display(0), _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.func = fi;
	fi->IncrRefCount();
}

Currency::Currency(StructData* in_data): type(TYPE_STRUCT), mask(MASK_NONE), out_name(NULL),
    _display(0), _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.sd = in_data;
}

Currency::Currency(OutputFormat* fmt) : type(TYPE_FORMAT), mask(MASK_NONE), out_name(NULL),
    _display(0), _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.format = fmt;
}

Currency::Currency(Currency* ptr) : type(TYPE_POINTER), out_name(NULL), mask(MASK_NONE),
    _display(0), _outputType (OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.cur_ptr = ptr;
}

Currency::Currency(OMLTree* tree) : type(TYPE_TREE), out_name(NULL), mask(MASK_NONE),
_display(0), _outputType(OUTPUT_TYPE_DEFAULT), message(NULL), classname(NULL), _is_utf8(false), _is_linear_range(false)
{
	data.tree = tree;
}

//------------------------------------------------------------------------------
// Constructor for swig bound objects
//------------------------------------------------------------------------------
Currency::Currency(void* obj, const std::string& name) 
    : type         (TYPE_BOUNDOBJECT)
    , out_name     (NULL)
    , _display     (0)
    , _outputType  (OUTPUT_TYPE_DEFAULT)
    , classname    (NULL)
    , mask         (MASK_NONE)
{
    assert(!name.empty());

    data.boundobj = obj;
	classname     = const_cast<std::string*>(vm.GetStringPointer(name));
}
//------------------------------------------------------------------------------
//! Copies given currency
//! \param[in] cur Given currency
//------------------------------------------------------------------------------
void Currency::Copy(const Currency& cur)
{
	hwMatrix*         old_matrix   = NULL;
	hwMatrixN*        old_matrix_n = NULL;
	hwMatrixS*        old_matrix_s = NULL;
	HML_CELLARRAY*    old_cells    = NULL;
	HML_ND_CELLARRAY* old_cells_n = NULL;
	StructData*       old_sd       = NULL;
	hwComplex*        old_complex  = NULL;
	bool              was_scalar   = false;
	FunctionInfo*     old_fi       = NULL;

	classname = NULL;

	if (type == TYPE_SCALAR)
		was_scalar = true;
	else if (type == TYPE_MATRIX)
		old_matrix = data.mtx;
	else if (type == TYPE_ND_MATRIX)
			old_matrix_n = data.mtxn;
	else if (type == TYPE_SPARSE)
		old_matrix_s = data.mtxs;
	else if (type == TYPE_CELLARRAY)
		old_cells = data.cells;
	else if (type == TYPE_ND_CELLARRAY)
		old_cells_n = data.cells_nd;
	else if ((type == TYPE_STRUCT) || (type == TYPE_OBJECT))
		old_sd = data.sd;
	else if (type == TYPE_COMPLEX)
		old_complex = data.complex;
	else if (type == TYPE_FUNCHANDLE)
		old_fi = data.func;

	type             = cur.type;
	mask             = cur.mask;
	out_name         = cur.out_name;
    _display         = NULL; // Do not copy display, it will be set later
    _outputType      = cur._outputType;
	_is_utf8         = cur._is_utf8;
	_is_linear_range = cur._is_linear_range;

	if (type == TYPE_SCALAR)
	{
		data.value   = cur.data.value;
	}
	else if (type == TYPE_MATRIX)
	{
		data.mtx = cur.data.mtx;

		if (!IgnoreCoW(data.mtx))
		{
			if (data.mtx)
				data.mtx->IncrRefCount();
		}
	}
	else if (type == TYPE_ND_MATRIX)
	{
		data.mtxn = cur.data.mtxn;

		if (data.mtxn)
			data.mtxn->IncrRefCount();
	}
	else if (type == TYPE_SPARSE)
	{
		data.mtxs = cur.data.mtxs;

		if (data.mtxs)
			data.mtxs->IncrRefCount();
	}
	else if (type == TYPE_COMPLEX)
	{
		data.complex = new hwComplex(cur.Real(), cur.Imag());
	}
	else if (type == TYPE_CELLARRAY)
	{
		data.cells = cur.data.cells;

		if (!IgnoreCoW(data.cells))
		{
			if (data.cells)
				data.cells->IncrRefCount();
		}
	}
	else if (type == TYPE_ND_CELLARRAY)
	{
		data.cells_nd = cur.data.cells_nd;

		if (data.cells_nd)
			data.cells_nd->IncrRefCount();
	}
	else if (type == TYPE_STRUCT)
	{
		data.sd = cur.data.sd;

		if (data.sd)
			data.sd->IncrRefCount();
	}
	else if (type == TYPE_OBJECT)
	{
		data.sd = cur.data.sd;

		if (data.sd)
			data.sd->IncrRefCount();

		classname = cur.classname;
	}
	else if (type == TYPE_FORMAT)
	{
        data.format = new OutputFormat(*cur.data.format);
	}
	else if (type == TYPE_POINTER)
	{
		data.cur_ptr = cur.data.cur_ptr;
	}
	else if (type == TYPE_FUNCHANDLE)
	{
		data.func	= cur.data.func;

		if (data.func)
			data.func->IncrRefCount();
	}
	else if (type == TYPE_ERROR)
	{
		message  = new std::string;
		*message = cur.Message();
	}
	else if (type == TYPE_TREE)
	{
		data.tree = cur.data.tree;
	}
    else if (type == TYPE_BOUNDOBJECT)
    {
        data.boundobj = cur.data.boundobj; //\todo: Add ref counting
        classname     = cur.classname;
    }
	
	if (type != TYPE_POINTER)
	{
		if (was_scalar)
			;
		else if (old_matrix)
			DeleteMatrix(old_matrix);
		else if (old_matrix_n)
			DeleteMatrixN(old_matrix_n);
		else if (old_matrix_s)
			DeleteMatrixS(old_matrix_s);
		else if (old_cells)
			DeleteCells(old_cells);
		else if (old_cells_n)
			DeleteCellsN(old_cells_n);
		else if (old_sd)
			DeleteStruct(old_sd);
		else if (old_complex)
			delete old_complex;
		else if (old_fi)
			DeleteFunctionInfo(old_fi);
	}
}

void Currency::DeleteMatrix(hwMatrix* matrix)
{
	if (matrix)
	{
		if (IgnoreCoW(matrix))
			return;

		if (!matrix->IsMatrixShared())
		{
			delete matrix;

			if (matrix == data.mtx)
				data.mtx = NULL;
		}
		else
		{
			matrix->DecrRefCount();
		}
	}
}

void Currency::DeleteMatrixN(hwMatrixN* matrix)
{
	if (matrix)
	{
		if (!matrix->IsMatrixShared())
		{
			delete matrix;

			if (matrix == data.mtxn)
				data.mtxn = NULL;
		}
		else
		{
			matrix->DecrRefCount();
		}
	}
}

void Currency::DeleteMatrixS(hwMatrixS* matrix)
{
	if (matrix)
	{
		if (!matrix->IsMatrixShared())
		{
			delete matrix;

			if (matrix == data.mtxs)
				data.mtxs = NULL;
		}
		else
		{
			matrix->DecrRefCount();
		}
	}
}

void Currency::DeleteCells(HML_CELLARRAY* cells)
{
	if (cells)
	{
		if (IgnoreCoW(cells))
			return;

		if (!cells->IsMatrixShared())
		{
			delete cells;

			if (cells == data.cells)
				data.cells = NULL;
		}
		else
		{
			cells->DecrRefCount();
		}
	}
}

void Currency::DeleteCellsN(HML_ND_CELLARRAY* cells)
{
	if (cells)
	{
		if (!cells->IsMatrixShared())
		{
			delete cells;

			if (cells == data.cells_nd)
				data.cells_nd = NULL;
		}
		else
		{
			cells->DecrRefCount();
		}
	}
}

void Currency::DeleteStruct(StructData* sd)
{
	if (sd)
	{
		if (sd->GetRefCount() == 1)
		{
			delete sd;

			if (sd == data.sd)
				data.sd = NULL;
		}
		else
		{
			sd->DecrRefCount();
		}
	}
}

void Currency::DeleteFunctionInfo(FunctionInfo* fi)
{
	if (fi)
	{
		fi->DecrRefCount();

		if (fi->GetRefCount() == 0)
		{
			delete fi;

			if (fi == data.func)
				data.func = NULL;
		}
	}
}

Currency::Currency(const Currency& cur)
{
	type = TYPE_SCALAR;
	Copy(cur);
}

Currency& Currency::operator=(const Currency& in)
{
	Copy(in);
	return *this;
}

#include <algorithm>

void Currency::Swap(Currency& cur)
{
	std::swap(type, cur.type);
	std::swap(data, cur.data);
	std::swap(mask, cur.mask);
	std::swap(out_name, cur.out_name);
	std::swap(_outputType, cur._outputType);
	std::swap(_display, cur._display);
	std::swap(classname, cur.classname);
	std::swap(message, cur.message);
	std::swap(_is_utf8, cur._is_utf8);
	std::swap(_is_linear_range, cur._is_linear_range);
}

Currency::Currency(Currency&& cur): type(TYPE_MATRIX), mask(MASK_DOUBLE), out_name(NULL)
    , _display(0)
    , _outputType (OUTPUT_TYPE_DEFAULT)
{
	data.mtx = NULL; 
	Swap(cur);
}

Currency& Currency::operator=(Currency&& in)
{
	Swap(in);
	return *this;
}

Currency::~Currency()
{	
	if (type == TYPE_MATRIX)
		DeleteMatrix(data.mtx);
	else if (type == TYPE_SCALAR)
		;
	else if (type == TYPE_ND_MATRIX)
		DeleteMatrixN(data.mtxn);
	else if (type == TYPE_SPARSE)
		DeleteMatrixS(data.mtxs);
	else if (type == TYPE_COMPLEX)
        delete data.complex;
	else if ((type == TYPE_STRUCT) || (type == TYPE_OBJECT))
		DeleteStruct(data.sd);
    else if (type == TYPE_FORMAT)
        delete data.format;
	else if (type == TYPE_CELLARRAY)
		DeleteCells(data.cells);
	else if (type == TYPE_ND_CELLARRAY)
		DeleteCellsN(data.cells_nd);
	else if (type == TYPE_ERROR)
		delete message;
	else if (type == TYPE_FUNCHANDLE)
		DeleteFunctionInfo(data.func);
}

void Currency::ReplaceComplex(hwComplex new_value)
{
	if (!data.complex)
		data.complex = new hwComplex;

	data.complex->Real(new_value.Real());
	data.complex->Imag(new_value.Imag());

	type = TYPE_COMPLEX;
}
//------------------------------------------------------------------------------
// Returns a string description of the currency type
//------------------------------------------------------------------------------
std::string Currency::GetTypeString() const
{
    std::string output;
	if (IsScalar())
	{
		if (IsLogical())
			output = "logical";
		else
			output = "number";
	}
	else if (IsString())
	{
		output = "string";
	}
	else if (IsComplex())
	{
		output = "complex";
	}
	else if (IsMatrix())
	{
		int rows = 0;
		int cols = 0;
		const hwMatrix* mtx = Matrix();

		if (mtx)
		{
			rows = mtx->M();
			cols = mtx->N();
		}

		char buffer[1024];
		sprintf_s(buffer, "matrix [%d x %d]", rows, cols);
		output = buffer;
	}
	else if (IsFunctionHandle())
	{
		output = "function handle";
	}
	else if (IsCellArray())
	{
		HML_CELLARRAY* mtx = CellArray();
		char buffer[1024];
		sprintf_s(buffer, "cell array [%d x %d]", mtx->M(), mtx->N());
		output = buffer;
	}
	else if (IsNDCellArray())
	{
		output = "cell array [";

		HML_ND_CELLARRAY* cells = CellArrayND();

		if (cells)
		{
			std::vector<int> dims = cells->Dimensions();

			size_t num_dims = dims.size();

			for (size_t j = 0; j < num_dims; j++)
			{
				output += std::to_string(static_cast<long long>(dims[j]));

				if (j != num_dims - 1)
					output += " x ";
			}
		}

		output += "]";
	}
	else if (IsStruct())
	{
		StructData* mtx = Struct();
        if (!mtx)
        {
            output = "struct [ ]";
        }
        else
        {
            char buffer[1024];
            sprintf_s(buffer, "struct [%d x %d]", mtx->M(), mtx->N());
            output = buffer;
        }
	}
    else if (IsNDMatrix()) 
    {
        std::ostringstream strstream;
        strstream << "matrix [";
        const hwMatrixN* mtxN = MatrixN();
        if( mtxN ) 
        {
            const std::vector<int>& dims = mtxN->Dimensions();
            size_t size = dims.size();
            for(int i = 0; i < size; ++i) 
            {
                strstream << dims[i];
                if( i < size-1 ) strstream << " x ";
            }
        }
        strstream << "]" << std::ends;
        output = strstream.str();
    }
    else if (IsBoundObject())
    {
        return GetClassname();
    }
    else if (IsSparse())
    {
        const hwMatrixS* mtx = MatrixS();
        if (!mtx)
        {
            return "sparse [ ]";
        }
        int m    = mtx->M();
        int n    = mtx->N();
        int nnz  = mtx->NNZ();

        output = "sparse [" + std::to_string(static_cast<long long>(m)) + " x "
            + std::to_string(static_cast<long long>(n)) + "], nnz = "
            + std::to_string(static_cast<long long>(nnz));
    }
	return output;
}
//------------------------------------------------------------------------------
//! Returns true if currency is type
//------------------------------------------------------------------------------
bool Currency::IsString()    const
{
	return (mask == MASK_STRING && data.mtx);
}
//------------------------------------------------------------------------------
// Returns false if this is a multi line string
//------------------------------------------------------------------------------
bool Currency::IsMultilineString() const
{
	bool ret_val = IsString();
	if (ret_val)
	{
		hwMatrix* mtx = data.mtx;
        if (mtx->Size() == 0 || mtx->M() == 1)  // Columns can be 0 sometimes
        {
            ret_val = false;
        }
	}
	return ret_val;
}
//------------------------------------------------------------------------------
//! Returns the string value for the currency
//------------------------------------------------------------------------------
std::string Currency::StringVal() const
{
    if (!data.mtx) return "";

	std::string st;

	st.reserve(data.mtx->Size());

	int rows = data.mtx->M();
	int cols = data.mtx->N();

	if (rows > 1) 
		st += "\n";  // Add newline if there are multiple rows

	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			if ((*data.mtx)(i, j) != 0x00)
            {
				st += static_cast<unsigned char>((*data.mtx)(i,j));
            }
			else
			{
				st += " ";
			}
		}

		if (i+1 < rows)
			st += "\n";
	}

 	return st;
}

bool Currency::IsVector() const
{
	bool ret_val = false;

	if (type == TYPE_MATRIX && data.mtx && mask != MASK_STRING)
	{
		if (data.mtx->Size() != 1)
		{
			if ((data.mtx->M() == 1) || (data.mtx->N() == 1))
				ret_val = true;
		}
	}

	return ret_val;
}

bool Currency::IsRealVector() const
{
	bool ret_val = false;

	if (type == TYPE_MATRIX && data.mtx && mask != MASK_STRING)
	{
		if ((data.mtx->M() == 1) || (data.mtx->N() == 1))
			ret_val = true;

		if (!data.mtx->IsReal())
			ret_val = false;
	}

	return ret_val;
}

std::vector<double> Currency::Vector() const
{ 
	std::vector<double> v;

	if (data.mtx->M() == 1)
	{
		for (int i=0; i < data.mtx->N(); i++)
			v.push_back ((*data.mtx)(0,i));
	}
	else if (data.mtx->N() == 1)
	{
		for (int i=0; i < data.mtx->M(); i++)
			v.push_back ((*data.mtx)(i,0));
	}
	return v; 
}

void Currency::MakeStruct()
{
	if (type != TYPE_STRUCT)
	{
		if (type == TYPE_CELLARRAY)
			DeleteCells(data.cells);
		else if (type == TYPE_MATRIX)
			DeleteMatrix(data.mtx);

		type = TYPE_STRUCT;
		mask = MASK_NONE;
		data.sd = new StructData;
	}
}

void Currency::MakeCell()
{
	if (type != TYPE_CELLARRAY)
	{
		type = TYPE_CELLARRAY;
		mask = MASK_NONE;
		data.cells = EvaluatorInterface::allocateCellArray(1, 1);
	}
}

bool Currency::IsScalar() const
{
	if (type == TYPE_SCALAR)
		return true;

	if ((type == TYPE_COMPLEX) && (data.complex->Imag() == 0.0))
		return true;

	if (type == TYPE_MATRIX && mask != MASK_STRING)
	{
		if (data.mtx && (data.mtx->Size() == 1) && (data.mtx->IsReal()))
			return true;
	}

	return false;
}

double Currency::Scalar() const
{
	if (type == TYPE_SCALAR)
		return data.value;

	if ((type == TYPE_COMPLEX) && (data.complex->Imag() == 0.0))
		return data.complex->Real();

	if ((type == TYPE_MATRIX) && data.mtx && (data.mtx->Size() == 1) && (data.mtx->IsReal()))
		return (*data.mtx)(0);

	return 0.0;
}

bool Currency::IsComplex() const
{
	if (type == TYPE_COMPLEX)
		return true;

	if (type == TYPE_MATRIX && mask != MASK_STRING)
	{
		if (data.mtx && (data.mtx->Size() == 1) && (!data.mtx->IsReal()))
			return true;
	}

	return false;
}

hwComplex Currency::Complex() const
{
	if (type == TYPE_COMPLEX)
		return *data.complex;

	if ((type == TYPE_MATRIX) && data.mtx && (data.mtx->Size() == 1) && (!data.mtx->IsReal()))
		return data.mtx->z(0);

	return hwComplex(0.0, 0.0);
}

bool Currency::IsInteger() const
{
	if (IsScalar())
	{
		double    temp  = Scalar();
		long long temp2;

		if (temp < 0)
			temp2 = ((long long)(temp - 1.0e-15));
		else
			temp2 = ((long long)(temp + 1.0e-15));

		double delta = temp - temp2;

    	if (abs(delta) < 1.0e-15)
			return true;
	}

	return false;
}

bool Currency::IsPositiveInteger() const
{
	if (IsScalar())
	{
		double temp = Scalar();

		if (temp > 0.0)
		{
			if (abs(temp - (int(temp+1.0e-15)) < 1.0e-15))
				return true;
		}
	}

	return false;
}

bool Currency::IsPositiveInteger64() const
{
	if (IsScalar())
	{
		double value1 = Scalar();

		if (value1 > 0.0)
		{
			int64_t value2 = static_cast<int64_t> (value1);

			if (value1 == static_cast<double> (value2))
				return true;
		}
	}

	return false;
}

bool Currency::IsPositiveIntegralVector() const
{
	if (IsLogical())
		return false;

	if (IsMatrix() && !IsString())
	{
		const hwMatrix* mtx = Matrix();

		if (!mtx->IsReal())
			return false;

		if (mtx->IsEmpty())
			return false;

		if (mtx->Size() == 1) // treat as a scalar instead
			return false;

		if ((mtx->M() != 1) && (mtx->N() != 1))
			return false;

		if (IsLinearRange())
		{
			double val = (*mtx)(0);

			if (val >= 1.0)
				return true;
			else
				return false;
		}

        for (int k=0; k<mtx->N(); k++)
        {
            for (int j=0; j<mtx->M(); j++)
            {
				double val = (*mtx)(j,k);

				if (val <= 0.0)
					return false;

				// making sure it's an integer (or close enough)
				if (abs(val - (int(val+1.0e-15)) > 1.0e-15))
					return false;				
			}
		}

		return true;
	}

	return false;
}

bool Currency::IsPositiveIntegralMatrix() const
{
	if (IsLogical())
		return false;

	if (IsMatrix() && !IsString())
	{
		const hwMatrix* mtx = Matrix();

		if (!mtx->IsReal())
			return false;

		if (mtx->IsEmpty())
			return false;

		for (int k = 0; k<mtx->Size(); k++)
		{
			double val = (*mtx)(k);

			// making sure it's an integer (or close enough)
			if (abs(val - (int(val + 1.0e-15)) > 1.0e-15))
				return false;
		}

		return true;
	}

	return false;
}

bool Currency::IsEmpty() const
{
	if (IsMatrixOrString())
	{
		if (data.mtx)
			return data.mtx->IsEmpty();
		else
			return true;
	}

	return false;
}

bool  Currency::IsMatrix() const 
{ 
	if (mask == MASK_STRING)
		return false;
	else
		return type == TYPE_MATRIX; 
}

bool  Currency::IsNDMatrix() const 
{ 
	// ND matrices can't be strings so don't bother checking
	return type == TYPE_ND_MATRIX; 
}

bool  Currency::IsSparse() const
{
	// Sparse matrices can't be strings so don't bother checking
	return type == TYPE_SPARSE;
}

bool  Currency::IsLinearRange() const
{
	if (type == TYPE_MATRIX)
		return _is_linear_range;

	return false;
}

void  Currency::SetLinearRange(bool range)
{
	if (type == TYPE_MATRIX)
		_is_linear_range = range;
}

bool Currency::IsCharacter() const 
{
	if (IsString())
	{
		if (data.mtx->Size() == 1)
			return true;
	}

	return false;
}

const hwMatrix* Currency::Matrix() const
{
	if (type == TYPE_MATRIX && !data.mtx)
		data.mtx = new hwMatrix();

	return data.mtx;
}

const hwMatrixN* Currency::MatrixN() const
{
	return data.mtxn;
}

const hwMatrixS* Currency::MatrixS() const
{
	return data.mtxs;
}

hwMatrix* Currency::GetWritableMatrix() 
{
	if (!data.mtx)
		data.mtx = ExprTreeEvaluator::allocateMatrix();

	_is_linear_range = false; // just to be certain

	return data.mtx;
}

hwMatrixN* Currency::GetWritableMatrixN() 
{
	if (!data.mtxn)
		data.mtxn = ExprTreeEvaluator::allocateMatrixN();

	return data.mtxn;
}

hwMatrixS* Currency::GetWritableMatrixS()
{
	if (!data.mtxs)
		data.mtxs = ExprTreeEvaluator::allocateMatrixS();

	return data.mtxs;
}

void Currency::ReplaceMatrix(hwMatrix* new_mtx)
{
	if (type == TYPE_MATRIX)
	{
		if (new_mtx != data.mtx)
		{
			DeleteMatrix(data.mtx);
			data.mtx = new_mtx;
		}
	}
	else if (type == TYPE_SCALAR)
	{
		type = TYPE_MATRIX;
		data.mtx = new_mtx;
	}
	else if (type == TYPE_COMPLEX)
	{
		type = TYPE_MATRIX;
		data.mtx = new_mtx;
	}
}

void Currency::ReplaceCellArray(HML_CELLARRAY* new_cells)
{
	if (type == TYPE_CELLARRAY)
	{
		if (new_cells != data.cells)
		{
			DeleteCells(data.cells);
			data.cells = new_cells;
		}
	}
}

void Currency::ReplaceStruct(StructData* new_sd)
{
	if ((type == TYPE_STRUCT) || (type == TYPE_OBJECT))
	{
		if (new_sd != data.sd)
		{
			DeleteStruct(data.sd);
			data.sd = new_sd;
		}
	}
}

bool Currency::IsLogical() const
{
	if (mask == MASK_LOGICAL)
		return true;
	else
		return false;
}

bool Currency::IsSyntaxError() const
{
	if (IsError() && data.value == -999.0)
		return true;

	return false;
}

bool Currency::IsCellList() const
{
    if (IsCellArray() && (mask == MASK_CELL_LIST))
		return true;

	return false;
}

const hwMatrix* Currency::ConvertToMatrix() const
{
	if (type == TYPE_SCALAR)
	{
		double old_val = data.value;
		data.mtx = ExprTreeEvaluator::allocateMatrix(1,1, hwMatrix::REAL);
		(*data.mtx)(0) = old_val;
	}
	else if (type == TYPE_COMPLEX)
	{
		hwComplex* old_complex = data.complex;
		data.mtx = ExprTreeEvaluator::allocateMatrix(1,1, hwMatrix::COMPLEX);
		data.mtx->z(0) = *old_complex;
	}
	else if (type == TYPE_MATRIX)
	{
		if (!data.mtx)
			data.mtx = ExprTreeEvaluator::allocateMatrix(0, 0, hwMatrix::REAL);
	}
	else if (IsCellList())
	{
		HML_CELLARRAY* cells   = CellArray();
		bool           is_real = true;

		int size = cells->Size();

		int rows_needed = 0;
		int cols_needed = 0;

		for (int j=0; j<size; j++)
		{
			Currency temp = (*cells)(j);

			if (temp.IsScalar())
			{
				cols_needed++;
			}
			else if (temp.IsComplex())
			{
				cols_needed++;
			}
			else if (temp.IsString())
			{
				const hwMatrix* string = temp.Matrix();
				cols_needed += string->N();

				mask = MASK_STRING;
			}
			else if (temp.IsMatrix())
			{
				const hwMatrix* inner_mtx = temp.Matrix();
				if (rows_needed == 0)
					rows_needed = inner_mtx->M();
				else if (rows_needed != inner_mtx->M())
					rows_needed = -1;

				cols_needed += inner_mtx->N();
			}
			else
			{
				mask = MASK_CELL_LIST;
				return NULL;
			}
		}

		if (!rows_needed)
			rows_needed = 1;

		hwMatrix* mtx = ExprTreeEvaluator::allocateMatrix(rows_needed, cols_needed, hwMatrix::REAL);
		int count = 0;

		int	row_offset = 0;
		int col_offset = 0;

		for (int j=0; j<size; j++)
		{
			Currency temp = (*cells)(j);

			if (temp.IsScalar())
			{
				if (mtx->IsReal())
					(*mtx)(count) = temp.Scalar();
				else
					mtx->z(count) = temp.Scalar();
				count++;
			}
			else if (temp.IsComplex())
			{
				mtx->MakeComplex();
				mtx->z(count) = temp.Complex();
				count++;
			}
			else if (temp.IsString())
			{
				const hwMatrix* string = temp.Matrix();
				
				for (int k=0; k<string->N(); k++)
				{
					(*mtx)(count) = (*string)(k);
					count++;
				}
			}
			else if (temp.IsMatrix() && (rows_needed != -1))
			{
				const hwMatrix* inner_mtx = temp.Matrix();

				mtx->WriteSubmatrix(row_offset, col_offset, *inner_mtx);

				col_offset += inner_mtx->N();
			}
		}

		if (!cells->IsMatrixShared())
			delete cells;
		else
			cells->DecrRefCount();

		data.mtx = mtx;
	}
	else
	{
		data.mtx = NULL;
	}

	type = TYPE_MATRIX;
	return data.mtx;
}

void Currency::ExpandMatrixPair(const hwMatrix& src1, const hwMatrix& src2,
								hwMatrix*& res1, hwMatrix*& res2)
{
	if (res1 != nullptr || res2 != nullptr)
		return;

	int m1 = src1.M();
	int n1 = src1.N();
	int m2 = src2.M();
	int n2 = src2.N();

	// reuse m and n values for the replication numbers
	if (m1 == m2)
	{
		m1 = 1;
		m2 = 1;
	}
	else if (m1 == 1)
	{
		m1 = m2;
		m2 = 1;
	}
	else if (m2 == 1)
	{
		m2 = m1;
		m1 = 1;
	}
	else
	{
		// trigger an error
		m1 = 1;
		m2 = 1;
	}

	if (n1 == n2)
	{
		n1 = 1;
		n2 = 1;
	}
	else if (n1 == 1)
	{
		n1 = n2;
		n2 = 1;
	}
	else if (n2 == 1)
	{
		n2 = n1;
		n1 = 1;
	}
	else
	{
		// trigger an error
		n1 = 1;
		n2 = 1;
	}

	res1 = ExprTreeEvaluator::allocateMatrix();
	res2 = ExprTreeEvaluator::allocateMatrix();

	hwMathStatus status;

	status = res1->Repmat(src1, m1, n1);

	if (!status.IsOk())
		throw hwMathException(status.GetMsgCode());

	status = res2->Repmat(src2, m2, n2);

	if (!status.IsOk())
		throw hwMathException(status.GetMsgCode());
}

void Currency::ExpandMatrixPair(const hwMatrixN& src1, const hwMatrixN& src2,
								hwMatrixN*& res1, hwMatrixN*& res2)
{
	if (res1 != nullptr || res2 != nullptr)
		return;

	const std::vector<int>& dims_src1 = src1.Dimensions();
	const std::vector<int>& dims_src2 = src2.Dimensions();
	std::vector<int> dims_rep1;
	std::vector<int> dims_rep2;

	size_t minDim = _min(dims_src1.size(), dims_src2.size());

	for (size_t i = 0; i < minDim; ++i)
	{
	    if (dims_src1[i] == dims_src2[i])
	    {
	    	dims_rep1.push_back(1);
	    	dims_rep2.push_back(1);
	    }
		else if (dims_src1[i] == 1)
		{
			dims_rep1.push_back(dims_src2[i]);
			dims_rep2.push_back(1);
		}
		else if (dims_src2[i] == 1)
		{
			dims_rep1.push_back(1);
			dims_rep2.push_back(dims_src1[i]);
		}
		else
		{
			// trigger an error
			dims_rep1.push_back(1);
			dims_rep2.push_back(1);
		}
	}

	for (size_t i = minDim; i < dims_src1.size(); ++i)
	{
		dims_rep1.push_back(1);
	}

	for (size_t i = minDim; i < dims_src2.size(); ++i)
	{
		dims_rep2.push_back(1);
	}

	res1 = ExprTreeEvaluator::allocateMatrixN();
	res2 = ExprTreeEvaluator::allocateMatrixN();

	res1->Repmat(src1, dims_rep1);
	res2->Repmat(src2, dims_rep2);
}

void Currency::ExpandMatrixPair(const hwMatrix& src1, const hwMatrixN& src2,
							    hwMatrixN*& res1, hwMatrixN*& res2)
{
	if (res1 != nullptr || res2 != nullptr)
		return;

	int m = src1.M();
	int n = src1.N();
	const std::vector<int>& dims_src2 = src2.Dimensions();
	std::vector<int> dims_rep1;
	std::vector<int> dims_rep2;

	if (m == dims_src2[0])
	{
		dims_rep1.push_back(1);
		dims_rep2.push_back(1);
	}
	else if (m == 1)
	{
		dims_rep1.push_back(dims_src2[0]);
		dims_rep2.push_back(1);
	}
	else if (dims_src2[0] == 1)
	{
		dims_rep1.push_back(1);
		dims_rep2.push_back(m);
	}
	else
	{
		// trigger an error
		dims_rep1.push_back(1);
		dims_rep2.push_back(1);
	}

	if (n == 1)
	{
		dims_rep1.push_back(dims_src2[1]);
		dims_rep2.push_back(1);
	}
	else if (dims_src2[1] == 1)
	{
		dims_rep1.push_back(1);
		dims_rep2.push_back(n);
	}
	else if (n == dims_src2[1])
	{
		dims_rep1.push_back(1);
		dims_rep2.push_back(1);
	}
	else
	{
		// trigger an error
		dims_rep1.push_back(1);
		dims_rep2.push_back(1);
	}

	for (size_t i = 2; i < dims_src2.size(); ++i)
	{
		dims_rep1.push_back(dims_src2[i]);
		dims_rep2.push_back(1);
	}

	res1 = ExprTreeEvaluator::allocateMatrixN();
	res2 = ExprTreeEvaluator::allocateMatrixN();

	res1->Repmat(src1, dims_rep1);
	res2->Repmat(src2, dims_rep2);
}

HML_CELLARRAY* Currency::ConvertToCellArray()
{
	bool ret_val = true;

	if (type == TYPE_SCALAR)
	{
		double old_val = data.value;
		data.cells = ExprTreeEvaluator::allocateCellArray(1,1);
		(*data.cells)(0) = old_val;
	}
	else if (type == TYPE_COMPLEX)
	{
		hwComplex* old_complex = data.complex;
		data.cells = ExprTreeEvaluator::allocateCellArray(1,1);
		(*data.cells)(0) = *old_complex;
	}
	else if (type == TYPE_MATRIX)
	{
		hwMatrix* old_mtx = data.mtx;
		data.cells = ExprTreeEvaluator::allocateCellArray(1,1);
		(*data.cells)(0) = old_mtx;

		if (mask == MASK_STRING)
		{
			(*data.cells)(0).SetMask(MASK_STRING);
			mask = MASK_NONE;
		}
	}
	else if (type == TYPE_CELLARRAY)
	{
	}
	else
	{
		data.cells = NULL;
	}

	type = TYPE_CELLARRAY;

	return data.cells;
}

void Currency::FlattenCellList()
{
	if (IsCellList())
	{
		HML_CELLARRAY* new_cells = NULL;
		HML_CELLARRAY* old_cells = CellArray();

		int new_size = 0;
		
		for (int j = 0; j < old_cells->Size(); ++j)
		{
			Currency loc_temp = (*old_cells)(j);

			if (loc_temp.IsCellArray())
			{
				HML_CELLARRAY* loc_cells = loc_temp.CellArray();

				if (loc_cells->M() == 1)
					new_size += loc_cells->N();
				else
					return; // failed conversion
			}
			else if (loc_temp.IsScalar() || loc_temp.IsString() || loc_temp.IsMatrix())
			{
				new_size += 1;
			}
			else if (loc_temp.IsStruct())
			{
				StructData* sd = loc_temp.Struct();
				
				if (sd->Size() == 1)
					new_size += 1;
				else
					return;
			}
			else
			{
				return;
			}
		}

		new_cells        = new HML_CELLARRAY(1, new_size, HML_CELLARRAY::REAL);
		int new_cell_cnt = 0;

		for (int j = 0; j < old_cells->Size(); ++j)
		{
			Currency loc_temp = (*old_cells)(j);

			if (loc_temp.IsCellArray())
			{
				HML_CELLARRAY* loc_cells = loc_temp.CellArray();

				for (int k = 0; k < loc_cells->Size(); ++k)
				{
					(*new_cells)(new_cell_cnt) = (*loc_cells)(k);
					++new_cell_cnt;
				}
			}
			else if (loc_temp.IsScalar() || loc_temp.IsString() || loc_temp.IsMatrix() || loc_temp.IsStruct())
			{
				(*new_cells)(new_cell_cnt) = loc_temp;
				++new_cell_cnt;
			}
		}

		ReplaceCellArray(new_cells);
	}
}

void Currency::ConvertToStruct()
{
	if (IsStruct())
	{
		return;
	}
	else if (IsCellArray())
	{
		HML_CELLARRAY* cells = CellArray();

		if (cells->Size() == 0)
		{
			delete cells;
			data.sd = new StructData;
			type = TYPE_STRUCT;
		}
		else if (IsCellList())
		{
			bool valid_switch = true;
			for (int j = 0; j < cells->Size(); j++)
			{
				if ((*cells)(j).IsStruct())
				{
				}
				else
				{
					valid_switch = false;
				}
			}

			if (valid_switch)
			{
				type = TYPE_STRUCT;
				StructData* temp_sd = new StructData;
				temp_sd->DimensionNew(1, cells->Size());

				for (int j = 0; j < cells->Size(); j++)
					temp_sd->SetElement(j+1, -1, (*cells)(j).Struct());

				delete cells;
				data.sd = temp_sd;
			}
		}
	}
	else if (IsMatrix())
	{
		hwMatrix* mtx = data.mtx;

	    if (mtx && mtx->Size() == 0)
			delete mtx;

		data.sd = new StructData;
		type = TYPE_STRUCT;
	}
}

OMLTree* Currency::ConvertToTree() const
{
	if (IsTree())
	{
		// nothing to be done
	}
	else if (IsString())
	{
		OMLTree* tree = new OMLTree(IDENT, StringVal(), NULL, 0, 0);
		
		data.mtx->DecrRefCount();

		if (data.mtx->GetRefCount() == 0)
			delete data.mtx;

		data.tree = tree;

		type = TYPE_TREE;
	}
	else if (IsScalar())
	{
		char buffer[256];

		sprintf(buffer, "%g", data.value);
		
		OMLTree* tree = new OMLTree(NUMBER, buffer, NULL, 0, 0);

		data.tree = tree;

		type = TYPE_TREE;
	}
	else
	{
		return NULL;
	}

	return Tree();
}

void Currency::SetOutputName(const std::string& name) const
{
	out_name = vm.GetStringPointer(name);
}

void Currency::SetOutputName(const std::string* name) const
{
	out_name = name;
}

void Currency::ClearOutputName()
{
   static const std::string* empty_str = NULL;
   
   if (!empty_str)
	   empty_str = vm.GetStringPointer("");

   out_name = empty_str;
}

const std::string* Currency::GetOutputNamePtr() const
{
	if (out_name)
		return out_name;
	else
		return vm.GetStringPointer("");
}

//------------------------------------------------------------------------------
// Gets output for printing
//------------------------------------------------------------------------------
std::string Currency::GetOutputString(const OutputFormat* fmt) const
{
	std::ostringstream os;

	const std::string* output_name = GetOutputNamePtr();

    if (!output_name->empty() && !IsDispOutput() && !IsError() && !IsPrintfOutput())
			os << *output_name << " = ";

	if (IsScalar())
	{
        return CurrencyDisplay::ScalarToString(fmt, Scalar(), os);
	}
	else if (IsComplex())
	{
        return CurrencyDisplay::ComplexToString(fmt, Complex(), os);
	}
	else if (IsString())
	{
        std::string val (StringVal());
        if (CurrencyDisplay::CanPaginate(val))
        {
            return (GetDisplay()->GetOutput(fmt, os));
        }
        else
        {
		    os << val;
        }
	}
	else if (IsMatrix() || IsCellArray() || IsStruct() || IsNDMatrix() ||
             IsObject() || IsNDCellArray() || IsSparse() || IsCellList())
	{
        return (GetDisplay()->GetOutput(fmt, os));
	}
    else if (IsFunctionHandle())
	{
		// This should never happen, but there's some weird UI cases that actually do trigger it
		if (data.func)
		{
			if (data.func->IsAnonymous())
			{
				std::vector<const std::string*> params = data.func->Parameters();

				os << "@(";

				for (int j = 0; j < params.size(); ++j)
				{
					os << *params[j];

					if (j != params.size() - 1)
						os << ", ";
				}

				os << ") " << data.func->Statements()->GetStringRepresentation();
			}
			else
			{
				os << "@" << data.func->FunctionName();
			}
		}
	}
	else if (IsError())
	{
		os << *message;
	}
    else if (IsBoundObject())
    {
        os << GetClassname();
    }
	else if (IsPointer())
	{
		os << "<pointer " << Pointer() << ">";
	}
    std::string output (os.str());
	return output;
}

std::string Currency::GetClassname() const
{
	if (IsObject() || IsBoundObject())
		return *classname;
	else if (IsPointer())
		return data.cur_ptr->GetClassname();
	else
		return "";
}


const std::string* StringManager::GetStringPointer(const std::string& var)
{
	StringStorage::iterator iter;

	iter = _strings.find(&var);

	if (iter == _strings.end())
	{
		std::string* new_str = new std::string(var);
		_strings.insert(new_str);
		return new_str;
	}
	else
	{
		return (*iter);
	}
}
//------------------------------------------------------------------------------
//! Creates, if needed and returns display
//------------------------------------------------------------------------------
CurrencyDisplay* Currency::GetDisplay() const
{
    if (_display)
    {
        return _display;
    }

    if (IsCellArray())
    {
        _display = new CellDisplay(*this);
    }
    else if (IsMatrix())
    {
        _display = new MatrixDisplay(*this);
    }
    else if (IsStruct() || IsObject())
    {
        _display = new StructDisplay(*this);
    }
    else if (IsNDMatrix())
    {
        _display = new MatrixNDisplay(*this);
    }
    else if (IsString())
    {
        std::string val (StringVal());
        if (CurrencyDisplay::CanPaginate(val))
        {
            _display = new StringDisplay(*this);
        }
    }
    else if (IsNDCellArray())
    {
        _display = new CellNDisplay(*this);
    }
    else if (IsSparse())
    {
        _display = new SparseDisplay(*this);
    }
    return _display;
}
//------------------------------------------------------------------------------
//! Deletes display
//------------------------------------------------------------------------------
void Currency::DeleteDisplay()
{
    delete _display;
    _display = 0;
}

void Currency::SetClass(const std::string& class_name)
{
	type        = TYPE_OBJECT;
	mask        = MASK_NONE;
	_outputType = OUTPUT_TYPE_DEFAULT;

	classname = (std::string*)vm.GetStringPointer(class_name);
}
//------------------------------------------------------------------------------
//! Utility to get just the values, without any header(s), if applicable
//------------------------------------------------------------------------------
std::string Currency::GetValues(const OutputFormat* fmt) const
{    
    if (IsScalar() || IsComplex())
    {
        return GetOutputString(fmt);
    }

    bool stringpaginate = (IsString() && CurrencyDisplay::CanPaginate(StringVal()));

	if (IsMatrix() || IsCellArray() || IsStruct() || IsNDMatrix() || 
        IsObject() || stringpaginate || IsSparse())
    {
        std::string out = GetDisplay()->GetValues(fmt);
        const_cast<Currency*>(this)->DeleteDisplay();
        return out;
    }

    return GetOutputString(fmt);
} 
//------------------------------------------------------------------------------
// Utility to get just the values without any headers or currency name
//------------------------------------------------------------------------------
std::string Currency::GetValuesForDisplay(const OutputFormat* fmt) const
{    
    OUTPUT_TYPE origtype = _outputType;

    const_cast<Currency*>(this)->_outputType = OUTPUT_TYPE_DISP;
    std::string out = GetOutputString(fmt);
    const_cast<Currency*>(this)->_outputType = origtype;

    return out;
} 

bool Currency::IgnoreCoW(void* ptr)
{
    if (!Currency::ignore_cow_pointers.empty())
    {
        if (Currency::ignore_cow_pointers.find(ptr) != Currency::ignore_cow_pointers.end())
            return true;
    }

	return false;
}

void Currency::CheckForUTF8Characters()
{
	if (!IsString())
		return;

	int length = data.mtx->Size();

	for (int i = 0; i < length; i++)
	{
		if ((*data.mtx)(i) >= 192)
			_is_utf8 = true;
	}
}

void Currency::SetMask(int new_mask) 
{ 
	mask = MaskType(new_mask);

	if (mask == MASK_STRING)
		CheckForUTF8Characters();
}

StringManager::~StringManager()
{
	StringStorage::iterator iter;

	for (iter = _strings.begin(); iter != _strings.end(); iter++)
		delete *iter;
}
 //------------------------------------------------------------------------------
 // Update output type to save to output log
 //------------------------------------------------------------------------------
 void Currency::SetOutputTypeToLog()
 {
     if (!IsLogOutput())
     {
         int val = _outputType | OUTPUT_TYPE_LOG;
         _outputType = static_cast<OUTPUT_TYPE>(val);
     }
 }
 //------------------------------------------------------------------------------
 // Sets output type to be from 'disp' command
 //------------------------------------------------------------------------------
 void Currency::UpdateOutputType(OUTPUT_TYPE type)
 {
     if (!IsLogOutput())
     {
         _outputType = type;
     }
     else
     {
         int val = type | OUTPUT_TYPE_LOG;
         _outputType = static_cast<OUTPUT_TYPE>(val);
     }
 }
//------------------------------------------------------------------------------
// Returns true if this output is the same as the given type
//------------------------------------------------------------------------------
 bool Currency::IsOutputType(OUTPUT_TYPE ref) const
 {
     if (_outputType == ref)
     {
         return true;
     }
     return (((int)_outputType & ref) != 0);  // Check for combination of outputs
 }
//------------------------------------------------------------------------------
// Returns true if currency is type
//------------------------------------------------------------------------------
 bool Currency::IsSingle() const
 {
     return (mask == MASK_SINGLE &&
            (IsScalar() || IsComplex() || type == TYPE_MATRIX));
 }
