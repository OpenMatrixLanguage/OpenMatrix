/**
* @file OMLInterface.cpp
* @date January 2017
* Copyright (C) 2017-2018 Altair Engineering, Inc.  
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

#include "OMLInterface.h"

#include "Currency.h"
#include "StructData.h"
#include "FunctionInfo.h"
#include "hwMatrixN.h"
#include "OML_Error.h"

void OMLInterfaceImpl::RegisterFunction(const char* name, ALT_FUNCPTR fp)
{
	_eval->RegisterBuiltInFunction(name, fp);
}

void OMLInterfaceImpl::ThrowError(const char* message)
{
	throw OML_Error(message);
}

int OMLInterfaceImpl::Nargout() const
{
	return _eval->GetNargoutValue();
}

OMLCurrencyList* OMLInterfaceImpl::CreateCurrencyList()
{
	return new OMLCurrencyListImpl();
}

const OMLCurrency* OMLInterfaceImpl::CallFunction(const OMLFunctionHandle* handle, OMLCurrencyList* inputs)
{
	std::vector<Currency> ins;

	for (int j=0; j<inputs->Size(); j++)
	{
		const OMLCurrency* temp = inputs->Get(j);
		OMLCurrencyImpl* impl = (OMLCurrencyImpl*)temp;

		ins.push_back(impl->GetCurrency());
	}

	OMLFunctionHandleImpl* fh = (OMLFunctionHandleImpl*)handle;
	Currency result = _eval->CallInternalFunction(fh->GetFunctionInfo(), ins);
	return new OMLCurrencyImpl(result);
}

bool OMLCurrencyImpl::IsScalar() const
{
	return _cur.IsScalar();
}

bool OMLCurrencyImpl::IsComplex() const
{
	return _cur.IsComplex();
}

bool OMLCurrencyImpl::IsString() const
{
	return _cur.IsString();
}

bool OMLCurrencyImpl::IsMatrix() const
{
	return _cur.IsMatrix();
}

bool OMLCurrencyImpl::IsNDMatrix() const
{
	return _cur.IsNDMatrix();
}

bool OMLCurrencyImpl::IsCellArray() const
{
	return _cur.IsCellArray();
}

bool OMLCurrencyImpl::IsStruct() const
{
	return _cur.IsStruct();
}

bool OMLCurrencyImpl::IsFunctionHandle() const
{
	return _cur.IsFunctionHandle();
}

double OMLCurrencyImpl::GetScalar() const
{
	return _cur.Scalar();
}

const char* OMLCurrencyImpl::GetString() const
{
	std::string temp = _cur.StringVal();
	size_t      temp_len = temp.size();

	char* ret = new char [temp_len+1];
	strcpy(ret, temp.c_str());
	ret[temp_len] = 0;
	return ret;
}

const OMLCellArray* OMLCurrencyImpl::GetCellArray() const
{
	return new OMLCellArrayImpl(_cur.CellArray());
}

const OMLStruct* OMLCurrencyImpl::GetStruct() const
{
	return new OMLStructImpl(_cur.Struct());
}

const OMLComplex* OMLCurrencyImpl::GetComplex() const
{
	hwComplex cplx = _cur.Complex();
	return new OMLComplexImpl(cplx.Real(), cplx.Imag());
}

const OMLMatrix* OMLCurrencyImpl::GetMatrix() const
{
	return new OMLMatrixImpl(_cur.Matrix());
}

const OMLNDMatrix* OMLCurrencyImpl::GetNDMatrix() const
{
	return new OMLNDMatrixImpl(_cur.MatrixN());
}

const OMLFunctionHandle* OMLCurrencyImpl::GetFunctionHandle() const
{
	return new OMLFunctionHandleImpl(_cur.FunctionHandle());
}

bool OMLMatrixImpl::IsReal() const
{
	return _mtx->IsReal();
}

int OMLMatrixImpl::GetRows() const
{
	return _mtx->M();
}

int OMLMatrixImpl::GetCols() const
{
	return _mtx->N();
}

const double* OMLMatrixImpl::GetRealData() const
{
	if (_mtx->IsReal())
		return _mtx->GetRealData();
	
	hwMatrix temp;
	_mtx->UnpackComplex(&temp, NULL);
	return temp.GetRealData();
}

const double* OMLMatrixImpl::GetImaginaryData() const
{
	if (_mtx->IsReal())
		return NULL;

	hwMatrix temp;
	_mtx->UnpackComplex(NULL, &temp);
	return temp.GetRealData();
}

OMLCurrency* OMLMatrixImpl::GetCurrency() const
{
	return (new OMLCurrencyImpl(_mtx));
}

hwMatrix* OMLMatrixImpl::GetMatrixPointer() const
{
	return _mtx;
}

bool OMLNDMatrixImpl::IsReal() const
{
	return _mtx->IsReal();
}

int OMLNDMatrixImpl::GetNumDimension() const
{
	std::vector<int> dims = _mtx->Dimensions();
	return (int)dims.size();
}

int OMLNDMatrixImpl::GetDimension(int num) const
{
	std::vector<int> dims = _mtx->Dimensions();
	
	return dims[num];
}

const double* OMLNDMatrixImpl::GetRealData() const
{
	return &((*_mtx)(0));
}

const double* OMLNDMatrixImpl::GetImaginaryData() const
{
	return NULL;
}

OMLCurrency* OMLNDMatrixImpl::GetCurrency() const
{
	return (new OMLCurrencyImpl(_mtx));
}

hwMatrixN* OMLNDMatrixImpl::GetMatrixPointer() const
{
	return _mtx;
}

OMLCurrency* OMLCellArrayImpl::GetValue(int index1) const
{
	return new OMLCurrencyImpl((*_cells)(index1));
}

OMLCurrency* OMLCellArrayImpl::GetValue(int index1, int index2) const
{
	return new OMLCurrencyImpl((*_cells)(index1, index2));
}

void OMLCellArrayImpl::SetValue(int index1, OMLCurrency* val)
{
	OMLCurrencyImpl* ci = (OMLCurrencyImpl*)val;
	(*_cells)(index1)   = ci->GetCurrency();
}

void OMLCellArrayImpl::SetValue(int index1, int index2, OMLCurrency* val)
{
	OMLCurrencyImpl* ci = (OMLCurrencyImpl*)val;
	(*_cells)(index1, index2)   = ci->GetCurrency();
}

int OMLCellArrayImpl::GetRows() const
{
	return _cells->M();
}

int OMLCellArrayImpl::GetCols() const
{
	return _cells->N();
}

OMLCurrency* OMLCellArrayImpl::GetCurrency() const
{
	return (new OMLCurrencyImpl(_cells));
}

HML_CELLARRAY* OMLCellArrayImpl::GetCells() const
{
	return _cells;
}

int OMLCurrencyListImpl::Size() const
{
	return _count;
}

const OMLCurrency* OMLCurrencyListImpl::Get(int idx) const
{
	if (idx < 0)
		return NULL;

	if (idx >= _count)
		return NULL;

	return _list[idx];
}

void OMLCurrencyListImpl::Expand()
{
	OMLCurrency** new_list = new OMLCurrency* [_count+1];

	if (_list)
	{
		memcpy(new_list, _list, sizeof(OMLCurrency*)*_count);
		delete [] _list;
	}
	
	_list = new_list;
	_count++;
}

void OMLCurrencyListImpl::AddScalar(double val)
{
	Expand();
	_list[_count-1] = new OMLCurrencyImpl(val);
}

void OMLCurrencyListImpl::AddString(const char* str)
{
	Expand();
	_list[_count-1] = new OMLCurrencyImpl(str);
}

void OMLCurrencyListImpl::AddMatrix(OMLMatrix* mtx)
{
	Expand();
	_list[_count-1] = mtx->GetCurrency();
}

void OMLCurrencyListImpl::AddMatrix(const hwMatrix* mtx)
{
	Expand();
	_list[_count-1] = new OMLCurrencyImpl((hwMatrix*)mtx);
}

void OMLCurrencyListImpl::AddNDMatrix(OMLNDMatrix* mtx)
{
	Expand();
	_list[_count-1] = mtx->GetCurrency();
}

void OMLCurrencyListImpl::AddNDMatrix(const hwMatrixN* mtx)
{
	Expand();
	_list[_count-1] = new OMLCurrencyImpl((hwMatrixN*)mtx);
}

void OMLCurrencyListImpl::AddComplex(OMLComplex* cplx)
{
	Expand();
	_list[_count-1] = cplx->GetCurrency();
}

void OMLCurrencyListImpl::AddCellArray(HML_CELLARRAY* cells)
{
	Expand();
	_list[_count-1] = new OMLCurrencyImpl(cells);
}

void OMLCurrencyListImpl::AddCellArray(OMLCellArray* array)
{
	Expand();
	_list[_count-1] = array->GetCurrency();
}

void OMLCurrencyListImpl::AddStruct(OMLStruct* in_struct)
{
	Expand();
	_list[_count-1] = in_struct->GetCurrency();
}

void OMLCurrencyListImpl::AddStruct(StructData* in_sd)
{
	Expand();
	_list[_count-1] = new OMLCurrencyImpl(new StructData(*in_sd));
}

void OMLCurrencyListImpl::AddFunctionHandle(FunctionInfo* fi)
{
	Expand();
	_list[_count-1] = new OMLCurrencyImpl(new FunctionInfo(*fi));
}

OMLCellArray* OMLCurrencyListImpl::CreateCellArray(int rows, int cols)
{
	HML_CELLARRAY* cells = new HML_CELLARRAY(rows, cols, HML_CELLARRAY::REAL);
	return new OMLCellArrayImpl(cells);
}

OMLStruct* OMLCurrencyListImpl::CreateStruct(int rows, int cols)
{
	StructData* sd = new StructData();
	sd->DimensionNew(rows, cols);
	return new OMLStructImpl(sd);
}

OMLMatrix* OMLCurrencyListImpl::CreateMatrix(int rows, int cols, double* data)
{
	hwMatrix* mtx = new hwMatrix(rows, cols, data, hwMatrix::REAL);
	return new OMLMatrixImpl(mtx);
}

OMLNDMatrix* OMLCurrencyListImpl::CreateNDMatrix(int num_dims, int* dims, double* data)
{
	std::vector<int> dim_vec;

	for (int j=0; j<num_dims; j++)
		dim_vec.push_back(dims[j]);

	hwMatrixN* mtx = new hwMatrixN(dim_vec, data, hwMatrixN::REAL);
	return new OMLNDMatrixImpl(mtx);
}

OMLNDMatrix* OMLCurrencyListImpl::CreateNDMatrix(int num_dims, int* dims, double* real, double* imag)
{
	std::vector<int> dim_vec;

	for (int j=0; j<num_dims; j++)
		dim_vec.push_back(dims[j]);

	// still need to handle imaginary data
	hwMatrixN* mtx = new hwMatrixN(dim_vec, real, hwMatrixN::REAL);
	return new OMLNDMatrixImpl(mtx);
}

OMLMatrix* OMLCurrencyListImpl::CreateMatrix(int rows, int cols, double* real, double* imag)
{
	// need to pass in the imaginary data too somehow
	hwMatrix* mtx = new hwMatrix(rows, cols, real, hwMatrix::COMPLEX);
	return new OMLMatrixImpl(mtx);
}

OMLCurrency* OMLCurrencyListImpl::CreateCurrencyFromDouble(double dbl)
{
	return new OMLCurrencyImpl(dbl);
}

OMLCurrency* OMLCurrencyListImpl::CreateCurrencyFromString(const char* str)
{
	return new OMLCurrencyImpl(str);
}

OMLComplex* OMLCurrencyListImpl::CreateComplex(double real, double imag)
{
	return NULL;
}

double OMLComplexImpl::GetReal() const
{
	return cplx.Real();
}

double OMLComplexImpl::GetImag() const
{
	return cplx.Imag();
}

OMLCurrency* OMLComplexImpl::GetCurrency() const
{
	return (new OMLCurrencyImpl(cplx));
}

OMLCurrency* OMLStructImpl::GetValue(int index1, const char* field) const
{
	Currency cur = _sd->GetValue(index1, -1, field);
	return new OMLCurrencyImpl(cur);
}

OMLCurrency* OMLStructImpl::GetValue(int index1, int index2, const char* field) const
{
	Currency cur = _sd->GetValue(index1, index2, field);
	return new OMLCurrencyImpl(cur);
}

int OMLStructImpl::GetRows() const
{
	return _sd->M();
}

int OMLStructImpl::GetCols() const
{
	return _sd->N();
}

void OMLStructImpl::SetValue(int index, const char* field, OMLCurrency* val)
{
	OMLCurrencyImpl* ci = (OMLCurrencyImpl*)val;
	_sd->SetValue(index, -1, field, ci->GetCurrency());
}

void OMLStructImpl::SetValue(int index1, int index2, const char* field, OMLCurrency* val)
{
	OMLCurrencyImpl* ci = (OMLCurrencyImpl*)val;
	_sd->SetValue(index1, index2, field, ci->GetCurrency());
}

OMLCurrency* OMLStructImpl::GetCurrency() const
{
	return (new OMLCurrencyImpl(_sd));
}

StructData* OMLStructImpl::GetStructData() const
{
	return _sd;
}