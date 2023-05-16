/**
* @file OMLInterface.h
* @date January 2017
* Copyright (C) 2017-2018 Altair Engineering, Inc.  
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

#pragma once
#ifndef __OMLInterface_h
#define __OMLInterface_h

#include "OMLInterfacePublic.h"

#include "StructData.h"
#include "Evaluator.h"

#include <mutex>

template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;
typedef hwTMatrix<Currency, void*> HML_CELLARRAY;

// NEVER create any of these Impl classes on the stack.  That will
// break the new, automatic garbage collection in ways that will be difficult
// to find via debugging

class OMLImplBase
{
public:
	OMLImplBase(EvaluatorInterface* in_eval);
	virtual ~OMLImplBase();

protected:
	EvaluatorInterface* _eval;
};

class OMLInterfaceImpl : public OMLInterface5, OMLImplBase
{
public:
	OMLInterfaceImpl(EvaluatorInterface* in_eval);
	~OMLInterfaceImpl();

	void RegisterFunction(const char*, ALT_FUNCPTR);
	void RegisterHiddenFunction(const char*, ALT_FUNCPTR);
	void RegisterFunctionWithMetadata(const char*, ALT_FUNCPTR, const char*, int, int);
	void RegisterFunctionWithMetadata(const char*, ALT_FUNCPTR, const char*, int, int, bool);

	void ThrowError(const char*);

	int  Nargout() const;

	const OMLCurrency* GetGlobalValue(const char*);

	const OMLCurrency* CallFunction(const OMLFunctionHandle* handle, OMLCurrencyList* inputs);
	const OMLCurrency* CallFunction(const char* name, OMLCurrencyList* inputs);

	OMLCurrencyList* CreateCurrencyList();
};

class OMLCurrencyImpl : public OMLCurrency4, OMLImplBase
{
public:
	OMLCurrencyImpl(EvaluatorInterface* in_eval, Currency in_cur);
	~OMLCurrencyImpl();

	bool IsScalar() const;
	bool IsComplex() const;
	bool IsString() const;
	bool IsMatrix() const;
	bool IsNDMatrix() const;
	bool IsCellArray() const;
	bool IsNDCellArray() const;
	bool IsSparseMatrix() const;
	bool IsStruct() const;
	bool IsFunctionHandle() const;
	bool IsLogical() const;

	Currency GetCurrency() const { return _cur; }

	double                   GetScalar() const;
	const char*              GetString() const;
	const OMLCellArray*      GetCellArray() const;
	const OMLMatrix*         GetMatrix() const;
	const OMLNDMatrix*       GetNDMatrix() const;
	const OMLComplex*        GetComplex() const;
	const OMLStruct*         GetStruct() const;
	const OMLFunctionHandle* GetFunctionHandle() const;
	bool                     GetLogical() const;

private:
	Currency _cur;
};

class OMLComplexImpl : public OMLComplex, OMLImplBase
{
public:
	OMLComplexImpl(EvaluatorInterface* in_eval, double real, double imag);
	~OMLComplexImpl();

	double GetReal() const;
	double GetImag() const;

	OMLCurrency* GetCurrency() const;

	static void GarbageCollect();

private:
	hwComplex      cplx;
};

class OMLMatrixImpl : public OMLMatrix, OMLImplBase
{
public:
	OMLMatrixImpl(EvaluatorInterface* in_eval, const hwMatrix* in_mtx);
	~OMLMatrixImpl();

	bool    IsReal() const;

	int     GetRows() const;
	int     GetCols() const;

	const double* GetRealData() const;
	const double* GetImaginaryData() const;

	OMLCurrency*   GetCurrency() const;
	hwMatrix*      GetMatrixPointer() const;

private:
	hwMatrix*  _mtx;
};

class OMLNDMatrixImpl : public OMLNDMatrix, OMLImplBase
{
public:
	OMLNDMatrixImpl(EvaluatorInterface* in_eval, const hwMatrixN* in_mtx);
	~OMLNDMatrixImpl();

	bool    IsReal() const;

	int     GetNumDimension() const;
	int     GetDimension(int) const;

	const double* GetRealData() const;
	const double* GetImaginaryData() const;

	OMLCurrency*   GetCurrency() const;
	hwMatrixN*     GetMatrixPointer() const;

private:
	hwMatrixN* _mtx;
};

class OMLSparseMatrixImpl : public OMLSparseMatrix, OMLImplBase
{
public:
	OMLSparseMatrixImpl(EvaluatorInterface* in_eval, const hwMatrixS* in_mtx);
	~OMLSparseMatrixImpl();

	bool    IsReal() const;

	int     GetRows() const;
	int     GetCols() const;

	const double* GetRealData() const;
	const double* GetImaginaryData() const;

	const int* GetRowVector() const;
	const int* GetColumnVector() const;

	OMLCurrency*   GetCurrency() const;
	hwMatrixS*      GetMatrixPointer() const;

private:
	hwMatrixS*  _mtxs;
};

class OMLCellArrayImpl : public OMLCellArray, OMLImplBase
{
public:
	OMLCellArrayImpl(EvaluatorInterface* in_eval, HML_CELLARRAY* in_cells);
	OMLCellArrayImpl(EvaluatorInterface* in_eval, HML_CELLARRAY* in_cells, bool temp);

	~OMLCellArrayImpl();

	OMLCurrency* GetValue(int index1) const;
	OMLCurrency* GetValue(int index1, int index2) const;

	int          GetRows() const;
	int          GetCols() const;

	void         SetValue(int index1, OMLCurrency* val);
	void         SetValue(int index1, int index2, OMLCurrency* val);

	OMLCurrency*   GetCurrency() const;
	HML_CELLARRAY* GetCells() const;

private:
	HML_CELLARRAY* _cells;
	bool           _is_temp;

};

class OMLNDCellArrayImpl : public OMLNDCellArray, OMLImplBase
{
public:
	OMLNDCellArrayImpl(EvaluatorInterface* in_eval, HML_ND_CELLARRAY* in_cells);
	~OMLNDCellArrayImpl();

	int     GetNumDimension() const;
	int     GetDimension(int) const;

	OMLCurrency* GetValue(int index1) const;
	void         SetValue(int index1, OMLCurrency* val);

	OMLCurrency*      GetCurrency() const;
	HML_ND_CELLARRAY* GetCells() const;

private:
	HML_ND_CELLARRAY* _cells;
};

class OMLStructImpl : public OMLStruct, OMLImplBase
{
public:
	OMLStructImpl(EvaluatorInterface* in_eval, StructData* in_sd);
	~OMLStructImpl();

	OMLCurrency* GetValue(int index1, const char* field) const;
	OMLCurrency* GetValue(int index1, int index2, const char* field) const;

	int          GetRows() const;
	int          GetCols() const;

	void         SetValue(int index, const char* field, OMLCurrency* val);
	void         SetValue(int index1, int index2, const char* field, OMLCurrency* val);

	OMLCurrency* GetCurrency() const;
	StructData*  GetStructData() const;

private:
	StructData* _sd;
};

class OMLFunctionHandleImpl : public OMLFunctionHandle, OMLImplBase
{
public:
	OMLFunctionHandleImpl(EvaluatorInterface* in_eval, const FunctionInfo* in_fh);
	~OMLFunctionHandleImpl();

	FunctionInfo* GetFunctionInfo() const;

private:
	FunctionInfo* _fi;
};

class OMLCurrencyListImpl : public OMLCurrencyList3, OMLImplBase
{
public:
	OMLCurrencyListImpl(EvaluatorInterface* in_eval);
	~OMLCurrencyListImpl();

	int Size() const;
	const OMLCurrency* Get(int idx) const;

	void AddScalar(double);
	void AddString(const char*);
	void AddLogical(bool);
	void AddCellArray(OMLCellArray*);
	void AddCellArray(HML_CELLARRAY* cells);
	void AddNDCellArray(HML_ND_CELLARRAY* cells);
	void AddNDCellArray(OMLNDCellArray*);
	void AddMatrix(OMLMatrix*);
	void AddMatrix(const hwMatrix*);
	void AddNDMatrix(OMLNDMatrix*);
	void AddNDMatrix(const hwMatrixN*);
	void AddSparseMatrix(OMLSparseMatrix*);
	void AddSparseMatrix(const hwMatrixS*);
	void AddComplex(OMLComplex*);
	void AddComplex(hwComplex);
	void AddStruct(OMLStruct*);
	void AddStruct(StructData*);
	void AddFunctionHandle(FunctionInfo*);

	double* AllocateData(int size);

	// I'd love for these to be static, but since there are no static virtual functions,
	// I have to either do this or play the factory game
	OMLCurrency*  CreateCurrencyFromDouble(double dbl);
	OMLCurrency*  CreateCurrencyFromString(const char* str);

	OMLCellArray* CreateCellArray(int rows, int cols);
	OMLMatrix*    CreateMatrix(int rows, int cols, double* data);
	OMLMatrix*    CreateMatrix(int rows, int cols, double* real, double* imag);
	OMLNDMatrix*  CreateNDMatrix(int num_dims, int* dims, double* real);
	OMLNDMatrix*  CreateNDMatrix(int num_dims, int* dims, double* real, double* imag);
	OMLComplex*   CreateComplex(double real, double imag);
	OMLStruct*    CreateStruct(int rows, int cols);


	OMLNDCellArray*  CreateNDCellArray(int num_dims, int* dims);
	OMLSparseMatrix* CreateSparseMatrix(int num_vals, int* ivec, int* jvec, double* vals, int rows, int cols);

	OMLCellArray* CreateTemporaryCellArray(int rows, int cols);

private:
	void Expand();

	OMLCurrency** _list;
	int           _count;
};

#endif
