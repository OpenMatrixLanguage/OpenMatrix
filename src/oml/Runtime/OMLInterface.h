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

#include "hwMatrix.h"
#include "hwComplex.h"
#include "StructData.h"
#include "Evaluator.h"

typedef hwTMatrix<Currency, void*> HML_CELLARRAY;

class OMLInterfaceImpl : public OMLInterface4
{
public:
	OMLInterfaceImpl(EvaluatorInterface* in_eval) { _eval = in_eval; }

	void RegisterFunction(const char*, ALT_FUNCPTR);
	void RegisterHiddenFunction(const char*, ALT_FUNCPTR);
	void RegisterFunctionWithMetadata(const char*, ALT_FUNCPTR, const char*, int, int);

	void ThrowError(const char*);

	int  Nargout() const;

	const OMLCurrency* CallFunction(const OMLFunctionHandle* handle, OMLCurrencyList* inputs);
	const OMLCurrency* CallFunction(const char* name, OMLCurrencyList* inputs);

	OMLCurrencyList* CreateCurrencyList();

private:
	EvaluatorInterface* _eval;
};

class OMLCurrencyImpl : public OMLCurrency3
{
public:
	OMLCurrencyImpl(Currency in_cur) { _cur = in_cur; cached_pointers.push_back(this); }

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

	Currency GetCurrency() const { return _cur; }

	double                   GetScalar() const;
	const char*              GetString() const;
	const OMLCellArray*      GetCellArray() const;
	const OMLMatrix*         GetMatrix() const;
	const OMLNDMatrix*       GetNDMatrix() const;
	const OMLComplex*        GetComplex() const;
	const OMLStruct*         GetStruct() const;
	const OMLFunctionHandle* GetFunctionHandle() const;
	
	static void GarbageCollect();

private:
	Currency _cur;

	static std::vector<OMLCurrencyImpl*> cached_pointers;
};

class OMLComplexImpl : public OMLComplex
{
public:
	OMLComplexImpl(double real, double imag) : cplx(real, imag) { cached_pointers.push_back(this); }

	double GetReal() const;
	double GetImag() const;

	OMLCurrency* GetCurrency() const;

	static void GarbageCollect();

private:
	hwComplex      cplx;

	static std::vector<OMLComplexImpl*> cached_pointers;
};

class OMLMatrixImpl : public OMLMatrix
{
public:
	OMLMatrixImpl(const hwMatrix* in_mtx)  { _mtx = (hwMatrix*)in_mtx; cached_pointers.push_back(this); }

		bool    IsReal() const;

		int     GetRows() const;
		int     GetCols() const;

		const double* GetRealData() const;
		const double* GetImaginaryData() const;

		OMLCurrency*   GetCurrency() const;
		hwMatrix*      GetMatrixPointer() const;

		static void GarbageCollect();

private:
		hwMatrix*  _mtx;
		
		static std::vector<OMLMatrixImpl*> cached_pointers;
};

class OMLNDMatrixImpl : public OMLNDMatrix
{
public:
		OMLNDMatrixImpl(const hwMatrixN* in_mtx) { _mtx = (hwMatrixN*)in_mtx; cached_pointers.push_back(this); } 

		bool    IsReal() const;

		int     GetNumDimension() const;
		int     GetDimension(int) const;

		const double* GetRealData() const;
		const double* GetImaginaryData() const;

		OMLCurrency*   GetCurrency() const;
		hwMatrixN*     GetMatrixPointer() const;

		static void GarbageCollect();

private:
		hwMatrixN* _mtx;

		static std::vector<OMLNDMatrixImpl*> cached_pointers;
};

class OMLSparseMatrixImpl : public OMLSparseMatrix
{
public:
	OMLSparseMatrixImpl(const hwMatrixS* in_mtx) { _mtxs = (hwMatrixS*)in_mtx; cached_pointers.push_back(this); }

	bool    IsReal() const;

	int     GetRows() const;
	int     GetCols() const;

	const double* GetRealData() const;
	const double* GetImaginaryData() const;

	const int* GetRowVector() const;
	const int* GetColumnVector() const;

	OMLCurrency*   GetCurrency() const;
	hwMatrixS*      GetMatrixPointer() const;

	static void GarbageCollect();

private:
	hwMatrixS*  _mtxs;

	static std::vector<OMLSparseMatrixImpl*> cached_pointers;
};

class OMLCellArrayImpl : public OMLCellArray
{
public:
	OMLCellArrayImpl(HML_CELLARRAY* in_cells) { _cells = in_cells; cached_pointers.push_back(this); } 

	OMLCurrency* GetValue(int index1) const;
	OMLCurrency* GetValue(int index1, int index2) const;

	int          GetRows() const;
	int          GetCols() const;

	void         SetValue(int index1, OMLCurrency* val);
	void         SetValue(int index1, int index2, OMLCurrency* val);

	OMLCurrency*   GetCurrency() const;
	HML_CELLARRAY* GetCells() const;

	static void GarbageCollect();

private:
	HML_CELLARRAY* _cells;

	static std::vector<OMLCellArrayImpl*> cached_pointers;
};

class OMLNDCellArrayImpl : public OMLNDCellArray
{
public:
	OMLNDCellArrayImpl(HML_ND_CELLARRAY* in_cells) { _cells = in_cells; cached_pointers.push_back(this); }

	int     GetNumDimension() const;
	int     GetDimension(int) const;

	OMLCurrency* GetValue(int index1) const;
	void         SetValue(int index1, OMLCurrency* val);

	OMLCurrency*      GetCurrency() const;
	HML_ND_CELLARRAY* GetCells() const;

	static void GarbageCollect();

private:
	HML_ND_CELLARRAY* _cells;

	static std::vector<OMLNDCellArrayImpl*> cached_pointers;
};

class OMLStructImpl : public OMLStruct
{
public:
	OMLStructImpl(StructData* in_sd) { _sd = in_sd; cached_pointers.push_back(this); } 

	OMLCurrency* GetValue(int index1, const char* field) const;
	OMLCurrency* GetValue(int index1, int index2, const char* field) const;

	int          GetRows() const;
	int          GetCols() const;

	void         SetValue(int index, const char* field, OMLCurrency* val);
	void         SetValue(int index1, int index2, const char* field, OMLCurrency* val);

	OMLCurrency* GetCurrency() const;
	StructData*  GetStructData() const;

	static void GarbageCollect();

private:
	StructData* _sd;

	static std::vector<OMLStructImpl*> cached_pointers;
};

class OMLFunctionHandleImpl : public OMLFunctionHandle
{
public:
		OMLFunctionHandleImpl(const FunctionInfo* in_fh) { _fi = (FunctionInfo*)in_fh; cached_pointers.push_back(this); } 

		FunctionInfo*      GetFunctionInfo() const { return _fi; }

		static void GarbageCollect();

private:
		FunctionInfo* _fi;

		static std::vector<OMLFunctionHandleImpl*> cached_pointers;
};

class OMLCurrencyListImpl : public OMLCurrencyList2
{
public:
	OMLCurrencyListImpl() { _list = nullptr; _count = 0; }
	~OMLCurrencyListImpl();

	int Size() const;
	const OMLCurrency* Get(int idx) const;

	void AddScalar(double);
	void AddString(const char*);
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

private:
	void Expand();

	OMLCurrency** _list;
	int           _count;
};

#endif
