/**
* @file OMLInterface.h
* @date January 2017
* Copyright (C) 2017-2024 Altair Engineering, Inc.  
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

	OMLImplBase() : _eval(nullptr) {}
};

class OMLInterfaceImpl : public OMLInterface5, OMLImplBase
{
public:
	OMLInterfaceImpl(EvaluatorInterface* in_eval);
	~OMLInterfaceImpl();

	void RegisterFunction(const char*, ALT_FUNCPTR) override;
	void RegisterHiddenFunction(const char*, ALT_FUNCPTR) override;
	void RegisterFunctionWithMetadata(const char*, ALT_FUNCPTR, const char*, int, int) override;
	void RegisterFunctionWithMetadata(const char*, ALT_FUNCPTR, const char*, int, int, bool) override;

	void ThrowError(const char*) override;

	int  Nargout() const override;

	const OMLCurrency* GetGlobalValue(const char*) override;

	const OMLCurrency* CallFunction(const OMLFunctionHandle*, OMLCurrencyList*) override;
	const OMLCurrency* CallFunction(const char*, OMLCurrencyList*) override;

	OMLCurrencyList* CreateCurrencyList() override;\
};

class OMLCurrencyImpl : public OMLCurrency4, OMLImplBase
{
public:
	OMLCurrencyImpl(EvaluatorInterface* in_eval, const Currency& in_cur);
	~OMLCurrencyImpl();

	bool IsScalar() const override;
	bool IsComplex() const override;
	bool IsString() const override;
	bool IsMatrix() const override;
	bool IsNDMatrix() const override;
	bool IsCellArray() const override;
	bool IsNDCellArray() const override;
	bool IsSparseMatrix() const override;
	bool IsStruct() const override;
	bool IsFunctionHandle() const override;
	bool IsLogical() const override;

	Currency GetCurrency() const { return _cur; }

	double                   GetScalar() const override;
	const char*              GetString() const override;
	const OMLCellArray*      GetCellArray() const override;
	const OMLMatrix*         GetMatrix() const override;
	const OMLNDMatrix*       GetNDMatrix() const override;
	const OMLComplex*        GetComplex() const override;
	const OMLStruct*         GetStruct() const override;
	const OMLFunctionHandle* GetFunctionHandle() const override;
	bool                     GetLogical() const override;

private:
	Currency _cur;
};

class OMLComplexImpl : public OMLComplex, OMLImplBase
{
public:
	OMLComplexImpl(EvaluatorInterface* in_eval, double real, double imag);
	~OMLComplexImpl();

	double GetReal() const override;
	double GetImag() const override;

	OMLCurrency* GetCurrency() const override;

	static void GarbageCollect();

private:
	hwComplex      cplx;
};

class OMLMatrixImpl : public OMLMatrix, OMLImplBase
{
public:
	OMLMatrixImpl(EvaluatorInterface* in_eval, const hwMatrix* in_mtx);
	~OMLMatrixImpl();

	bool    IsReal() const override;

	int     GetRows() const override;
	int     GetCols() const override;

	const double* GetRealData() const override;
	const double* GetImaginaryData() const override;

	OMLCurrency*   GetCurrency() const override;
	hwMatrix*      GetMatrixPointer() const;

private:
	hwMatrix*  _mtx;
};

class OMLNDMatrixImpl : public OMLNDMatrix, OMLImplBase
{
public:
	OMLNDMatrixImpl(EvaluatorInterface* in_eval, const hwMatrixN* in_mtx);
	~OMLNDMatrixImpl();

	bool    IsReal() const override;

	int     GetNumDimension() const override;
	int     GetDimension(int) const override;

	const double* GetRealData() const override;
	const double* GetImaginaryData() const override;

	OMLCurrency*   GetCurrency() const override;
	hwMatrixN*     GetMatrixPointer() const;

private:
	hwMatrixN* _mtx;
};

class OMLSparseMatrixImpl : public OMLSparseMatrix, OMLImplBase
{
public:
	OMLSparseMatrixImpl(EvaluatorInterface* in_eval, const hwMatrixS* in_mtx);
	~OMLSparseMatrixImpl();

	bool    IsReal() const override;

	int     GetRows() const override;
	int     GetCols() const override;

	const double* GetRealData() const override;
	const double* GetImaginaryData() const override;

	const int* GetRowVector() const override;
	const int* GetColumnVector() const override;

	OMLCurrency*   GetCurrency() const override;
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

	OMLCurrency* GetValue(int index1) const override;
	OMLCurrency* GetValue(int index1, int index2) const override;

	int          GetRows() const override;
	int          GetCols() const override;

	void         SetValue(int index1, OMLCurrency* val) override;
	void         SetValue(int index1, int index2, OMLCurrency* val) override;

	OMLCurrency*   GetCurrency() const override;
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

	int     GetNumDimension() const override;
	int     GetDimension(int) const override;

	OMLCurrency* GetValue(int index1) const override;
	void         SetValue(int index1, OMLCurrency* val) override;

	OMLCurrency*      GetCurrency() const override;
	HML_ND_CELLARRAY* GetCells() const;

private:
	HML_ND_CELLARRAY* _cells;
};

class OMLStructImpl : public OMLStruct, OMLImplBase
{
public:
	OMLStructImpl(EvaluatorInterface* in_eval, StructData* in_sd);
	~OMLStructImpl();

	OMLCurrency* GetValue(int index1, const char* field) const override;
	OMLCurrency* GetValue(int index1, int index2, const char* field) const override;

	int          GetRows() const override;
	int          GetCols() const override;

	void         SetValue(int index, const char* field, OMLCurrency* val) override;
	void         SetValue(int index1, int index2, const char* field, OMLCurrency* val) override;

	OMLCurrency* GetCurrency() const override;
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


class OMLASTImpl : public OMLAST
{
public:
	OMLASTImpl(EvaluatorInterface* in_eval, OMLTree* in_tree);
	~OMLASTImpl();

	virtual void AddChild(OMLAST* child_tree) override;
	OMLTree* GetTree();

private:
	OMLTree* _tree;

	OMLASTImpl() : _tree (nullptr) {}
};


class OMLCurrencyListImpl : public OMLCurrencyList4, OMLImplBase
{
public:
	OMLCurrencyListImpl(EvaluatorInterface* in_eval);
	~OMLCurrencyListImpl();

	int Size() const override;
	const OMLCurrency* Get(int idx) const override;

	void AddScalar(double) override;
	void AddString(const char*) override;
	void AddLogical(bool) override;
	void AddCellArray(OMLCellArray*) override;
	void AddCellArray(HML_CELLARRAY* cells);
	void AddNDCellArray(HML_ND_CELLARRAY* cells);
	void AddNDCellArray(OMLNDCellArray*) override;
	void AddMatrix(OMLMatrix*) override;
	void AddMatrix(const hwMatrix*);
	void AddNDMatrix(OMLNDMatrix*);
	void AddNDMatrix(const hwMatrixN*);
	void AddSparseMatrix(OMLSparseMatrix*) override;
	void AddSparseMatrix(const hwMatrixS*);
	void AddComplex(OMLComplex*) override;
	void AddComplex(hwComplex);
	void AddStruct(OMLStruct*) override;
	void AddStruct(const StructData*);
	void AddFunctionHandle(FunctionInfo*);

	double* AllocateData(int size) override;

	// I'd love for these to be static, but since there are no static virtual functions,
	// I have to either do this or play the factory game
	OMLCurrency*  CreateCurrencyFromDouble(double dbl) override;
	OMLCurrency*  CreateCurrencyFromString(const char* str) override;

	OMLCellArray* CreateCellArray(int rows, int cols) override;
	OMLMatrix*    CreateMatrix(int rows, int cols, double* data) override;
	OMLMatrix*    CreateMatrix(int rows, int cols, double* real, double* imag) override;
	OMLNDMatrix*  CreateNDMatrix(int num_dims, int* dims, double* real) override;
	OMLNDMatrix*  CreateNDMatrix(int num_dims, int* dims, double* real, double* imag) override;
	OMLComplex*   CreateComplex(double real, double imag) override;
	OMLStruct*    CreateStruct(int rows, int cols) override;

	OMLNDCellArray*  CreateNDCellArray(int num_dims, int* dims);
	OMLSparseMatrix* CreateSparseMatrix(int num_vals, int* ivec, int* jvec, double* vals, int rows, int cols);

	OMLCellArray* CreateTemporaryCellArray(int rows, int cols);

	OMLAST*       CreateAST(int type, const char* label);

private:
	void Expand();

	OMLCurrency** _list;
	int           _count;
};


#endif
