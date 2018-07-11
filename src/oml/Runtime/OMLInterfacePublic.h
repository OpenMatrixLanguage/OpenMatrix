/**
* @file OMLInterfacePublic.cpp
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

#pragma once
#ifndef __OMLInterfacePublic_h
#define __OMLInterfacePublic_h

class OMLCurrency;
class OMLCellArray;
class OMLMatrix;
class OMLNDMatrix;
class OMLComplex;
class OMLStruct;
class OMLCurrencyList;
class OMLFunctionHandle;
class OMLInterface;

typedef bool (*ALT_FUNCPTR) (OMLInterface*, const OMLCurrencyList* ins, OMLCurrencyList* outs);

class OMLInterface
{
public:
	virtual void RegisterFunction(const char*, ALT_FUNCPTR) = 0;

	virtual void ThrowError(const char*) = 0;

	virtual int  Nargout() const = 0;

	virtual ~OMLInterface() {};
};

class OMLInterface2 : public OMLInterface
{
public:
	virtual const OMLCurrency* CallFunction(const OMLFunctionHandle* handle, OMLCurrencyList* inputs) = 0;

	virtual OMLCurrencyList* CreateCurrencyList() = 0;
};

class OMLCurrency
{
public:
	virtual bool IsScalar() const = 0;
	virtual bool IsComplex() const = 0;
	virtual bool IsString() const = 0;
	virtual bool IsMatrix() const = 0;
	virtual bool IsNDMatrix() const = 0;
	virtual bool IsCellArray() const = 0;
	virtual bool IsStruct() const = 0;

	virtual double              GetScalar() const = 0;
	virtual const char*         GetString() const = 0;
	virtual const OMLCellArray* GetCellArray() const = 0;
	virtual const OMLMatrix*    GetMatrix() const = 0;
	virtual const OMLNDMatrix*  GetNDMatrix() const = 0;
	virtual const OMLComplex*   GetComplex() const = 0;
	virtual const OMLStruct*    GetStruct() const = 0;

	virtual ~OMLCurrency() {};
};

class OMLCurrency2 : public OMLCurrency
{
public:
	virtual bool IsFunctionHandle() const = 0;

	virtual const OMLFunctionHandle* GetFunctionHandle() const = 0;
};

class OMLComplex
{
public:
	virtual double GetReal() const = 0;
	virtual double GetImag() const = 0;

	virtual OMLCurrency* GetCurrency() const = 0;

	virtual ~OMLComplex() {};
};

class OMLMatrix
{
public:
	virtual bool    IsReal() const = 0;

	virtual int     GetRows() const = 0;
	virtual int     GetCols() const = 0;

	virtual const double* GetRealData() const = 0;
	virtual const double* GetImaginaryData() const = 0;

	virtual OMLCurrency* GetCurrency() const = 0;

	virtual ~OMLMatrix() {};
};

class OMLNDMatrix
{
public:
	virtual bool    IsReal() const = 0;

	virtual int     GetNumDimension() const = 0;
	virtual int     GetDimension(int) const = 0;

	virtual const double* GetRealData() const = 0;
	virtual const double* GetImaginaryData() const = 0;

	virtual OMLCurrency* GetCurrency() const = 0;

	virtual ~OMLNDMatrix() {};
};

class OMLCellArray
{
public:
	virtual OMLCurrency* GetValue(int index1) const = 0;
	virtual OMLCurrency* GetValue(int index1, int index2) const = 0;

	virtual int          GetRows() const = 0;
	virtual int          GetCols() const = 0;

	virtual void         SetValue(int index1, OMLCurrency* val) = 0;
	virtual void         SetValue(int index1, int index2, OMLCurrency* val) = 0;

	virtual OMLCurrency* GetCurrency() const = 0;

	virtual ~OMLCellArray() {};
};

class OMLStruct
{
public:
	virtual OMLCurrency* GetValue(int index1, const char* field) const = 0;
	virtual OMLCurrency* GetValue(int index1, int index2, const char* field) const = 0;

	virtual int          GetRows() const = 0;
	virtual int          GetCols() const = 0;

	virtual void         SetValue(int index, const char* field, OMLCurrency* val) = 0;
	virtual void         SetValue(int index1, int index2, const char* field, OMLCurrency* val) = 0;

	virtual OMLCurrency* GetCurrency() const = 0;

	virtual ~OMLStruct() {};
};

class OMLFunctionHandle
{
};

class OMLCurrencyList
{
public:
	virtual int Size() const = 0;
	virtual const OMLCurrency* Get(int idx) const = 0;

	virtual void AddScalar(double) = 0;
	virtual void AddString(const char*) = 0;
	virtual void AddCellArray(OMLCellArray*) = 0;
	virtual void AddStruct(OMLStruct*) = 0;
	virtual void AddMatrix(OMLMatrix*) = 0;
	virtual void AddNDMatrix(OMLNDMatrix*) = 0;
	virtual void AddComplex(OMLComplex*) = 0;

	virtual double* AllocateData(int size) = 0;

	// These two are needed to push doubles/strings into cells/structs
	virtual OMLCurrency*  CreateCurrencyFromDouble(double dbl) = 0;
	virtual OMLCurrency*  CreateCurrencyFromString(const char* str) = 0;

	virtual OMLCellArray* CreateCellArray(int rows, int cols) = 0;
	virtual OMLStruct*    CreateStruct(int rows, int cols) = 0;
	virtual OMLMatrix*    CreateMatrix(int rows, int cols, double* data) = 0;
	virtual OMLMatrix*    CreateMatrix(int rows, int cols, double* real, double* imag) = 0;
	virtual OMLNDMatrix*  CreateNDMatrix(int num_dims, int* dims, double* real) = 0;
	virtual OMLNDMatrix*  CreateNDMatrix(int num_dims, int* dims, double* real, double* imag) = 0;
	virtual OMLComplex*   CreateComplex(double real, double imag) = 0;

	virtual ~OMLCurrencyList() {};
};

#endif