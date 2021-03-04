/**
* @file Currency.h
* @date August 2013
* Copyright (C) 2013-2018 Altair Engineering, Inc.  
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

#ifndef __Currency_h
#define __Currency_h

#include "OMLDll.h"

#include <string>
#include <vector>
#include <set>

#include "hwComplex.h"

template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

template <typename T1, typename T2> class hwTMatrixN;
typedef hwTMatrixN<double, hwTComplex<double> > hwMatrixN;

template <typename T1, typename T2> class hwTMatrixS;
typedef hwTMatrixS<double, hwTComplex<double> > hwMatrixS;

class OMLDLL_DECLS Currency;
class OMLDLL_DECLS CurrencyDisplay;
class OutputFormat;
class FunctionInfo;
class StructData;

typedef hwTMatrix<Currency, void*> HML_CELLARRAY;
typedef hwTMatrixN<Currency, void*> HML_ND_CELLARRAY;

typedef Currency (*EXTPTR) (const std::string&);


class StringManager
{
public:
	StringManager() {}
	~StringManager();

	struct StringManagerLTPred
    {
        bool operator () (const std::string* lhs, const std::string* rhs) const
		{        
            return (*lhs < *rhs);
        }
    };
	typedef std::set<const std::string*, StringManagerLTPred> StringStorage;

	const std::string* GetStringPointer(const std::string&);

private:
	StringStorage _strings;
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//! Currency implementation
//------------------------------------------------------------------------------
class OMLDLL_DECLS Currency
{
public:
	Currency(double val);
	Currency(int val);
	Currency(size_t val);
	Currency(double val, int in_type); // used for break, error, etc.
	Currency(double val, std::string error_msg);
	Currency(const std::string& str);
	Currency(const char* str);
	Currency(const std::vector<double>& data);
	Currency(bool logical_val);
	Currency(hwMatrix* data);
	Currency(hwMatrixN* data);
	Currency(hwMatrixS* data);
	Currency(const hwComplex& cplx);
	Currency(); // Microsoft STL forces this 
	Currency(const Currency& cur);
	Currency(HML_CELLARRAY* cells);
	Currency(HML_ND_CELLARRAY* cells);
	Currency(FunctionInfo* fi);
	Currency(StructData* sd);
    Currency(OutputFormat* fmt);
    Currency(Currency* ptr);	
    //! Constructor for swig bound objects
    //! \param[in] obj  Pointer to bound object
    //! \param[in] name Class name associated with the bound object
    Currency( void*              obj, 
              const std::string& name);            

	~Currency();
	Currency& operator= (const Currency&);
 
	Currency(Currency&& cur);
	Currency& operator= (Currency&&);
	void Swap(Currency&);

    //! Gets output as a string for printing
    //! \param[in] fmt Output format
	std::string GetOutputString( const OutputFormat* fmt) const;
	std::string GetTypeString() const;

	int GetType() const { return type; }

    //! Returns true if currency is of the given type
    bool  IsScalar()    const;
	bool  IsString()    const;
	bool  IsUTF8String()    const { return _is_utf8; }
	void  CheckForUTF8Characters();
	bool  IsMultilineString() const;
	bool  IsVector()    const;
	bool  IsRealVector() const;
	bool  IsMatrix()    const;
	bool  IsNDMatrix()    const;
	bool  IsSparse()    const;
	bool  IsComplex()   const;
	bool  IsMatrixOrString() const { return type == TYPE_MATRIX; }
	bool  IsColon()     const      { return type == TYPE_COLON; }
	bool  IsBreak()      const     { return type == TYPE_BREAK; }
	bool  IsContinue()      const  { return type == TYPE_CONTINUE; }
	bool  IsReturn()      const    { return type == TYPE_RETURN; }
	bool  IsError()     const      { return type == TYPE_ERROR; }
	bool  IsSyntaxError() const;
	bool  IsCellArray() const      { return type == TYPE_CELLARRAY; }
	bool  IsNDCellArray() const    { return type == TYPE_ND_CELLARRAY; }
	bool  IsFunctionHandle() const { return type == TYPE_FUNCHANDLE; }
	bool  IsStruct() const         { return type == TYPE_STRUCT; }
	bool  IsObject() const         { return type == TYPE_OBJECT; }
	bool  IsNothing() const        { return type == TYPE_NOTHING; }
    bool  IsFormat() const         { return type == TYPE_FORMAT; }
    bool  IsBreakpoint() const     { return type == TYPE_BREAKPOINT; }
	bool  IsEmpty() const;
	bool  IsInteger() const;
	bool  IsPositiveInteger() const;
	bool  IsPositiveIntegralVector() const;
	bool  IsPositiveIntegralMatrix() const;
	bool  IsLinearRange() const;
	bool  IsLogical() const;
	bool  IsCharacter() const;
	bool  IsCellList() const;
	bool  IsPointer()       const { return type == TYPE_POINTER; }
    bool  IsBoundObject()   const { return type == TYPE_BOUNDOBJECT; }

	void  MakeStruct();
	void  MakeCell();

	void  ReplaceMatrix(hwMatrix* new_matrix);
	void  ReplaceCellArray(HML_CELLARRAY* new_cells);
	void  ReplaceStruct(StructData* new_sd);

	double              Scalar()    const;
	std::string         StringVal() const; 
	std::string         Message() const        { return *message; } 
	std::vector<double> Vector() const;
	const hwMatrix*     Matrix() const;
	hwMatrix*           GetWritableMatrix();
	const hwMatrixN*    MatrixN() const;
	const hwMatrixS*    MatrixS() const;
	hwMatrixN*          GetWritableMatrixN();
	hwMatrixS*          GetWritableMatrixS();
	hwComplex           Complex() const;
	double              Real() const           { return data.complex->Real(); }
	double              Imag() const           { return data.complex->Imag(); }
	HML_CELLARRAY*      CellArray() const      { return data.cells; }
	HML_ND_CELLARRAY*   CellArrayND() const    { return data.cells_nd; }
	FunctionInfo*       FunctionHandle() const { return data.func; }
	StructData*         Struct() const         { return data.sd; }
    OutputFormat*       Format() const         { return data.format; }
	Currency*           Pointer() const        { return data.cur_ptr; }
    void*               BoundObject() const    { return data.boundobj; }

	void                SetMask(int new_mask);
	int                 GetMask() const        { return mask; }
	void                SetOutputName(const std::string& name) const; // the const is a mistake and needs to be fixed -- JDS
	void                SetOutputName(const std::string* name) const; // the const is a mistake and needs to be fixed -- JDS
	void                ClearOutputName();
	const std::string*  GetOutputNamePtr() const;
	std::string         GetOutputName() const { return *GetOutputNamePtr(); }

	void                SetClass(const std::string& name);
	std::string         GetClassname() const;

	void                ReplaceScalar(double new_value) { data.value = new_value; type = TYPE_SCALAR; }
	void                ReplaceComplex(hwComplex new_value);

	const hwMatrix*     ConvertToMatrix() const;
	const hwMatrix*     ExpandMatrix(const hwMatrix*) const;
	HML_CELLARRAY*      ConvertToCellArray();

	void                FlattenCellList();

	void                ConvertToStruct();

	void                SetLinearRange(bool);

	enum CurrencyType { TYPE_SCALAR, TYPE_STRING, TYPE_MATRIX, TYPE_COLON, TYPE_COMPLEX, TYPE_CELLARRAY, TYPE_ERROR, TYPE_BREAK, TYPE_RETURN, TYPE_FUNCHANDLE, TYPE_STRUCT, TYPE_NOTHING, TYPE_FORMAT, TYPE_BREAKPOINT, TYPE_POINTER, TYPE_CONTINUE, TYPE_ND_MATRIX, TYPE_OBJECT, TYPE_BOUNDOBJECT, TYPE_ND_CELLARRAY, TYPE_SPARSE };
	enum MaskType { MASK_NONE, MASK_DOUBLE, MASK_STRING, MASK_LOGICAL, MASK_CELL_LIST, MASK_EXPLICIT_COMPLEX };

	static StringManager vm;
	static StringManager pm;

	static std::set<void*> ignore_cow_pointers;
	static bool IgnoreCoW(void*);

    //! Creates, if needed and returns display
    CurrencyDisplay* GetDisplay() const;
    //! Sets matrix display
    //! \param[in] display Given display
    void SetDisplay( CurrencyDisplay* display) { _display = display; }
    //! Deletes display
    void DeleteDisplay();

    //! Returns true if evaluator is in experimental mode
    static bool GetExperimental() { return _experimental; }
    //! Sets the experimental mode flag (-ex)
    //! \param[in] val Sets to true if the evaluator is in experimental mode
    static void SetExperimental( bool val) { _experimental = val; }

    //! True if output is from 'disp' command
    bool IsDispOutput() const { return (_outputType == OUTPUT_TYPE_DISP); }
    //! Sets output type to be from 'disp' command
    void DispOutput() { _outputType = OUTPUT_TYPE_DISP; }
    //! 
    //! Resets to default output
    //!
    void ResetOutputType() { _outputType = OUTPUT_TYPE_DEFAULT; }

    //! True if output is from 'printf'/'fprintf' commands
    bool IsPrintfOutput() const { return (_outputType == OUTPUT_TYPE_PRINTF); }
    //! Sets output type to be from printf/fprintf commands
    void PrintfOutput() { _outputType = OUTPUT_TYPE_PRINTF; }

    //! Utility to get just the values, without any header(s), if applicable
    //! \param[in] fmt Format
    std::string GetValues( const OutputFormat* fmt) const;
    //!
    //! Utility to get just the values without any headers or currency name
    //! \param fmt Format
    //!
    std::string GetValuesForDisplay(const OutputFormat* fmt) const;

private:
	Currency(const hwMatrix* data); // stubbed
	Currency(const hwMatrixS* data); // stubbed
	Currency(const hwMatrixN* data); // stubbed

	void  Copy(const Currency&);
	void  DeleteMatrix(hwMatrix*);
	void  DeleteMatrixN(hwMatrixN*);
	void  DeleteMatrixS(hwMatrixS*);
	void  DeleteCells(HML_CELLARRAY*);
	void  DeleteCellsN(HML_ND_CELLARRAY*);
	void  DeleteStruct(StructData*);
	void  DeleteFunctionInfo(FunctionInfo*);

	mutable CurrencyType type;
	mutable MaskType mask;

	union DataStorage
	{
		double               value;
		hwComplex*           complex;
		hwMatrix*            mtx;
		hwMatrixN*           mtxn;
		hwMatrixS*			 mtxs;
		HML_CELLARRAY*       cells;
		HML_ND_CELLARRAY*    cells_nd;
		FunctionInfo*        func;
		StructData*          sd;
		OutputFormat*        format;
		Currency*            cur_ptr;   
        void*                boundobj;        //! Bound object pointer
	} mutable data;

    enum OUTPUT_TYPE                          //! Currency output type
    {
        OUTPUT_TYPE_DEFAULT,                  //! Default
        OUTPUT_TYPE_DISP,                     //! disp output
        OUTPUT_TYPE_PRINTF                    //! printf/fprintf output
    };

    OUTPUT_TYPE                _outputType;   //! Type of output
	std::string*               message;
	std::string*               classname;
	mutable const std::string* out_name;
    mutable CurrencyDisplay*   _display;        //! Displays for output
	bool                       _is_utf8;
	bool                       _is_linear_range;
    
    static bool _experimental;                  //! True if experimental mode is active
};

#endif
