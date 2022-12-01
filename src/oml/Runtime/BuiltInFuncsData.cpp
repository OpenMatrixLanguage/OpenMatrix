/**
* @file BuiltInFuncsData.cpp
* @date June 2016
* Copyright (C) 2016-2022 Altair Engineering, Inc.  
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

// Begin defines/includes
#include "BuiltInFuncsData.h"
#include "CellND.cc"
#include <cassert>
#include <climits>
#include <memory>

#include "BuiltInFuncsUtils.h"
#include "Evaluator.h"
#include "MatrixNDisplay.h"
#include "OML_Error.h"
#include "StructData.h"

#include "hwComplex.h"
#include "hwMatrix.h"

typedef hwTMatrix< double,  hwTComplex<double> > hwMatrix;

//------------------------------------------------------------------------------
// Returns true if the first complex number is less than the second
//------------------------------------------------------------------------------
bool IsGreaterThan(const std::pair<double, int>& p1,
                   const std::pair<double, int>& p2)
{
#ifndef OS_WIN
    if (IsInf_T(p1.first) || IsNaN_T(p1.first))
    {
        return true;
    }
    else if (IsInf_T(p2.first) || IsNaN_T(p2.first))
    {
        return false;
    }
    else if (IsNegInf_T(p1.first))
    {
        return false;
    }
    else if (IsNegInf_T(p2.first))
    {
        return true;
    }
#endif
    return (p1.first > p2.first);
}
//------------------------------------------------------------------------------
// Returns true if the first complex number is less than the second
//------------------------------------------------------------------------------
bool IsLessThan(const std::pair<double, int>& p1,
                const std::pair<double, int>& p2)
{
#ifndef OS_WIN
    if (IsInf_T(p1.first) || IsNaN_T(p1.first))
    {
        return false;
    }
    else if (IsInf_T(p2.first) || IsNaN_T(p2.first))
    {
        return true;
    }
    else if (IsNegInf_T(p1.first))
    {
        return true;
    }
    else if (IsNegInf_T(p2.first))
    {
        return false;
    }
#endif
    return (p1.first < p2.first);
}
//------------------------------------------------------------------------------
// Returns true if the first complex number is less than the second
//------------------------------------------------------------------------------
bool IsLessThanComplex(const std::pair<hwComplex, int>& p1, 
                       const std::pair<hwComplex, int>& p2)
{
    double x = p1.first.Real();
    double y = p2.first.Real();
    if (x < y || x > y)
    {
        return (x < y);
    }

    return (p1.first.Imag() < p2.first.Imag());
}
//------------------------------------------------------------------------------
// Returns true if the first complex number is less than the second
//------------------------------------------------------------------------------
bool IsGreaterThanComplex(const std::pair<hwComplex, int>& p1,
                          const std::pair<hwComplex, int>& p2)
{
    double x = p1.first.Real();
    double y = p2.first.Real();
    if (x < y || x > y)
    {
        return (x > y);
    }
    return (p1.first.Imag() > p2.first.Imag());
}
//------------------------------------------------------------------------------
// Sets fields recursively. First currency is the input, last currency is value
//------------------------------------------------------------------------------
bool BuiltInFuncsData::Setfield(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.empty() ? 0 : inputs.size();
    if (nargin < 3) 
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency in = inputs[0];      
    if (!in.IsStruct() && !in.IsCellArray() && !in.IsMatrixOrString() &&
        !in.Scalar()   && !in.IsComplex())
        throw OML_Error(OML_ERR_CELLMTXSTRUCT, 1, OML_VAR_TYPE);
        
    Currency*           curr = &in;        // Current currency being processed
    std::string         currfield;         // Current field name for structs
    std::pair<int, int> curridx (0, -1);   // Current index
    bool                hasidx = false;    // True if index has been set explicitly

    BuiltInFuncsData funcs;

    // Store prev info as matrices need to be set in prev currency
    Currency*           prev = NULL;       // Previous currency 
    std::string         prevfield;         // Previous field name
    std::pair<int, int> previdx (0, -1);   // Previous index

    // Loop between input(first) and value(last) to get nested indicies
    std::vector<Currency>::const_iterator itr = inputs.begin() + 1; 
    for (int i = 1; itr != inputs.end() - 1 && curr; ++itr, ++i)
    {
        Currency tmp = *itr;
        // Cell arrays that have index info and strings have struct field names
        if (!tmp.IsCellArray() && !tmp.IsString())
            throw OML_Error(OML_ERR_CELLSTRING, i + 1, OML_VAR_TYPE);

        bool hasmoredata = (i < nargin - 2);  // Has more data to process              

        if (tmp.IsString())                   // Process field name - structs
        {
            if (!curr->IsStruct())
                throw OML_Error(OML_ERR_STRUCT, i + 1, OML_VAR_TYPE);

            StructData* sd = curr->Struct();
            if (!hasidx && hasmoredata && (sd->M() > 1 || sd->N() > 1))
                throw OML_Error(OML_ERR_INVALIDSTRUCTINDEX, i + 2, OML_VAR_TYPE);

            currfield = tmp.StringVal();
            if (!hasmoredata)
                continue;

            if (sd->Size() == 1 && !hasidx)  // Create an index
                curridx = std::pair<int, int>(0, -1);

            prev      = curr; 
            prevfield = currfield;
            previdx   = curridx;
            curr      = funcs.GetStructElement(curr->Struct(), currfield, curridx);
            if (!curr)
                throw OML_Error(OML_ERR_INVALIDSTRUCTINDEX, i + 1, OML_VAR_VALUE);

            // Reset as we are done processing field
            hasidx    = false; 
            currfield = "";
            continue;
        }

        // Processing a field index
        curridx = funcs.GetFieldIndex(tmp, i);
        
        bool isScalar = curr->IsScalar();

        if (curr->IsStruct()) 
        {
            if (!currfield.empty() && !hasidx && hasmoredata)
                throw OML_Error(OML_ERR_INVALIDSTRUCTINDEX, i + 1, OML_VAR_TYPE);
            hasidx = true;
            continue; // Need to get a field name to process structs
        }
        else if (curr->IsScalar() || curr->IsComplex() && (curridx.first > 1 || curridx.second > 1))
        {
            *curr = funcs.GetMatrixElement(eval, *curr, curridx);
            break; // No more nested indices to process for matrices
        }
        else if (!curr->IsCellArray()) 
        {
            if (hasmoredata) // This is not the end of the inputs
                throw OML_Error(OML_ERR_NUMARGIN);

            break;           // Only structs/cells can have nested fields
        }

        if (hasmoredata)
        {
            prev      = curr; 
            prevfield = currfield;
            previdx   = curridx;
            curr      = funcs.GetCellElement(eval, curr->CellArray(), curridx);
        }
    }
                          
    if (!curr)
    {
        outputs.push_back(in);
        return true;
    }

    // Set the value
    Currency value = inputs.back();

    if (curr->IsStruct())
        funcs.SetStructElement(value, currfield, curridx, hasidx, *curr);
    else if (curr->IsScalar() || curr->IsComplex())
    {
        if (!value.IsScalar() && !value.IsComplex())
            throw OML_Error(OML_ERR_SCALAR_COMPLEX, (int)nargin, OML_VAR_TYPE);
        *curr = value;
    }
    else if (curr->IsCellArray())
    {
        HML_CELLARRAY* cell = curr->CellArray();
        if (!cell)
            throw OML_Error(OML_ERR_INVALIDINDEX, (int)nargin-1, OML_VAR_VALUE);

        if (curridx.second > 0)
        {
            if (curridx.first - 1 < cell->M() && curridx.second - 1 < cell->N())
                (*curr->CellArray())(curridx.first - 1, curridx.second - 1) = value;
            else
                throw OML_Error(OML_ERR_INVALIDINDEX, (int)nargin-1, OML_VAR_VALUE);
        }
        else if (curridx.first - 1 < cell->Size())         
            (*curr->CellArray())(curridx.first - 1) = value;
        else
            throw OML_Error(OML_ERR_INVALIDINDEX, (int)nargin-1, OML_VAR_VALUE);
    }
    else if (curr->IsMatrixOrString())
    {
        Currency mtx = funcs.SetMatrixElement(eval, *curr, value, curridx, (int)nargin);
        if (prev)
            funcs.SetMatrixParent(mtx, previdx, prevfield, prev);
        else
            *curr = mtx;
    }
    else
        *curr = value;

    outputs.push_back(in);
    return true;
}
//------------------------------------------------------------------------------
// Helper method for setfield which grows struct and gets requested element
//------------------------------------------------------------------------------
Currency* BuiltInFuncsData::GetStructElement(StructData*                in,
                                             const std::string&         field,
                                             const std::pair<int, int>& index) const 
{
    assert(in);
    if (field.empty()) return NULL;

    int index1 = (index.first  > 0) ? index.first  - 1 : 0;
    int index2 = (index.second > 0) ? index.second - 1 : -1;

    // Check if struct needs to grow
    if (in->M() < index.first || (in->N() < index.second && index.second > 0)) 
        in->SetValue(index1, index2, field, EvaluatorInterface::allocateStruct());
    
    return const_cast<Currency*>(in->GetPointer(index1, index2, field));
}
//------------------------------------------------------------------------------
// Helper method for Setfield to get indices from cell array
//------------------------------------------------------------------------------
std::pair<int, int> BuiltInFuncsData::GetFieldIndex(const Currency& in, int idx) const
{
    std::pair<int, int> fieldindex(0, -1);

    assert(in.IsCellArray());
    HML_CELLARRAY* cell = in.CellArray();
    assert(cell);

    int cellsize = cell->Size();
    if (cellsize == 0 || cellsize > 2)
        throw OML_Error(OML_ERR_INVALIDSTRUCTINDEX, idx + 1, OML_VAR_TYPE);

    if (!(*cell)(0).IsPositiveInteger())
        throw OML_Error(OML_ERR_INVALID_INDEX, idx + 1, OML_VAR_INDEX);
    fieldindex.first = static_cast<int>((*cell)(0).Scalar());

    if (cellsize > 1)
    {
        if (!(*cell)(1).IsPositiveInteger())
            throw OML_Error(OML_ERR_INVALID_INDEX, idx + 1, OML_VAR_INDEX);
        fieldindex.second = static_cast<int>((*cell)(1).Scalar());
    }
    return fieldindex;
}
//------------------------------------------------------------------------------
// Helper method for setfield which returns matrix after setting an element
//------------------------------------------------------------------------------
Currency BuiltInFuncsData::SetMatrixElement(EvaluatorInterface         eval,
                                            const Currency&            in,
                                            const Currency&            value,
                                            const std::pair<int, int>& index,
                                            int                        argidx) const
{
    assert (in.IsMatrixOrString());

    bool isScalar = value.IsScalar();
    if (!isScalar && !value.IsComplex())
        throw OML_Error(OML_ERR_SCALAR_COMPLEX, argidx, OML_VAR_TYPE);

    if (in.IsString())
    {
        if (!isScalar)
            throw OML_Error(OML_ERR_SCALAR, argidx, OML_VAR_TYPE);
        if (value.Scalar() < 0 || value.Scalar() > UCHAR_MAX)
            throw OML_Error(OML_ERR_SCALAROUTOFRANGE, argidx, OML_VAR_TYPE);
    }

    hwMatrix* mtx = EvaluatorInterface::allocateMatrix(in.Matrix());
    int oldm = mtx->M();
    int oldn = mtx->N();

    int m = index.first;
    int n = (index.second >= 1) ? index.second : 1;

    if (index.second < 1 && oldm == 1)
    {
        m = 1;
        n = index.first;
    }

    if (oldm < m || oldn < n)
    {
        int newm = std::max(m, oldm);
        int newn = std::max(n, oldn);

        hwMathStatus status = mtx->Resize(newm, newn, true);
        BuiltInFuncsUtils::CheckMathStatus(eval, status);
    }

    // Use SetElement as it will take care of flipping the matrix type from 
    // real to complex, based on existing data and new value added
    if (isScalar)
        mtx->SetElement(m - 1, n - 1, value.Scalar());
    else
        mtx->SetElement(m - 1, n - 1, value.Complex());

    Currency out(mtx);
    out.SetMask(in.GetMask());

    return out;
}
//------------------------------------------------------------------------------
// Helper method for setfield which grows cell and gets requested element
//------------------------------------------------------------------------------
Currency* BuiltInFuncsData::GetCellElement(EvaluatorInterface         eval,
                                           HML_CELLARRAY*             cellin,
                                           const std::pair<int, int>& index) const
{
    HML_CELLARRAY* cell(cellin);
    assert(cell);

    int m = index.first;
    int n = index.second;

    if (n < 0)  // Only one index is used
    {
        if (cell->Size() < m)
        {
            hwMathStatus status = cell->Resize(m, true);
            BuiltInFuncsUtils::CheckMathStatus(eval, status);
            (*cell)(m - 1) = EvaluatorInterface::allocateCellArray();
        }
        return &(*cell)(m - 1);
    }
    
    if (cell->M() < m || cell->N() < n)
    {
        hwMathStatus status = cell->Resize(m, n, true);
        BuiltInFuncsUtils::CheckMathStatus(eval, status);
        (*cell)(m - 1, n - 1) = EvaluatorInterface::allocateCellArray();
    }
    return &(*cell)(m - 1, n - 1);
}
//------------------------------------------------------------------------------
// Helper method for setfield which sets a struct value
//------------------------------------------------------------------------------
void BuiltInFuncsData::SetStructElement(const Currency&            value,
                                        const std::string&         field,
                                        const std::pair<int, int>& index,
                                        bool                       hasindex,
                                        Currency&                  cur) const
{
    assert(cur.Struct());

    StructData* sd = cur.Struct();
    assert(sd);

	if (sd->GetRefCount() != 1)
	{
		StructData* temp = new StructData(*sd);
		cur.ReplaceStruct(temp);
	}

    bool        isStructArray = (sd->M() > 1 || sd->N() > 1);
    std::string err (HW_ERROR_UNSUPOP);

    // If an array of structs is being set, field index needs to be specified
    if (!hasindex && isStructArray)
        throw OML_Error(err);

    bool isTargetStruct = value.IsStruct();
    if (field.empty() && !isTargetStruct)
        throw OML_Error(err + "; field name must specified to set a value in a struct");

    if (isTargetStruct)
    {
        StructData* out = value.Struct();
        if ((value.Struct()->M() > 1 || value.Struct()->N() > 1))
            throw OML_Error(err + "; cannot set struct array as an element");

        std::map< std::string, int> outnames = out->GetFieldNames();
        std::map< std::string, int> innames  = sd->GetFieldNames();

        if (outnames.size() != innames.size())
            throw OML_Error(err + "; struct field names must match");

        for (std::map< std::string, int>::const_iterator itr = innames.begin();
                itr != innames.end(); ++itr)
        {
            std::string name = itr->first;
            if (outnames.find(name) == outnames.end())
                throw OML_Error(err + "; struct field names must match");
        }
        for (std::map< std::string, int>::const_iterator itr = outnames.begin();
                itr != outnames.end(); ++itr)
        {
            std::string name = itr->first;
            if (innames.find(name) == innames.end())
                throw OML_Error(err + "; struct field names must match");
        }
    }

    int m = (index.first  > 0) ? index.first - 1 : 0;
    int n = (index.second > 0) ? index.second - 1 : -1;

    cur.Struct()->SetValue(m, n, field, value);
}
//------------------------------------------------------------------------------
// Helper method for setfield which grows matrix and gets requested element
//------------------------------------------------------------------------------
Currency BuiltInFuncsData::GetMatrixElement(EvaluatorInterface         eval,
                                            const Currency&            in,
                                            const std::pair<int, int>& index) const
{
    assert(in.IsScalar() || in.IsComplex());

    // Scalar/complex now needs to become a mtx, preferably with a single row
    int m = (index.second < 0) ? 1 : index.first;
    int n = (index.second < 0) ? index.first : index.second;

    hwMatrix*          mat  = EvaluatorInterface::allocateMatrix(m, n, in.IsScalar());

    mat->SetElements(0.0);

    if (in.IsScalar())
        mat->SetElement(0, 0, in.Scalar());
    else
        mat->SetElement(0, 0, in.Complex());

    return mat;
}
//------------------------------------------------------------------------------
// Helper method for setfield which sets matrix parent
//------------------------------------------------------------------------------
void BuiltInFuncsData::SetMatrixParent(const Currency&            mtx,
                                       const std::pair<int, int>& index,
                                       const std::string&         field,
                                       Currency*&                 parent) const
{
    assert(parent);
    if (!parent->IsCellArray() && !parent->IsStruct()) 
        return;

    int m = index.first > 0 ? index.first - 1 : 0;

    if (parent->IsStruct())
    {
        int n = (index.second > 0) ? index.second - 1 : -1;
        parent->Struct()->SetValue(m, n, field, mtx);
        return;
    }
	HML_CELLARRAY* parent_cells = parent->CellArray();
	if (parent_cells->GetRefCount() != 1)
	{
		Currency* parent_cur = parent;

		HML_CELLARRAY* cells = EvaluatorInterface::allocateCellArray(parent_cells);
		parent_cur->ReplaceCellArray(cells);
    }
    if (index.second < 0) 
        (*parent->CellArray())(m) = mtx;
    else
        (*parent->CellArray())(m, index.second - 1) = mtx;
}
//------------------------------------------------------------------------------
// Returns true after converting matrix to cell array [mat2cell command]
//------------------------------------------------------------------------------
bool BuiltInFuncsData::Mat2Cell(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
{
    size_t nargin = (!inputs.empty()) ? inputs.size() : 0;
    if (nargin < 1) 
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency cur = inputs[0];
    if (!cur.IsMatrixOrString() && !cur.IsScalar())
        throw OML_Error(OML_ERR_MTXSTRING, 1, OML_VAR_TYPE);

    const hwMatrix* data = cur.ConvertToMatrix();
    assert(data);
    if (!data)
    {
        outputs.push_back(EvaluatorInterface::allocateCellArray());
        return true;
    }

    BuiltInFuncsData funcs;
    // Get row dimensions for sub matrices
    int dataM = data->M();
    std::vector<int> v1;
    if (nargin > 1)
        v1 = funcs.GetDimensions(inputs[1], 2, dataM);
    else   
        v1.push_back(dataM);

    // Get column dimensions for submatrices
    int dataN = data->N();
    std::vector<int> v2;
    if (nargin > 2)
        v2 = funcs.GetDimensions(inputs[2], 3, dataN);
    else
        v2.push_back(dataN);

    // Create output cell   
    int            nrows = static_cast<int>(v1.size());
    int            ncols = static_cast<int>(v2.size());
    HML_CELLARRAY* cell  = EvaluatorInterface::allocateCellArray(nrows, ncols);
    assert(cell);

    int datarowidx = 0;
    int i          = 0;
    for (std::vector<int>::const_iterator itr1 = v1.begin(); 
         itr1 != v1.end(); ++itr1, ++i)
    {
        int rows       = *itr1;
        int j          = 0;
        int datacolidx = 0;

        for (std::vector<int>::const_iterator itr2 = v2.begin(); 
             itr2 != v2.end(); ++itr2, ++j)
        {
            int       cols = *itr2;
            hwMatrix* mtx  = EvaluatorInterface::allocateMatrix(rows, cols, data->IsReal());

            for (int n = 0; n < cols && datacolidx+n < dataN; ++n)
            {
				if (data->IsReal())
				{
					for (int m = 0; m < rows && datarowidx + m < dataM; ++m)
						(*mtx)(m,n) =  (*data)(datarowidx + m, datacolidx+n);
				}
				else
				{
					for (int m = 0; m < rows && datarowidx + m < dataM; ++m)
						mtx->z(m,n) =  data->z(datarowidx + m, datacolidx+n);
				}
            }
            Currency tmp(mtx);
            tmp.SetMask(cur.GetMask());
            (*cell)(i, j) = tmp;

            datacolidx += cols;
        }
        datarowidx += rows;
    }

    outputs.push_back(cell);
    return true;
}
//------------------------------------------------------------------------------
// Gets dimensions of sub-matrices from given vector
//------------------------------------------------------------------------------
std::vector<int> BuiltInFuncsData::GetDimensions(const Currency& cur,
                                                 int             index,
                                                 int             ref) const
{
    std::vector<int> out;
    bool             isint = cur.IsPositiveInteger();
    if (!isint && !cur.IsVector())
        throw OML_Error(OML_ERR_POSINTEGER_VEC, index, OML_VAR_TYPE);

    std::string msg = "Error: invalid input in argument " +
                      std::to_string(static_cast<long long>(index)) + "; ";
    int dims = 0;
    if (isint)
    {
        dims = static_cast<int>(cur.Scalar());
        if (dims != ref)
        {
            throw OML_Error(msg + "value must match matrix dimension ["
                + std::to_string(static_cast<long long>(dims)) + " != "
                + std::to_string(static_cast<long long>(ref)) + "]");
        }
        out.push_back(dims);
        return out;
    }


    std::vector<double> in = cur.Vector();
    if (in.empty())
    {
        out.push_back(ref);
        return out;
    }
    out.reserve(in.size());

    int i   = 1;
    for (std::vector<double>::const_iterator itr = in.begin(); 
            itr != in.end(); ++itr, ++i)
    {
        Currency tmp (*itr);
        bool     isInt = tmp.IsInteger();
        int      val   = isInt ? static_cast<int>(tmp.Scalar()) : 0;

        if (!isInt || (isInt && val < 0))
            throw OML_Error(msg + "value must be a nonnegative integer in index "
                + std::to_string(static_cast<long long>(i)));
        
        dims += val;
        out.push_back(val);
    }
    if (dims != ref)
    {
        throw OML_Error(msg + "sum of elements must match matrix dimension ["
            + std::to_string(static_cast<long long>(dims)) + " != "
            + std::to_string(static_cast<long long>(ref)) + "]");
    }

    return out;
}
//------------------------------------------------------------------------------
// Returns true if input is a row vector [isrow]
//------------------------------------------------------------------------------
bool BuiltInFuncsData::IsRow(EvaluatorInterface           eval,
                             const std::vector<Currency>& inputs,
                             std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    const Currency& cur   = inputs[0];
    bool            isrow = false;

    if (cur.IsScalar() || cur.IsComplex())
    {
        isrow = true;
    }
    else if (cur.IsMatrix())
    {
        const hwMatrix* mtx = cur.Matrix();
        if (mtx && mtx->M() == 1 && mtx->N() >= 0)
        {
            isrow = true;
        }
    }
    else if (cur.IsString() && !cur.IsMultilineString())
    {
        isrow = true;
    }
	else if (cur.IsStruct())
	{
		StructData* sd = cur.Struct();

		if (sd && sd->M() == 1 && sd->N() >= 0)
			isrow = true;
	}
	else if (cur.IsCellArray())
	{
		HML_CELLARRAY* cells = cur.CellArray();

		if (cells && cells->M() == 1 && cells->N() >= 0)
			isrow = true;
	}

    outputs.push_back(isrow);
    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is a column vector [iscol]
//------------------------------------------------------------------------------
bool BuiltInFuncsData::IsColumn(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    const Currency& cur   = inputs[0];
    bool            iscol = false;

    if (cur.IsScalar() || cur.IsComplex())
    {
        iscol = true;
    }
    else if (cur.IsMatrixOrString())
    {
        const hwMatrix* mtx = cur.Matrix();
        if (mtx && mtx->M() >= 0 && mtx->N() == 1)
        {
            iscol = true;
        }
    }
	else if (cur.IsStruct())
	{
		StructData* sd = cur.Struct();

		if (sd && sd->N() == 1 && sd->M() >= 0)
			iscol = true;
	}
	else if (cur.IsCellArray())
	{
		HML_CELLARRAY* cells = cur.CellArray();

		if (cells && cells->N() == 1 && cells->M() >= 0)
			iscol = true;
	}

    outputs.push_back(iscol);
    return true;
}

bool Num2CellNDHelper(const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (!inputs.size())
		throw OML_Error(OML_ERR_NUMARGIN);

	int dim = -1;

	if (inputs.size() == 2)
	{
		if (!inputs[1].IsPositiveInteger())
		{
			throw OML_Error(OML_ERR_POSINTEGER, 2);
		}

		dim = static_cast<int>(inputs[1].Scalar());
	}

	const hwMatrixN* mat  = inputs[0].MatrixN();  
	std::vector<int> dims = mat->Dimensions();

	if (dim == -1)
	{
		HML_ND_CELLARRAY* cells = new HML_ND_CELLARRAY(dims, HML_ND_CELLARRAY::REAL);

		for (int j = 0; j < mat->Size(); j++)
		{
			if (mat->IsReal())
				(*cells)(j) = (*mat)(j);
			else
				(*cells)(j) = mat->z(j);
		}

		outputs.push_back(cells);
	}
	else
	{
		std::vector<int> new_dims;
		int              interior_size;

		for (int j = 0; j < dims.size(); ++j)
		{
			if (dim == (j + 1))
			{
				new_dims.push_back(1);
				interior_size = dims[j];
			}
			else
			{
				new_dims.push_back(dims[j]);
			}
		}

		HML_ND_CELLARRAY* cells = new HML_ND_CELLARRAY(new_dims, HML_ND_CELLARRAY::REAL);

		for (int j = 0; j < cells->Size(); ++j)
		{
			std::vector<int> lhs_dims = cells->IndexVector(j);
			std::vector<hwSliceArg> lhs_slices;
			std::vector<hwSliceArg> rhs_slices;

			for (int k = 0; k < lhs_dims.size(); ++k)
			{
				if (dim == (k + 1))
				{
					rhs_slices.push_back(hwSliceArg());
					lhs_slices.push_back(0);
				}
				else
				{
					rhs_slices.push_back(lhs_dims[k]);
					lhs_slices.push_back(lhs_dims[k]);
				}
			}

			hwMatrixN* temp_slice = EvaluatorInterface::allocateMatrixN();
			mat->SliceRHS(rhs_slices, *temp_slice);		
			cells->SliceLHS(lhs_slices, Currency(temp_slice));
		}

		outputs.push_back(cells);
	}

	return true;
}
//------------------------------------------------------------------------------
// Returns true after converting number/matrix to cell array [num2cell command]
//------------------------------------------------------------------------------
bool BuiltInFuncsData::Num2Cell(EvaluatorInterface           eval,
	                            const std::vector<Currency>& inputs,
	                            std::vector<Currency>&       outputs)
{
	size_t nargin = (!inputs.empty()) ? inputs.size() : 0;
	if (nargin < 1)
	{
		throw OML_Error(OML_ERR_NUMARGIN);
	}

    if (inputs[0].IsCellArray())
    {
        outputs.push_back(inputs[0]);
        return true;
    }

	if (inputs[0].IsScalar())
	{
		Currency temp = inputs[0];
		temp.ConvertToCellArray();
		outputs.push_back(temp);
		return true;
	}

	if (inputs[0].IsNDMatrix())
		return Num2CellNDHelper(inputs, outputs);

	if (!inputs[0].IsMatrixOrString())
	{
		throw OML_Error(OML_ERR_STRSCALARCOMPLEXMTX, 1);
	}
    const hwMatrix* mtx = inputs[0].Matrix();
    if (!mtx || mtx->Size() == 0)
    {
        outputs.push_back(EvaluatorInterface::allocateCellArray());
        return true;
    }


	int dim = -1;
	if (nargin > 1)
	{
		if (!inputs[1].IsPositiveInteger())
		{
			throw OML_Error(OML_ERR_POSINTEGER, 2);
		}
		dim = static_cast<int>(inputs[1].Scalar());
		if (dim > 2)
		{
			dim = -1;  // Ignore dimensions greater than 2
		}
	}

    bool isReal = mtx->IsReal();
    bool setStringMask = (inputs[0].IsString());

    int rows  = mtx->M();
    int cols  = mtx->N();
    int cellm = (dim == 1) ? 1 : rows;
    int celln = (dim == 2) ? 1 : cols;

    std::unique_ptr<HML_CELLARRAY> cell(
        EvaluatorInterface::allocateCellArray(cellm, celln));

	if (dim < 1)  // No restrictions on dimensions
	{
		for (int i = 0; i < rows; ++i)
		{
			for (int j = 0; j < cols; ++j)
			{
				if (isReal)
				{
					if (!setStringMask)
					{
						(*cell)(i, j) = (*mtx)(i, j);
					}
					else
					{
						std::string strval;
						strval += char((*mtx)(i, j));
						(*cell)(i, j) = strval;  // Setting mask does not work
					}
				}
				else
				{
					(*cell)(i, j) = mtx->z(i, j);
				}
			}
		}

		outputs.push_back(cell.release());
		return true;
	}

	// Restricted dimensions
    int m = (dim == 1) ? rows : 1;
    int n = (dim == 1) ? 1 : cols;

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            std::unique_ptr<hwMatrix> val = nullptr;

            int valindex  = (dim == 1) ? i : j;
            int cellindex = (dim == 1) ? j : i;

            if (valindex == 0)    // Create all the sub matrices
            {
                val.reset(EvaluatorInterface::allocateMatrix(m, n, isReal));
            }
            else
            {
                val.reset(EvaluatorInterface::allocateMatrix(
                    (*cell)(cellindex).Matrix()));
            }
            assert(val);
            if (isReal)
            {
                (*val)(valindex) = (*mtx)(i, j);
            }
            else
            {
                val->z(valindex) = mtx->z(i, j);
            }

            (*cell)(cellindex) = val.release();
            if (setStringMask)
            {
                Currency element = (*cell)(cellindex);
                element.SetMask(Currency::MASK_STRING);
                (*cell)(cellindex) = element;
            }

        }
    }
    outputs.push_back(cell.release());
    return true;
}
//------------------------------------------------------------------------------
// Returns true after sorting rows in the matrix [sortrows]
//------------------------------------------------------------------------------
bool BuiltInFuncsData::Sortrows(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsMatrixOrString())
    {
        throw OML_Error(OML_ERR_MTXSTRING, 1);
    }

    int nargout = eval.GetNargoutValue();

    const hwMatrix* data = inputs[0].Matrix();
    if (!data || data->Size() == 0 || data->IsEmpty())
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
        if (nargout > 1)
        {
            outputs.push_back(EvaluatorInterface::allocateMatrix());
        }
        return true;
    }

    int m = data->M();
    int n = data->N();

    std::unique_ptr<hwMatrix> dims(EvaluatorInterface::allocateMatrix(1, 1, true));
    (*dims)(0) = 1;  // Defaults to first column, 1-based index

    if (inputs.size() > 1)
    {
        if (inputs[1].IsInteger())
        {
            (*dims)(0) = static_cast<int>(inputs[1].Scalar());
        }
        else if (inputs[1].IsVector())
        {
            dims.reset(EvaluatorInterface::allocateMatrix(inputs[1].Matrix()));
            if (!dims->IsReal())
            {
                throw OML_Error(OML_ERR_REALVECTOR, 2);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARVECTOR, 2);
        }
    }

    const hwMatrix* m1 = inputs[0].Matrix();

    // Go through the dims and sort on the columns mentioned
    std::unique_ptr<hwMatrix> mtx(EvaluatorInterface::allocateMatrix(m1));
    std::unique_ptr<hwMatrix> idx = nullptr;
    if (nargout > 1)
    {
        idx.reset(EvaluatorInterface::allocateMatrix(m, 1, true));
    }

    BuiltInFuncsUtils utils;
    bool isreal = mtx->IsReal();
    int dimssize = dims->Size();

    for (int i = dimssize - 1; i >= 0; --i)
    {
        double val = (*dims)(i);
        if (IsNaN_T(val) || IsNegInf_T(val) || IsInf_T(val))
        {
            throw OML_Error(OML_ERR_FINITE, 2);
        }

        int col = static_cast<int>(val);
        if (col == 0 || col > n)
        {
            throw OML_Error(OML_ERR_INVALIDINDEX, 2);
        }

        bool descend = (col < 0);
        col = abs(col) - 1;

        std::unique_ptr<hwMatrix> out(EvaluatorInterface::allocateMatrix(m, n,
                                      mtx->IsReal()));
        if (isreal)
        {
            std::vector< std::pair<double, int> > tmp;
            tmp.reserve(m);
            for (int j = 0; j < m; ++j)
            {
                tmp.push_back(std::pair<double, int>((*mtx)(j, col), j));
            }
            if (!descend)
            {
                std::stable_sort(tmp.begin(), tmp.end(), &IsLessThan);
            }
            else
            {
                std::stable_sort(tmp.begin(), tmp.end(), &IsGreaterThan);
            }
            for (int j = 0; j < m; ++j)
            {
                int sortedidx = tmp[j].second;

                std::unique_ptr<hwMatrix> row (
                    EvaluatorInterface::allocateMatrix(1, n, mtx->IsReal()));
                utils.CheckMathStatus(eval, mtx->ReadRow(sortedidx, *row));
                utils.CheckMathStatus(eval, out->WriteRow(j, *row));
                if (idx && i == 0)
                {
                    for (int k = 0; k < m; ++k)
                    {
                        std::unique_ptr<hwMatrix> r1(
                            EvaluatorInterface::allocateMatrix(1, n, true));
                        utils.CheckMathStatus(eval, m1->ReadRow(k, *r1));
                        if (row->IsEqual(*r1))
                        {
                            (*idx)(j) = k + 1;
                        }
                    }
                }
            }
        }
        else
        {
            std::vector< std::pair< hwComplex, int> > tmp;
            tmp.reserve(m);
            for (int j = 0; j < m; ++j)
            {
                tmp.push_back(std::pair< hwComplex, int>(mtx->z(j, col), j));
            }
            if (!descend)
            {
                std::sort(tmp.begin(), tmp.end(), &IsLessThanComplex);
            }
            else
            {
                std::sort(tmp.begin(), tmp.end(), &IsGreaterThanComplex);
            }
            for (int j = 0; j < m; ++j)
            {
                int sortedidx = tmp[j].second;

                std::unique_ptr<hwMatrix> row(
                    EvaluatorInterface::allocateMatrix(1, n, mtx->IsReal()));
                utils.CheckMathStatus(eval, mtx->ReadRow(sortedidx, *row));
                utils.CheckMathStatus(eval, out->WriteRow(j, *row));

                if (idx && i == 0)
                {
                    for (int k = 0; k < m; ++k)
                    {
                        std::unique_ptr<hwMatrix> r1(
                            EvaluatorInterface::allocateMatrix(1, n, m1->IsReal()));
                        utils.CheckMathStatus(eval, m1->ReadRow(k, *r1));
                        if (row->IsEqual(*r1))
                        {
                            (*idx)(j) = k + 1;
                        }
                    }
                }
            }
        }
        mtx.reset(out.release());
    }

    Currency cur = mtx.release();
    cur.SetMask(inputs[0].GetMask());
    outputs.push_back(cur);

    if (nargout > 1)
    {
        outputs.push_back(idx.release());
    }
    return true;
}
