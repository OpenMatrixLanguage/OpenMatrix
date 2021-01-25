/**
* @file BuiltInFuncsData.h
* @date June 2016
* Copyright (C) 2016-2019 Altair Engineering, Inc.  
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

#ifndef __BUILTINFUNCSDATA__
#define __BUILTINFUNCSDATA__

// Begin defines/includes
#include <utility>

#include "EvaluatorInt.h"

// End defines/includes

//------------------------------------------------------------------------------
//!
//! \brief Class for built-in functions implementing data structures commands
//!
//------------------------------------------------------------------------------
class OMLDLL_DECLS BuiltInFuncsData
{
public:
    //!
    //! Destructor
    //!
    ~BuiltInFuncsData() {}
    //!
    //! Returns true after converting matrix to cell array [mat2cell command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Mat2Cell( EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs);
	//!
	//! Returns true after converting number/matrix to cell array [num2cell]
	//! \param eval    Evaluator interface
	//! \param inputs  Vector of inputs
	//! \param outputs Vector of outputs
	//!
	static bool Num2Cell(EvaluatorInterface           eval,
		                 const std::vector<Currency>& inputs,
		                 std::vector<Currency>&       outputs);
    //!
    //! Sets fields recursively. First currency is the input, last currency is value
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Setfield( EvaluatorInterface           eval,
                          const std::vector<Currency>& inputs,
                          std::vector<Currency>&       outputs);
    //!
    //! Returns true if input is a row vector [isrow]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool IsRow(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns true if input is a column vector [iscolumn]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool IsColumn(EvaluatorInterface           eval,
                         const std::vector<Currency>& inputs,
                         std::vector<Currency>&       outputs);
    //!
    //! Returns true after sorting rows in the matrix [sortrows]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Sortrows(EvaluatorInterface           eval,
                         const std::vector<Currency>& inputs,
                         std::vector<Currency>&       outputs);

private:
    //!
    //! Constructor
    //!
    BuiltInFuncsData() {}
    //!
    //! Helper method for Setfield to get indices from cell array
    //! \param in       Input currency  
    //! \param inputIdx Input currency index - for error messages
    //!
    std::pair<int, int> GetFieldIndex(const Currency& in,
                                      int             idx) const;
    //!
    //! Helper method for setfield which grows cell and gets requested element
    //! \param eval  Evaluator
    //! \param in    Cell
    //! \param index Index to the cell element
    //!
    Currency* GetCellElement( EvaluatorInterface         eval,
                              HML_CELLARRAY*             cell,
                              const std::pair<int, int>& index) const;
    //!
    //! Helper method for setfield which grows matrix and gets requested element
    //! \param eval  Evaluator
    //! \param in    Given input currency
    //! \param index Index to the matrix element
    //!
    Currency GetMatrixElement( EvaluatorInterface         eval,
                               const Currency&            in,
                               const std::pair<int, int>& index) const;
    //!
    //! Helper method for setfield which returns matrix after setting an element
    //! \param eval   Evaluator interface
    //! \param in     Matrix to set
    //! \param value  Element value
    //! \param index  Element index
    //! \param argidx Argument index - for error handling
    //!
    Currency SetMatrixElement(EvaluatorInterface         eval,
                              const Currency&            in,
                              const Currency&            value,
                              const std::pair<int, int>& index,
                              int                        argidx) const;
    //!
    //! Helper method for setfield which sets matrix parent
    //! \param mtx    Matrix
    //! \param index  Index in parent
    //! \param field  Field name in parent, if applicable
    //! \param parent Matrix parent
    //!
    void SetMatrixParent(const Currency&            mtx,
                         const std::pair<int, int>& index,
                         const std::string&         field,
                         Currency*&                 parent) const;
    //!
    //! Helper method for setfield which grows struct and gets requested element
    //! \param in        Struct
    //! \param fieldname Field name
    //! \param index     Index to the struct element
    //!
    Currency* GetStructElement( StructData*                in,
                                const std::string&         field,
                                const std::pair<int, int>& index) const;
    //!
    //! Helper method for setfield which sets a struct value
    //! \param value    Value to set
    //! \param field    Field name
    //! \param index    Field index
    //! \param hasindex True if field index was explicitly set
    //! \param cur      Currency to set
    //!
    void SetStructElement(const Currency&            value,
                          const std::string&         field,
                          const std::pair<int, int>& index,
                          bool                       hasindex,
                          Currency&                  cur) const;
    //!
    //! Gets dimensions of sub-matrices from given vector
    //! \param cur   Given currency
    //! \param index Input index for error messages
    //! \param ref   True if row dimensions are being set
    //!
    std::vector<int> GetDimensions(const Currency& cur,
                                   int             index,
                                   int             ref) const;
};

#endif // __BUILTINFUNCSDATA__


