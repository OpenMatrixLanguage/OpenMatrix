/**
* @file BuiltInFuncsData.h
* @date June 2016
* Copyright (C) 2016-2024 Altair Engineering, Inc.  
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
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //!
    static bool Mat2Cell(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
	//!
	//! Returns true after converting number/matrix to cell array [num2cell]
	//! \param Evaluator interface
	//! \param Vector of inputs
	//! \param Vector of outputs
	//!
	static bool Num2Cell(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Sets fields recursively. First currency is the input, last currency is value
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //!
    static bool Setfield(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Returns true if input is a row vector [isrow]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //!
    static bool IsRow(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Returns true if input is a column vector [iscolumn]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //!
    static bool IsColumn(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Sorts rows in a matrix [sortrows]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //!
    static bool Sortrows(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Outputs an array indicating which fields are empty in a struct [fieldempty]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //!
    static bool FieldEmpty(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Outputs a cell which contains the contents of the given field [fields2cell]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //!
    static bool Fields2Cell(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Converts cell array to struct of specified fields [cell2fields]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //!
    static bool Cell2Fields(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Concatenates structures along a specified dimension [structcat]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //!
    static bool Structcat(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);

private:
    //!
    //! Constructor
    //!
    BuiltInFuncsData() {}
    //!
    //! Helper method for Setfield to get indices from cell array
    //! \param Input currency  
    //! \param Input currency index for error messages
    //!
    std::pair<int, int> GetFieldIndex(const Currency&, int) const;
    //!
    //! Helper method for setfield which grows cell and gets requested element
    //! \param Evaluator
    //! \param Cell
    //! \param Index to the cell element
    //!
    Currency* GetCellElement(EvaluatorInterface, HML_CELLARRAY*, const std::pair<int, int>&) const;
    //!
    //! Helper method for setfield which grows matrix and gets requested element
    //! \param Evaluator
    //! \param Given input currency
    //! \param Index to the matrix element
    //!
    Currency GetMatrixElement(EvaluatorInterface, const Currency&, const std::pair<int, int>&) const;
    //!
    //! Helper method for setfield which returns matrix after setting an element
    //! \param Evaluator interface
    //! \param Matrix to set
    //! \param Element value
    //! \param Element index
    //! \param Argument index - for error handling
    //!
    Currency SetMatrixElement(EvaluatorInterface, const Currency&, const Currency&, const std::pair<int, int>&, int) const;
    //!
    //! Helper method for setfield which sets matrix parent
    //! \param Matrix
    //! \param Index in parent
    //! \param Field name in parent, if applicable
    //! \param Matrix parent
    //!
    void SetMatrixParent(const Currency&, const std::pair<int, int>&, const std::string&, Currency*&) const;
    //!
    //! Helper method for setfield which grows struct and gets requested element
    //! \param Struct
    //! \param Field name
    //! \param Index to the struct element
    //!
    Currency* GetStructElement(StructData*, const std::string&, const std::pair<int, int>&) const;
    //!
    //! Helper method for setfield which sets a struct value
    //! \param Value to set
    //! \param Field name
    //! \param Field index
    //! \param True if field index was explicitly set
    //! \param Currency to set
    //!
    void SetStructElement(const Currency&, const std::string&, const std::pair<int, int>&, bool, Currency&) const;
    //!
    //! Gets dimensions of sub-matrices from given vector
    //! \param Given currency
    //! \param Input index for error messages
    //! \param True if row dimensions are being set
    //!
    std::vector<int> GetDimensions(const Currency&, int, int) const;
    //!
    //! Returns 1 if struct field at given index is empty, returns 0 otherwise
    //! \param Struct
    //! \param Field
    //! \param Row
    //! \param Column
    //! 
    static double IsStructFieldEmpty(const StructData*, const std::string&, int, int);
};

#endif // __BUILTINFUNCSDATA__


