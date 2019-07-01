/**
* @file CellNDisplay.h
* @date January 2019
* Copyright (C) 2019 Altair Engineering, Inc.
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
#ifndef __CELLNDISPLAY_H__
#define __CELLNDISPLAY_H__

// Begin defines/includes
#include "Currency.h"
#include "CurrencyDisplay.h"
#include "hwMatrixN.h"

class Interpreter;
class OutputFormat;
// End defines/includes

//------------------------------------------------------------------------------
//! \class CellNDisplay
//! \brief Displays ND cell arrays
//------------------------------------------------------------------------------
class HML2DLL_DECLS CellNDisplay : public CurrencyDisplay
{
public:
	friend class Currency; //!< Only currency is allowed to construct

	//!
	//! Destructor 
	//!
	virtual ~CellNDisplay() {}

	//!
	//! Gets output - called from Currency::GetOutputString
	//! \param fmt Output format
	//! \param os  Output stream
	//!
	std::string GetOutput(const OutputFormat* fmt,
		                  std::ostringstream& os) const;
	//!
	//! Gets number of rows and cols in given currency
	//! \param rows Number of rows
	//! \param cols Number of columns
	//!
	virtual void GetCurrencySize(int& rows,
		                         int& cols) const;
	//!
	//! Returns true if parent is ND cell array
	//!
	virtual bool IsNDCellDisplay() const { return true; }

	//!
	//! Helper to slice a given ND cell array to 2D cells
	//! \param cell    Given cell array
	//! \param curs   Vector of currencies of 2D component cell arrays
	//! \param labels Vector of slice labels
	//!
	static void GetSlices(HML_ND_CELLARRAY*         cell,
		                  std::vector<Currency>&    curs,
		                  std::vector<std::string>& labels);

	//!
	//! Gets values as a string
	//! \param fmt Output format
	//!
	virtual std::string GetValues(const OutputFormat* fmt) const;

private:
	//!
	//! Constructor - Only currency is allowed to construct
	//! \param cur Currency associated with this display
	//!
	CellNDisplay(const Currency& cur);

	CellNDisplay();                                    // Stubbed out 
	CellNDisplay(const CellNDisplay& src);             // Stubbed out 
	CellNDisplay& operator=(const CellNDisplay& src);  // Stubbed out

	//!
	//! Gets output with no pagination
	//! \param fmt Format
	//!
	std::string GetOutputNoPagination(const OutputFormat* fmt) const;
	//!
	//! Helper to slice a given ND cell array to 2D cells
	//! \param cell   Given cell
	//! \param slices Slices
	//! \param curs   Vector of currencies of 2D component cells
	//! \param labels Vector of slice labels
	//
	static void GetSlicesHelper(HML_ND_CELLARRAY*         cell,
								std::vector<hwSliceArg>   slices,
								std::vector<Currency>&    curs,
								std::vector<std::string>& sliceLabels);
};
#endif


