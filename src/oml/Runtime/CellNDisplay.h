/**
* @file CellNDisplay.h
* @date January 2019
* Copyright (C) 2019-2024 Altair Engineering, Inc.
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
class OMLDLL_DECLS CellNDisplay : public CurrencyDisplay
{
public:
	friend class Currency; //!< Only currency is allowed to construct

	//!
	//! Destructor 
	//!
	virtual ~CellNDisplay() {}

	//!
	//! Gets output - called from Currency::GetOutputString
	//! \param Output format
	//! \param Output stream
	//!
	virtual std::string GetOutput(const OutputFormat*, std::ostringstream&) const override; // cppcheck-suppress missingOverride
	//!
	//! Gets number of rows and cols in given currency
	//! \param Number of rows
	//! \param Number of columns
	//!
	virtual void GetCurrencySize(int&, int&) const override;
	//!
	//! Returns true if parent is ND cell array
	//!
	virtual bool IsNDCellDisplay() const override { return true; }

	//!
	//! Helper to slice a given ND cell array to 2D cells
	//! \param Given cell array
	//! \param Vector of currencies of 2D component cell arrays
	//! \param Vector of slice labels
	//!
	static void GetSlices(HML_ND_CELLARRAY*, std::vector<Currency>&, std::vector<std::string>&);

protected:
    //!
    //! Sets data for forward pagination
    //!
    virtual void SetForwardDisplayData() override;
    //!
    //! Sets data for back pagination
    //!
    virtual void SetBackDisplayData() override;

    //!
    //! Gets number of rows that can be fit
    //!
    virtual int GetNumRowsToFit() const override;

    //!
    //! Gets values as a string
    //! \param Output format
    //!
    virtual std::string GetValues(const OutputFormat*) const override;
    //!
    //! True if rows are being processed during pagination
    //!
    virtual bool IsPaginatingRows() const override;
	//!
	//! True if paginating
	//!
	virtual bool IsPaginating() const { return false; }

private:
    mutable bool _childPaginating;     //!< True if child currency is paginating
	mutable int  _lastRowPrinted;      //!< Last row that was printed in interactive pagination
	//!
	//! Constructor - Only currency is allowed to construct
	//! \param Currency associated with this display
	//!
	CellNDisplay(const Currency&);

	CellNDisplay();                                    // Stubbed out  // cppcheck-suppress noExplicitConstructor
	CellNDisplay(const CellNDisplay& src);             // Stubbed out 
	CellNDisplay& operator=(const CellNDisplay& src);  // Stubbed out

	//!
	//! Gets output with no pagination
	//! \param Format
	//!
	std::string GetOutputNoPagination(const OutputFormat*) const;
    //!
    //! Gets outputfor forward/down pagination
    //! \param Format
    //!
    std::string GetOutputForwardPagination(const OutputFormat*) const;
    //!
    //! Gets outputfor back/up pagination
    //! \param Format
    //!
    std::string GetOutputBackPagination(const OutputFormat*) const;
    //!
    //! Updates the number of rows to fit
    //!
    void UpdateNumLinesPrinted();
    //!
    //! True if this ND matrix has only empty slices
    //! \param Slices
    //!
    bool HasOnlyEmptySlices(const std::vector<Currency>&) const;
    //!
    //! Gets display string if matrix ND has only empty slices
    //! 
    std::string GetOutputEmpty() const;

	//!
	//! Helper to slice a given ND cell array to 2D cells
	//! \param Given cell
	//! \param Slices
	//! \param Vector of currencies of 2D component cells
	//! \param Vector of slice labels
	//
	static void GetSlicesHelper(HML_ND_CELLARRAY*, std::vector<hwSliceArg>, std::vector<Currency>&, std::vector<std::string>&);
};
#endif


