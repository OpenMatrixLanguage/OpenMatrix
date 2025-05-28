/**
* @file CellDisplay.h
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

#ifndef __CELLDISPLAY_H__
#define __CELLDISPLAY_H__

// Begin defines/includes
#include "CurrencyDisplay.h"

class Interpreter;
class OutputFormat;
// End defines/includes

//------------------------------------------------------------------------------
//! \class CellDisplay
//! \brief Displays cell arrays
//------------------------------------------------------------------------------
class OMLDLL_DECLS CellDisplay : public CurrencyDisplay
{
public:
    friend class Currency; // Only currency is allowed to construct

    //!
    //! Destructor 
    //!
    virtual ~CellDisplay() {}

    //!
    //! Gets output - called from Currency::GetOutputString
    //! \param Output format
    //! \param Output stream
    //!
    virtual std::string GetOutput(const OutputFormat*, std::ostringstream&) const override;
    //!
    //! Returns true if end of pagination message needs to be printed
    //! \param Additional message that needs to be printed
    //!
    virtual bool GetPaginationEndMsg(std::string&) const override;
    //!
    //! Gets number of rows and cols in given currency
    //! \param Given cell array/matrix
    //! \param Number of rows
    //! \param Number of columns
    //!
    virtual void GetCurrencySize(int&, int&) const override;

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
    //! Gets values as a string
    //! \param Output format
    //!
    virtual std::string GetValues(const OutputFormat*) const override;

private:
    //!
    //! Constructor - Only currency is allowed to construct
    //! \param Currency associated with this display
    //!
    CellDisplay(const Currency&); // cppcheck-suppress noExplicitConstructor

    CellDisplay();                                // Stubbed out
    CellDisplay(const CellDisplay&) ;             // Stubbed out 
    CellDisplay& operator=(const CellDisplay&);   // Stubbed out

    //!
    //! Gets output with no pagination
    //! \param Format
    //! \param Output stream
    //!
    std::string GetOutputNoPagination(const OutputFormat*, std::ostringstream& os) const;
    //!
    //! Gets outputfor forward/down pagination
    //! \param Format
    //! \param Output stream
    //!
    std::string GetOutputForwardPagination(const OutputFormat*, std::ostringstream&) const;
    //!
    //! Gets outputfor forward/down pagination for cell list
    //! \param Format
    //!
    std::string GetCellListOutputForwardPagination(const OutputFormat*) const;
    //!
    //! Gets outputfor back/up pagination
    //! \param Format
    //!
    std::string GetOutputBackPagination(const OutputFormat*) const;
    //!
    //! Gets output for back/up pagination for cell list
    //! \param Format
    //!
    std::string GetCellListOutputBackPagination(const OutputFormat*) const;

    //!
    //! Updates the number of rows to fit
    //!
    void UpdateNumLinesPrinted();
};
#endif

// End of file:

