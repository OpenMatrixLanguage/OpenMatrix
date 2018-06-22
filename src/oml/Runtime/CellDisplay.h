/**
* @file CellDisplay.h
* @date June 2016
* Copyright (C) 2016-2018 Altair Engineering, Inc.  
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
class HML2DLL_DECLS CellDisplay : public CurrencyDisplay
{
public:
    friend class Currency; // Only currency is allowed to construct

    //!
    //! Destructor 
    //!
    virtual ~CellDisplay() {}

    //!
    //! Gets output - called from Currency::GetOutputString
    //! \param fmt Output format
    //! \param os  Output stream
    //!
    std::string GetOutput(const OutputFormat* fmt,
                          std::ostringstream& os) const;
    //!
    //! Returns true if end of pagination message needs to be printed
    //! \param msg Additional message that needs to be printed
    //!
    virtual bool GetPaginationEndMsg(std::string& msg) const;
    //!
    //! Gets number of rows and cols in given currency
    //! \param cur   Given cell array/matrix
    //! \param rows Number of rows
    //! \param cols Number of columns
    //!
    virtual void GetCurrencySize(int& rows, 
                                 int& cols) const;
protected:
    //!
    //! Sets data for forward pagination
    //!
    virtual void SetForwardDisplayData();
    //!
    //! Sets data for back pagination
    //!
    virtual void SetBackDisplayData();
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
    CellDisplay(const Currency& cur);

    CellDisplay();                                    // Stubbed out 
    CellDisplay(const CellDisplay& src) ;             // Stubbed out 
    CellDisplay& operator=(const CellDisplay& src);   // Stubbed out

    //!
    //! Gets output with no pagination
    //! \param fmt Format
    //! \param os  Output stream
    //!
    std::string GetOutputNoPagination(const OutputFormat* fmt,
                                      std::ostringstream& os) const;
    //!
    //! Gets outputfor forward/down pagination
    //! \param fmt Format
    //! \param os  Output stream
    //!
    std::string GetOutputForwardPagination(const OutputFormat* fmt,
                                           std::ostringstream& os) const;
    //!
    //! Gets outputfor back/up pagination
    //! \param fmt Format
    //!
    std::string GetOutputBackPagination(const OutputFormat* fmt) const;
    //!
    //! Updates the number of rows to fit
    //!
    void UpdateNumLinesPrinted();
};
#endif

// End of file:

