/**
* @file SparseDisplay.h
* @date June, 2019
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
#ifndef __SPARSEDISPLAY_H__
#define __SPARSEDISPLAY_H__

// Begin defines/includes
#include "Hml2Dll.h"

#include "CurrencyDisplay.h"
#include "DisplayFormatVars.h"

class OutputFormat;

// End defines/includes

//------------------------------------------------------------------------------
//!
//! \class SparseDisplay
//! \brief Displays formatted sparse output
//!
//------------------------------------------------------------------------------
class HML2DLL_DECLS SparseDisplay : public CurrencyDisplay
{
public:
    friend class Currency; //! Only currency is allowed to construct

    //!
    //! Destructor 
    //!
    virtual ~SparseDisplay() {}

    //!
    //! Gets number of rows and cols in given currency
    //! \param rows Number of rows
    //! \param cols Number of columns
    //!
    virtual void GetCurrencySize(int& rows,
                                 int& cols) const;
    //!
    //! Gets values as a string
    //! \param fmt Output format
    //!
    virtual std::string GetValues(const OutputFormat* fmt) const;

    //!
    //! Initialize
    //! \param fmt    Output format   
    //! \param interp Interpreter
    //! \param parent Parent display
    //!
    virtual void Initialize(const OutputFormat* fmt,
                            Interpreter*        interp,
                            CurrencyDisplay*    parent = 0);
protected:
    //!
    //! Sets data for forward pagination
    //!
    virtual void SetForwardDisplayData();
    //!
    //! Sets data for back pagination
    //!
    virtual void SetBackDisplayData();
    ////!
    ////! Gets values as a string
    ////! \param fmt Output format
    ////!

private:
    mutable int _realwidth;                //!< Widest real value width
    mutable int _imagwidth;                //!< Widest imag value width
    mutable DisplayFormatVars _formatvars; //!< Format variables

    //!
    //! Constructor - Only currency is allowed to construct
    //! \param cur Currency associated with this display
    //!
    SparseDisplay(const Currency& cur);

    // Stubbed out default, copy constructors and assignment operator
    SparseDisplay();
    SparseDisplay(const SparseDisplay& src);
    SparseDisplay& operator=(const SparseDisplay& src);

    //!
    //! Gets output
    //! \param fmt Output format
    //! \param os  Output stream
    //!
    virtual std::string GetOutput(const OutputFormat* fmt,
                                  std::ostringstream& os) const;
    //!
    //! Gets output
    //! \param index      Index
    //! \param isreal     True if this is a real mtx
    //! \param isrealdata True if this is a complex mtx with no imaginary parts
    //! \param output     Output
    //!
    void GetOutput(int          index,
                   bool         isreal,
                   bool         isrealdata,
                   std::string& output) const;
    //!
    //! Gets data with no pagination
    //! \param fmt              Format
    //! \param os  Output stream
    //!
    std::string GetOutputNoPagination(const OutputFormat* fmt) const;
    //!
    //! Gets outputfor forward/down pagination
    //! \param fmt Format
    //! \param os  Output stream
    //!
    std::string GetOutputForwardPagination(const OutputFormat* fmt) const;
    //!
    //! Gets outputfor back/up pagination
    //! \param fmt Format
    //!
    std::string GetOutputBackPagination(const OutputFormat* fmt) const;


    //!
    //! Returns true if display was paginating
    //!
    virtual bool IsPaginating() const;
    //!
    //! Returns true if display was paginating
    //!
    virtual bool WasPaginating() const;
    //!
    //! Gets pagination info for printing
    //! \param rows Number of rows
    //!
    std::string GetPaginationHeader(int rows) const;

    //!
    //! Sets matrix format
    //! \param fmt    Format
    //! \param interp Interpreter
    //!
    virtual void SetFormat(const OutputFormat* fmt,
                           Interpreter*        interp);
    //!
    //! Scan matrix and set the width of different columns
    //! \param interp Interpreter
    //!
    virtual void SetWidth(Interpreter* interp);
};


#endif

// End of file:


