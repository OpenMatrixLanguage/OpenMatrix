/**
* @file StringDisplay.h
* @date February 2018
* Copyright (C) 2018 Altair Engineering, Inc.  
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

#ifndef __STRINGDISPLAY_H__
#define __STRINGDISPLAY_H__

// Begin defines/includes
#include "CurrencyDisplay.h"

class Interpreter;
class OutputFormat;

// End defines/includes

//------------------------------------------------------------------------------
//!
//! \class StringDisplay
//! \brief Displays formatted multiline string output
//!
//------------------------------------------------------------------------------
class HML2DLL_DECLS StringDisplay : public CurrencyDisplay
{
public:
    friend class Currency; //! Only currency is allowed to construct

    //!
    //! Destructor 
    //!
    virtual ~StringDisplay() {}

    //!
    //! Gets output - called from Currency::GetOutputString
    //! \param fmt Output format
    //! \param os  Output stream
    std::string GetOutput( const OutputFormat* fmt,
                           std::ostringstream& os) const;

    //!
    //! Gets number of rows and cols in given currency -  not used
    //! \param cur  Given cell array/matrix
    //! \param rows Number of rows
    //! \param cols Number of columns
    //!
    virtual void GetCurrencySize(int& rows, 
                                 int& cols) const {}
    //!
    //! Gets values as a string
    //! \param fmt Output format
    //!
    virtual std::string GetValues(const OutputFormat* fmt) const;

    //!
    //! Returns true if this is a chained display, which needs to be deleted
    //! along with given currency display
    //! \param display Given display
    //!
    virtual bool IsChainedDisplay(CurrencyDisplay* display) const;
    //!
    //! Gets chained display
    //!
    virtual long GetChainedDisplay() const { return _chaineddisplay; }
    //!
    //! Sets chained display
    //! \param id Chained display id
    //!
    virtual void SetChainedDisplay(long id) { _chaineddisplay = id; }

protected:
    //!
    //! Sets data for forward pagination
    //!
    virtual void SetForwardDisplayData();
    //!
    //! Sets data for back pagination
    //!
    virtual void SetBackDisplayData();

private:
    long int _chaineddisplay;  //!< Chained display id

    //!
    //! Constructor - Only currency is allowed to construct
    //! \param cur Currency associated with this display
    //!
    StringDisplay(const Currency& cur);
    
    StringDisplay();                                      // Stubbed out 
    StringDisplay(            const StringDisplay& src) ; // Stubbed out 
    StringDisplay& operator=( const StringDisplay& src);  // Stubbed out
    
    //!
	//! Gets matrix data with no pagination - using defaults
	//! \param fmt Format
    //!
	std::string GetOutputNoPagination(const OutputFormat* fmt,
                                      std::ostringstream& os) const;
    //!
	//! Gets matrix data with forward pagination
	//! \param fmt Format
    //!
	std::string GetOutputForwardPagination(const OutputFormat* fmt,
                                           std::ostringstream& os) const;
    //!
	//! Gets matrix data with back pagination
	//! \param fmt Format
    //!
	std::string GetOutputBackPagination(const OutputFormat* fmt,
                                        std::ostringstream& os) const; 
};

#endif

