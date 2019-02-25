/**
* @file SignalHandler.h
* @date June 2016
* Copyright (C) 2016-2018 Altair Engineering, Inc.  
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

#ifndef __SIGNALHANDLER_H__
#define __SIGNALHANDLER_H__

// Begin defines/includes
#include "Runtime/SignalHandlerBase.h"

class WrapperBase;
// End defines/includes

//------------------------------------------------------------------------------
//! Class for handling actions from oml core in client
//------------------------------------------------------------------------------
class SignalHandler : public SignalHandlerBase
{
public:
    SignalHandler() : _src (NULL) {}                  //! Constructor
    virtual ~SignalHandler() {}                       //! Destructor    

    //! Creates clone
    virtual SignalHandlerBase* CreateClone();
    //! Returns true of this is a clone
    virtual bool IsClone() const { return _src ? true : false; }
    
    //! Sets the console wrapper
    void SetWrapper( WrapperBase* wrapper) { _wrapper = wrapper; }

    //! Returns class info
    virtual const char* ClassName() const { return "SignalHandler"; }

    //! Clear results
    virtual void OnClearResultsHandler();
    //! Print result
    //! \param[in] cur Result to print
    virtual void OnPrintResultHandler( const Currency& cur);
    //! Start pause
    //! \param[in] msg  User message to display
    //! \param[in] wait True if waiting for a keystroke input from user
    virtual void OnPauseStartHandler( const std::string& msg,
                                      bool               wait);
    //! Prints prompt for save on exit
    virtual void OnSaveOnExitHandler();
    //! Adds nested display for pagination
    //! \param[in] display Nested display to add
    virtual void OnAddDisplayHandler( CurrencyDisplay* display);
    //! Get user input
    //! \param[in]  prompt Prompt to display to user
    //! \param[in]  type   Type, if specified
    //! \param[out] input  Input from user
    virtual void OnUserInputHandler( const std::string& prompt,
                                     const std::string& type,
                                     std::string&       input);

private:
    const SignalHandler*       _src;            //! Source, if this is a clone                 
    WrapperBase* _wrapper; //! Console wrapper

    SignalHandler( const SignalHandler* src);   //! Copy constructor
};

#endif
// End of file:
