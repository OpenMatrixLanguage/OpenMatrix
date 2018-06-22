/**
* @file WrapperBase.h
* @date January 2015
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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

#ifndef __WRAPPER_BASE__
#define __WRAPPER_BASE__

// Begin defines/includes
#include "Runtime/Currency.h"

class CurrencyDisplay;
class Interpreter;

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
//! Base class for handling interpretor signals and methods for console and atom client
//-------------------------------------------------------------------------------------
class WrapperBase
{
public:
    //! Constructor
    //! \param[in] interpretor
    WrapperBase(Interpreter* interp);
    //! Destructor
    virtual ~WrapperBase() = 0;
    //! Prints currency to console
    //! \param[in] cur Currency to print
    virtual void HandleOnPrintResult( const Currency& cur) {}
    //! Clears results and pagination related to results
    virtual void HandleOnClearResults() {}
    //! Handles application exit
    virtual void HandleOnSaveOnExit() {}
    //! Displays a prompt and gets user input
	//! \param[in]  prompt    Prompt to display to user
    //! \param[in]  type      Type, if specified
	//! \param[out] userInput Input from user
    virtual void HandleOnGetUserInput( const std::string& prompt,
                               const std::string& type,
                               std::string&       userInput) {}
    //! Initiates a user defined pause
    //! \param[in] msg  User message to display
    //! \param[in] wait True if waiting for a keystroke input from user
	virtual void HandleOnPauseStart( const std::string& msg, 
                                     bool               wait) {}
    //! Handles nested displays during pagination
    //! \param[in] display Display to be added
    virtual void HandleOnAddDisplay( CurrencyDisplay* display) {}

protected:
    
    // Protected default constructor
    WrapperBase() : _interp(0) {}

    //! Connects oml signal handler
    void ConnectOmlSignalHandler();
    //! Disconnect oml signal handler
    void DisconnectOmlSignalHandler();

    Interpreter* _interp;  //! Interpreter
};

#endif
