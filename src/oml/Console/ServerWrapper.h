/**
* @file ServerWrapper.h
* @date June 2017
* Copyright (C) 2016-2017 Altair Engineering, Inc.  
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

#ifndef __SERVER_WRAPPER__
#define __SERVER_WRAPPER__

#include "WrapperBase.h"

class Interpreter;
class OmlServer;
// End defines/includes

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//! Class for handling interpretor signals and methods
//------------------------------------------------------------------------------
class ServerWrapper : public WrapperBase
{
public:
    //! Constructor
    //! \param[in] interpretor
    //! \param[in] OmlServer
    ServerWrapper(Interpreter* interp, OmlServer* server);

    //! Destructor
    ~ServerWrapper();

    //! Prints currency to console
    //! \param[in] cur Currency to print
	void HandleOnPrintResult( const Currency& cur);

    //! Displays a prompt and gets user input
	//! \param[in]  prompt    Prompt to display to user
    //! \param[in]  type      Type, if specified
	//! \param[out] userInput Input from user
    void HandleOnGetUserInput( const std::string& prompt,
                               const std::string& type,
		                       std::string&       userInput);
    //! Initiates a user defined pause
    //! \param[in] msg  User message to display
    //! \param[in] wait True if waiting for a keystroke input from user
	void HandleOnPauseStart( const std::string& msg, 
                             bool               wait);

private:

    // Private default constructor
    ServerWrapper() {}

    //! Prints to stdout
    //! \param[in] msg Message to print
    void PrintToStdout( const std::string& msg);

	OmlServer* _omlServer;
};


#endif

// End of file:
