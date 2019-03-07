/**
* @file ServerWrapper.cpp
* @date June 2017
* Copyright (C) 2017-2018 Altair Engineering, Inc.  
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

// Begin defines/includes

#include "ServerWrapper.h"

#include "SignalHandler.h"

#include <cassert>
#include <sstream>

#include "Runtime/Interpreter.h"
#include "Runtime/OutputFormat.h"

#include "OmlServer.h"

static OmlServer* s_server = nullptr;

bool server_clc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if(s_server)
    {
        std::string status = "{\"type\":\"control\",\"data\":\"clearcmdwindow\"}";
        s_server->Send(status);
    }
	return true;
}
//------------------------------------------------------------------------------
//! Constructor
//! \param[in] interpretor
//! \param[in] socket descriptor
//------------------------------------------------------------------------------
ServerWrapper::ServerWrapper(Interpreter* interp, OmlServer* server)
    : WrapperBase (interp)
    , _omlServer(server)
{
    s_server = server;
    interp->RegisterBuiltInFunction("clc", server_clc, FunctionMetaData(0, 1, "CoreMinimalInterpreter"));
}
//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
ServerWrapper::~ServerWrapper()
{
}
//------------------------------------------------------------------------------
//! Slot called when printing result to console
//! \param[in] cur Currency to print
//------------------------------------------------------------------------------
void ServerWrapper::HandleOnPrintResult(const Currency& cur)
{	
    assert(_interp);
    const OutputFormat* format = _interp->GetOutputFormat();
    PrintToStdout(cur.GetOutputString(format));
}
//------------------------------------------------------------------------------
//! Slot which displays a prompt and gets user input
//! \param[in]  prompt    Prompt to display to user
//! \param[in]  type      Type, if specified
//! \param[out] userInput Input from user
//------------------------------------------------------------------------------
void ServerWrapper::HandleOnGetUserInput(const std::string& prompt,
                                                     const std::string& type,
													 std::string&       userInput)
{
    fflush(stdout);
	std::cout << prompt;
	std::string waiting = "{\"type\":\"userinput\",\"data\":1}";
    _omlServer->Send(waiting);
    std::getline(std::cin, userInput);
    fflush(stdout);
	if (type == "s" || type == "S")
	{
		userInput += "'";
		userInput.insert(0, "'");
	}
	waiting = "{\"type\":\"userinput\",\"data\":0}";
    _omlServer->Send(waiting);
}
//------------------------------------------------------------------------------
//! Slot called when initiating a user defined pause. Any character typed will
//! break out of the pause
//! \param[in] msg  User message to display
//! \param[in] wait True if waiting for a keystroke input from user
//------------------------------------------------------------------------------
void ServerWrapper::HandleOnPauseStart(const std::string& msg, 
                                                   bool               wait)
{
    if (!wait)
    {
        PrintToStdout(msg);
        return;
    }

    PrintToStdout(msg);

    while(1)
	{
		if(_omlServer)
		{
			fflush(stdout);
			std::string waiting = "{\"type\":\"pause\",\"data\":1}";
            _omlServer->Send(waiting);

			std::string userInput;
			std::getline(std::cin, userInput);
			waiting = "{\"type\":\"pause\",\"data\":0}";
            _omlServer->Send(waiting);
		}
		break;
	}
}
//------------------------------------------------------------------------------
//! Prints message at standard o/p
//! \param[in] msg Message to print
//------------------------------------------------------------------------------
void ServerWrapper::PrintToStdout(const std::string& msg)
{
    std::cout << msg << std::endl;
    fflush(stdout);
}
