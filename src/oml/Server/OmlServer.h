/**
* @file OmlServer.h
* @date April 2023
* Copyright (C) 2023 Altair Engineering, Inc.  
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

#ifndef OML_SERVER_H__
#define OML_SERVER_H__

// Begin defines/includes
#include <ServerDll.h>

#include <string>

// End defines/includes

class Interpreter;
class SignalHandlerBase;

struct ServerArgs
{
    ServerArgs() : interp(nullptr), handler(nullptr) {}
    std::string port;
    Interpreter* interp;
    SignalHandlerBase* handler;
};

extern "C" {

    SERVERDLL_DECLS void RunOmlServer(ServerArgs* args, bool inWorkerThread = true);
}

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
//! Class for providing server capability to communicate with OML lanaguage over socket
//-------------------------------------------------------------------------------------
class SERVERDLL_DECLS OmlServer
{
public:
    //! Constructor
    //! \param[in] interpreter
    OmlServer(Interpreter* interp, SignalHandlerBase* handler);
    //! Destructor
    ~OmlServer();

    //! Start OML server at given port
    //! \param[in] port port no in string format
    void Start(const std::string& port);
    //! Send data to server
    //! \param[in] data in string format<json>
    void Send(const std::string &data);

private:
    //! class to encapsulate internal details of server
    class Internal;
    Internal* _internal;
};


#endif // OML_SERVER_H__
