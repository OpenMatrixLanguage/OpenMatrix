/**
* @file OmlServer.cpp
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

// Begin defines/includes

#include "OmlServer.h"

#include "ServerWrapper.h"

#include "Runtime/Interpreter.h"
#include "Runtime/SignalHandlerBase.h"

#ifdef OS_WIN
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <winsock2.h>
#include <ws2tcpip.h>
#define RW_OFF SD_BOTH
#else
#include <unistd.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <pthread.h>
#define SOCKET int
#define SOCKET_ERROR 0xFFFFFFFF
#define INVALID_SOCKET 0xFFFFFFFF
#define RW_OFF SHUT_RDWR 
#endif

#include <cstring>
#include <numeric>
#include <functional>
#include <memory>

#define BUFLEN 2048

// End defines/includes

#ifdef OS_WIN
DWORD WINAPI RunServer(LPVOID args)
#else
void* RunServer(void* args)
#endif
{
    ServerArgs* servArgs = (ServerArgs*)(args);
    if (servArgs)
    {
        std::unique_ptr<OmlServer> omlserver(new OmlServer(servArgs->interp, servArgs->handler));
        omlserver->Start(servArgs->port);
    }

    return 0;
}

void RunOmlServer(ServerArgs* args, bool inWorkerThread)
{
    if (inWorkerThread)
    {
#ifdef OS_WIN

        HANDLE  hThread;
        DWORD   dwThreadId;
        hThread = CreateThread(
            NULL,        // default security attributes
            0,           // use default stack size  
            RunServer,   // thread function name
            args,        // argument to thread function 
            0,           // use default creation flags 
            &dwThreadId  // returns the thread identifier 
        );
        CloseHandle(hThread);
#else
        pthread_t tid;
        pthread_create(&tid, NULL, &RunServer, args);
#endif
    }
    else
    {
        RunServer(args);
    }
}

namespace // anonymous
{
    //! Interpreter states
    enum InterpState
    {
        READY   = 1,
        BUSY    = 2,
        UNKNOWN = -1
    };
    //! A brief command structure to pass arguments to worker thread
    struct CommandArgs 
    {
        CommandArgs() : interp(nullptr) {}
	    std::string mode;
	    std::string command;
        Interpreter* interp;
        std::function<void(InterpState)> sendinterpstate_func;
    };

//! A C method to run in worker thread
#ifdef OS_WIN
DWORD WINAPI RunScript(LPVOID arg)
#else
void* RunScript(void *arg)
#endif
{
    CommandArgs *pArgs = (CommandArgs*)arg;
    if(pArgs && pArgs->interp)
    {
        if ("evalfile" == pArgs->mode)
        {
           pArgs->interp->DoFile(pArgs->command);
        }
        else if ("evalstring" == pArgs->mode)
        {
            pArgs->interp->DoString(pArgs->command);
        }

        if(pArgs->sendinterpstate_func)
            pArgs->sendinterpstate_func(READY);

        delete pArgs; 
    }
    return 0;
}

} // end of anonymous namespace

//----------------------------------------------------------------
//! A brief class to encapsulate server communication details
//----------------------------------------------------------------
class OmlServer::Internal
{
public:
    //! Constructor
    //! \param[in] interpreter
    //! \param[in] omlserrver handle
    Internal(Interpreter* interp, OmlServer* omlserver, SignalHandlerBase* handler):
      _sock(INVALID_SOCKET), _interp(interp), _omlServer(omlserver), _handler(handler)
    {
        InitMutex();
    }
    //! Destructor
    ~Internal()
    {
        DeleteMutex();
    }
    //! start a socket server communication at a given port
    //! port port no in string format
    bool Start(const std::string& port)
    {
#ifdef OS_WIN
        // Initialize Winsock
        WSADATA wsaData;
        if (WSAStartup(MAKEWORD(2,2), &wsaData) != 0)
            return false;
#endif
        struct addrinfo *result = NULL;
        struct addrinfo hints;
        memset(&hints, 0, sizeof(hints));

        hints.ai_family = AF_INET;
        hints.ai_socktype = SOCK_STREAM;
        hints.ai_protocol = IPPROTO_TCP;
        hints.ai_flags = AI_PASSIVE;

        int gotInfo = getaddrinfo(NULL, port.c_str(), &hints, &result);
        if ( gotInfo != 0 ) 
        {
            CleanUp();
            return false;
        }

        SOCKET sockfd = socket(result->ai_family, result->ai_socktype, result->ai_protocol);
        if (sockfd == INVALID_SOCKET) 
        {
            freeaddrinfo(result);
            CloseSocket(sockfd);
            CleanUp();
            return false;
        }

        // Setup listening socket
        int binded = bind(sockfd, result->ai_addr, (int)result->ai_addrlen);
        if (binded == SOCKET_ERROR) 
        {
            freeaddrinfo(result);
            CloseSocket(sockfd);
            CleanUp();
            return false;
        }

        freeaddrinfo(result);

        int listened = listen(sockfd, SOMAXCONN);
        if (listened == SOCKET_ERROR) 
        {
            CloseSocket(sockfd);
            CleanUp();
            return false;
        }

        while((_sock = accept(sockfd,0,0))==INVALID_SOCKET);

        SendInterpreterState(READY);

        std::unique_ptr<ServerWrapper> Wrapper(new ServerWrapper(_interp, _omlServer, _handler));

        char recvbuf[BUFLEN];
        int recvbuflen = BUFLEN;
        memset(recvbuf, 0, sizeof(recvbuf));
        while( recv(_sock, recvbuf, recvbuflen, 0) > 0)
        {
            std::string cmd(recvbuf);
                        
            size_t splitIndex = cmd.find_first_of(",");
            int cmdNumber = -1;
			std::string command = "";
            if(splitIndex > 0)
            {
                command = cmd.substr(0,splitIndex);
            }
            if("stop" == command) break;

			else if("evalstring" == command || "evalfile" == command)
            {
                CommandArgs *cmdArgs = new CommandArgs;
                cmdArgs->mode = command;
                cmdArgs->command = cmd.substr(splitIndex+1);
                cmdArgs->interp = _interp;
                using std::placeholders::_1;
                cmdArgs->sendinterpstate_func = std::bind(&OmlServer::Internal::SendInterpreterState,this, _1);

                SendInterpreterState(BUSY);

                RunInWorkerThread(cmdArgs);
            }
            else if("interrupt" == command)
            {
                _interp->TriggerInterrupt();
            }
			else if("ispartialexpression" ==  command)
			{
				bool isPartial = _interp->IsPartialExpression(cmd.substr(splitIndex+1));

				std::string status;
				if(_interp->IsPartialExpression(cmd.substr(splitIndex+1)))
				{
					status = "{\"type\":\"partialexpression\",\"data\":true}";
				}
				else 
				{
					status = "{\"type\":\"partialexpression\",\"data\":false}";
				}

				Send(status);
			}
            else if ("signature" == command)
            {
                std::vector<Currency> inputs;
                std::string funcName = cmd.substr(splitIndex + 1);
                inputs.emplace_back(funcName.c_str());
                Currency cur = _interp->CallFunction("getsyntax", inputs);
                if (cur.IsString())
                {
                    std::string result_str(cur.StringVal());
                    if (!result_str.empty())
                    {
                        result_str = result_str.substr(0, result_str.find_last_of("\n"));
                        std::string signatures = "{\"type\":\"signatures\",\"data\":\"" + result_str + "\"}";
                        Send(signatures);
                        std::string status = "{\"type\":\"status\",\"data\":\"ready\"}";
                        Send(status);
                    }
                }
            }
            memset(recvbuf, 0, sizeof(recvbuf));
        }
        shutdown(_sock, RW_OFF);
        CloseSocket(_sock);
        CleanUp();
        return true;
    }
    //! Send data to socket
    void Send(const std::string &data)
    {
        if( Lock() == 0)
        {
            if(_sock != INVALID_SOCKET)
            {
                send(_sock, data.c_str(), (int)data.length(), 0);
            }

            Unlock();
        }
    }

private:
    //! Create a new worker thread and run 'RunScript' method in new thread
    void RunInWorkerThread(CommandArgs *args)
    {
#ifdef OS_WIN
        
        HANDLE  hThread;
        DWORD   dwThreadId;
        hThread = CreateThread( 
                                NULL,        // default security attributes
                                0,           // use default stack size  
                                RunScript,   // thread function name
                                args,        // argument to thread function 
                                0,           // use default creation flags 
                                &dwThreadId  // returns the thread identifier 
                          );
        CloseHandle(hThread);
#else
        pthread_t tid;
        pthread_create(&tid, NULL, &RunScript, args);
#endif
    }
    //! Send interpreter state notification to socket
    void SendInterpreterState(InterpState state)
    {
        std::string status;
		if(state == BUSY)
		{
			status = "{\"type\":\"status\",\"data\":\"busy\"}";
		}
		else if(state == READY)
		{
			status = "{\"type\":\"status\",\"data\":\"ready\"}";
			std::vector<std::string> keywordList = _interp->GetKeywords();
			string keywords;
			keywords = "{\"type\":\"keywordslist\",\"data\":\""+std::accumulate(keywordList.begin(), keywordList.end(), keywords, 
											[] (const string& key1, const string& key2) -> string { 
												return key1.empty() ? key2 : key1 + " " + key2; })+"\"}";
			Send(keywords);

			std::vector<std::string> funcList = _interp->GetFunctionNames();
			string functions;
			functions = "{\"type\":\"functionslist\",\"data\":\""+std::accumulate(funcList.begin(), funcList.end(), functions, 
											[] (const string& func1, const string& func2) -> string { 
												return func1.empty() ? func2 : func1 + " " + func2; })+"\"}";

			Send(functions);
		}
        else
        {
            //TODO:
        }
		Send(status);
    }
    //! clean up post connection is closed
    void CleanUp()
    {
#ifdef OS_WIN
        WSACleanup();
#endif
    }
    //! close given socket for further communication
    void CloseSocket(SOCKET sock)
    {
#ifdef OS_WIN
        closesocket(sock);
#else
        close(sock);
#endif
    }
    //! Initialize mutex for thread synchronization
    void InitMutex()
    {
#ifdef OS_WIN
        _mutex = CreateMutex(NULL, FALSE, NULL);
#else
        pthread_mutex_init(&_mutex, nullptr);
#endif
    }
    //! Delete mutex
    void DeleteMutex()
    {
#ifdef OS_WIN
        CloseHandle(_mutex);
#else
        pthread_mutex_destroy(&_mutex);
#endif
    }
    //! Acquire mutex
    int Lock()
    {
#ifdef OS_WIN
        DWORD dwWaitResult = WaitForSingleObject(_mutex,INFINITE);
        return (int)dwWaitResult;
#else
        return pthread_mutex_lock(&_mutex);
#endif
    }
    //! Release mutex
    int Unlock()
    {
#ifdef OS_WIN
        return ReleaseMutex(_mutex);
#else
        return pthread_mutex_unlock(&_mutex);
#endif

    }

    Interpreter*       _interp;
    OmlServer*         _omlServer;
    SignalHandlerBase* _handler;
    SOCKET             _sock;
#ifdef OS_WIN
    HANDLE          _mutex;
#else
    pthread_mutex_t _mutex;
#endif
};

//-------------------------------------------------------------
//! Constructor
//! \param[in] interpreter
OmlServer::OmlServer(Interpreter* interp, SignalHandlerBase* handler)
{
    _internal = new Internal(interp, this, handler);
}
//! Destructor
OmlServer::~OmlServer()
{
    delete _internal;
}
//! Start OML server at given port
//! \param[in] port port no in string format
void OmlServer::Start(const std::string& port)
{
    _internal->Start(port);
}
//! Send data to server
//! \param[in] data in string format<json>
void OmlServer::Send(const std::string &data)
{
    _internal->Send(data);
}
