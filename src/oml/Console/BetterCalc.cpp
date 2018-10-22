/**
* @file BetterCalc.cpp
* @date June 2014
* Copyright (C) 2014-2018 Altair Engineering, Inc.  
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

#include "../Runtime/Hml2Dll.h"
#include "../Runtime/Currency.h"
#include "../Runtime/CurrencyDisplay.h"
#include "../Runtime/Interpreter.h"
#include "../Runtime/EvaluatorDebug.h"
#include "../Runtime/BuiltInFuncsUtils.h"
#include "../Runtime/StructData.h"	// hml2 for now 

#include <cassert>
#include <clocale>
#include <fstream>
#include <memory>

#include "BetterCalc.h"
#include "ConsoleWrapper.h"
#include "SignalHandler.h"
#include "OmlServer.h"

// Shows new prompt on console
// \param interpWrapper Interpreter wrapper
// \param license       License
void CallNewConsolePrompting(ConsoleWrapper* wrapper);

// global variables
Interpreter*    interp  = NULL;
ConsoleWrapper* wrapper = nullptr;
std::string     dummyFilename;
bool            ClearCommandString = false;
#define DEFAULT_PORT "9099"

// User interrupt related variables, decls
static int g_numUserInterrupts = 0;  //! Number of control-C interrrupts pressed
#ifdef OS_WIN
#   include <Windows.h>
// Called when control key is down
// \param controlType Control type
BOOL OnControlKeyDown(DWORD controlType);
#else
#   include <signal.h>
#   include <time.h>
#   include <stdio.h>
#   include <unistd.h>
// Called when control key is down
// \param controlType Control type
void OnControlKeyDown(int controlType);
#endif

//------------------------------------------------------------------------------
// Entry point for console
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    interp = new Interpreter;
    char* cdir = getenv("OML_APPDIR");
    if (cdir)
    {
        interp->SetApplicationDir(std::string(cdir));
    }

    bool paginateVal = CurrencyDisplay::GetPaginate();
    CurrencyDisplay::SetPaginate(0);  // Disable pagination

    // Flags for command line arguments

	bool continueAfterScript = false; // Needed if commands/files are in cmd line
    bool continueRequired    = false; // Needed if commands/files are in cmd line	  
	bool toolbox_load        = false; // True if toolboxes need to be loaded
#if !defined(_DEBUG)
	toolbox_load = true;              // always load the toolboxes in release
#endif 

    std::string port;
    bool nextIsPort = false, isServer = false;

    // File to execute, if specified
	std::string scriptPath;

    // First pass through arguments.
    std::vector<std::string> argsv;
    argsv.reserve(argc);
    assert(argc > 0);
    argsv.push_back(std::string(argv[0]));

    std::vector<std::string> argsToProcess;
    argsToProcess.reserve(argc);

	for (int ix = 1; ix < argc; ++ix)
	{
		std::string arg (argv[ix]);
        argsv.push_back(arg);

		// convert command line argument to lower case for comparison
		std::string lower_str = arg;
		std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);

        std::string cmd;
        if (!lower_str.empty() && (lower_str[0] == '-' || lower_str[0] == '/'))
        {
            cmd = (lower_str.size() != 1) ? lower_str.substr(1) : "";
        }

		if (cmd == "toolbox")
		{
			toolbox_load = true;  // always load the toolboxes
		}
		else if (cmd == "notoolbox")
		{
			toolbox_load = false;  // provide capability to not load toolboxes
		}
		else if(cmd == "x")
		{
            interp->SetExperimental(true);
		}
        else if (cmd == "continue")
        {
			continueAfterScript = true;
        }
		else if (cmd == "version")
		{
            std::cout << GetVersion(interp->GetApplicationDir()) << std::endl;
			continueRequired = true;
		}
		else if (cmd == "help")
		{
			OML_help();
		    delete interp;
			return 0;
		}
        else if (cmd == "port")
		{
			isServer = true;
			nextIsPort = true;
		}
        else if (nextIsPort)
		{
			nextIsPort = false;
			port = arg;
		}
        else
        {
            argsToProcess.push_back(arg);
        }
	}

    // Wrapper for interpreter signals and methods
    ConsoleWrapper* wrapper = nullptr; 
    if (!isServer)
    {
        wrapper = new ConsoleWrapper (interp);
        assert(wrapper);

        wrapper->SetArgv(argsv);
    }
    RegisterBuiltInFuncs();
    
    // Set the interrupt handler after the toolboxes are loaded as this 
    // interferes with fortran libraries loaded in toolboxes and how they handle
    // control C
    SetUserInterruptHandler();

	bool e_argument_flag = false;  // -e command line argument
	bool f_argument_flag = false;  // -f command line argument

	// second argument parse - evaluate all other arguments
    for (std::vector<std::string>::const_iterator itr = argsToProcess.begin();
         itr != argsToProcess.end(); ++itr)
	{
		std::string arg (*itr);
        if (arg.empty())
        {
            continue;
        }

		// convert command line argument to lower case for comparison
		std::string lower_str = arg;
		std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);

        char ch = arg[0];
        if (ch == '-')
		{
			e_argument_flag = false;
            f_argument_flag = false;  // clear -e and -f 
		}

		if (lower_str == "-e")
		{
			// process arguments following as oml commands until the next argument which begins with a dash '-'
			e_argument_flag = true;
			continue;
		}
		else if(e_argument_flag)
		{
			// process the argument as a oml command
			interp->DoString(arg);
			continueRequired = true;
		}
		else if (lower_str == "-f")
		{
			// process arguments following as script files until the next argument which begins with a dash '-'
			f_argument_flag = true;
			continue;
		}
		else if(f_argument_flag)
		{
			// process the argument as a script file
			interp->DoFile(arg);
			continueRequired = true;
		}
		else if (lower_str.compare(0,10,"/filename=") == 0)
        {
			dummyFilename = arg.substr(10, arg.length()-10);
        }
		else if(lower_str.compare(0,5,"/ansi") == 0)
		{
#if OS_WIN
			SetConsoleCP(1252);
			SetConsoleOutputCP(1252);
#endif
		}
		else if(lower_str.compare(0,5,"/utf8") == 0)
		{
#if OS_WIN
			SetConsoleCP(65001);
			SetConsoleOutputCP(65001);
#endif 
		}
		else if (lower_str == "/quiet" || lower_str == "-quiet")
        {
            if (wrapper)
            {
                wrapper->SetQuietMode(true);
            }
        }

        else
			scriptPath = arg; 
	}

    if (isServer)
	{
        if (port.empty()) 
        {
            port = DEFAULT_PORT;
        }
        std::unique_ptr<OmlServer> omlserver(new OmlServer(interp));
        omlserver->Start(port);
    }
    else
    {
	    std::setlocale(LC_ALL, "");
	    std::setlocale(LC_NUMERIC, "C");

        if(!scriptPath.empty())
	    {
		    Currency output = interp->DoFile(scriptPath);	
		    continueRequired = true;
	    }
	
	    if(continueRequired && !continueAfterScript)
	    {
		    delete wrapper;
		    delete interp;
		    return 0;
	    }

        // Start of interactive mode
        CurrencyDisplay::SetPaginate(paginateVal);  // Reset pagination value
        PrintBanner();  

        if (wrapper)  
        {
            CallNewConsolePrompting(wrapper);
        }
    }

    delete wrapper;
	delete interp;

	return 0;
}
//------------------------------------------------------------------------------
// Shows new prompt on console
//------------------------------------------------------------------------------
void CallNewConsolePrompting(ConsoleWrapper* wrapper)
{
    assert(interp);
    assert(wrapper);
    wrapper->SetPaginationVariables();

    int  hearbeatcounter = 1;
    bool quietmode       = wrapper->GetQuietMode();
	while (1)
    {
        bool breakloop = HandleUserInterrupt(wrapper);
        if (breakloop)
            break;

        bool isPaginating = wrapper->IsPaginating();
        if (isPaginating)
        {
            wrapper->Paginate();
            continue;
        }

        wrapper->PrintNewPrompt();
        
        std::string strCommand = GetInputCommand(wrapper);
        Currency    output = interp->DoString(strCommand);

        if (quietmode && output.IsError()) 
        {
            break;
        }
    }
}
//------------------------------------------------------------------------------
// Gets input command
//------------------------------------------------------------------------------
std::string GetInputCommand(ConsoleWrapper* wrapper)
{
    assert(wrapper);

    std::string strCommand;
    bool        quietmode = wrapper->GetQuietMode();
    while (1)
    {
		if (ClearCommandString)
		{
			strCommand.clear();
			ClearCommandString = false;
			std::cin.clear();
		}

		// if running in quietMode and redirected file has ended
		// break out of the input loop.
		if (quietmode && std::cin.eof()) break;

		if (!strCommand.empty())
            wrapper->PrintToStdout("??>");

        std::string partialCommand;
		std::getline(std::cin, partialCommand);

        if (!partialCommand.empty())
            g_numUserInterrupts = 0; // Reset control-C interrupt attempts

        strCommand += partialCommand;

        if (!strCommand.empty())
            strCommand += "\r\n";

        bool isPartialExpression = interp->IsPartialExpression(strCommand);
        if (!isPartialExpression) break;
			
	}
    if (!strCommand.empty())
        g_numUserInterrupts = 0; // Reset control-C interrupt attempts

    strCommand.append("\r\n");
    return strCommand;
}
#ifdef OS_WIN
//------------------------------------------------------------------------------
// Called when control key is down
//------------------------------------------------------------------------------
BOOL OnControlKeyDown(DWORD controlType)
{
    ClearCommandString = true;

    // This routine may be called in a separate thread. So, just set interrupt
    // flag instead of calling exit handler which will destroy interpreter
    if (controlType == CTRL_C_EVENT)
        SetUserInterrupt();

    return true;
}
#else
//------------------------------------------------------------------------------
// Called when control key is down
//------------------------------------------------------------------------------
void OnControlKeyDown(int controlType)
{
    ClearCommandString = true;

    // This routine may be called in a separate thread. So, just set interrupt
    // flag instead of calling exit handler which will destroy interpreter
    if (controlType == SIGINT)
        SetUserInterrupt();
}
#endif

//------------------------------------------------------------------------------
// Prints banner
//------------------------------------------------------------------------------
void PrintBanner()
{
#ifndef _DEBUG
    std::string line ("//----------------------------------------------------------------------");

    std::cout << line << std::endl;
	std::cout << GetVersion(interp->GetApplicationDir()) << std::endl;
	std::cout << "(c) Altair Engineering, Inc. and Contributors. (2007-2018)"  << std::endl;
	std::cout << line << std::endl;
#endif
}
//------------------------------------------------------------------------------
// Sets user interrupt
//------------------------------------------------------------------------------
void SetUserInterrupt()
{
    g_numUserInterrupts++;

    // First interruption will cause interpreter to interrupt. Subsequent 
    // interrupts will cause application to quit
    if (g_numUserInterrupts == 1 && interp)
        interp->TriggerInterrupt();
}
//------------------------------------------------------------------------------
// Handles user interrupt returns true if application needs to quit
//------------------------------------------------------------------------------
bool HandleUserInterrupt(ConsoleWrapper* omlWrapper)
{
    if (g_numUserInterrupts <= 0) 
    {
        return false;
    }

    assert(omlWrapper);
    if (!omlWrapper) 
    {
        return true;
    }

    std::cout << std::endl;
    
    if (g_numUserInterrupts == 1)
    {
        omlWrapper->HandleOnClearResults();
        return false;  // Only an interrupt was requested
    }
    
    omlWrapper->PrintNewPrompt();
    omlWrapper->PrintToStdout("Quitting the application...");
    omlWrapper->HandleOnSaveOnExit();

    return true;
}
//------------------------------------------------------------------------------
// Sets user interrupt handler - hooks to control C
//------------------------------------------------------------------------------
void SetUserInterruptHandler()
{
#ifdef OS_WIN
	SetConsoleCtrlHandler((PHANDLER_ROUTINE) OnControlKeyDown, TRUE);
#else
    signal(SIGINT, OnControlKeyDown);
#endif
}
//------------------------------------------------------------------------------
// Registers built in functions
//------------------------------------------------------------------------------
void RegisterBuiltInFuncs()
{
    assert(interp);

    if (wrapper)
    {
	    interp->RegisterBuiltInFunction("clc", hml_clc,    
            FunctionMetaData(0, 1, "CoreMinimalInterpreter")); 
    }
    interp->RegisterBuiltInFunction("version",            OmlVersion, 
        FunctionMetaData(0, 1, "CoreMinimalInterpreter"));
	interp->RegisterBuiltInFunction("getargc",            OmlGetArgC,
        FunctionMetaData(1, 0, "CoreMinimalInterpreter"));
	interp->RegisterBuiltInFunction("getargv",            OmlGetArgV, 
        FunctionMetaData(1, 1, "CoreMinimalInterpreter"));
}