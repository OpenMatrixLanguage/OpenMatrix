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
#include "../Runtime/Interpreter.h"
#include "../Runtime/EvaluatorDebug.h"
#include "../Runtime/StructData.h"	// hml2 for now 

#include <cassert>
#include <clocale>
#include <fstream>
#include <memory>
#include "BetterCalc.h"
#include "ConsoleWrapper.h"
#include "SignalHandler.h"

#include "OmlServer.h"

//! Shows new prompt on console
//! \param[in] wrapper Interpreter wrapper
void CallNewConsolePrompting( ConsoleWrapper* wrapper);
// global variables for now
Interpreter* interp = NULL;
std::string dummyFilename;

bool ClearCommandString = false;
#define DEFAULT_PORT "9099"

// User interrupt related variables, decls
static int g_numUserInterrupts = 0;  //! Number of control-C interrrupts pressed
#ifdef OS_WIN
#   include <Windows.h>
//! Called when control key is down
//! \param[in] controlType Control type
BOOL OnControlKeyDown( DWORD controlType);
#else
#   include <signal.h>
#   include <time.h>
#   include <stdio.h>
#   include <unistd.h>
//! Called when control key is down
//! \param[in] controlType Control type
void OnControlKeyDown( int controlType);
#endif

class MyDebugListener : public EvaluatorDebugInterface
{
public:
	MyDebugListener(void* debuggee):_debuggee(debuggee) {}

	// we'll need to push the dsi info into the debuggee first, then call ProcessWhatever (dependent on type of event)
	virtual void PreFunctionCall(DebugStateInfo dsi) { std::cout << "Call " << dsi.function_name << std::endl;}
	virtual void PreStatement(DebugStateInfo dsi) { std::cout << "Execute " << dsi.filename << " line " << dsi.line_number << std::endl;}
	virtual void PreFunctionReturn(DebugStateInfo dsi) { std::cout << "Return " << dsi.function_name << std::endl;}

private:
	void* _debuggee; // this should be an actual hwdbgHMLDebuggee pointer instead, but I can't do that here yet
};


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//! Entry point for console/batch
int main(int argc, char* argv[])
{
    // Flags for command line arguments
    std::string port;
    bool nextIsPort = false, isServer = false;
    // -continue flag is required if commands/files are specified in the command line.
	bool continueAfterScript = false; 
    // -continue flag is only required if commands/files are specified in the command line.
	bool continueRequired = false;    

    // True when toolboxes are loaded
	bool toolbox_loaded = false;
    // True if toolboxes need to be loaded, not loaded by default in debug
	bool toolbox_load   = false; // No toolbox functions for debug.
#if !defined(_DEBUG)
	toolbox_load = true; // always load the toolboxes in release
#endif 

    // File to execute, if specified
	std::string scriptPath;

    // First pass through arguments.
	for (int ix = 1; ix < argc; ++ix)
	{
		std::string arg (argv[ix]);
		// convert command line argument to lower case for comparison
		std::string lower_str = arg;
		std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);

		if(lower_str == "/toolbox" || lower_str == "-toolbox")
		{
			toolbox_load = true;  // always load the toolboxes
		}
		else if(lower_str == "/notoolbox" || lower_str == "-notoolbox")
		{
			toolbox_load = false;  // provide capability to not load toolboxes
		}
		else if(lower_str.compare("/x") == 0 || lower_str.compare("-x") == 0)
		{
            interp->SetExperimental(true);
		}
        else if (lower_str.compare("-port") == 0)
		{
			isServer = true;
			nextIsPort = true;
		}
        else if (nextIsPort)
		{
			nextIsPort = false;
			port = arg;
		}
	}


	interp = new Interpreter;
    
    char* cdir = getenv("OML_APPDIR");
    if (cdir)
        interp->SetApplicationDir(std::string(cdir));
    
    interp->RegisterBuiltInFunction("version", OmlVersion, FunctionMetaData(0, 1, "CoreMinimalInterpreter"));
    // Wrapper for interpreter signals and methods
    ConsoleWrapper* wrapper = nullptr; 
    if(! isServer )
    {
        interp->RegisterBuiltInFunction("clc",	   hml_clc,    FunctionMetaData(0, 1, "CoreMinimalInterpreter")); // FunctionMetaData(OML_NO_NARG, OML_NO_NARG));
        wrapper = new ConsoleWrapper (interp);
        assert(wrapper);
    }

	/*if(toolbox_load && !toolbox_loaded)
	{
        std::string mathtoolbox = interp->GetApplicationDir();
        if (!mathtoolbox.empty())
            mathtoolbox += "/";
        mathtoolbox += "/plugins/scriptview/mathtoolbox/init_console.oml";
        interp->DoFile(mathtoolbox);
		toolbox_loaded = true;
	}*/

    // Set the interrupt handler after the toolboxes are loaded as this 
    // interferes with fortran libraries loaded in toolboxes and how they handle
    // control C
    SetUserInterruptHandler();

	bool e_argument_flag = false;  // -e command line argument
	bool f_argument_flag = false;  // -f command line argument

	// second argument parse - evaluate all other arguments
	for (int ix = 1; ix < argc; ++ix)
	{
		std::locale loc;
		std::string arg (argv[ix]);
		// convert command line argument to lower case for comparison
		std::string lower_str = arg;
		std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);

        if(*argv[ix] == '-')
		{
			e_argument_flag = f_argument_flag = false;  // clear -e and -f 
		}

		if (lower_str == "/toolbox"   || lower_str == "-toolbox"   || 
            lower_str == "/notoolbox" || lower_str == "-notoolbox" ||
            lower_str == "/x"         || lower_str == "-x")
		{
			continue;  // These options have been processed
		}
		else if (lower_str == "-version")
		{
            std::cout << GetVersion(interp->GetApplicationDir()) << std::endl;
			continueRequired = true;
			continue;  
		}
		else if (lower_str == "-e")
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
		else if (lower_str == "-help")
		{
			OML_help();
			return 0;
		}
		else if (lower_str == "/quiet" || lower_str == "-quiet")
        {
            if(wrapper) wrapper->SetQuietMode(true);
        }
				
        else if(lower_str == "/continue" || lower_str == "-continue")
			continueAfterScript=true;

		else if(lower_str.compare(0,10,"/filename=") == 0)
			dummyFilename = arg.substr(10, arg.length()-10);

		else if(lower_str.compare(0,5,"/ansi") == 0)
		{
#if OS_WIN
			// For testing console interface in Windows.  
			// I don't know if there are Linux equivalents.
			// set code page to 1252
			
			// interp->DoString("system('chcp 1252');");
			SetConsoleCP(1252);
			SetConsoleOutputCP(1252);
 
			// other windows code pages:
			// chcp 437 US
			// chcp 65001 unicode
			// chcp 850 latin 1
			// chcp 1252 ansi
#endif
		}
		else if(lower_str.compare(0,5,"/utf8") == 0)
		{
#if OS_WIN
			// For testing console interface in Windows.  
			// I don't know if there are Linux equivalents.
			SetConsoleCP(65001);
			SetConsoleOutputCP(65001);
#endif 
		}
        else
			scriptPath = arg; 
	}

    if (isServer)
	{

        if ("" == port) 
            port = DEFAULT_PORT;

        std::unique_ptr<OmlServer> omlserver(new OmlServer(interp));
        omlserver->Start(port);
    }
    else
    {
	    std::setlocale(LC_ALL, "");
	    std::setlocale(LC_NUMERIC, "C");

        if(scriptPath.compare("") != 0)
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
        PrintBanner();  

        if(wrapper)  CallNewConsolePrompting(wrapper);
    }

    delete wrapper;
	delete interp;

	return 0;
}
//------------------------------------------------------------------------------
//! Shows new prompt on console
//! \param[in] wrapper Interpreter wrapper
//------------------------------------------------------------------------------
void CallNewConsolePrompting(ConsoleWrapper* wrapper)
{
    assert(interp);
    assert(wrapper);
    wrapper->SetEnablePagination(true);

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

        if (quietmode && output.IsError()) break;

    }
}
//------------------------------------------------------------------------------
//! Gets input command
//! \param[in] wrapper Oml interp wrapper
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
//! Called when control key is down
//! \param[in] controlType Control event type
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
//! Called when control key is down
//! \param[in] controlType Control event type
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
//! Prints banner
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
//! Sets user interrupt
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
//! Handles user interrupt returns true if application needs to quit
//! \param[in] omlWrapper Oml interpreter wrapper
//------------------------------------------------------------------------------
bool HandleUserInterrupt(ConsoleWrapper* omlWrapper)
{
    if (g_numUserInterrupts <= 0) return false;

    assert(omlWrapper);
    if (!omlWrapper) return true;

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
//! Sets user interrupt handler - hooks to control C
//------------------------------------------------------------------------------
void SetUserInterruptHandler()
{
#ifdef OS_WIN
	SetConsoleCtrlHandler((PHANDLER_ROUTINE) OnControlKeyDown, TRUE);
#else
    signal(SIGINT, OnControlKeyDown);
#endif
}
