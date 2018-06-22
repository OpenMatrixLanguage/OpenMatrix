/**
* @file BetterCalcFuncs.cpp
* @date June 2015
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

// Begin defines/includes

#include "../Runtime/BuiltInFuncsCore.h"
#include "../Runtime/Interpreter.h"
#include "../Runtime/StructData.h"
#include "../Runtime/ErrorInfo.h"
#include "../Runtime/OML_Error.h"
#include "BetterCalc.h"

#if OS_WIN
#include <Windows.h>
#endif

extern Interpreter* interp;
#define OML_PRODUCT "OpenMatrix "
#define OML_VERSION "1.0"
// End defines/includes


bool hml_keylist(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    // std::cout << "hml_keylist running..." << std::endl;
	std::vector<std::string> funclist = interp->GetKeywords();
	std::string funcstring;
	for (std::vector<std::string>::iterator it = funclist.begin() ; it != funclist.end(); ++it)
      std::cout << ' ' << *it;
	std::cout << "" << std::endl;
	// std::cout << "hml_keylist completed." << std::endl;
	return true;
}

#if OS_WIN
void WindowsClear() {
    COORD topLeft  = { 0, 0 };
    HANDLE console = GetStdHandle(STD_OUTPUT_HANDLE);
    CONSOLE_SCREEN_BUFFER_INFO screen;
    DWORD written;

	GetConsoleScreenBufferInfo(console, &screen);
    FillConsoleOutputCharacterA(
        console, ' ', screen.dwSize.X * screen.dwSize.Y, topLeft, &written
    );
    FillConsoleOutputAttribute(
        console, FOREGROUND_GREEN | FOREGROUND_RED | FOREGROUND_BLUE,
        screen.dwSize.X * screen.dwSize.Y, topLeft, &written
    );
    SetConsoleCursorPosition(console, topLeft);
}

// The following does not compile for me.  I'm not aware of any problems with the
// above code on windows.  But I'll keep this code for now.
#if 1
void WindowsScreenFilling()
{
	CONSOLE_SCREEN_BUFFER_INFO Info;
	HANDLE ConsoleHandle;
	ConsoleHandle = CreateConsoleScreenBuffer(GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, CONSOLE_TEXTMODE_BUFFER, NULL);
		// dwDesiredAccess
		// dwShareMode
		// SECURITY_ATTRIBUTES *lpSecurityAttributes,
		// dwFlags this is fixed
		// Reserved
	// wrong? GetConsoleWindowInfo(ConsoleHandle, &Info) ;
	GetConsoleScreenBufferInfo(ConsoleHandle, &Info) ;
	SHORT Width = Info.srWindow.Right - Info.srWindow.Left + 1 ;
	for ( SHORT N = Info.srWindow.Top ; N <= Info.srWindow.Bottom ; ++N ) {
		DWORD Chars ;
		COORD Pos = { Info.srWindow.Left, N } ;
		FillConsoleOutputCharacter(ConsoleHandle, ' ', Width, Pos, &Chars) ;
		// FillConsoleOutputAttribute(ConsoleHandle, attr, Width, Pos, &Chars) ;
		FillConsoleOutputAttribute(ConsoleHandle, NULL, Width, Pos, &Chars) ;
	}
	COORD TopLeft = { Info.srWindow.Left, Info.srWindow.Top } ;
	SetConsoleCursorPosition(ConsoleHandle, TopLeft) ;
}

void WindowsScreenScrolling()
{
	CONSOLE_SCREEN_BUFFER_INFO Info;
	HANDLE ConsoleHandle;
	ConsoleHandle = CreateConsoleScreenBuffer(
			GENERIC_READ | GENERIC_WRITE, // dwDesiredAccess
			FILE_SHARE_READ | FILE_SHARE_WRITE, // dwShareMode
			NULL, // SECURITY_ATTRIBUTES *lpSecurityAttributes,
			CONSOLE_TEXTMODE_BUFFER, // dwFlags this is fixed
			NULL // Reserved
			);
	// Wrong? GetConsoleWindowInfo(ConsoleHandle, &Info) ;
	GetConsoleScreenBufferInfo(ConsoleHandle, &Info) ;
	CHAR_INFO space ;
	space.Char.AsciiChar = ' ' ;
	space.Attributes = NULL; // attr ;
	SHORT Height = Info.srWindow.Bottom - Info.srWindow.Top + 1 ;
	COORD Origin = { Info.srWindow.Left, Info.srWindow.Top - Height } ;
	ScrollConsoleScreenBuffer(ConsoleHandle, &Info.srWindow, NULL, Origin, &space) ;
	COORD TopLeft = { Info.srWindow.Left, Info.srWindow.Top } ;
	SetConsoleCursorPosition(ConsoleHandle, TopLeft) ;
}
#endif
#endif

bool hml_clc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
#if OS_WIN
	// windows clear with system call
	system("cls");
	return true;
#else
	// try this for Linux build.
	system("clear");
	return true;
#endif
#if 0
	// generic clear by writing multiple rows of blank lines.
	int screen_size = 50;
	int screen_width = 50;
	std::string blanks;
	blanks.assign(screen_width, ' ');
	std::string command = "disp('";
	command.append(blanks);
	command.append("')");
	for(int i=0; i < screen_size; i++)
		// std::cout << blanks << std::endl;
    	interp->DoString(command);
	return true;
#endif 
}

bool OML_help()
{
	// This could be added as a OML function if needed.   I did not implement as such because:
	// The requested argument to invoke this help list was -help, not a valid function name.
	// A help function already exists in OML.DLL in the BuildInFuncs.cpp module.
	std::cout << OML_PRODUCT << " Console Help:" << std::endl;
	std::cout << "-e \"any valid OML command\" --> execute the OML command in " << OML_PRODUCT << " and then quit." << std::endl;
	std::cout << "-f foo.oml --> load and execute the OML script called foo.oml in " << OML_PRODUCT << " and then quit" << std::endl;
	std::cout << "-help      --> give the list of supported optional arguments" << std::endl;
	std::cout << "-version   --> give the version of " << OML_PRODUCT << std::endl;
	return true;
}
//------------------------------------------------------------------------------
//! Gets the version string of the application
//------------------------------------------------------------------------------
std::string GetVersion(const std::string appdir)
{
    std::string version (OML_PRODUCT);
    version += " Version ";
    version += OML_VERSION;
#if 0
    std::string vfile = appdir;
    if (!vfile.empty())
        vfile += "/";
    vfile += "config/hypermath/version.xml";

    version += BuiltInFuncsCore::GetBuildNumber(vfile);
#endif 
    return version;
}
//------------------------------------------------------------------------------
//! Returns true if successful in getting the version string
//! \param[in]  eval    Evaluator interface
//! \param[in]  inputs  Vector of inputs
//! \param[out] outputs Vector of outputs
//------------------------------------------------------------------------------
bool OmlVersion(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs)
{
	std::string version = GetVersion(eval.GetApplicationDir());
	Currency out (version);
    out.DispOutput();

    outputs.push_back(out);
    return true;
}

