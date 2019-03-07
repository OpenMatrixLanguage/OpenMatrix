/**
* @file ConsoleWrapper.cpp
* @date June 2015
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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

#include "ConsoleWrapper.h"

#include "SignalHandler.h"

#include <algorithm>
#include <cassert>
#include <sstream>

#include "Runtime/BuiltInFuncsUtils.h"
#include "Runtime/CurrencyDisplay.h"
#include "Runtime/Interpreter.h"
#include "Runtime/OutputFormat.h"

#ifdef OS_WIN
#   include <Windows.h>
#else
#   include <sys/ioctl.h>
#   include <termios.h>
#   include <unistd.h>
#endif

//# define CONSOLEWRAPPER_DBG 1  // Uncomment to print debug info
#ifdef CONSOLERWRAPPER_DBG
#    define CONSOLEWRAPPER_PRINT(m) { std::cout << m; }
#    define CONSOLEWRAPPER_PRINT_PAGINATEINFO(rows, cols) {  \
            std::cout << "Rows (" << rows << ") Cols (" << cols << ")"; \
            std::cout << std::endl; }
#else
#    define CONSOLEWRAPPER_PRINT(m) 0
#    define CONSOLEWRAPPER_PRINT_PAGINATEINFO(rows, cols) 0
#endif

int maxcols = 0;
// End defines/includes

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
ConsoleWrapper::ConsoleWrapper(Interpreter* interp)
    : WrapperBase(interp)
    , _quietMode        (false)
    , _appendOutput     (false)
    , _addnewline       (false)
{
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
ConsoleWrapper::~ConsoleWrapper()
{
}
//------------------------------------------------------------------------------
// Slot called when printing result to console
//------------------------------------------------------------------------------
void ConsoleWrapper::HandleOnPrintResult(const Currency& cur)
{	
    bool couldPrint = PrintResult(cur);
    if (!couldPrint)
        _resultsToPrint.push_back(cur);  // Print later, once pagination is done
}
//------------------------------------------------------------------------------
// Gets visible rows and columns from command window screen size
//------------------------------------------------------------------------------
void ConsoleWrapper::GetCommandWindowInfo() const
{
    int rows = 0;
    int cols = 0;

    maxcols = 0;
#ifdef OS_WIN  // Windows specific code
    HANDLE hConsoleOutput = GetStdHandle(STD_OUTPUT_HANDLE);
    CONSOLE_SCREEN_BUFFER_INFO srBuffInfo;
    BOOL result = GetConsoleScreenBufferInfo(hConsoleOutput, &srBuffInfo);
    if (result)
    {
        // Lines and columns use the screen size and not the buffer size
        rows = abs(srBuffInfo.srWindow.Top   - srBuffInfo.srWindow.Bottom);
        cols = abs(srBuffInfo.srWindow.Right - srBuffInfo.srWindow.Left);

    }
#else
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    rows = w.ws_row;
    cols = w.ws_col;
#endif    

    maxcols = cols;
    CurrencyDisplay::SetMaxCols(cols);
	CurrencyDisplay::SetMaxRows(rows);
    CONSOLEWRAPPER_PRINT_PAGINATEINFO(rows, cols);
}
//------------------------------------------------------------------------------
// Print currency. Return false if this currency could not be printed and 
// needs to be processed later
//------------------------------------------------------------------------------
bool ConsoleWrapper::PrintResult(const Currency& cur)
{
    assert(_interp);
    const OutputFormat* format = _interp->GetOutputFormat();

    // Check the pagination environment variable and get command window info
    SetPaginationVariables();  

    bool canPaginate = CurrencyDisplay::CanPaginate(cur);
    if (!CurrencyDisplay::IsValidDisplaySize()) // No window size available so just print
    {
        PrintToConsole(cur.GetOutputString(format), cur.IsPrintfOutput());
        return true;
    }

    if (!_displayStack.empty()) return false;  // Can't print yet as display stack is not empty

    if (!canPaginate)
    {               
        // There is no print in progress, so don't need to worry
        PrintToConsole(cur.GetOutputString(format), cur.IsPrintfOutput());
        return true;
    }

    // At this point, the currency that can paginate. Cache if possible
    CacheDisplay(cur);
    bool wasPaginating = false;
    
    if (_interp->IsInterrupt())  // Check if there has been an interrupt
    {
        PrintPaginationMessage();
        HandleOnClearResults();
        return true;  // Done with printing
    }

    CurrencyDisplay* display = GetCurrentDisplay();
    assert(display);
    PrintToConsole(cur.GetOutputString(format), cur.IsPrintfOutput());

    CurrencyDisplay* topdisplay = GetCurrentDisplay();
    assert(topdisplay);

    if (topdisplay != display)
    {
        // Nested pagination
        Currency newcur = const_cast<Currency&>(topdisplay->GetCurrency());
        newcur.SetDisplay(topdisplay);
        if (topdisplay)
        {
            wasPaginating = (topdisplay->IsPaginatingRows() || topdisplay->IsPaginatingCols());
        }
        else
        {
            wasPaginating = false;
        }
    }

    if (topdisplay->IsPaginatingCols() || topdisplay->IsPaginatingRows())
    {
        PrintPaginationMessage(topdisplay->CanPaginateColumns());
        return true;               // Done printing for pagination
    }

    EndPagination(wasPaginating);           // Done paginating

    if (!_displayStack.empty())
    {
        display = _displayStack.top();
        assert(display);
        if (display->IsPaginatingCols() || display->IsPaginatingRows())
            PrintPaginationMessage(display->CanPaginateColumns());
        else
            EndPagination(true);  // Done paginating
        return true;
    }

    return true;
}
//------------------------------------------------------------------------------
// Prints message for continuing/skipping pagination
//------------------------------------------------------------------------------
void ConsoleWrapper::PrintPaginationMessage(bool showColControl)
{
#ifdef OS_WIN
    HANDLE hConsoleOutput = GetStdHandle(STD_OUTPUT_HANDLE);
	if (hConsoleOutput == INVALID_HANDLE_VALUE) return;

	CONSOLE_SCREEN_BUFFER_INFO csbiInfo; 
	if (!GetConsoleScreenBufferInfo(hConsoleOutput, &csbiInfo)) return;

	WORD oldColorAttrs = csbiInfo.wAttributes; // Cache color info

    // Gray background with white letters
	SetConsoleTextAttribute(hConsoleOutput, BACKGROUND_INTENSITY |
        FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
    std::cout << "-- (f)orward, (b)ack, (q)uit, (e)xit, (\x18), (\x19)";
    if (showColControl)
        std::cout << ", (\x1a), (\x1b)";
	SetConsoleTextAttribute(hConsoleOutput, oldColorAttrs); // Reset color info
#else

    std::cout << "-- (f)orward, (b)ack, (q)uit, (e)xit";
#endif
}
//------------------------------------------------------------------------------
// Prints end of pagination message
//------------------------------------------------------------------------------
void ConsoleWrapper::PrintPaginationMessage()
{
    CurrencyDisplay* display = GetCurrentDisplay();
    if (!display) return;

    std::string strmsg;
    display->GetPaginationEndMsg(strmsg);
    if (!strmsg.empty())
        std::cout << strmsg << std::endl;
}
//------------------------------------------------------------------------------
// Prints string to console
//------------------------------------------------------------------------------
void ConsoleWrapper::PrintToConsole(const std::string& result,
                                    bool               isprintoutput)
{	
    if (result.empty()) return;

    if (!_appendOutput && _addnewline)
        std::cout << std::endl;
    else if (!isprintoutput && _appendOutput)
        std::cout << std::endl;

    bool hasTrailingNewline = (result[result.length()-1] == '\n');
    if (!isprintoutput)
    {
        std::cout << result;
        if (!hasTrailingNewline)
            std::cout << std::endl;

        _appendOutput = false;
        _addnewline   = false;
        return;
    }
    
    if (!hasTrailingNewline)
        std::cout << result;
    else
    {
        std::string msg = result.substr(0, result.length()-1);
        std::cout << msg;
        _addnewline = true;
    }

    _appendOutput = !hasTrailingNewline;   
}
//------------------------------------------------------------------------------
// Clears results and pagination related to data
//------------------------------------------------------------------------------
void ConsoleWrapper::HandleOnClearResults()
{
    CurrencyDisplay::ClearLineCount();
    _appendOutput = false;

    _resultsToPrint.clear();

    while (!_displayStack.empty())
    {
        EndPagination(false); 
    }
}
//------------------------------------------------------------------------------
// Processes pagination
//------------------------------------------------------------------------------
void ConsoleWrapper::ProcessPagination()
{
    CurrencyDisplay* display = GetCurrentDisplay();
    assert(display);

    if (display->GetMode() == CurrencyDisplay::DISPLAYMODE_EXIT)
    {
        std::cout << "\r" << std::endl;

        while (!_displayStack.empty())
        {
            EndPagination(false);
        }
    }
    else
    {
        if (display->GetMode() != CurrencyDisplay::DISPLAYMODE_QUIT)
        {
	        // Print to console
	        GetCommandWindowInfo();
            display->SetModeData();

            const Currency& curBeingPrinted = display->GetCurrency();
            const_cast<Currency&>(curBeingPrinted).SetDisplay(display);

            std::string out = curBeingPrinted.GetOutputString(_interp->GetOutputFormat());
            if (display->GetDeleteLine())
            {
                std::cout << "\r";
                display->SetDeleteLine(false);
            }
            else
            {
                std::cout << std::endl;
            }
            bool isPrintOutput = curBeingPrinted.IsPrintfOutput();

            // Reget display as there could be nested pagination
            CurrencyDisplay* topdisplay = GetCurrentDisplay();
            assert(topdisplay);
            if (topdisplay != display)
            {
                Currency newcur = const_cast<Currency&>(topdisplay->GetCurrency());
                newcur.SetDisplay(topdisplay);
                isPrintOutput = newcur.IsPrintfOutput();
            }

            PrintToConsole(out, isPrintOutput);
	        if (topdisplay->IsPaginatingCols() || topdisplay->IsPaginatingRows()) 
            {
                PrintPaginationMessage(topdisplay->CanPaginateColumns());
                return; // Still printing
            }
        }
        else
        {
            std::cout << "\r" << std::endl;
            display->SetDeleteLine(false);
        }

        EndPagination(true);  // Done with pagination

        if (!_displayStack.empty())
        {
            display = _displayStack.top();
            if (display && display->IsPaginatingCols() || display->IsPaginatingRows())
            {
                if (!display->CanPrintRows())
                {
                    PrintPaginationMessage(display->CanPaginateColumns());
                    return;  // Nothing more to print
                }
                // Print the next rows
                display->SetModeData();
                const Currency& curBeingPrinted = display->GetCurrency();
                const_cast<Currency&>(curBeingPrinted).SetDisplay(display);
            
                std::string out = curBeingPrinted.GetOutputString(_interp->GetOutputFormat());
                bool        isprintoutput = curBeingPrinted.IsPrintfOutput();
                PrintToConsole(out, isprintoutput);

                CurrencyDisplay* topdisplay = GetCurrentDisplay();
                assert(topdisplay);
                if (topdisplay != display)
                {
                    // Nested pagination
                    Currency newcur = const_cast<Currency&>(topdisplay->GetCurrency());
                    newcur.SetDisplay(topdisplay);
                    display = topdisplay;
                }

                if (topdisplay->IsPaginatingCols() || topdisplay->IsPaginatingRows())
                {
                    PrintPaginationMessage(topdisplay->CanPaginateColumns());
                    return;
                }
                
            }
            if (display && !(display->IsPaginatingCols() || display->IsPaginatingRows()))
            {
                EndPagination(true);
            }
        }
    }

    // Process all the other results
    for (std::vector<Currency>::iterator itr = _resultsToPrint.begin();
         itr != _resultsToPrint.end();)
    {
        if (IsPaginating()) break;

        EndPagination(true);

        bool couldPrint = PrintResult(*itr);
        if (!couldPrint)
            ++itr;
        else
            itr = _resultsToPrint.erase(itr);
    }
}
//------------------------------------------------------------------------------
// Paginates matrix and returns true if user has not pressed quit 
//------------------------------------------------------------------------------
void ConsoleWrapper::Paginate()
{
    assert(_interp);
    assert(!_displayStack.empty());

#ifdef OS_WIN
	HANDLE consoleInput = GetStdHandle(STD_INPUT_HANDLE); // Console input buffer
    if (consoleInput == INVALID_HANDLE_VALUE)             // Something is wrong
    {
        CONSOLEWRAPPER_PRINT("Error getting console input\n");
        HandleOnClearResults();
        return;
    }

    // Go in a loop till we either get the pagination controls or a quit
    while(1)
    {   
        DWORD        charToRead = 1;  // Num chars to read
        DWORD        numEvents  = 0;  // Num events read
        INPUT_RECORD input;	          // Input buffer data
   
		// Stop from echoing.
		FlushConsoleInputBuffer(consoleInput);
	    BOOL read = ReadConsoleInput(consoleInput, &input, charToRead, &numEvents);
        if (!read)
		{
			CONSOLEWRAPPER_PRINT("Error reading console input\n");
			HandleOnClearResults();
			return;
		}

        if (_interp->IsInterrupt()) return;

		if (input.EventType != KEY_EVENT || !input.Event.KeyEvent.bKeyDown)  
			continue; // Not a key event with key down

        CurrencyDisplay* display = _displayStack.top();
        assert(display);
        bool canPaginateCols = display->CanPaginateColumns();

        WORD key = input.Event.KeyEvent.wVirtualKeyCode;

		if (key == 0x46 || key == 0x66)      // (f)orward pagination
			display->SetMode(CurrencyDisplay::DISPLAYMODE_FORWARD);
		
		else if (key == 0x51 || key == 0x71)      // (q)uit
        {
            DeleteChainedDisplays(display);
			display->SetMode(CurrencyDisplay::DISPLAYMODE_QUIT);
        }

		else if (key == 0x42 || key == 0x62)      // (b)ack pagination
			display->SetMode(CurrencyDisplay::DISPLAYMODE_BACK);

		else if (key == VK_RIGHT)                // right pagination
        {
            if (canPaginateCols)
            {
		        display->SetMode(CurrencyDisplay::DISPLAYMODE_RIGHT);
            }
            else
            {
                display->SetMode(CurrencyDisplay::DISPLAYMODE_FORWARD);
            }
        }
		else if (key == VK_LEFT)                // left pagination
        {
            if (canPaginateCols)
            {
		        display->SetMode(CurrencyDisplay::DISPLAYMODE_LEFT);
            }
            else
            {
                display->SetMode(CurrencyDisplay::DISPLAYMODE_BACK);
            }
            break;
        }

		else if (key == VK_UP)                // Up pagination
			display->SetMode(CurrencyDisplay::DISPLAYMODE_UP);

		else if (key == VK_DOWN)                // Down pagination
			display->SetMode(CurrencyDisplay::DISPLAYMODE_DOWN);

	    else if (key == 0x45 || key == 0x65)    // (e)xit
        {
            DeleteChainedDisplays(display);
		    display->SetMode(CurrencyDisplay::DISPLAYMODE_EXIT);
        }

        // Delete the last line has the pagination message
        std::string deletedline ("\r");
        for (int i = 1; i < maxcols - 1; ++i)
        {
            deletedline += " ";
        }
        std::cout << deletedline;
		
		ProcessPagination();
        bool isPaginating = IsPaginating();
        if (!isPaginating) 
        {
            return;    // Done paginating
        }
    }
#else
    // Get the terminal settings so there is no echo of input during 
    // pagination and user does not have to press the enter key

    termios savedt;
    tcgetattr(STDIN_FILENO, &savedt);
    termios tmpt = savedt;
    tmpt.c_lflag &= ~(ICANON | ECHO);  // Disable echo and waiting for EOL
    tcsetattr(STDIN_FILENO, TCSANOW, &tmpt);

    while (1)
    {        
        if (_interp->IsInterrupt()) break;

        CurrencyDisplay* display = _displayStack.top();
        assert(display);
        if (!display) break;

        int key=getchar();
        if (key == EOF)
            break;

		if (key == 0x46 || key == 0x66)           // (f)orward pagination
			display->SetMode(CurrencyDisplay::DISPLAYMODE_FORWARD);
		
		else if (key == 0x51 || key == 0x71)      // (q)uit
			display->SetMode(CurrencyDisplay::DISPLAYMODE_QUIT);

		else if (key == 0x42 || key == 0x62)      // (b)ack pagination
			display->SetMode(CurrencyDisplay::DISPLAYMODE_BACK);

	    else if (key == 0x45 || key == 0x65)    // (e)xit
		    display->SetMode(CurrencyDisplay::DISPLAYMODE_EXIT);

		else continue;
		
		ProcessPagination();

        bool isPaginating = IsPaginating();

        if (!isPaginating) break;    // Done paginating
    }
    tcsetattr(STDIN_FILENO, TCSANOW, &savedt);
#endif
}
//------------------------------------------------------------------------------
// Slot called when interpreter is deleted
//------------------------------------------------------------------------------
void ConsoleWrapper::HandleOnSaveOnExit()
{	
    DisconnectOmlSignalHandler();

    delete _interp;
    _interp = NULL;
    exit(EXIT_SUCCESS);
	
}
//------------------------------------------------------------------------------
// Slot which displays a prompt and gets user input
//------------------------------------------------------------------------------
void ConsoleWrapper::HandleOnGetUserInput(const std::string& prompt,
                                          const std::string& type,
										  std::string&       userInput)
{
    bool waspaginating = false;
    while (1)
    {
        if (IsPaginating())
        {
            waspaginating = true;
            Paginate();
            continue;
        }
        break;  // Nothing pending at this point
    }
    if (waspaginating)
        PrintToStdout(">>>");   

    fflush(stdout);
    std::cout << prompt;
    std::getline(std::cin, userInput);

    if (type == "s" || type == "S") 
    {
        userInput += "'";
        userInput.insert(0, "'");
    }
}
//------------------------------------------------------------------------------
// Sets pagination variables
//------------------------------------------------------------------------------
void ConsoleWrapper::SetPaginationVariables()
{
    // Environment variable takes precedence
    if (BuiltInFuncsUtils::IsPaginationEnvEnabled())
    {
        GetCommandWindowInfo();
    }
    else
    {
        CurrencyDisplay::SetMaxCols(0);
        CurrencyDisplay::SetMaxRows(0);
    }
}
//------------------------------------------------------------------------------
// Initializes and caches display for the given currency
//------------------------------------------------------------------------------
void ConsoleWrapper::CacheDisplay(const Currency& cur)
{
    if (!_displayStack.empty()) return;

    CurrencyDisplay* display = cur.GetDisplay();
    assert(display);

    assert(_interp);

    display->Initialize(_interp->GetOutputFormat(), _interp);

    _displayStack.push(display);  // This becomes the new top display
}
//------------------------------------------------------------------------------
// Clears display
//------------------------------------------------------------------------------
void ConsoleWrapper::ClearDisplay()
{
    if (_displayStack.empty()) return;
    CurrencyDisplay* display = GetCurrentDisplay();
    if (display)
    {
        CurrencyDisplay::DeleteDisplay(display);
        _displayStack.pop();
    }
    if (_displayStack.empty())
    {
        CurrencyDisplay::ClearLineCount();
    }
}
//------------------------------------------------------------------------------
// Gets current display
//------------------------------------------------------------------------------
CurrencyDisplay* ConsoleWrapper::GetCurrentDisplay() const
{
    if (_displayStack.empty()) return 0;

    CurrencyDisplay* display = _displayStack.top();
    assert(display);
    
    return display;
}
//------------------------------------------------------------------------------
// Returns true if pagination is in process
//------------------------------------------------------------------------------
bool ConsoleWrapper::IsPaginating()
{
    if (_displayStack.empty()) return false;

    CurrencyDisplay* display = GetCurrentDisplay();
    assert(display);

    if (display->IsPaginatingCols() || display->IsPaginatingRows()) return true;

    return false;
}
//------------------------------------------------------------------------------
// Slot for when a new nested display needs to be added
//------------------------------------------------------------------------------
void ConsoleWrapper::HandleOnAddDisplay(CurrencyDisplay* display)
{
    if (!display) return;

    assert(_interp);

    display->Initialize(_interp->GetOutputFormat(), _interp, GetCurrentDisplay());

    _displayStack.push(display);
}
//------------------------------------------------------------------------------
// Slot called when initiating a user defined pause. Any character typed will
// break out of the pause
//------------------------------------------------------------------------------
void ConsoleWrapper::HandleOnPauseStart(const std::string& msg, bool wait)
{
    bool waspaginating = false;
    while (1)
    {
        if (IsPaginating())
        {
            waspaginating = true;
            Paginate();
            continue;
        }
        break;  // Nothing pending at this point
    }
    PrintInfoToOmlWindow(msg);
    
    if (!wait)
    {
        return;
    }

#ifdef OS_WIN
	HANDLE consoleInput = GetStdHandle(STD_INPUT_HANDLE); // Console input buffer
    if (consoleInput == INVALID_HANDLE_VALUE)             // Something is wrong
    {
        CONSOLEWRAPPER_PRINT("Error getting console input\n");
        HandleOnClearResults();
        return;
    }

    // Go in a loop till we either get the pagination controls or a quit
    while(1)
    {        
        DWORD        charToRead = 1;  // Num chars to read
        DWORD        numEvents  = 0;  // Num events read
        INPUT_RECORD input;	          // Input buffer data
   
		// Stop from echoing.
		FlushConsoleInputBuffer(consoleInput);
	    BOOL read = ReadConsoleInput(consoleInput, &input, charToRead, &numEvents);
        if (!read)
		{
			CONSOLEWRAPPER_PRINT("Error reading console input\n");
			HandleOnClearResults();
			return;
		}

        if (_interp->IsInterrupt()) return;

		if (input.EventType != KEY_EVENT || !input.Event.KeyEvent.bKeyDown)  
			continue; // Not a key event with key down

        break;
    }
#else

	fflush(stdout);
    std::string userInput;
    std::getline(std::cin, userInput);
#endif
}
//------------------------------------------------------------------------------
// Cleans up after pagination like clearing display
//------------------------------------------------------------------------------
void ConsoleWrapper::EndPagination(bool printMsg)
{
    if (_displayStack.empty()) return;

    CurrencyDisplay* display = GetCurrentDisplay();
    if (display)
    {
        if (printMsg)
            PrintPaginationMessage();

        CurrencyDisplay::DeleteDisplay(display);
        _displayStack.pop();
    }
}
//------------------------------------------------------------------------------
// Prints info message (in different color than results) to command window
//------------------------------------------------------------------------------
void ConsoleWrapper::PrintInfoToOmlWindow(const std::string& msg)
{
#ifdef OS_WIN
    HANDLE hConsoleOutput = GetStdHandle(STD_OUTPUT_HANDLE);
    WORD   oldColorAttrs;
    bool   setColorAttrs = false;
	if (hConsoleOutput != INVALID_HANDLE_VALUE)
    {
	    CONSOLE_SCREEN_BUFFER_INFO csbiInfo; 
	    if (GetConsoleScreenBufferInfo(hConsoleOutput, &csbiInfo))
        {
            setColorAttrs = true;
    	    oldColorAttrs = csbiInfo.wAttributes; // Cache color info

            // Gray background with black letters
            SetConsoleTextAttribute(hConsoleOutput, BACKGROUND_INTENSITY | 0);
        }
    }
    std::cout << msg;

    if (setColorAttrs)
	    SetConsoleTextAttribute(hConsoleOutput, oldColorAttrs); // Reset color info
    std::cout << std::endl;
#else
    std::cout << msg << std::endl;
#endif
}
//------------------------------------------------------------------------------
// Prints info message, silent in quiet mode
//------------------------------------------------------------------------------
void ConsoleWrapper::PrintToStdout(const std::string& msg)
{
    if (!_quietMode)
        std::cout << msg;
}
//------------------------------------------------------------------------------
// Prints new prompt, resets append flags
//------------------------------------------------------------------------------
void ConsoleWrapper::PrintNewPrompt()
{
    if (_appendOutput)
    {
        std::cout << std::endl;    // Add newline, in case there was a printf
        _appendOutput = false;     // Reset printf output
    }
    else if (_addnewline)
    {
        std::cout << std::endl;    // Printf with newline
        _addnewline = false;
    }

    PrintToStdout(">>>");        // New prompt
}
//------------------------------------------------------------------------------
// Gets argv at the given index
//------------------------------------------------------------------------------
std::string ConsoleWrapper::GetArgv(int idx) const
{
    if (idx < 0 || _argv.empty() || idx >= static_cast<int>(_argv.size()))
        return "";

    return _argv[idx];
}
//------------------------------------------------------------------------------
// Gets argc
//------------------------------------------------------------------------------
int ConsoleWrapper::GetArgc() const
{
    return (_argv.empty()) ? 0 : static_cast<int>(_argv.size());
}
//------------------------------------------------------------------------------
// Deletes chained displays, if any
//------------------------------------------------------------------------------
void ConsoleWrapper::DeleteChainedDisplays(CurrencyDisplay* display)
{
    if (!display || _resultsToPrint.empty())
    {
        return;
    }

    for (std::vector<Currency>::iterator itr = _resultsToPrint.begin();
         itr != _resultsToPrint.end();)
    {
        Currency cur = *itr;
        CurrencyDisplay* curdisp = cur.GetDisplay();
        if (display->IsChainedDisplay(curdisp))
        {
            CurrencyDisplay::DeleteDisplay(curdisp);
            cur.SetDisplay(nullptr);
            itr = _resultsToPrint.erase(itr);
        }
        else
        {
            ++itr;
        }
    }
}