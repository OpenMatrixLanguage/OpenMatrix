/**
* @file ConsoleWrapper.cpp
* @date June 2015
* Copyright (C) 2015-2022 Altair Engineering, Inc.  
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

#ifndef OS_WIN
#   include <sys/ioctl.h>
#   include <termios.h>
#   include <unistd.h>
#endif

#include "Runtime/BuiltInFuncsUtils.h"
#include "Runtime/BuiltInFuncsCore.h"
#include "Runtime/CurrencyDisplay.h"
#include "Runtime/Interpreter.h"
#include "Runtime/OutputFormat.h"

//# define CONSOLEWRAPPER_DBG 1  // Uncomment to print debug info
#ifdef CONSOLERWRAPPER_DBG
#    define CONSOLEWRAPPER_PRINT(m) { std::cout << m; }
#else
#    define CONSOLEWRAPPER_PRINT(m) 0
#endif

int maxcols = 0;
#ifndef OS_WIN
enum ARROWKEY
{
    ARROWKEY_UP = 65,
    ARROWKEY_LEFT = 68,
    ARROWKEY_RIGHT = 67,
    ARROWKEY_DOWN = 66
};
#endif

const size_t g_buffsize = 256 * 1024;
static char  g_buf[g_buffsize + 1];
// Utility to flush output buffer
void FlushStdout()
{
    fflush(stdout);
    memset(g_buf, 0, sizeof(g_buf));
}

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
ConsoleWrapper::ConsoleWrapper(Interpreter* interp)
    : WrapperBase(interp)
    , _appendOutput     (false)
    , _addnewline       (false)
    , _getWindowSize    (false)
    , _childHandler     (nullptr)
{
    memset(g_buf, 0, sizeof(g_buf));
    std::cout.rdbuf()->pubsetbuf(g_buf, g_buffsize);
    std::cin.tie(nullptr);
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
ConsoleWrapper::~ConsoleWrapper()
{
    BuiltInFuncsUtils::CloseOutputLog();

    std::cout.rdbuf()->pubsetbuf(nullptr, 0);
}
//------------------------------------------------------------------------------
// Slot called when printing result to console
//------------------------------------------------------------------------------
void ConsoleWrapper::HandleOnPrintResult(const Currency& cur)
{	
    bool islog = cur.IsLogOutput();
    bool closelog = false;
    if (islog && !BuiltInFuncsUtils::IsOutputLogOpen() && _displayStack.empty())
    {
        BuiltInFuncsUtils::OpenOutputLogForAppend(); // This currency is coming from a child evaluator
        closelog = true;
    }

    // Force to paginate on for large matrices so that output is displayed
    // faster and the output string is not too big
    if (cur.IsMatrix())
    {
        const hwMatrix* mtx = cur.Matrix();
        if (mtx && mtx->Size() >= CurrencyDisplay::GetSkipFormat())
        {
            CurrencyDisplay* disp = cur.GetDisplay();
            if (disp)
            {
                disp->SetLocalPaginate();
            }
        }
    }

    bool couldPrint = PrintResult(cur);
    if (closelog)
    {
        BuiltInFuncsUtils::CloseOutputLog();
    }

    if (!couldPrint)
        _resultsToPrint.push_back(cur);  // Print later, once pagination is done
}
//------------------------------------------------------------------------------
// Gets visible rows and columns from command window screen size
//------------------------------------------------------------------------------
void ConsoleWrapper::GetCommandWindowInfo()
{
    if (!_getWindowSize)
    {
        return;
    }
    int rows = 0;
    int cols = 0;

    maxcols = 0;
#ifdef OS_WIN  // Windows specific code
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(_outhandle, &csbi);
    // Lines and columns use the screen size and not the buffer size
    rows = csbi.srWindow.Bottom - csbi.srWindow.Top;
    cols = csbi.srWindow.Right  - csbi.srWindow.Left;
#else
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    rows = w.ws_row;
    cols = w.ws_col;
    if (rows > 2)
    {
        rows--;  // Reduce one row to account for pagination message
    }
#endif    

    maxcols = cols;
    CurrencyDisplay::SetMaxCols(cols);
	CurrencyDisplay::SetMaxRows(rows);
    _getWindowSize = false;
}
//------------------------------------------------------------------------------
// Print currency. Return false if this currency could not be printed and 
// needs to be processed later
//------------------------------------------------------------------------------
bool ConsoleWrapper::PrintResult(const Currency& cur)
{
    assert(_interp);
    const OutputFormat* format = _interp->GetOutputFormat();

    if (!_displayStack.empty())
    {
        return false;  // Can't print yet as display stack is not empty
    }

    bool canPaginate = CurrencyDisplay::CanPaginate(cur);
    bool forceflush = (cur.IsPrintfOutput() || cur.IsDispOutput() || cur.IsError());
    if (!canPaginate)
    {
        // There is no print in progress, so don't need to worry
        PrintToConsole(cur.GetOutputString(format), cur.IsPrintfOutput(), forceflush);
        return true;
    }
    GetCommandWindowInfo();
    if (!CurrencyDisplay::IsValidDisplaySize()) // No window size available so just print
    {
        PrintToConsole(cur.GetOutputString(format), cur.IsPrintfOutput(), forceflush);
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
    PrintToConsole(cur.GetOutputString(format), cur.IsPrintfOutput(), forceflush);

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
    if (CurrencyDisplay::IsPaginateInteractive())
    {
        std::string msg("-- (f)orward, (b)ack, (q)uit, (e)xit");

        msg += ", (up), (down)";
        if (showColControl)
        {
            msg += ", (left), (right)";
        }

       FlushStdout();

#ifdef OS_WIN
        // Gray background with white letters
        SetConsoleTextAttribute(_outhandle, BACKGROUND_INTENSITY |
            FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);

        Print(msg, true);
        FlushStdout();

        SetConsoleTextAttribute(_outhandle, _defaultWinAttr); // Reset color info
#else
        Print(msg, true);
#endif
    }
}
//------------------------------------------------------------------------------
// Prints end of pagination message
//------------------------------------------------------------------------------
void ConsoleWrapper::PrintPaginationMessage()
{
    CurrencyDisplay* display = GetCurrentDisplay();
    if (display)
    {
        std::string strmsg;
        display->GetPaginationEndMsg(strmsg);
        if (!strmsg.empty())
        {
            Print(strmsg + '\n', true);
        }
    }
}
//------------------------------------------------------------------------------
// Prints string to console
//------------------------------------------------------------------------------
void ConsoleWrapper::PrintToConsole(const std::string& result,
                                    bool               isprintoutput,
                                    bool               forceflush)
{	
    if (!result.empty())
    {
        std::string msg;
        if ((!_appendOutput && _addnewline) || (!isprintoutput && _appendOutput))
        {
            msg += '\n';
        }

        bool hasTrailingNewline = (result.back() == '\n');
        if (!isprintoutput)
        {
            msg += result;
            if (!hasTrailingNewline)
            {
                msg += '\n';
            }
            _appendOutput = false;
            _addnewline   = false;
            PrintToStdout(msg);
            if (forceflush)
            {
                FlushStdout();
            }

            return;
        }

        if (!hasTrailingNewline)
        {
            msg += result;
        }
        else
        {
            std::string tmp(result);
            tmp.pop_back();
            msg += tmp;
            _addnewline = true;
        }
        PrintToStdout(msg);
        if (forceflush)
        {
            FlushStdout();
        }

        _appendOutput = !hasTrailingNewline;
    }
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
    if (!display)
    {
        return;
    }
    else if (_interp->IsInterrupt())
    {
        HandleOnClearResults();   // Quit all printing
        return;
    }

    if (display->GetMode() == CurrencyDisplay::DISPLAYMODE_EXIT)
    {
        if (CurrencyDisplay::IsPaginateInteractive())
        {
            std::cout << '\r' << '\n';
            BuiltInFuncsUtils::SaveToOutputLog("\n");
        }
        HandleOnClearResults();   // Quit all printing
    }
    else
    {
        bool printmsg = true;
        if (display->GetMode() != CurrencyDisplay::DISPLAYMODE_QUIT)
        {
	        // Print to console
            display->SetModeData();

            const Currency& curBeingPrinted = display->GetCurrency();
            const_cast<Currency&>(curBeingPrinted).SetDisplay(display);

            std::string out = curBeingPrinted.GetOutputString(_interp->GetOutputFormat());
            if (CurrencyDisplay::IsPaginateInteractive())
            {
                if (display->GetDeleteLine())
                {
                    std::cout << '\r';
                    display->SetDeleteLine(false);
                }
                else
                {
                    std::cout << '\n';
                    BuiltInFuncsUtils::SaveToOutputLog("\n");
                }
            }
            bool isPrintOutput = curBeingPrinted.IsPrintfOutput();
            bool forceflush = (curBeingPrinted.IsDispOutput() || 
                               curBeingPrinted.IsError()      || 
                               curBeingPrinted.IsPrintfOutput());
            // Reget display as there could be nested pagination
            CurrencyDisplay* topdisplay = GetCurrentDisplay();
            assert(topdisplay);
            if (topdisplay != display)
            {
                Currency newcur = const_cast<Currency&>(topdisplay->GetCurrency());
                newcur.SetDisplay(topdisplay);
            }

            PrintToConsole(out, isPrintOutput, forceflush);
	        if (topdisplay->IsPaginatingCols() || topdisplay->IsPaginatingRows()) 
            {
                PrintPaginationMessage(topdisplay->CanPaginateColumns());
                return; // Still printing
            }
        }
        else
        {
            if (CurrencyDisplay::IsPaginateInteractive())
            {
                std::cout << '\r' << '\n';
                display->SetDeleteLine(false);
            }
            printmsg = false;
        }
        EndPagination(printmsg);  // Done with pagination

        if (!_displayStack.empty())
        {
            display = _displayStack.top();
            if (display && display->IsPaginatingCols() || display->IsPaginatingRows())
            {
                if (!display->CanPrintRows() || CurrencyDisplay::IsPaginateInteractive())
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
                bool forceflush = (curBeingPrinted.IsDispOutput() ||
                    curBeingPrinted.IsError() ||
                    curBeingPrinted.IsPrintfOutput());

                PrintToConsole(out, isprintoutput, forceflush);

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
                EndPagination(printmsg);
            }
            if (!_displayStack.empty())
            {
                if (CurrencyDisplay::IsPaginateInteractive())
                {
                    display = GetCurrentDisplay();
                    if (display && (
                        (display->IsPaginatingRows() || display->IsPaginatingCols()) ||
                        display->IsNDCellDisplay()))
                    {
                        PrintPaginationMessage(display->CanPaginateColumns());
                        return;               // Done printing for pagination
                    }
                }
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
    // Go in a loop till we either get the pagination controls or a quit
    while(1)
    {   
        CurrencyDisplay* display = _displayStack.top();
        bool canPaginateCols = display ? display->CanPaginateColumns() : false;

        if (CurrencyDisplay::IsPaginateInteractive())
        {
            DWORD        charToRead = 1;  // Num chars to read
            DWORD        numEvents = 0;  // Num events read
            INPUT_RECORD input;	          // Input buffer data

            // Stop from echoing.
            FlushConsoleInputBuffer(_inhandle);
            BOOL read = ReadConsoleInput(_inhandle, &input, charToRead, &numEvents);
            if (!read)
            {
                CONSOLEWRAPPER_PRINT("Error reading console input\n");
                HandleOnClearResults();
                return;
            }

            if (_interp->IsInterrupt()) return;

            if (input.EventType != KEY_EVENT || !input.Event.KeyEvent.bKeyDown)
                continue; // Not a key event with key down

            WORD key = input.Event.KeyEvent.wVirtualKeyCode;

            if (key == 0x46 || key == 0x66 || key == 0x0D)      // (f)orward pagination
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

            else
            {
                // Treat any other key as a quit command
                DeleteChainedDisplays(display);
                display->SetMode(CurrencyDisplay::DISPLAYMODE_QUIT);
            }

            // Delete the last line has the pagination message
            std::string deletedline("\r");
            for (int i = 1; i < maxcols - 1; ++i)
            {
                deletedline += " ";
            }
            std::cout << deletedline;
        }
        else
        {
            display->SetMode(CurrencyDisplay::DISPLAYMODE_FORWARD);
        }
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

    struct termios savedt;
    tcgetattr(STDIN_FILENO, &savedt);
    struct termios tmpt = savedt;
    tmpt.c_lflag &= ~(ICANON | ECHO);  // Disable echo and waiting for EOL
    tcsetattr(STDIN_FILENO, TCSANOW, &tmpt);

    while (1)
    {        
        if (_interp->IsInterrupt()) break;

        CurrencyDisplay* display = _displayStack.top();
        assert(display);
        if (!display) break;

        bool canPaginateCols = display ? display->CanPaginateColumns() : false;

        if (CurrencyDisplay::IsPaginateInteractive())
        {
            int key = getchar();
            if (key == EOF)
                break;

            if (key == 0x46 || key == 0x66)           // (f)orward pagination
                display->SetMode(CurrencyDisplay::DISPLAYMODE_FORWARD);

            else if (key == 0x51 || key == 0x71)      // (q)uit
            {
                DeleteChainedDisplays(display);
                display->SetMode(CurrencyDisplay::DISPLAYMODE_QUIT);
            }

            else if (key == 0x42 || key == 0x62)      // (b)ack pagination
                display->SetMode(CurrencyDisplay::DISPLAYMODE_BACK);

            else if (key == 0x45 || key == 0x65)    // (e)xit
            {
                DeleteChainedDisplays(display);
                display->SetMode(CurrencyDisplay::DISPLAYMODE_EXIT);
            }
            else if (key == 27) // ANSI escape sequence for arrow keys
            {
                key = getchar();
                if (key == 91)
                {
                    key = getchar();
                }

                switch (key)
                {
                    case ARROWKEY_UP:
                        display->SetMode(CurrencyDisplay::DISPLAYMODE_UP);
                        break;

                    case ARROWKEY_DOWN:
                        display->SetMode(CurrencyDisplay::DISPLAYMODE_DOWN);
                        break;

                    case ARROWKEY_LEFT:
                        if (canPaginateCols)
                        {
                            display->SetMode(CurrencyDisplay::DISPLAYMODE_LEFT);
                        }
                        else
                        {
                            display->SetMode(CurrencyDisplay::DISPLAYMODE_BACK);
                        }
                        break;

                    case ARROWKEY_RIGHT:
                        if (canPaginateCols)
                        {
                            display->SetMode(CurrencyDisplay::DISPLAYMODE_RIGHT);
                        }
                        else
                        {
                            display->SetMode(CurrencyDisplay::DISPLAYMODE_FORWARD);
                        }
                        break;

                    default: // Treat any other key as quit
                        DeleteChainedDisplays(display);
                        display->SetMode(CurrencyDisplay::DISPLAYMODE_QUIT);
                        break;
                }
            }
            else
            {
                // Treat any other key as quit
                DeleteChainedDisplays(display);
                display->SetMode(CurrencyDisplay::DISPLAYMODE_QUIT);
            }
        }
        else
        {
            display->SetMode(CurrencyDisplay::DISPLAYMODE_FORWARD);
        }
		
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
void ConsoleWrapper::HandleOnSaveOnExit(int returnCode)
{	
    DisconnectOmlSignalHandler();

    // Do not delete the interpreter. If the exit command is in a function, this
    // causes a crash as the function is not complete and the memory scope is still 
    // open. This is not a solution, just a band-aid fix. We should ideally be
    // setting a flag here and stop further execution of commands
    //delete _interp;
    //_interp = NULL;

    exit(returnCode);
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
    {
        std::cout << ">>>";
        FlushStdout();
    }

    std::cout << prompt;
    FlushStdout();
    std::getline(std::cin, userInput);

    if (type == "s" || type == "S") 
    {
        userInput += "'";
        userInput.insert(0, "'");
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
    //PrintInfoToOmlWindow(msg); // Printed earlier
    
    if (!wait)
    {
        return;
    }

#ifdef OS_WIN
    // Go in a loop till we either get the pagination controls or a quit
    while(1)
    {        
        DWORD        charToRead = 1;  // Num chars to read
        DWORD        numEvents  = 0;  // Num events read
        INPUT_RECORD input;	          // Input buffer data
   
		// Stop from echoing.
		FlushConsoleInputBuffer(_inhandle);
	    BOOL read = ReadConsoleInput(_inhandle, &input, charToRead, &numEvents);
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
    // Gray background with black letters
    SetConsoleTextAttribute(_outhandle, BACKGROUND_INTENSITY | 0);

    Print(msg + '\n', true);
    SetConsoleTextAttribute(_outhandle, _defaultWinAttr); // Reset color info
#else
    Print(msg + '\n', true);
#endif
}
//------------------------------------------------------------------------------
// Prints info message
//------------------------------------------------------------------------------
void ConsoleWrapper::PrintToStdout(const std::string& msg)
{
    size_t len = msg.size();
    if (len < g_buffsize)
    {
        Print(msg, false);
        return;
    }

    size_t numprinted = 0;
    while (numprinted < len)
    {
        FlushStdout(); // Flush before printing large output
        if (numprinted + g_buffsize > len)
        {
            Print(msg.substr(numprinted), false);
            break;
        }
        else
        {
            Print(msg.substr(numprinted, g_buffsize), false);
            numprinted += g_buffsize;
        }
    }
    FlushStdout();
}
//------------------------------------------------------------------------------
// Prints new prompt, resets append flags
//------------------------------------------------------------------------------
void ConsoleWrapper::PrintNewPrompt()
{
    std::cout << std::flush;

    std::string msg;
    if (_appendOutput)
    {
        msg = '\n';    // Add newline, in case there was a printf
        _appendOutput = false;     // Reset printf output
    }
    else if (_addnewline)
    {
        msg = '\n';    // Printf with newline
        _addnewline = false;
    }
    BuiltInFuncsUtils::SaveToOutputLog(msg);
    msg += ">>>";
    std::cout << msg << std::flush;        // New prompt
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
//------------------------------------------------------------------------------
// Initializes command window screen buffer for interactive mode
//------------------------------------------------------------------------------
void ConsoleWrapper::InitCommandWindowInfo()
{
#ifdef OS_WIN
    _inhandle  = GetStdHandle(STD_INPUT_HANDLE); 
    _outhandle = GetStdHandle(STD_OUTPUT_HANDLE);

    CONSOLE_SCREEN_BUFFER_INFO csbiInfo;
    GetConsoleScreenBufferInfo(_outhandle, &csbiInfo);
    _defaultWinAttr = csbiInfo.wAttributes; // Cache color info

#endif
}
//------------------------------------------------------------------------------
// Prints to stdout
//------------------------------------------------------------------------------
void ConsoleWrapper::Print(const std::string& msg, bool forceFlush)
{
    if (!_childHandler)  // Main interp, print to stdout
    {
        std::cout << msg;
        if (forceFlush)
        {
            FlushStdout();
        }
        BuiltInFuncsUtils::SaveToOutputLog(msg);
    }
}
