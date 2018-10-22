/**
* @file ConsoleWrapper.h
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

#ifndef __CONSOLE_WRAPPER__
#define __CONSOLE_WRAPPER__

// Begin defines/includes

#include "WrapperBase.h"

#include "Runtime/Currency.h"

#include <stack>
#include <string>
#include <vector>

class CurrencyDisplay;
class Interpreter;
// End defines/includes

//------------------------------------------------------------------------------
//! 
//! \class ConsoleWrapper
//! \brief Handles interpretor signals and methods for console
//!
//------------------------------------------------------------------------------
class ConsoleWrapper : public WrapperBase
{
public:
    //!
    //! Constructor
    //! \param[in] interpretor
    //!
    ConsoleWrapper(Interpreter* interp);
    //!
    //! Destructor
    //!
    virtual ~ConsoleWrapper();

    //!
    //! Returns true if pagination is in process
    //!
    bool IsPaginating();
    //!
    //! Paginates matrix
    //!
    void Paginate();
    //!
    //! Sets pagination variables
    //!
    void SetPaginationVariables();

    //!
    //! True if quiet mode is enabled
    //!
    bool GetQuietMode() const { return _quietMode; }
    //!
    //! Sets quiet mode
    //! \param val True if quiet mode is enabled
    //!
    void SetQuietMode(bool val) { _quietMode = val; }

    //!
    //! Prints new prompt and resets append flags
    //!
    void PrintNewPrompt();
    //!
    //! Prints to stdout if not in quiet mode
    //! \param msg Message to print
    //!
    void PrintToStdout(const std::string& msg);

    //!
    //! Prints currency to console
    //! \param cur Currency to print
    //!
	void HandleOnPrintResult(const Currency& cur);
    //!
    //! Clears results and pagination related to results
    //!
    void HandleOnClearResults();
    //!
    //! Handles application exit
    //!
	void HandleOnSaveOnExit();
    //!
    //! Displays a prompt and gets user input
	//! \param prompt    Prompt to display to user
    //! \param type      Type, if specified
	//! \param userInput Input from user
    //!
    void HandleOnGetUserInput(const std::string& prompt,
                              const std::string& type,
		                      std::string&       userInput);
    //!
    //! Initiates a user defined pause
    //! \param msg  User message to display
    //! \param wait True if waiting for a keystroke input from user
    //!
	void HandleOnPauseStart(const std::string& msg, 
                            bool               wait);
    //!
    //! Handles nested displays during pagination
    //! \param display Display to be added
    //!
    void HandleOnAddDisplay(CurrencyDisplay* display);

    //!
    //! Gets argc
    //!
    int GetArgc() const;
    //!
    //! Gets argv at the given index
    //! \param idx 1-based index
    //!
    std::string GetArgv(int idx) const;
    //!
    //! Sets argv
    //! \param args Args
    //!
    void SetArgv(const std::vector<std::string>& args) { _argv = args; }

private:
    bool _appendOutput;                           //! True if output is appended
    bool _quietMode;                              //! True if in quiet mode
    bool _addnewline;                             //! True if new line should be added

    std::vector<Currency>        _resultsToPrint; //! Results to print
    std::stack<CurrencyDisplay*> _displayStack;   //! Stack of currenct displays

    std::vector<std::string> _argv;               //!< Arg v

    // Private default constructor
    ConsoleWrapper() : WrapperBase(nullptr) {}

    //! True if pagination environment is enabled
    bool IsPaginationEnvEnabled() const;

    //!
    //! Cleans up after pagination like clearing display
    //! \param printMsg True if end of pagination message needs to be printed
    //!
    void EndPagination(bool printMsg = true);
    //!
    //! Processes pagination
    //!
    void ProcessPagination();
    //!
    //! Gets visible rows and columns from command window screen size
    //!
    void GetCommandWindowInfo() const;
    //!
    //! Print currency. Return false if this currency could not be printed and 
    //! needs to be processed later
    //! \param cur Currency (result) to print
    //!
    bool PrintResult(const Currency& cur);
    //!
    //! Prints string to console
    //! \param result        Result
    //! \param isprintoutput True if this is a result from printf/fprintf
    //!
    void PrintToConsole(const std::string& result,
                        bool               isprintOutput);
    //!
    //! Prints message for continuing/skipping pagination
    //! \param colPaginate True if column control keys are dispalyed
    //!
    void PrintPaginationMessage(bool showColControl);
    //!
    //! Prints end of pagination message
    //!
    void PrintPaginationMessage();
    //!
    //! Prints info message (in different color than results) to command window
    //! \param msg  Message to print
    //!
    void PrintInfoToOmlWindow(const std::string& msg);
    //!
    //! Initializes and caches display given the currency
    //! \param cur Given currency
    //!
    void CacheDisplay(const Currency& cur);
    //!
    //! Clears display
    //!
    void ClearDisplay();
    //!
    //! Gets current display
    //!
    CurrencyDisplay* GetCurrentDisplay() const;
    //!
    //! Deletes chained displays, if any
    //! \param display Given display
    //!
    void DeleteChainedDisplays(CurrencyDisplay* display);
};


#endif

// End of file:
