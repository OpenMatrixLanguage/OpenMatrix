/**
* @file SignalHandler.cpp
* @date June 2016
* Copyright (C) 2016-2020 Altair Engineering, Inc.  
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
#include "SignalHandler.h"

#include <cassert>

#include "WrapperBase.h"

// End defines/includes

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
SignalHandler::SignalHandler(const SignalHandler* src) 
    : _src (src)
    , _wrapper (NULL)
    , _batchmode(true)
{
    // Set the console wrapper for handling actions from oml core in client
    if (_src)
    {
        _wrapper = _src->_wrapper;
        _batchmode = _src->_batchmode;
	}
}
//------------------------------------------------------------------------------
// Clone signal handler
//------------------------------------------------------------------------------
SignalHandlerBase* SignalHandler::CreateClone()
{
    return (new SignalHandler(this));
}
//------------------------------------------------------------------------------
// Emits signal to clear results
//------------------------------------------------------------------------------
void SignalHandler::OnClearResultsHandler()
{
    // Clear results only once, in the signal handler which is not cloned
    if (!_src && _wrapper)  
        _wrapper->HandleOnClearResults();
}
//------------------------------------------------------------------------------
// Prints prompt for save on exit
//------------------------------------------------------------------------------
void SignalHandler::OnSaveOnExitHandler()
{
    if (_wrapper)  
        _wrapper->HandleOnSaveOnExit();
}
//------------------------------------------------------------------------------
// Start pause
//------------------------------------------------------------------------------
void SignalHandler::OnPauseStartHandler(const std::string& msg, bool wait)
{
    if (_wrapper)  
        _wrapper->HandleOnPauseStart(msg, wait);
}
//------------------------------------------------------------------------------
// Get user input
//------------------------------------------------------------------------------
void SignalHandler::OnUserInputHandler(const std::string& prompt,
                                       const std::string& type,
                                       std::string&       input)
{
    if (_wrapper)  
        _wrapper->HandleOnGetUserInput(prompt, type, input);
}
//------------------------------------------------------------------------------
// Print result
//------------------------------------------------------------------------------
void SignalHandler::OnPrintResultHandler(const Currency& cur)
{
    // Printing is done only by the source console wrapper
    if (!_src && _wrapper)  
        _wrapper->HandleOnPrintResult(cur);
}
//------------------------------------------------------------------------------
// Add nested display
//------------------------------------------------------------------------------
void SignalHandler::OnAddDisplayHandler(CurrencyDisplay* display)
{
     if (!_src && _wrapper)  
        _wrapper->HandleOnAddDisplay(display);
}
// End of file:
