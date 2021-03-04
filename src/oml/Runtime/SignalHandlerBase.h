/**
* @file SignalHandlerBase.h
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

#ifndef __SIGNALHANDLERBASE_H__
#define __SIGNALHANDLERBASE_H__

// Begin defines/includes
#include "OMLDll.h"
#include "Currency.h"
#include "EvaluatorInt.h"
#include "OML_Error.h"

#include <string>

class CurrencyDisplay;
class FunctionInfo;
class Interpreter;
// End defines/includes

//------------------------------------------------------------------------------
//! Base class for handling signals in evaluator
//------------------------------------------------------------------------------
class OMLDLL_DECLS SignalHandlerBase
{
public:
    //!
    //! Destructor
    //!
    virtual ~SignalHandlerBase() {}

    //!
    //! Creates clone
    //!
    virtual SignalHandlerBase* CreateClone() { return 0; }
    //!
    //! Deletes clone
    //!
    virtual void DeleteClone() {}
    //!
    //! Returns true of this is a clone
    //!
    virtual bool IsClone() const { return false; }
    //!
    //! Copies signals for clones, needed only for GUI
    //!
    virtual void CopySignals() {}
    //!
    //! Copies signals for clones that are not related to printing, needed only for GUI
    //! This method should be called in cases where printing is handled separately
    //!
    virtual void CopyNonPrintSignals() {}
    //!
    //! Connects print signal
    //!
    virtual void ConnectPrintSignal() {}
    //!
    //! Disconnect print signal
    //!
    virtual void DisconnectPrintSignal() {}

    //!
    //! Returns true if worker thread is being executed
    //!
    virtual bool IsInWorkerThread() { return true; }
    //!
    //! Returns true if application can execute in GUI thread
    //!
    virtual bool CanExecuteInGuiThread() { return false; }
    //!
    //! Returns result after running given function in GUI thread
    //! \param         Function to run
    //! \param         Evaluator
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    virtual bool RunInGuiThread( FUNCPTR                      fptr,
                                 EvaluatorInterface&          eval, 
                                 const std::vector<Currency>& inputs, 
                                 std::vector<Currency>&       outputs) { return true; }
    //!
    //! Throws oml error - needed for thread safe throw in GUI
    //! \param e Oml error
    //!
    virtual void ThrowThreadSafeError( const OML_Error& e) { throw e; }

    //!
    //! Returns class info
    //!
    virtual const char* ClassName() const { return "SignalHandlerBase"; }

    // Handlers which emit the signals in derived classes from core to client
    // All implementations need to be in the derived classes

    //!
    //! Emits signal to clear results
    //!
    virtual void OnClearResultsHandler() {}
    //!
    //! Emits signal to print result
    //! \param cur Result to print
    //!
    virtual void OnPrintResultHandler( const Currency& cur) {}
    //!
    //! Emits signal to start pause
    //! \param msg  User message to display
    //! \param wait True if waiting for a keystroke input from user
    //!
    virtual void OnPauseStartHandler( const std::string& msg, 
                                      bool               wait) {}
    //!
    //! Emits signal to end pause
    //!
    virtual void OnPauseEndHandler() {}

    //!
    //! Change current working directory
    //! \param dir Fully qualified path of the new directory
    //!
    virtual void OnChangeDirHandler( const std::string& dir) {}
    //!
    //! Refreshes directories in client
    //!
    virtual void OnRefreshDirsHandler() {}

    //!
    //! Add nested display
    //! \param display Nested display to add
    //!
    virtual void OnAddDisplayHandler(CurrencyDisplay* display) {}
    //!
    //! Prompts save on exit in client
    //! \param returnCode Code to exit the application with
    //!
    virtual void OnSaveOnExitHandler(int returnCode) {}
    //!
    //! Update function list in language
    //!
    virtual void OnUpdateFuncListHandler() {}
    //!
    //! Emits signal to get user input from client
    //! \param  prompt Prompt to display to user
    //! \param  type   Type, if specified
    //! \param input   Input from user
    //!
    virtual void OnUserInputHandler( const std::string& prompt,
                                     const std::string& type,
                                     std::string&       input) {}
    //!
    //! Gets current file name that is being updated in the (GUI) editor
    //!
    virtual std::string GetEditorFileName() const { return ""; }
    //!
    //! Sets current file name that is being updated in the (GUI) editor
    //! \param name Name of the file being edited
    //!
    virtual void SetEditorFileName( const std::string& name) {}
    //!
    //! Execute UI call back
    //! \param funinfo handler to function
    //! \param fname   Function name
    //! \param inputs  Inputs to the function
    //!
    virtual void ExecuteCallBack(FunctionInfo*                finfo,
								 const std::string&           fname,
                                 const std::vector<Currency>& inputs) {}
    //!
    //! Set the current interpreter
    //! \param interp handle to current interpreter
    //!
    virtual void SetCurrentInterp(Interpreter* interp) {}
    //!
    //! Returns the current interpreter
    //!
    virtual Interpreter* GetCurrentInterp() const { return NULL; }

    //!
    //! Cancel printing if this variable is being printed
    //! \param cur Variable name
    //!
    virtual void CancelPrinting(const std::string& name) {}

    //!
    //! Helper method which returns true if application is in Batch mode
    //!
    virtual bool IsInBatchMode() const { return false; }
    //!
    //! Helper method which returns true if application is in Console mode
    //!
    virtual bool IsInConsoleMode() const { return false; }
    //!
    //! Helper method which returns true if application is in GUI mode
    //!
    virtual bool IsInGuiMode() const { return false; }
    //!
    //! Returns true if in Console, non-interactive mode
    //!
    virtual bool IsInConsoleBatchMode() const { return false; }
    //!
    //! Shows status
    //! \param Message to display in status bar, if applicable
    //!
    virtual void ShowStatus(const std::string& msg) {}
    //!
    //! Hides status bar, if applicable
    //!
    virtual void HideStatus() {}
protected:
    //!
    //! Constructor
    //!
    SignalHandlerBase() {} 
private:

};

#endif
// End of file:

