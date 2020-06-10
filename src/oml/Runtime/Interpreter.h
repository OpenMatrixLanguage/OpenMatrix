/**
* @file Interpreter.h
* @date February 2015
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

#ifndef __Interpreter_h
#define __Interpreter_h

// Begin defines/includes
#include "Hml2Dll.h"

#include "Currency.h"
#include "EvaluatorDebug.h"
#include "EvaluatorInt.h"

class ExprTreeEvaluator;
class InterpreterImpl;
class OutputFormat;
class SignalHandlerBase;
class FunctionInfo;
// End defines/includes

//------------------------------------------------------------------------------
//! oml interpreter
//------------------------------------------------------------------------------
class HML2DLL_DECLS Interpreter
{
public:
    Interpreter(ExprTreeEvaluator* source);        //! Constructor
    Interpreter(EvaluatorInterface& source);       //! Constructor
    Interpreter();                                  //! Default constructor
    virtual ~Interpreter();                         //! Destructor

    //! Adds toolbox
    void AddToolbox(const std::string& filepath); 

    //! Clear
    void Clear(const std::string& varname);     
    //! Resets function search cache
    void ResetFuncSearchCache();

    //! Gets result after processing file
	Currency DoFile(const std::string& filename);   

    //! Gets result after processing string  
	Currency DoString(const std::string& instring, bool store_suppressed=true);    

    //! Gets result and clears, if needed
	std::vector<Currency> GetOutputCurrencyList(bool clear = false); 

    //! Gets value
	const Currency& GetValue(const std::string& varname, 
                              int                offset = 0) const;  
    //! Sets value
	bool SetValue(const std::string& varname,            
                   Currency           value); 

    //! Gets a global value
	const Currency&  GetGlobalValue(const std::string& varname);

    //! Sets a global value
	bool SetGlobalValue(const std::string& varname, Currency value);  

    //! Gets function names
    std::vector<std::string> GetFunctionNames() const;               
    //! Gets keywords
	std::vector<std::string> GetKeywords() const;            
    //! Gets variable names
    std::vector<std::string> GetVariableNames(int offset = 0) const; 

    //! Gets call stack 
    std::vector<DebugStateInfo> GetCallStack() const;  
    //! Gets stack depth
	int GetStackDepth() const;        
    //! Gets stack info
    void GetStackInfo(int          level,     
                       std::string& filename, 
                       int&         linenumber);

    //! Returns true if name exists
	bool Exists(const std::string& varname) const; 
    //! Returns true if name is global
	bool IsGlobalVariable(const std::string& varname) const; 
    //! True if expression is partial
	bool IsPartialExpression(const std::string& expr) const;   
    //! True if interrupt has been requested by evaluator
    bool IsInterrupt() const;

    //!
    //! True if execution has been paused by evaluator
    //!
    bool IsPaused() const;
    //!
    //! True if a pause request is pending in the evaluator
    //!
    bool IsPauseRequestPending() const;

          
    //! Returns true if successful in registering a built in function
    bool RegisterBuiltInFunction(const std::string& funcname, 
                                  FUNCPTR            fp, 
                                  int                nargin, 
                                  int                nargout);
	//! Returns true if successful in registering a built in function
    bool RegisterBuiltInFunction(const std::string& funcname, 
                                  FUNCPTR            fp, 
								  FunctionMetaData   fmd);
    //! Registers debug interface
	void RegisterDebugListener(EvaluatorDebugInterface* edi);

	void RegisterExternalVariableFunction(EXTPTR ep);
    //!
    //! Returns true if successful in locking a built-in function
    //! \param funcname Built in function to lock
    //!
    bool LockBuiltInFunction(const std::string& funcname);

    //! Rename
    int RenameVariable(const std::string& oldname,     
                        const std::string& newname);

    //! Triggers interrupt  
    void TriggerInterrupt();    

    //! Triggers pause  
    void TriggerPause(bool state);   

    //! Function callback
	Currency CallFunction(const std::string&           funcname, 
                           const std::vector<Currency>& inputs);
    //! Function callback
	Currency CallFunction(FunctionInfo*                  finfo, 
                           const std::vector<Currency>& inputs);

	//! Function callback
	Currency CallFunctionInCurrentScope(FunctionInfo* finfo,
		const std::vector<Currency>& inputs);
	
    //! Function callback
	void CallFunction(const std::string&           funcname, 
                       const std::vector<Currency>& inputs, 
                       std::vector<Currency>&       outputs, 
                       int                          nargout = 1);

    const OutputFormat* GetOutputFormat() const;

    //! Returns true if quitting
    bool IsQuit() const;
    //! Sets the quit flag
    //! \param[in] val True if the evaluator is quitting
    void SetQuit( bool val);

    //! Returns true if evaluator is in experimental mode
    bool GetExperimental() const;
    //! Sets the experimental mode flag (-x)
    //! \param[in] val Sets to true if the evaluator is in experimental mode
    void SetExperimental(bool val);

	//! Search User Docs
	std::string GetHelpModule(std::string funcName);
	std::string GetHelpDirectory(std::string funcName);
	std::string GetHelpFilename(std::string funcName);

    //! Gets signal handler
    SignalHandlerBase* GetSignalHandler() const;
    //! Sets the signal handler
    //! \param[in] signalHandler Signal handler instance
    void SetSignalHandler(SignalHandlerBase* signalHandler);
    //! Gets list of functions for a file
    //! \param[in] fileName a valid file path
	std::vector<std::string> FunctionListInFile(const std::string& fileName);
    //! Gets no of args for a function in a file
    //! \param[in] functionName name of the function
    //! \param[in] fileName a valid file path
	int ArgumentsForFunction(const std::string& functionName, const std::string& fileName);

    FunctionInfo* FunctionInfoFromName(const std::string& funcname);

    //! Gets evaluator
    ExprTreeEvaluator* GetEvaluator() const;            

    //! Gets the application directory
    std::string GetApplicationDir() const;
    //! Sets the application directory
    void SetApplicationDir( const std::string& dir);

	int         GetSyntaxErrorLine() const;
	std::string GetSyntaxErrorFile() const;

	void* GetUserData()            { return _user_data; }
	void  SetUserData(void* _data) { _user_data = _data; }

private:
    friend class InterpreterImpl; //! Access private/protected members
    InterpreterImpl* _impl;       //! Internal implementation

    //! Copy constructor - disallow
	Interpreter(const Interpreter&);  
    //! Assignment op - disallow
	Interpreter& operator=(const Interpreter&);

	void* _user_data;
};
#endif
// End of file:

