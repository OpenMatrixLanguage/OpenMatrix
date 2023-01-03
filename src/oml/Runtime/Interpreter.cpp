/**
* @file Interpreter.cpp
* @date February 2015
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
#include "ANTLRoverride.h"    // keep this first

#include "Interpreter.h"
#include "SignalHandlerBase.h"

#include <algorithm>
#include <cassert>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>

#include "ErrorInfo.h"
#include "OML_Error.h"
#include "Evaluator.h"
#include "FunctionInfo.h"
#include "ANTLRData.h"
#include "BuiltInFuncs.h"
#include "BuiltInFuncsUtils.h"
#include "OutputFormat.h"
#include "OMLTree.h"
#include "hwMathException.h"

#ifdef OS_WIN
#include <Windows.h>
#include <WinBase.h>
#include <Shlwapi.h>
#define S_ISDIR(mode) (mode & S_IFDIR)
#endif
std::map< void *, ErrorInfo > errMap;
// End defines/includes

#if 0
#    define DEBUG_PRINT(s) { std::cout << s << std::endl; }
#else
#    define DEBUG_PRINT(s)
#endif

//------------------------------------------------------------------------------
//! Interpreter implementation
//------------------------------------------------------------------------------
class InterpreterImpl
{
public:
    //! Constructor
    InterpreterImpl(Interpreter*   intp,
                    ExprTreeEvaluator* src) 
                    : _interp (intp), _eval (src) { assert(_interp); }
    //! Constructor
    InterpreterImpl(Interpreter* intp) : _interp (intp) { assert(_interp); }
    //! Destructor
    ~InterpreterImpl()
    {
        SignalHandlerBase* s = _eval.GetSignalHandler();
        if (s && s->IsClone())
        {
            s->DeleteClone();
            delete s;
            _eval.SetSignalHandler(NULL);
        }
    }

    //! Gets evaluator
    ExprTreeEvaluator* GetEvaluator() const { return const_cast<ExprTreeEvaluator*>(&_eval); }

    //! Adds tool box
    void AddToolbox(const std::string& filepath) { DoFile(filepath); }  

	
    //! Clears name
    void Clear(const std::string& name) { _eval.Clear(name); }     
    //! Renames variable
    int RenameVariable(const std::string& oldname, 
                        const std::string& newname);
    //! Trigger input
    void TriggerInterrupt() { _eval.SetInterrupt(true); }

    //! Trigger input
    void TriggerPause(bool state) { _eval.SetPauseRequestPending(state); }

    //! Reads file and returns results
	Currency DoFile(const std::string& filename);

    //! Reads string and returns results
	Currency DoString(const std::string& input, bool store_suppressed);

    //! Gets output and clears buffer if needed
    std::vector<Currency> GetOutputCurrencyList(bool clear);

    //! Gets function names
    std::vector<std::string> GetFunctionNames() const { return _eval.GetFunctionNames(); }
    //! Gets keywords
    std::vector<std::string> GetKeywords() const { return _eval.GetKeywords(); }
    //! Gets variable names
    std::vector<std::string> GetVariableNames(int offset) const { return _eval.GetVariableNames(offset); }

    //! Gets value
	const Currency& GetValue(const std::string& varname, 
                              int                offset) const { return _eval.GetValue(varname, offset); }
    //! Returns true if successful in setting value
	bool SetValue(const std::string& varname, 
                   Currency           value) { return _eval.SetValue(varname, value); }

    //! Gets global value
    const Currency& GetGlobalValue(const std::string& varname) const { return _eval.GetGlobalValue(varname); }

    //! Returns true if successful in setting value
	bool SetGlobalValue(const std::string& varname, Currency value) { return _eval.SetGlobalValue(varname, value); }

    // Gets call stack
    std::vector<DebugStateInfo> GetCallStack() const { return _eval.GetCallStack(); }
    //! Get stack depth
	int GetStackDepth() const { return _eval.GetStackDepth(); }
    //! Gets stack info
    void GetStackInfo(int          level, 
                       std::string& filename, 
                       int&         linenumber) { return _eval.GetStackInfo(level, filename, linenumber); }

    //! Returns true if name exists
    bool Exists(const std::string& name) const { return _eval.Contains(name); }
    //! Returns true if name is global
    bool IsGlobalVariable(const std::string& name) const { return _eval.IsGlobal(name); }
    //! True if expression is partial
	bool IsPartialExpression(const std::string& expr);

    //! Register built in functions
	bool RegisterBuiltInFunction(const std::string& funcname, 
                                  FUNCPTR            fp, 
                                  int                nargin, 
                                  int                nargout) { return _eval.RegisterBuiltInFunction(funcname, fp, nargin, nargout); }

	    //! Register built in functions
	bool RegisterBuiltInFunction(const std::string& funcname, 
                                  FUNCPTR            fp, 
                                  FunctionMetaData   fmd) { return _eval.RegisterBuiltInFunction(funcname, fp, fmd); }

    //! Register debug interface
	void RegisterDebugListener(EvaluatorDebugInterface* edi) { _eval.RegisterDebugListener(edi); }
		
	void RegisterExternalVariableFunction(EXTPTR ep) {_eval.RegisterExternalVariableFunction(ep); }

    // Function callback
	Currency CallFunction(const std::string&, const std::vector<Currency>&, bool = true);

	//! Function callback
	Currency CallFunctionInCurrentScope(FunctionInfo* finfo,
		const std::vector<Currency>& inputs);

    Currency CallFunction(FunctionInfo*                  finfo, 
                           const std::vector<Currency>& inputs);
    //! Function callback
	void CallFunction(const std::string&           funcname, 
                       const std::vector<Currency>& inputs, 
                       std::vector<Currency>&       outputs, 
                       int                          nargout) { _eval.CallFunction(funcname, inputs, outputs, nargout); }

    const OutputFormat* GetOutputFormat() const { return _eval.GetOutputFormat(); }

    //! Returns true if quitting
    bool IsQuit() const { return _eval.IsQuit(); }
    //! Sets the quit flag
    //! \param[in] val True if the evaluator is quitting
    void SetQuit( bool val) { _eval.SetQuit(val); }

    //! Returns true if evaluator is in experimental mode
    bool GetExperimental() const { return _eval.GetExperimental(); }
    //! Sets the experimental mode flag (-ex)
    //! \param[in] val Sets to true if the evaluator is in experimental mode
    void SetExperimental(bool val) { _eval.SetExperimental(val); }

	std::string GetHelpModule(const std::string&, const bool&);
	std::string GetHelpDirectory(const std::string& funcName);
	std::string GetHelpFilename(const std::string&);

	int ArgumentsForFunction(const std::string& functionName);

    FunctionInfo* FunctionInfoFromName(const std::string& funcname);

	void ClearUserFunctions();

    bool LockBuiltInFunction(const std::string& name) { return _eval.LockBuiltInFunction(name); }

    void LogInput(const std::string&); // Logs input  in diary
    void LogOutput(const Currency&);   // Logs output in diary

	bool IsInPaths(const std::string& path);

	Currency GetProfileData() const;
	void     ClearProfileData();
	void     Profile(bool on);

	int         syntax_error_line;
	std::string syntax_error_file;

private:
	ExprTreeEvaluator _eval;    //! Evaluator
    Interpreter*      _interp;  //! Interpreter
	std::string       AST;

    InterpreterImpl() : _interp(0) {} // Private constructor
	//! Stubbed out copy constructor
	InterpreterImpl(const InterpreterImpl&);
    //! Stubbed out assignment operator
	InterpreterImpl& operator=(const InterpreterImpl&);

    //! Evaluates input and returns currency
	Currency CommonEvaluate(pANTLR3_INPUT_STREAM input);
    //! True if this is a valid name
	bool IsValidVarname(const std::string& name) const;

    //! Slot for handling push results
    void OnPushResultsHandler();
};
//------------------------------------------------------------------------------
//! Gets output as a list and clears results buffer if needed
std::vector<Currency> InterpreterImpl::GetOutputCurrencyList(bool clear)
{
	std::vector<Currency> results (_eval.GetOutputResults());

	if (clear)
    {
		_eval.ClearOutputResults();
        _eval.OnClearResults();  // Clear display results on the client 
    }

	return results;
}	
//------------------------------------------------------------------------------
//! Renames variable
int InterpreterImpl::RenameVariable(const std::string& oldname, const std::string& newname)
{
	if (!IsValidVarname(newname))
		return -2;

	return _eval.RenameVariable(oldname, newname);
}
//------------------------------------------------------------------------------
//! Reads file and returns results
Currency InterpreterImpl::DoFile(const std::string& filename)
{
	_eval.SetDebugInfo(nullptr, 0);
	size_t      dot_index = filename.rfind('.');
	std::string extension = filename.substr(dot_index + 1);

	if (_eval.IsExtensionEncrypted(extension))
	{
		try
		{
			_eval.Mark();
			_eval.RunEncryptedFile(extension, filename);
			_eval.Unmark();
		}
		catch (const OML_Error& e)
		{
			std::string error_str = e.GetFormatMessage() ?
				_eval.FormatErrorMessage(e.GetErrorMessage()) : e.GetErrorMessage();
			_eval.SetLastErrorMessage(e.GetErrorMessage());
			Currency cur(-1.0, error_str);
			_eval.PushResult(cur);
			_eval.Restore();
			return cur;
		}
		return _eval.GetLastResult();
	}

	pANTLR3_INPUT_STREAM input = ANTLRData::InputFromFilename(filename);

	if (!input)
	{
        std::string err ("Error: Input file not found: " + filename);
        Currency cur(-1.0, err);
        _eval.SetLastErrorMessage(err);
        _eval.PushResult(cur);
		return cur;
	}

	_eval.SetInterrupt(false);
	_eval.Mark();
	_eval.ClearFunctionsFromFile(filename);

	try 
	{
		CommonEvaluate(input);
        _eval.Unmark();
		return _eval.GetLastResult();
	}
	catch (OML_Error& error)
	{
		std::string error_str = error.GetFormatMessage() ? 
            _eval.FormatErrorMessage(error.GetErrorMessage()) : error.GetErrorMessage();
        _eval.SetLastErrorMessage(error.GetErrorMessage());
		Currency cur(-1.0, error_str);

		_eval.PushResult(cur);
		_eval.Restore();

		return cur;
	}
}

//------------------------------------------------------------------------------
//! Reads string and returns results
Currency InterpreterImpl::DoString(const std::string& instring, bool store_suppressed)
{
	pANTLR3_INPUT_STREAM input = ANTLRData::InputFromExpression(instring, "dummy");

	_eval.SetInterrupt(false);
    _eval.Mark();   

	std::string dbg_file = _eval.GetCurrentDebugFile();
	int         dbg_line = _eval.GetCurrentDebugLine();

	const std::string* dbg_file_ptr = Currency::pm.GetStringPointer(dbg_file);

	std::string instringcopy (instring);
#if 0
    // Commenting this for now, should not be needed
	instringcopy.erase(std::remove(instringcopy.begin(), instringcopy.end(), '\n'), instringcopy.end());
#endif

    LogInput(instring);
	try
	{
		_eval.StoreSuppressedResults(store_suppressed);
		Currency temp = CommonEvaluate(input);

        _eval.Unmark();                         
		_eval.StoreSuppressedResults(false);
		_eval.SetDebugInfo(dbg_file_ptr, dbg_line);
		_eval.SetInterrupt(false);

		if (temp.IsReturn())
			return temp;

		Currency result = _eval.GetLastResult();
        LogOutput(result);
		return result;
	}
	catch (OML_Error& error)
	{
		std::string error_str = error.GetFormatMessage() ? 
            _eval.FormatErrorMessage(error.GetErrorMessage()) : error.GetErrorMessage();
		Currency cur(Currency::TYPE_ERROR, error_str);
		_eval.SetLastErrorMessage(error.GetErrorMessage());
		_eval.PushResult(cur);
        _eval.Restore();       
		_eval.SetInterrupt(false);
        LogOutput(cur);
		return cur;
	}

	
}
//------------------------------------------------------------------------------
//! Evaluates input and returns result
Currency InterpreterImpl::CommonEvaluate(pANTLR3_INPUT_STREAM input)
{
#if OS_WIN
#	if _MSC_VER < 1900
	_set_output_format(_TWO_DIGIT_EXPONENT);
#	endif
#endif

	static bool bread_crumb = false;

    // Clear all results first
   _eval.ClearOutputResults();      // clear any outputs this may have generated
   _eval.OnClearResults();     // Emit signals for any GUI/Batch to handle

	ANTLRData ad(input, true);
	pANTLR3_COMMON_TOKEN_STREAM tokens = ad.GetTokens();
                                         	
//	ad.DumpTokenStream();
	ANTLRData::PreprocessTokenStream(tokens);
//	ad.DumpTokenStream();

	ad.CreateParser(tokens);
	pExprCppTreeParser parser = ad.GetParser();

	syntax_error_line = 0;
	syntax_error_file = "";

	ExprCppTreeParser_prog_return r = parser->prog(parser);
	
	if (parser->pParser->rec->getNumberOfSyntaxErrors(parser->pParser->rec) > 0)
	{		
		std::map<void*, ErrorInfo>::iterator iter = errMap.find(parser->pParser);

		if (iter == errMap.end())
        {
			_eval.PushResult(Currency(-999.0, "Parse Error"));
            _eval.SetLastErrorMessage("Parse Error");
        }
		else
        {
			ErrorInfo ei = iter->second;

			syntax_error_line = ei.GetLine();
			syntax_error_file = ei.GetFilename();

			_eval.PushResult(Currency(-999.0, ei.GetErrStr()));
            _eval.SetLastErrorMessage(ei.GetErrStr());
        }
		return _eval.GetLastResult();
	}

	pANTLR3_BASE_TREE tree = r.tree;

    _eval.PreprocessAST(tree, tokens);
	OMLTree* oml_tree = OMLTree::ConvertTree(tree);

	std::string error;
	Currency    rr;
    bool        formaterr = true;
	try
	{
		if (input && input->fileName)
			_eval.SetScriptName((char*)input->fileName->chars);

		rr = _eval.RunTree(oml_tree);

	}
	catch (OML_Error& e)
	{
		error = e.GetErrorMessage();
        formaterr = e.GetFormatMessage();
		_eval.ErrorCleanup();
	}
    catch (hwMathException& e)
    {
        //throw OML_Error(e.Status().GetMessageString());
        throw OML_Error(e.what());
    }	
    catch (std::bad_alloc& e)
	{
		OML_Error err(HW_MATH_ERR_ALLOCFAILED);
		error = err.GetErrorMessage();
		_eval.ErrorCleanup();
	}
	if (error.size())
		throw OML_Error(error, formaterr);

	_eval.ClearTemporaryExternalVariables();

	delete oml_tree;
	return rr;
}
//------------------------------------------------------------------------------
//! True if this is a valid name
bool InterpreterImpl::IsValidVarname(const std::string& name) const
{
	pANTLR3_INPUT_STREAM input = ANTLRData::InputFromExpression(name, "dummy");
	ANTLRData ad(input, true);
	pANTLR3_COMMON_TOKEN_STREAM tokens = ad.GetTokens();

	pANTLR3_VECTOR       vec      = tokens->getTokens(tokens);
	pANTLR3_COMMON_TOKEN tok      = (pANTLR3_COMMON_TOKEN)vec->get(vec, 0);
	char*                tok_text = (char*)tok->getText(tok)->chars;

	bool                 is_valid = (name == tok_text);

	return is_valid;
}
//------------------------------------------------------------------------------
// Calls the function of the given name
//------------------------------------------------------------------------------
Currency InterpreterImpl::CallFunction(const std::string&           funcname, 
	                                   const std::vector<Currency>& inputs,
	                                   bool                         printError)
{
    GetOutputCurrencyList(true); // clear any outputs this may have generated

    _eval.Mark();

	std::string dbg_file = (printError) ? _eval.GetCurrentDebugFile() : "";
	int         dbg_line = (printError) ? _eval.GetCurrentDebugLine() : 0;
	 
	const std::string* dbg_file_ptr = Currency::pm.GetStringPointer(dbg_file);

    
    // Need a try-catch block to close scopes properly in case of any error
    try
    {
        Currency ret = _eval.CallFunction(funcname, inputs);
        GetOutputCurrencyList(true); // clear any outputs this may have generated
        return ret;
    }
	catch (OML_Error& error)
	{
		std::string error_str = (error.GetFormatMessage() && printError) ? 
            _eval.FormatErrorMessage(error.GetErrorMessage()) : error.GetErrorMessage();
        _eval.SetLastErrorMessage(error.GetErrorMessage());
		_eval.SetDebugInfo(dbg_file_ptr, dbg_line);

		Currency cur(-1.0, error_str);

		if (printError)
		{
			_eval.PushResult(cur);
		}
		_eval.Restore();
		return cur;
	}

    return Currency();  // keep compiler happy
}
//------------------------------------------------------------------------------
//! Function callback
Currency InterpreterImpl::CallFunction(FunctionInfo* finfo, const std::vector<Currency>& inputs)
{
 //   return _eval.CallInternalFunction(finfo, inputs);
    GetOutputCurrencyList(true); // clear any outputs this may have generated

    _eval.Mark();
    
    // Need a try-catch block to close scopes properly in case of any error
    try
    {
        Currency ret = _eval.CallInternalFunction(finfo, inputs);
        GetOutputCurrencyList(true); // clear any outputs this may have generated
        return ret;
    }
	catch (OML_Error& error)
	{
		std::string error_str = error.GetFormatMessage() ? 
            _eval.FormatErrorMessage(error.GetErrorMessage()) : error.GetErrorMessage();
        _eval.SetLastErrorMessage(error.GetErrorMessage());
		Currency cur(-1.0, error_str);

		_eval.PushResult(cur);
		_eval.Restore();

		return cur;
	}

    return Currency(); // keep compiler happy
}

//------------------------------------------------------------------------------
//! Function callback
Currency InterpreterImpl::CallFunctionInCurrentScope(FunctionInfo* finfo, const std::vector<Currency>& inputs)
{
	//   return _eval.CallInternalFunction(finfo, inputs);
	GetOutputCurrencyList(true); // clear any outputs this may have generated

	_eval.Mark();

	MemoryScopeManager* msm       = _eval.GetMSM();
	MemoryScope*        cur_scope = msm->GetCurrentScope();
    MemoryScope         old_scope;
    if (cur_scope)
    {
        old_scope = *cur_scope; // Done so any local variables declared in callbacks don't get propogated
    }
	FunctionInfo*       cur_fi    = cur_scope->GetFunctionInfo();

	// Need a try-catch block to close scopes properly in case of any error
	try
	{
		if (!cur_fi)
			cur_scope->SetFunctionInfo(finfo);

		Currency ret = _eval.CallInternalFunction(finfo, inputs, false, "", Currency());
		GetOutputCurrencyList(true); // clear any outputs this may have generated

		if (!cur_fi)
			cur_scope->SetFunctionInfo(NULL);

        if (cur_scope)
        {
            *cur_scope = old_scope;
        }
		return ret;
	}
	catch (OML_Error & error)
	{
		std::string error_str = error.GetFormatMessage() ?
			_eval.FormatErrorMessage(error.GetErrorMessage()) : error.GetErrorMessage();
		_eval.SetLastErrorMessage(error.GetErrorMessage());
		Currency cur(-1.0, error_str);

		_eval.PushResult(cur);
		_eval.Restore();

		if (!cur_fi)
			cur_scope->SetFunctionInfo(NULL);

        if (cur_scope)
        {
            *cur_scope = old_scope;
        }

		return cur;
	}

	return Currency(); // keep compiler happy
}

//------------------------------------------------------------------------------
//! True if expression is partial
//------------------------------------------------------------------------------
bool InterpreterImpl::IsPartialExpression(const std::string& expr)
{
 	pANTLR3_INPUT_STREAM input = ANTLRData::InputFromExpression(expr, "dummy");
	ANTLRData ad(input, true);   

    pANTLR3_COMMON_TOKEN_STREAM tokens = ad.GetTokens();
               
    ANTLRData::PreprocessTokenStream(tokens);

	ad.CreateParser(tokens);
	pExprCppTreeParser parser = ad.GetParser();

    ExprCppTreeParser_prog_return r = parser->prog(parser);

    bool is_partial = false;

    // This isn't what I really wanted to do.  Unfortunately, since the parsing didn't 
    // complete properly, there is no tree - so I can't use it to see what the current statement is.
    // Thus, I'm walking backwards through tokens to see if I hit a "magic" one.  It's not
    // bulletproof, but it should work in the majority of cases.
    if (parser->pParser->rec->getNumberOfSyntaxErrors(parser->pParser->rec) > 0)
    {
        pANTLR3_VECTOR       vec      = tokens->getTokens(tokens);
        pANTLR3_COMMON_TOKEN orig_tok = (pANTLR3_COMMON_TOKEN)parser->pParser->rec->state->exception->token;

        if ((parser->pParser->rec->state->exception->expecting == END) || (orig_tok->getType(orig_tok) == EOF))
        {
			int start_token_index = vec->size(vec);

			if (parser->pParser->rec->state->exception->index > 1)
				start_token_index = (int)parser->pParser->rec->state->exception->index;

            for (int j=start_token_index-1; j>=0; j--)
            {
                pANTLR3_COMMON_TOKEN tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, j);
                int                  tok_type = tok->getType(tok);

                // not sure if we need to deal with the line continuation (...) here - I don't see why we would
                if ((tok_type == FOR) || (tok_type == WHILE) || (tok_type == IF) || (tok_type == SWITCH) || (tok_type == TRY))
                {
                    is_partial = true;
                    break;
                }
            }
        }
        else if ((parser->pParser->rec->state->exception->expecting == RBRACKET)  || (parser->pParser->rec->state->exception->expecting == RCURLY))
        {
            // Unfortunately, this code never gets hit any more for RBRACKET.  Apparently the parser goes down the multi-return function path
            // instead of the matrix path when it gives up.
            is_partial = true;
        }
        else if (parser->pParser->rec->state->exception->expecting == 0) 
        {
			int orig_tok_type = orig_tok->getType(orig_tok);

            if ((orig_tok_type == FOR) || (orig_tok_type == FUNCTION) || (orig_tok_type == IF) || (orig_tok_type == SWITCH) || (orig_tok_type == WHILE))
            {
                is_partial = true;
            }
            else
            {
                int num_tokens = vec->size(vec);

				if (num_tokens >= 2)
				{
					pANTLR3_COMMON_TOKEN last_tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, num_tokens - 1);
					int                  last_tok_type = last_tok->getType(last_tok);
					pANTLR3_COMMON_TOKEN last_tok_2 = (pANTLR3_COMMON_TOKEN)vec->get(vec, num_tokens - 2);
					int                  last_tok_type_2 = last_tok->getType(last_tok_2);

					if ((last_tok_type == NEWLINE) && (last_tok_type_2 == CONT))
					{
						is_partial = true;
					}
				}

                // This code is unpretty, but it's the best I could come up with to detect
                // matrices broken across lines (i.e. the CR is the row end token)
                // I'm not sure if I have to loop backwards or not.  The loop code was there originally,
                // but for the latest bug I fixed, the exception->index token was always LBRACKET.
                for (int j=(int)parser->pParser->rec->state->exception->index; j>=0; j--)
                {
                    pANTLR3_COMMON_TOKEN tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, j);
                    int                  tok_type = tok->getType(tok);

                    if ((tok_type == LBRACKET) || (tok_type == LCURLY))
                    {
						int exit_token = RBRACKET;

						if (tok_type == LCURLY)
							exit_token = RCURLY;

						bool fake_partial = false;

						for (unsigned int k=(unsigned int)parser->pParser->rec->state->exception->index; k<vec->size(vec); k++)
						{
							pANTLR3_COMMON_TOKEN tok2      = (pANTLR3_COMMON_TOKEN)vec->get(vec, k);
							int                  tok_type2 = tok2->getType(tok2);

							if (tok_type2 == exit_token)
							{
								fake_partial = true;
								break;
							}
							else if ((tok_type2 == NEWLINE) || (tok_type2 == DUMMY))
							{
								break;
							}
						}

                        is_partial = !fake_partial;
                        break;
                    }
                    else if ((tok_type == RBRACKET) || (tok_type == RCURLY))
                    {
                        break;
                    }

                }
            }
        }
    }
    else
    { 
        pANTLR3_BASE_TREE    tree   = (pANTLR3_BASE_TREE)r.tree;
        pANTLR3_BASE_TREE    child  = (pANTLR3_BASE_TREE)tree->getChild(tree, 0); 

		if (child)
		{
			pANTLR3_COMMON_TOKEN tok    = child->getToken(child);

			if (tok->getType(tok) == FUNC_DEF)
			{
				int                  num_func_children = child->getChildCount(child);
				pANTLR3_BASE_TREE    end_child         = (pANTLR3_BASE_TREE)child->getChild(child, num_func_children-1);
				pANTLR3_COMMON_TOKEN end_tok           = end_child->getToken(end_child);

				if (end_tok->getType(end_tok) == EOF)
					is_partial = true;
				else
					is_partial = false;
			}
			else if (tok->getType(tok) == PARSE_DUMMY)
			{
				is_partial = true;
			}
		}
    }

    return is_partial;
}

int InterpreterImpl::ArgumentsForFunction(const std::string& functionName)
{
	return _eval.NarginFor(functionName);
}
//------------------------------------------------------------------------------
FunctionInfo* InterpreterImpl::FunctionInfoFromName(const std::string& funcName)
{
    FunctionInfo *fi   = nullptr;
    FUNCPTR       fptr = nullptr;
	ALT_FUNCPTR   aptr = nullptr; // ignore this for now

    _eval.FindFunctionByName(funcName, &fi, &fptr, &aptr);

	if (!fi && fptr)
		return new FunctionInfo(funcName, fptr);

	if (!fi && aptr)
		return new FunctionInfo(funcName, aptr);

    return fi;
}

void InterpreterImpl::ClearUserFunctions()
{
	_eval.ClearFunctions();
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//! Default constructor
Interpreter::Interpreter()
{
    _impl = new InterpreterImpl(this);
	_user_data = NULL;
}
//------------------------------------------------------------------------------
//! Constructor
//! \param[in] source Evaluator source
Interpreter::Interpreter(ExprTreeEvaluator* source)
{
    _impl = new InterpreterImpl(this, source);
	_user_data = NULL;
}
//------------------------------------------------------------------------------
//! Constructor
//! \param[in] source Evaluator source
Interpreter::Interpreter(EvaluatorInterface& source)
{
    _impl = new InterpreterImpl(this, source.eval);
	_user_data = NULL;
}
//------------------------------------------------------------------------------
//! Destructor
Interpreter::~Interpreter()
{
    delete _impl;
}
//------------------------------------------------------------------------------
//! Gets evaluator
//! \param[in] source Evaluator source
ExprTreeEvaluator* Interpreter::GetEvaluator() const
{
    return _impl->GetEvaluator();
}
//------------------------------------------------------------------------------
//! Adds tool box
void Interpreter::AddToolbox(const std::string& filepath)
{
	_impl->AddToolbox(filepath);
}
//------------------------------------------------------------------------------
//! Searches for User Docs
#ifdef OS_WIN
std::vector<std::string> files;
static void wildcardFileMatch_win(const char* target)
{
    WIN32_FIND_DATA fileData;

    HANDLE currentHandle = FindFirstFile(target, &fileData);
	char *fileName = fileData.cFileName;
	
    if (currentHandle != INVALID_HANDLE_VALUE)
	{
		do
		{
			if((fileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0 && strcmp(fileData.cFileName, ".") != 0 && strcmp(fileData.cFileName, "..") != 0)
			{
				std::string newStr = std::string(target);
				if (newStr.find("*.*") != std::string::npos)
					newStr.resize(newStr.size()-3);
				if (newStr.find("*..*") != std::string::npos)
					newStr.resize(newStr.size()-4);
				newStr = newStr + std::string(fileData.cFileName) + "\\*.html";
				wildcardFileMatch_win(newStr.c_str());
			}
			else
			{
				std::string newStr = std::string(target);
				newStr.resize(newStr.size()-6);
				 newStr = newStr + "\\" + fileData.cFileName;
				files.push_back(newStr);
			}
		}
		while (FindNextFile(currentHandle, &fileData));
	}
}
#endif
std::string Interpreter::GetHelpModule (std::string funcName, bool loadLibraryBrowser)
{
	return _impl->GetHelpModule(funcName, loadLibraryBrowser);
}

std::string Interpreter::GetHelpDirectory(std::string funcName)
{
	return _impl->GetHelpDirectory(funcName);
}

std::string Interpreter::GetHelpFilename (std::string funcName){
	return _impl->GetHelpFilename(funcName);
}
//----------------------------------------------------------------------------------------------
//! Processes a file and returns currency data
Currency Interpreter::DoFile(const std::string& filename)
{
	return _impl->DoFile(filename);
}
//---------------------------------------------------------------------------------------------
//! Processes a string and returns currency data
Currency Interpreter::DoString(const std::string& input, bool store_suppressed)
{
	return _impl->DoString(input, store_suppressed);
}
//------------------------------------------------------------------------------
//! Returns true if given name exists
//! \param[in] varname Given name
bool Interpreter::Exists(const std::string& varname) const
{
	return _impl->Exists(varname);
}
//------------------------------------------------------------------------------
//! Trigger input
void Interpreter::TriggerInterrupt()
{
	_impl->TriggerInterrupt();
}
//------------------------------------------------------------------------------
//! Triggers pause
void Interpreter::TriggerPause(bool state) 
{ 
    _impl->TriggerPause(state); 
}
//------------------------------------------------------------------------------
//! Gets keywords
std::vector<std::string> Interpreter::GetKeywords() const
{
	return _impl->GetKeywords();
}
//------------------------------------------------------------------------------
//! Gets function names
std::vector<std::string> Interpreter::GetFunctionNames() const
{
	return _impl->GetFunctionNames();
}
//------------------------------------------------------------------------------
//! Gets variable names
std::vector<std::string> Interpreter::GetVariableNames(int offset) const
{
	return _impl->GetVariableNames(offset);
}
//------------------------------------------------------------------------------
//! Gets value
const Currency& Interpreter::GetValue(const std::string& varname, int offset) const
{
	return _impl->GetValue(varname, offset);
}
//------------------------------------------------------------------------------
//! Returns true if successful in setting value
bool Interpreter::SetValue(const std::string& varname, Currency value)
{
	return _impl->SetValue(varname, value);
}
//------------------------------------------------------------------------------
//! Gets a global value
const Currency&  Interpreter::GetGlobalValue(const std::string& varname)
{
    return _impl->GetGlobalValue(varname);
}
//------------------------------------------------------------------------------
//! Returns true if successful in setting global value
bool Interpreter::SetGlobalValue(const std::string& varname, Currency value)
{
	return _impl->SetGlobalValue(varname, value);
}
//------------------------------------------------------------------------------
//! Gets output as a list
std::vector<Currency> Interpreter::GetOutputCurrencyList(bool clear)
{
	return _impl->GetOutputCurrencyList(clear);
}
//------------------------------------------------------------------------------
//! Gets stack info
void Interpreter::GetStackInfo(int level, std::string& filename, int& linenumber)
{
	return _impl->GetStackInfo(level, filename, linenumber);
}

//------------------------------------------------------------------------------
//! Clears varname
void Interpreter::Clear(const std::string& varname)
{
	return _impl->Clear(varname);
}
//------------------------------------------------------------------------------
//! Renames variable
int Interpreter::RenameVariable(const std::string& oldname, 
                                    const std::string& newname)
{
	return _impl->RenameVariable(oldname, newname);
}
//------------------------------------------------------------------------------
//! True if expression is partial
bool Interpreter::IsPartialExpression(const std::string& expr) const
{
	return _impl->IsPartialExpression(expr);
}
//------------------------------------------------------------------------------
//! Returns true if name is global
bool Interpreter::IsGlobalVariable(const std::string& varname) const
{
	return _impl->IsGlobalVariable(varname);
}
//------------------------------------------------------------------------------
//! Returns true if successful in registering a built in function
bool Interpreter::RegisterBuiltInFunction(const std::string& funcname, 
                                              FUNCPTR            fp, 
                                              int                nargin, 
                                              int                nargout)
{
	return _impl->RegisterBuiltInFunction(funcname, fp, nargin, nargout);
}
//! Returns true if successful in registering a built in function
bool Interpreter::RegisterBuiltInFunction(const std::string& funcname, 
                                              FUNCPTR            fp, 
											  FunctionMetaData   fmd)
{
	return _impl->RegisterBuiltInFunction(funcname, fp, fmd);
}
//------------------------------------------------------------------------------
//! Registers debug interface
void Interpreter::RegisterDebugListener(EvaluatorDebugInterface* edi)
{
	return _impl->RegisterDebugListener(edi);
}
//------------------------------------------------------------------------------
//! Gets stack depth
int Interpreter::GetStackDepth() const
{
	return _impl->GetStackDepth();
}
//------------------------------------------------------------------------------
//! Gets call stack 
std::vector<DebugStateInfo> Interpreter::GetCallStack() const 
{ 
	return _impl->GetCallStack(); 
}
//------------------------------------------------------------------------------
// Function callback
Currency Interpreter::CallFunction(const std::string&           funcname,  
                                   const std::vector<Currency>& inputs,
	                               bool                         printError)
{
	return _impl->CallFunction(funcname, inputs, printError);
}

//------------------------------------------------------------------------------
//! Function callback
Currency Interpreter::CallFunctionInCurrentScope(FunctionInfo* finfo,
	const std::vector<Currency>& inputs)
{
	return _impl->CallFunctionInCurrentScope(finfo, inputs);
}
//------------------------------------------------------------------------------
//! Function callback
Currency Interpreter::CallFunction(FunctionInfo*                  finfo,  
                                       const std::vector<Currency>& inputs)
{
	return _impl->CallFunction(finfo, inputs);
}
//------------------------------------------------------------------------------
//! Function callback
void Interpreter::CallFunction(const std::string&           funcname, 
                               const std::vector<Currency>& inputs, 
                               std::vector<Currency>&       outputs, 
                               int                          nargout)
{
	_impl->CallFunction(funcname, inputs, outputs, nargout);
}
//------------------------------------------------------------------------------
//! Gets output format
const OutputFormat* Interpreter::GetOutputFormat() const
{
    return _impl->GetOutputFormat();
}
//------------------------------------------------------------------------------
//! Returns true if quitting
//------------------------------------------------------------------------------
bool Interpreter::IsQuit() const
{
    return _impl->IsQuit();
}
//------------------------------------------------------------------------------
//! Returns true if quitting
//------------------------------------------------------------------------------
void Interpreter::RegisterExternalVariableFunction(EXTPTR ep)
{
    return _impl->RegisterExternalVariableFunction(ep);
}

std::string InterpreterImpl::GetHelpModule(const std::string& funcName, const bool& loadLibraryBrowser)
{
	if (_eval.IsKeyword(funcName))
		return "CoreMinimalInterpreter";
	else if (_eval.IsOperator(funcName))
		return "CoreMinimalInterpreter";
	else
		return _eval.GetHelpModule(funcName, loadLibraryBrowser);
}

std::string InterpreterImpl::GetHelpDirectory(const std::string& funcName)
{
	if (_eval.IsKeyword(funcName))
		return "CoreMinimalInterpreter";
	else if (_eval.IsOperator(funcName))
		return "CoreMinimalInterpreter";
	else
		return _eval.GetHelpDirectory(funcName);
}

std::string InterpreterImpl::GetHelpFilename(const std::string& funcName)
{
	std::string module = GetHelpDirectory(funcName);

	if (_eval.IsOperator(funcName))
	{
		if (funcName == "+")
			module += "/plus_sign.htm";
		else if (funcName == "-")
			module += "/minus_sign.htm";
		else if (funcName == "*")
			module += "/mtimes_sign.htm";
		else if (funcName == "/")
			module += "/mrdivide_sign.htm";
		else if (funcName == "\\")
			module += "/mldivide_sign.htm";
		else if (funcName == ".*")
			module += "/times_sign.htm";
		else if (funcName == "./")
			module += "/rdivide_sign.htm";
		else if (funcName == ".\\")
			module += "/ldivide_sign.htm";
		else if (funcName == "~")
			module += "/not_sign.htm";
		else if (funcName == "^")
			module += "/mpower_sign.htm";
		else if (funcName == ".^")
			module += "/power_sign.htm";
		else if (funcName == "==")
			module += "/eq_sign.htm";
		else if (funcName == "~=")
			module += "/ne_sign.htm";
		else if (funcName == "<")
			module += "/lt_sign.htm";
		else if (funcName == ">")
			module += "/gt_sign.htm";
		else if (funcName == "<=")
			module += "/le_sign.htm";
		else if (funcName == ">=")
			module += "/ge_sign.htm";
		else if (funcName == "&")
			module += "/and_sign.htm";
		else if (funcName == "&&")
			module += "/and_shortcircuit_sign.htm";
		else if (funcName == "|")
			module += "/or_sign.htm";
		else if (funcName == "||")
			module += "/or_shortcircuit_sign.htm";
		else if (funcName == ":")
			module += "/range.htm";
	}
	else
	{
		// module += "/";
		module += DIRECTORY_DELIM;
		module += funcName;
		module += ".htm";
	}
	// if the module include "/help/" the file must exist
	if (module.find(DIRECTORY_DELIM + std::string("help") + DIRECTORY_DELIM) != std::string::npos)
	{
		if (!BuiltInFuncsUtils::FileExists(module))
		{
			module = DIRECTORY_DELIM;
			module += funcName;
			module += ".htm";
		}
	}
	return module;
}
//------------------------------------------------------------------------------
//! True if evaluator is interuppted
//------------------------------------------------------------------------------
bool Interpreter::IsInterrupt() const
{
    ExprTreeEvaluator* eval = _impl->GetEvaluator();
    assert(eval);
    if (!eval) false;
    return eval->IsInterrupt();
}
//------------------------------------------------------------------------------
// True if execution has been paused by evaluator
//------------------------------------------------------------------------------
bool Interpreter::IsPaused() const
{
    ExprTreeEvaluator* eval = _impl->GetEvaluator();
    assert(eval);
    if (eval) 
    {
        return eval->IsPaused();
    }
    return false;
}
//------------------------------------------------------------------------------
//! //! Resets function search cache, called after cd in file browser
//------------------------------------------------------------------------------
void Interpreter::ResetFuncSearchCache()
{
    ExprTreeEvaluator* eval = _impl->GetEvaluator();
    assert(eval);
    if (eval)
        eval->ResetFuncSearchCache();
}
//------------------------------------------------------------------------------
//! Gets signal handler
//------------------------------------------------------------------------------
SignalHandlerBase* Interpreter::GetSignalHandler() const
{
    ExprTreeEvaluator* eval = _impl->GetEvaluator();
    assert(eval);
    if (!eval) return NULL;

    return eval->GetSignalHandler();
}
//------------------------------------------------------------------------------
//! Sets signal handler
//! \param[in] signalHandler Signal handler instance
//------------------------------------------------------------------------------
void Interpreter::SetSignalHandler(SignalHandlerBase* signalHandler)
{
    ExprTreeEvaluator* eval = _impl->GetEvaluator();
    assert(eval);

    if (eval)
        eval->SetSignalHandler(signalHandler);
}
//------------------------------------------------------------------------------
//! Returns true if evaluator is in experimental mode
//------------------------------------------------------------------------------
bool Interpreter::GetExperimental() const
{
    return Currency::GetExperimental(); 
}
//------------------------------------------------------------------------------
//! Sets the experimental mode flag (-x)
//! \param[in] val Sets to true if the evaluator is in experimental mode
//------------------------------------------------------------------------------
void Interpreter::SetExperimental(bool val)
{
     Currency::SetExperimental(val); 
}

//------------------------------------------------------------------------------
//! Returns true if evaluator is in verbose mode
//------------------------------------------------------------------------------

void CheckForFunctionsInTree(pANTLR3_BASE_TREE tree, std::vector<std::string>& results)
{
	int child_count = tree->getChildCount(tree);

	for (int j=0; j<child_count; j++)
	{
		pANTLR3_BASE_TREE child = (pANTLR3_BASE_TREE)tree->getChild(tree, j);

		pANTLR3_COMMON_TOKEN tok = child->getToken(child);
		int tok_type = tok->getType(tok);

		if (tok_type == FUNC_DEF)
		{
			pANTLR3_BASE_TREE func_child = (pANTLR3_BASE_TREE)child->getChild(child, 0);
			pANTLR3_STRING    str = func_child->getText(func_child);
			results.push_back((char*)str->chars);
		}
		else
		{
			CheckForFunctionsInTree(child, results);
		}
	}
}

std::vector<std::string> Interpreter::FunctionListInFile(const std::string& fileName)
{
	std::vector<std::string> results;

	pANTLR3_INPUT_STREAM input = ANTLRData::InputFromFilename(fileName);

	if (input)
	{
		ANTLRData ad(input, true);
		pANTLR3_COMMON_TOKEN_STREAM tokens = ad.GetTokens();
                                         	
		ANTLRData::PreprocessTokenStream(tokens);

		ad.CreateParser(tokens);
		pExprCppTreeParser parser = ad.GetParser();

		ExprCppTreeParser_prog_return r = parser->prog(parser);
	
		if (parser->pParser->rec->getNumberOfSyntaxErrors(parser->pParser->rec) == 0)
		{
			pANTLR3_BASE_TREE tree = r.tree;

			GetEvaluator()->PreprocessAST(tree, tokens);

			CheckForFunctionsInTree(tree, results);
		}
	}

	return results;
}
//------------------------------------------------------------------------------
//! Gets the application directory
//------------------------------------------------------------------------------
std::string Interpreter::GetApplicationDir() const 
{ 
    return GetEvaluator()->GetApplicationDir(); 
}
//------------------------------------------------------------------------------
//! Sets the application directory
//------------------------------------------------------------------------------
void Interpreter::SetApplicationDir(const std::string& dir)
{ 
    GetEvaluator()->SetApplicationDir(dir);
}
//------------------------------------------------------------------------------
//! Sets the quit flag
//------------------------------------------------------------------------------
void Interpreter::SetQuit(bool val)
{ 
    _impl->SetQuit(val);
}

int Interpreter::GetSyntaxErrorLine() const
{
	return _impl->syntax_error_line;
}

std::string Interpreter::GetSyntaxErrorFile() const
{
	return _impl->syntax_error_file;
}

void Interpreter::ClearUserFunctions()
{
	_impl->ClearUserFunctions();
}

//------------------------------------------------------------------------------
//! search for a given function and return no. of input args the fun takes
//------------------------------------------------------------------------------
int ArgumentsForFunctionInTree(pANTLR3_BASE_TREE tree, const std::string &funcname)
{
    int args = -1;

	int child_count = tree->getChildCount(tree);

	for (int j=0; j<child_count; j++)
	{
		pANTLR3_BASE_TREE child = (pANTLR3_BASE_TREE)tree->getChild(tree, j);

		pANTLR3_COMMON_TOKEN tok = child->getToken(child);
		int tok_type = tok->getType(tok);

		if (tok_type == FUNC_DEF)
		{
			pANTLR3_BASE_TREE func_child = (pANTLR3_BASE_TREE)child->getChild(child, 0);
			pANTLR3_STRING    str = func_child->getText(func_child);
            char* func = (char*)str->chars;
            if(func && (strcmp(funcname.c_str(), func) == 0))
            {
                int num_func_children = child->getChildCount(child);
                if (num_func_children >= 2) // func + args
	            {
                    pANTLR3_BASE_TREE func_args = (pANTLR3_BASE_TREE)child->getChild(child, 2);
                    return func_args->getChildCount(func_args);
                }
                return 0;
            }
		}
		else
		{
			args = ArgumentsForFunctionInTree(child, funcname);
            if(args != -1) break;
		}
	}

    return args;
}

//------------------------------------------------------------------------------
//! Gets no of args for a functions in a file
//------------------------------------------------------------------------------
int  Interpreter::ArgumentsForFunction(const std::string& functionName, const std::string& fileName)
{
    int args = -1;

    pANTLR3_INPUT_STREAM input = ANTLRData::InputFromFilename(fileName);

	if (input)
	{
		ANTLRData ad(input, true);
		pANTLR3_COMMON_TOKEN_STREAM tokens = ad.GetTokens();
                                         	
		ANTLRData::PreprocessTokenStream(tokens);

		ad.CreateParser(tokens);
		pExprCppTreeParser parser = ad.GetParser();

		ExprCppTreeParser_prog_return r = parser->prog(parser);
	
		if (parser->pParser->rec->getNumberOfSyntaxErrors(parser->pParser->rec) == 0)
		{
			pANTLR3_BASE_TREE tree = r.tree;

			GetEvaluator()->PreprocessAST(tree, tokens);

			args = ArgumentsForFunctionInTree(tree, functionName);
		}
	}

    return args;
}

//------------------------------------------------------------------------------
//! Gets function info from function name
//------------------------------------------------------------------------------
FunctionInfo* Interpreter::FunctionInfoFromName(const std::string& funcname)
{
    return _impl->FunctionInfoFromName(funcname);
}
//------------------------------------------------------------------------------
// True if a pause request is pending in the evaluator
//------------------------------------------------------------------------------
bool Interpreter::IsPauseRequestPending() const
{
    ExprTreeEvaluator* eval = _impl->GetEvaluator();
    assert(eval);
    if (eval) 
    {
        return eval->IsPauseRequestPending();
    }
    return false;
}
//------------------------------------------------------------------------------
// Returns true if successful in locking a built-in function
//------------------------------------------------------------------------------
bool Interpreter::LockBuiltInFunction(const std::string& name)
{
    return _impl->LockBuiltInFunction(name);
}
//------------------------------------------------------------------------------
// Returns true if path is in the list of paths
//------------------------------------------------------------------------------
bool Interpreter::IsInPaths(const std::string& path)
{
	return _impl->IsInPaths(path);
}
//------------------------------------------------------------------------------
// Returns profile data
//------------------------------------------------------------------------------
Currency Interpreter::GetProfileData() const
{
	return _impl->GetProfileData();
}
//------------------------------------------------------------------------------
// Clears profile data
//------------------------------------------------------------------------------
void Interpreter::ClearProfileData()
{
	_impl->ClearProfileData();
}
//------------------------------------------------------------------------------
// Turn on the profile mode
//------------------------------------------------------------------------------
void Interpreter::Profile(bool on)
{
	_impl->Profile(on);
}
//------------------------------------------------------------------------------
// Logs input in diary, if open
//------------------------------------------------------------------------------
void InterpreterImpl::LogInput(const std::string& in)
{
    std::ofstream& ofs = _eval.GetDiary();
    if (ofs.is_open())
    {
        std::string tmp(in);
        BuiltInFuncsUtils::StripTrailingNewline(tmp);
        if (!tmp.empty())
        {
            ofs << "> " << tmp << std::endl;
        }
    }
}
//------------------------------------------------------------------------------
// Logs output in diary, if open
//------------------------------------------------------------------------------
void InterpreterImpl::LogOutput(const Currency& cur)
{
    std::ofstream& ofs = _eval.GetDiary();
    if (ofs.is_open())
    {
        std::string out(cur.GetOutputString(_eval.GetOutputFormat()));
        if (!out.empty())
        {
            ofs << out << std::endl;
        }
    }
}
bool InterpreterImpl::IsInPaths(const std::string& path)
{
	return _eval.IsInPaths(path);
}

Currency InterpreterImpl::GetProfileData() const
{
	return _eval.GetProfileData();
}

void InterpreterImpl::ClearProfileData()
{
	_eval.ClearProfileData();
}

void InterpreterImpl::Profile(bool on)
{
	_eval.Profile(on);
}

// End of file: