/**
* @file EvaluatorInt.cpp
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

// Begin defines/includes
#include "EvaluatorInt.h"

#include <cassert>

#include "ANTLRData.h"
#include "BoundClassInfo.h"
#include "Evaluator.h"
#include "FunctionInfo.h"
#include "SignalHandlerBase.h"

// End defines/includes

EvaluatorInterface::EvaluatorInterface(const EvaluatorInterface& other, bool base_context, bool pop_nargs) : delete_eval(true)
{
    eval = other.eval->MakeContextCopy(base_context, pop_nargs);
}

EvaluatorInterface::~EvaluatorInterface()
{
    if (delete_eval)
        delete eval;
}

bool EvaluatorInterface::RegisterBuiltInFunction(const std::string& func_name, FUNCPTR fp, int nargin, int nargout)
{
    return eval->RegisterBuiltInFunction(func_name, fp, nargin, nargout);
}

bool EvaluatorInterface::RegisterBuiltInFunction(const std::string& func_name, FUNCPTR fp, const FunctionMetaData& md)
{
    return eval->RegisterBuiltInFunction(func_name, fp, md);
}

bool EvaluatorInterface::RegisterBuiltInFunction(const std::string& func_name, ALT_FUNCPTR fp)
{
    return eval->RegisterBuiltInFunction(func_name, fp);
}

std::vector<std::string> EvaluatorInterface::GetBuiltinFunctionNames() const
{
    return eval->GetBuiltinFunctionNames();
}

std::vector<std::string> EvaluatorInterface::GetKeywords() const
{
    return eval->GetKeywords();
}

const OutputFormat* EvaluatorInterface::GetOutputFormat() const
{
    return eval->GetOutputFormat();
}

const Currency& EvaluatorInterface::GetValue(std::string varname) const
{
    return eval->GetValue(varname);
}

const Currency& EvaluatorInterface::GetGlobalValue(std::string varname) const
{
    return eval->GetGlobalValue(varname);
}

bool EvaluatorInterface::FindFunctionByName(const std::string& func_name, FunctionInfo** fi, FUNCPTR* fptr)
{
    return eval->FindFunctionByName(func_name, fi, fptr);
}

bool EvaluatorInterface::SetValue(std::string varname, const Currency& value)
{
    return eval->SetValue(varname, value);
}

bool EvaluatorInterface::IsUserFunction(const std::string& func_name)
{
    return eval->IsUserFunction(func_name);
}

bool EvaluatorInterface::IsStdFunction(const std::string& func_name)
{
    return eval->IsStdFunction(func_name);
}

FUNCPTR EvaluatorInterface::GetStdFunction(const std::string& func_name) const
{
    return eval->GetStdFunction(func_name);
}

std::FILE* EvaluatorInterface::GetFile(int i)
{
    return eval->GetFile(i);
}

bool EvaluatorInterface::FindFileInPath(const std::string& file_plus_ext, std::string& filepath) const
{
    return eval->FindFileInPath(file_plus_ext, filepath);
}

bool EvaluatorInterface::IsGlobal(const std::string& varname) const
{
    return eval->IsGlobal(varname);
}

std::vector<std::string> EvaluatorInterface::GetVariableNames() const
{
    return eval->GetVariableNames(0);
}

Currency EvaluatorInterface::CallFunction(const std::string& func_name, const std::vector<Currency>& params)
{
    return eval->CallFunction(func_name, params);
}

std::vector<Currency> EvaluatorInterface::DoMultiReturnFunctionCall(FunctionInfo* fi, std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<std::string>* out_vars)
{
    return eval->DoMultiReturnFunctionCall(fi, param_values, num_ins, num_rets, suppress_output, out_vars);
}

std::vector<Currency> EvaluatorInterface::DoMultiReturnFunctionCall(FUNCPTR fptr, std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<std::string>* out_vars)
{
    return eval->DoMultiReturnFunctionCall(fptr, "", param_values, num_ins, num_rets, suppress_output, out_vars);
}

HML_CELLARRAY* EvaluatorInterface::CreateVararginCell(const std::vector<Currency>& params, int start_index)
{
    return eval->CreateVararginCell(params, start_index);
}

bool EvaluatorInterface::RemovePath(std::string &pathname)
{
    return eval->RemovePath(pathname);
}

void EvaluatorInterface::ClearPath()
{
    eval->ClearPath();
}

void EvaluatorInterface::PushResult(Currency res, bool to_output)
{
    eval->PushResult(res, to_output);
}

std::string EvaluatorInterface::GetFileMode (int i)
{
    return eval->GetFileMode(i);
}

bool EvaluatorInterface::CloseAllFiles()
{
    return eval->CloseAllFiles();
}

bool EvaluatorInterface::CloseFile (int i)
{
    return eval->CloseFile(i);
}

std::vector<int> EvaluatorInterface::GetFileIndices(int start)
{
    return eval->GetFileIndices(start);
}

int EvaluatorInterface::AddFile(std::FILE *newfile, const std::string &fname, const std::string &fmode)
{
    return eval->AddFile(newfile, fname, fmode);
}

std::string EvaluatorInterface::GetFileName (int i)
{
    return eval->GetFileName(i);
}

void EvaluatorInterface::PrintResult(Currency res)
{
    eval->PrintResult(res);
}

void EvaluatorInterface::Clear(const std::regex& varname)
{
    eval->Clear(varname);
}

void EvaluatorInterface::ClearFromGlobals(const std::regex& varname)
{
    eval->ClearFromGlobals(varname);
}

void EvaluatorInterface::ClearFromFunctions(const std::regex& varname)
{
    eval->ClearFromFunctions(varname);
}

void EvaluatorInterface::ClearFromVariables(const std::regex& varname)
{
    eval->ClearFromVariables(varname);
}

void EvaluatorInterface::Clear(const std::string& varname)
{
    eval->Clear(varname);
}

void EvaluatorInterface::ClearFromGlobals(const std::string& varname)
{
    eval->ClearFromGlobals(varname);
}

void EvaluatorInterface::ClearFromFunctions(const std::string& varname)
{
    eval->ClearFromFunctions(varname);
}

void EvaluatorInterface::ClearFromVariables(const std::string& varname)
{
    eval->ClearFromVariables(varname);
}

void EvaluatorInterface::ClearFunctions()
{
    eval->ClearFunctions();
}

void EvaluatorInterface::ClearVariables()
{
    eval->ClearVariables();
}

void EvaluatorInterface::ClearGlobals()
{
    eval->ClearGlobals();
}

Currency EvaluatorInterface::CallInternalFunction(FunctionInfo*fi, const std::vector<Currency>& param_values)
{
    return eval->CallInternalFunction(fi, param_values);
}

void EvaluatorInterface::ResetFuncSearchCache()
{
    eval->ResetFuncSearchCache();
}

bool EvaluatorInterface::Contains(const std::string& varname) const
{
    return eval->Contains(varname);
}

int EvaluatorInterface::GetContextEndValue() const
{
    return eval->GetContextEndValue();
}

int EvaluatorInterface::GetNargoutValue() const
{
    return eval->GetNargoutValue();
}

int EvaluatorInterface::GetNarginValue() const
{
    return eval->GetNarginValue();
}

int EvaluatorInterface::GetNumFiles()
{
    return eval->GetNumFiles();
}

void EvaluatorInterface::AddPath(std::string pathname, bool end)
{
    eval->AddPath(pathname, end);
}

void EvaluatorInterface::AddPath2(const std::string& pathname, const std::vector<std::string> funcs)
{
	eval->AddPath2(pathname, funcs);
}

const std::vector<std::string>& EvaluatorInterface::GetPaths() const
{
    return eval->GetPaths();
}

Currency EvaluatorInterface::EqualityOperator(const Currency& lhs, const Currency&rhs, int op)
{
    return eval->EqualityOperator(lhs, rhs, op);
}

Currency EvaluatorInterface::InequalityOperator(const Currency& lhs, const Currency&rhs, int op)
{
	Currency ret;

	switch (op)
	{
		case LTHAN:
			ret = eval->LessThanOperator(lhs, rhs);
		break;

		case GTHAN:
			ret = eval->GreaterThanOperator(lhs, rhs);
		break;

		case GEQ:
			ret = eval->GreaterEqualOperator(lhs, rhs);
		break;

		case LEQ:
			ret = eval->LessEqualOperator(lhs, rhs);
		break;
	}

	ret.SetMask(Currency::MASK_LOGICAL);
	return ret;
}

Currency EvaluatorInterface::BinaryOperator(const Currency& lhs, const Currency&rhs, int op)
{
    return eval->BinaryOperator(lhs, rhs, op);
}

Currency EvaluatorInterface::LogicalOperator(const Currency& lhs, const Currency&rhs, int op)
{
    return eval->LogicalOperator(lhs, rhs, op);
}

Currency EvaluatorInterface::UnaryOperator(const Currency& operand, int op)
{
    return eval->UnaryOperator(operand, op);
}

void EvaluatorInterface::RestorePath()
{
    eval->RestorePath();
}

std::string EvaluatorInterface::GetLastErrorMessage()
{
    return ExprTreeEvaluator::GetLastErrorMessage();
}

void EvaluatorInterface::SetLastErrorMessage(const std::string& msg)
{
    ExprTreeEvaluator::SetLastErrorMessage(msg);
}

#undef FormatMessage
std::string EvaluatorInterface::FormatMessage(const std::string& base_message)
{
    return eval->FormatErrorMessage(base_message);
}

std::string EvaluatorInterface::GetLastWarning()
{
    return ExprTreeEvaluator::GetLastWarning();
}

void EvaluatorInterface::SetLastWarning(const std::string& msg)
{
    ExprTreeEvaluator::SetLastWarning(msg);
}

hwMatrix* EvaluatorInterface::allocateMatrix()
{
    return ExprTreeEvaluator::allocateMatrix();
}

hwMatrix* EvaluatorInterface::allocateMatrix(const hwMatrix* mtx)
{
    return ExprTreeEvaluator::allocateMatrix(mtx);
}

hwMatrix* EvaluatorInterface::allocateMatrix(int m, int n, void* data, hwMatrix::DataType type)
{
    return ExprTreeEvaluator::allocateMatrix(m, n, data, type);
}

hwMatrix* EvaluatorInterface::allocateMatrix(int m, int n, hwMatrix::DataType type)
{
    return ExprTreeEvaluator::allocateMatrix(m, n, type);
}

hwMatrix* EvaluatorInterface::allocateMatrix(int m, int n, double value)
{
    return ExprTreeEvaluator::allocateMatrix(m, n, value);
}

hwMatrix* EvaluatorInterface::allocateMatrix(int m, int n, hwComplex& value)
{
    return ExprTreeEvaluator::allocateMatrix(m, n, value);
}

const hwMatrix* EvaluatorInterface::allocateColumn(const hwMatrix* mtx, int col)
{
    return ExprTreeEvaluator::allocateColumn(mtx, col);
}

hwMatrix* EvaluatorInterface::allocateMatrix(int m, int n, hwComplex&& value)
{
    return ExprTreeEvaluator::allocateMatrix(m, n, value);
}

hwMatrixI* EvaluatorInterface::allocateMatrix(int m, int n, hwMatrixI::DataType type)
{
    return ExprTreeEvaluator::allocateMatrix(m, n, type);
}

hwMatrixI* EvaluatorInterface::allocateMatrix(int m, int n, int val)
{
    return ExprTreeEvaluator::allocateMatrix(m, n, val);
}
HML_CELLARRAY* EvaluatorInterface::allocateCellArray()
{
    return ExprTreeEvaluator::allocateCellArray();
}

hwMatrixN* EvaluatorInterface::allocateMatrixN()
{
    return ExprTreeEvaluator::allocateMatrixN();
}

hwMatrixN* EvaluatorInterface::allocateMatrixN(const hwMatrixN* mtx)
{
    return ExprTreeEvaluator::allocateMatrixN(mtx);
}

HML_CELLARRAY* EvaluatorInterface::allocateCellArray(int m, int n)
{
    return ExprTreeEvaluator::allocateCellArray(m, n);
}

HML_CELLARRAY* EvaluatorInterface::allocateCellArray(const HML_CELLARRAY* cell)
{
    return ExprTreeEvaluator::allocateCellArray(cell);
}

StructData* EvaluatorInterface::allocateStruct(const StructData* strct)
{
    return ExprTreeEvaluator::allocateStruct(strct);
}

StructData* EvaluatorInterface::allocateStruct()
{
    return ExprTreeEvaluator::allocateStruct();
}

bool EvaluatorInterface::isUsedForEvalin() const
{
    return eval->IsUsedForEvalin();
}

void EvaluatorInterface::transferResultsFrom(const EvaluatorInterface& other)
{
    std::vector<Currency> results = other.eval->GetOutputResults();
    for (size_t i = 0; i < results.size(); ++i)
        eval->PushResult(results[i]);
}

bool EvaluatorInterface::HasBuiltin(const std::string& func_name) const
{
    return eval->HasBuiltin(func_name);
}

int EvaluatorInterface::NargoutFor(const std::string& func_name) const
{
    return eval->NargoutFor(func_name);
}

int EvaluatorInterface::NarginFor(const std::string& func_name) const
{
    return eval->NarginFor(func_name);
}

Currency EvaluatorInterface::VariableIndex(const Currency& target, const std::vector<Currency>& params)
{
    return eval->VariableIndex(target, params);
}

Currency EvaluatorInterface::CellIndex(const Currency& target, const std::vector<Currency>& params)
{
    return eval->CellValueHelper(target, params);
}

Currency EvaluatorInterface::Assignment(Currency& target, const std::vector<Currency>& indices, const Currency& value)
{
    eval->AssignHelper(target, indices, value);
	return target;
}

Currency EvaluatorInterface::CellAssignment(Currency& target, const std::vector<Currency>& params, const Currency& value)
{
    eval->CellAssignmentHelper(target, params, value);
	return target;
}

void EvaluatorInterface::LockCurrent()
{    
    eval->LockCurrent();
}    
     
void EvaluatorInterface::UnlockCurrent()
{    
    eval->UnlockCurrent();
}    
     
void EvaluatorInterface::Unlock(const std::string& fname)
{    
    eval->Unlock(fname);
}    
     
bool EvaluatorInterface::IsLocked(const std::string& fname) const
{
    return eval->IsLocked(fname);
}
     
bool EvaluatorInterface::IsCurrentLocked() const
{
    return eval->IsCurrentLocked();
}

std::string EvaluatorInterface::GetCurrentFilename() const
{
    return eval->GetCurrentFilename();
}

void EvaluatorInterface::SetDiary(bool val)
{
	eval->SetDiary(val);
}

bool EvaluatorInterface::IsDiaryOpen()
{
	return eval->IsDiaryOpen();
}

void EvaluatorInterface::SetDiary(std::string filename)
{
	eval->SetDiary(filename);
}

//------------------------------------------------------------------------------
//! Prompts to save before exiting
//------------------------------------------------------------------------------
void EvaluatorInterface::OnSaveOnExit()
{
    assert(eval);
    eval->OnSaveOnExit();
}
//------------------------------------------------------------------------------
//! True if evaluator is quitting
//------------------------------------------------------------------------------
bool EvaluatorInterface::IsQuit()
{
    assert(eval);
    return eval->IsQuit();
}
//------------------------------------------------------------------------------
//! Sets the quit flag
//! \param[in] val Sets to true if the evaluator is quitting
//------------------------------------------------------------------------------
void EvaluatorInterface::SetQuit(bool val)
{
    assert(eval);
    eval->SetQuit(val);
}

FunctionInfo* EvaluatorInterface::FunctionInfoFromString(const std::string& str) const
{
	return eval->FunctionInfoFromString(str);
}
//------------------------------------------------------------------------------
//! Start pause
//! \param[in] msg  User message to display
//! \param[in] wait True if waiting for a keystroke input from user
//------------------------------------------------------------------------------
void EvaluatorInterface::OnPauseStart(const std::string& msg, bool wait)
{
    assert(eval);
    eval->OnPauseStart(msg, wait);
}
//------------------------------------------------------------------------------
//! End pause
//------------------------------------------------------------------------------
void EvaluatorInterface::OnPauseEnd()
{
    assert(eval);
    eval->OnPauseEnd();
}
//------------------------------------------------------------------------------
//! Gets user input
//! \param[in]  prompt Prompt to display to the user
//! \param[in]  type   Type, if specified
//! \param[out] input  Input from user
//------------------------------------------------------------------------------
void EvaluatorInterface::OnUserInput(const std::string& prompt,
                                     const std::string& type,
		                             std::string&       input)
{
    assert(eval);
    eval->OnUserInput(prompt, type, input);
}
//------------------------------------------------------------------------------
//! Refreshes directories in client
//------------------------------------------------------------------------------
void EvaluatorInterface::OnRefreshDirs()
{
    assert(eval);
    eval->OnRefreshDirs();
}
//------------------------------------------------------------------------------
//! Change current working directory in client
//! \param[in] dir Fully qualified path of the new directory
//------------------------------------------------------------------------------
void EvaluatorInterface::OnChangeDir(const std::string& dir)
{
    assert(eval);
    eval->OnChangeDir(dir);
}

// fenv functions
int EvaluatorInterface::GetBaseEnvHandle()
{
	return eval->GetBaseEnvHandle();
}

int EvaluatorInterface::GetCurrentEnvHandle()
{
	return eval->GetCurrentEnvHandle();
}

int EvaluatorInterface::GetNewEnvHandle()
{
	return eval->GetNewEnvHandle();
}

Currency EvaluatorInterface::GetEnvValue(int handle, std::string varname)
{
	return eval->GetEnvValue(handle, varname);
}

void EvaluatorInterface::SetEnvValue(int handle, std::string varname, const Currency& new_val)
{
	eval->SetEnvValue(handle, varname, new_val);
}

void EvaluatorInterface::RemoveEnvValue(int handle, std::string varname)
{
	eval->RemoveEnvValue(handle, varname);
}

void EvaluatorInterface::ImportEnv(int handle1, int handle2)
{
	eval->ImportEnv(handle1, handle2);
}
//------------------------------------------------------------------------------
//! Returns true if evaluator is in experimental mode
//------------------------------------------------------------------------------
bool EvaluatorInterface::GetExperimental() const
{
    assert(eval);
    return eval->GetExperimental();
}
//------------------------------------------------------------------------------
//! Sets the experimental mode flag (-ex)
//! \param[in] val Sets to true if the evaluator is in experimental mode
//------------------------------------------------------------------------------
void EvaluatorInterface::SetExperimental(bool val)
{
    assert(eval);
    eval->SetExperimental(val);
}

void EvaluatorInterface::CacheLineInfomation()
{
	eval->CacheLineInfomation();
}

void EvaluatorInterface::UncacheLineInfomation()
{
	eval->UncacheLineInfomation();
}
//------------------------------------------------------------------------------
//! Returns true if successful in registering (swig) bound class
//! \param[in] name Bound class name
//! \param[in] info Bound class info
//------------------------------------------------------------------------------
bool EvaluatorInterface::RegisterBoundClass(const std::string& name,
                                            BoundClassInfo*    info)
{
    assert(eval);
    return eval->RegisterBoundClass(name, info);
}
//------------------------------------------------------------------------------
//! Gets signal handler
//------------------------------------------------------------------------------
SignalHandlerBase* EvaluatorInterface::GetSignalHandler() const
{
    assert(eval);
    return eval->GetSignalHandler();
}
//------------------------------------------------------------------------------
//! Returns a vector of function names
//------------------------------------------------------------------------------
std::vector<std::string> EvaluatorInterface::GetFunctionNames() const
{
    assert(eval);
    return eval->GetFunctionNames();
}
//------------------------------------------------------------------------------
//! Returns a count of functions
//------------------------------------------------------------------------------
int EvaluatorInterface::GetFunctionCount() const
{
    std::vector<std::string> funcnames (GetFunctionNames());

    return (funcnames.empty() ? 0 : static_cast<int>(funcnames.size()));
}
//------------------------------------------------------------------------------
//! Suspend function list updates in language
//------------------------------------------------------------------------------
void EvaluatorInterface::SuspendFuncListUpdate()
{
    eval->SuspendFuncListUpdate();
}
//------------------------------------------------------------------------------
//! Unsuspends function list updates, call OnUpdateFuncList to refresh
//------------------------------------------------------------------------------
void EvaluatorInterface::UnsuspendFuncListUpdate()
{
    eval->UnsuspendFuncListUpdate();
}
//------------------------------------------------------------------------------
//! Updates function list in language
//------------------------------------------------------------------------------
void EvaluatorInterface::OnUpdateFuncList()
{
    eval->OnUpdateFuncList();
}
//------------------------------------------------------------------------------
//! True if an interrupt has been requested
//------------------------------------------------------------------------------
bool EvaluatorInterface::IsInterrupt() const
{
    return eval->IsInterrupt();
}
//------------------------------------------------------------------------------
//! True if an interrupt has been requested
//------------------------------------------------------------------------------
void EvaluatorInterface::RegisterOMLDecryptor(const std::string& extension, ENCRPTR ptr)
{
    eval->RegisterOMLDecryptor(extension, ptr);
}
//------------------------------------------------------------------------------
//! Gets the application directory
//------------------------------------------------------------------------------
std::string EvaluatorInterface::GetApplicationDir() const
{
    return eval->GetApplicationDir();
}
//------------------------------------------------------------------------------
//! Sets the application directory
//------------------------------------------------------------------------------
void EvaluatorInterface::SetApplicationDir(const std::string& dir)
{
    eval->SetApplicationDir(dir);
}

void EvaluatorInterface::RunEncryptedFile(const std::string& file, const std::string& extension)
{
	eval->RunEncryptedFile(extension, file);
}

bool EvaluatorInterface::IsExtensionEncrypted(const std::string& extension)
{
	return eval->IsExtensionEncrypted(extension);
}

void EvaluatorInterface::RunPrecompiledFile(const std::string& filename)
{
	eval->RunPrecompiledFile(filename);
}

void EvaluatorInterface::RunFile(const std::string& filename)
{
	eval->RunFile(filename);
}

//------------------------------------------------------------------------------
//! Returns true if the given string is a keyword
//------------------------------------------------------------------------------
bool EvaluatorInterface::IsKeyword(const std::string& in) const
{
    return eval->IsKeyword(in);
}
//------------------------------------------------------------------------------
//! Returns true if the given string is an operator
//------------------------------------------------------------------------------
bool EvaluatorInterface::IsOperator(const std::string& in) const
{
    return eval->IsOperator(in);
}
//------------------------------------------------------------------------------
//! Returns true if the given string is a valid OML expression
//------------------------------------------------------------------------------	
std::string EvaluatorInterface::IsValidString(const std::string& in) const
{
	return eval->IsValidString(in);
}

void EvaluatorInterface::WritePFile(const std::string& infile, const std::string& outfile)
{
	eval->WritePFile(infile, outfile);
}

Currency EvaluatorInterface::Analyze(const std::string& infile)
{
	return eval->Analyze(infile);
}