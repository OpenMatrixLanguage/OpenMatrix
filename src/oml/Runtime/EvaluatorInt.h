/**
* @file EvaluatorInt.h
* @date June 2014
* Copyright (C) 2014-2021 Altair Engineering, Inc.  
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

#ifndef __ExprTreeEvaluatorInterface_h
#define __ExprTreeEvaluatorInterface_h
#include "OMLDll.h"
#include "Currency.h"
#include "FunctionMetaData.h"
#include "OMLInterfacePublic.h"
#include <string>
#include <vector>
#include <regex>
#define OML_NO_NARG 0xFFFFFF00

class BoundClassInfo;
class ExprTreeEvaluator;
class EvaluatorInterface;
class FunctionInfo;
class OutputFormat;
class SignalHandlerBase;
class OMLImplBase;

typedef bool (*FUNCPTR) (EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
typedef bool (*ALT_FUNCPTR) (OMLInterface*, const OMLCurrencyList* ins, OMLCurrencyList* outs);
typedef void* (*ENCRPTR) (const std::string& filename);

template <typename T> class hwTComplex;
typedef hwTComplex<double> hwComplex;

template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;
typedef hwTMatrix<int, hwTComplex<int> > hwMatrixI;

template <typename T1, typename T2> class hwTMatrixN;
typedef hwTMatrixN<double, hwTComplex<double> > hwMatrixN;

class OMLDLL_DECLS EvaluatorInterface
{
    friend class Interpreter;
public:
    EvaluatorInterface(ExprTreeEvaluator* ptr) : eval(ptr), delete_eval(false) {}
    EvaluatorInterface(const EvaluatorInterface& other, bool base_context, bool pop_nargs = true);
    EvaluatorInterface(const EvaluatorInterface& other) : eval(other.eval), delete_eval(false) {}
    ~EvaluatorInterface();

    bool RegisterBuiltInFunction(const std::string& func_name, ALT_FUNCPTR fp, bool thread_safe= true);
	bool RegisterBuiltInFunction(const std::string& func_name, ALT_FUNCPTR fp, const FunctionMetaData& md, bool thread_safe = true);
    bool RegisterBuiltInFunction(const std::string& func_name, FUNCPTR fp, int nargin, int nargout, bool thread_safe = true);
    bool RegisterBuiltInFunction(const std::string& func_name, FUNCPTR fp, const FunctionMetaData& md, bool thread_safe = true);
    std::vector<std::string> GetBuiltinFunctionNames() const;

    //! Returns true if successful in registering (swig) bound class
    //! \param[in] name Bound class name
    //! \param[in] info Bound class info
    bool RegisterBoundClass( const std::string& name,
                             BoundClassInfo*    info);

    std::vector<std::string> GetKeywords() const;
    const OutputFormat* GetOutputFormat() const;
    bool FindFunctionByName(const std::string& func_name, FunctionInfo** fi, FUNCPTR* fptr, ALT_FUNCPTR* aptr);
    
    const Currency& GetValue(std::string varname) const;
	const Currency& GetGlobalValue(std::string varname) const;
    bool SetValue(std::string varname, const Currency& value);
    
    bool IsUserFunction(const std::string& func_name);
    bool IsStdFunction(const std::string& func_name);
    FUNCPTR GetStdFunction(const std::string& func_name) const;
    
    bool IsGlobal(const std::string& varname) const;
    bool Contains(const std::string& varname) const;
    std::vector<std::string> GetVariableNames() const;

	bool IsA(const Currency& target, const std::string& classname) const;

    Currency CallFunction(const std::string& func_name, const std::vector<Currency>& params);
    Currency CallInternalFunction(FunctionInfo* fi, const std::vector<Currency>& param_values);
    Currency CallStaticClassMethod(const std::string& class_name, const std::string& mehtod_name, const std::vector<Currency>& params);
    std::vector<Currency> DoMultiReturnFunctionCall(FUNCPTR fptr, std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<std::string>* out_vars = nullptr);
	std::vector<Currency> DoMultiReturnFunctionCall(ALT_FUNCPTR aptr, std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<std::string>* out_vars = nullptr);
    std::vector<Currency> DoMultiReturnFunctionCall(FunctionInfo* fi, std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<std::string>* out_vars = nullptr);
    HML_CELLARRAY* CreateVararginCell(const std::vector<Currency>& params, int start_index);

	void WritePFile(const std::string& infile, const std::string& outfile);

   	Currency Analyze(const std::string& infile);
   	Currency GetMetadata(const std::string& infile);

	std::vector<std::string> GetProperties(const std::string& classname);
	std::vector<std::string> GetMethods(const std::string& classname, bool public_only = true);

    void ClearPath();
    bool RemovePath(std::string &pathname);
	bool RemoveHiddenPath(std::string& pathname);
    bool RemoveHiddenRootPath(std::string& pathname);
    void AddPath(std::string pathname, bool end);
	void AddHiddenPath(std::string pathname);
	void AddPath2(const std::string& pathname, const std::vector<std::string> funcs);
	void RegisterLibraryAlias(const std::string& path, const std::string& alias);
	std::string GetLibraryAlias(const std::string& path);
    std::vector<std::string> GetPaths() const;
    void ResetFuncSearchCache();
    void TreatAsBuiltin(const std::string& path);
    
    std::FILE* GetFile(int i);
    bool FindFileInPath(const std::string& file_plus_ext, std::string& filepath) const;
    int AddFile(std::FILE *newfile, const std::string& fname, const std::string& fmode);
    std::string GetFileMode (int i);
    bool CloseAllFiles();
    bool CloseFile (int i);
    std::vector<int> GetFileIndices(int start);
    std::string GetFileName (int i);
    
    void PushResult(Currency res, bool to_output = true);
    void PrintResult(Currency res);

    void Clear             (const std::regex& varname);
    void ClearFromGlobals  (const std::regex& varname);
    void ClearFromFunctions(const std::regex& varname);
    void ClearFromVariables(const std::regex& varname);
    void Clear             (const std::string& varname);
    void ClearFromGlobals  (const std::string& varname);
    void ClearFromFunctions(const std::string& varname);
    void ClearFromVariables(const std::string& varname);
    void ClearFromVariablesExcept(const std::vector<std::string>& varname,
        const std::vector<std::regex>& varwildname,
        const std::vector<std::string>& exceptname,
        const std::vector<std::regex>& exceptwildname);
    void ClearFromGlobalsExcept(const std::vector<std::string>& varname,
        const std::vector<std::regex>& varwildname,
        const std::vector<std::string>& exceptname,
        const std::vector<std::regex>& exceptwildname);
    void ClearFunctions();
    void ClearGlobals();
    void ClearVariables();
    void ClearClassesAndObjects();

	bool RemoveLibrary(const std::string& lib_name);

    int GetContextEndValue() const;
	int GetNargoutValue() const;
	int GetNarginValue() const;
    int	GetNumFiles();

    void  Mark();
    void  Unmark();
    void  Restore();

    Currency EqualityOperator  (const Currency& lhs, const Currency&rhs, int op);
    Currency InequalityOperator(const Currency& lhs, const Currency&rhs, int op);
    Currency LogicalOperator   (const Currency& lhs, const Currency&rhs, int op);
    Currency BinaryOperator    (const Currency& lhs, const Currency&rhs, int op);
    Currency UnaryOperator     (const Currency& op1, int op);

    void RestorePath();
    void RefreshPathCache();

    static std::string GetLastErrorMessage();
    static void SetLastErrorMessage(const std::string& msg);
    static std::string GetLastWarning();
    static void SetLastWarning(const std::string& msg);
	std::string FormatMessage(const std::string& base_message);
    
    static hwMatrix* allocateMatrix();
    static hwMatrix* allocateMatrix(const hwMatrix*);
    static hwMatrix* allocateMatrix(int m, int n, void* data, bool isReal);
    static hwMatrix* allocateMatrix(int m, int n, bool isReal);
    static hwMatrix* allocateMatrix(int m, int n, double value);
    static hwMatrix* allocateMatrix(int m, int n, const hwComplex&  value);
    static hwMatrix* allocateMatrix(int m, int n, const hwComplex&& value);

	// I would like to eliminate these two, but the statistics toolbox uses them
    static hwMatrixI* allocateMatrixI(int m, int n, bool isReal);
    static hwMatrixI* allocateMatrixI(int m, int n, int val);

    static const hwMatrix* allocateColumn(const hwMatrix* mtx, int col);

	static hwMatrixN* allocateMatrixN();
    static hwMatrixN* allocateMatrixN(const std::vector<int>& dims, bool isReal);
    static hwMatrixN* allocateMatrixN(const hwMatrixN*);

    static HML_CELLARRAY* allocateCellArray();
    static HML_CELLARRAY* allocateCellArray(int m, int n);
    static HML_CELLARRAY* allocateCellArray(const HML_CELLARRAY*);

    static HML_ND_CELLARRAY* allocateNDCellArray();
    static HML_ND_CELLARRAY* allocateNDCellArray(std::vector<int> dims);
    static HML_ND_CELLARRAY* allocateNDCellArray(const HML_ND_CELLARRAY*);

    static StructData* allocateStruct(const StructData*);
    static StructData* allocateStruct();

	std::string   IsValidString(const std::string&) const;
	FunctionInfo* FunctionInfoFromString(const std::string&) const;

    bool isUsedForEvalin() const;

    void transferResultsFrom(const EvaluatorInterface& other);
    
    bool HasBuiltin(const std::string& func_name) const;
    int NargoutFor(const std::string& func_name) const;
    int NarginFor(const std::string& func_name) const;

    Currency VariableIndex(const Currency& target, const std::vector<Currency>& params);
    Currency     CellIndex(const Currency& target, const std::vector<Currency>& params);

    Currency     Assignment(Currency& target, const std::vector<Currency>& indices, const Currency& value);
    Currency CellAssignment(Currency& target, const std::vector<Currency>& params, const Currency& value);

    void LockCurrent();
    void UnlockCurrent();
    void Unlock(const std::string& fname);
    bool IsCurrentLocked() const;
    bool IsLocked(const std::string& fname) const;
    bool LockBuiltInFunction(const std::string& fname, bool hide = true);

    std::string GetCurrentFilename() const;
	int         GetCurrentLinenumber() const;

	void CacheLineInfomation();
	void UncacheLineInfomation();

	void DisableWarning(const std::string& id);
	void EnableWarning(const std::string& id);
	bool IsWarningDisabled(const std::string& id);

	void SetDLLContext(const std::string& dll_name);
	void SetDLLHierarchy(const std::string& dll_hierarchy);

    //! Returns a vector of function names
    std::vector<std::string> GetFunctionNames() const;
    //! Returns a count of functions
    int GetFunctionCount() const;

	void EncryptOMLFile(const std::string& in_file, const std::string& out_file);
	bool IsExtensionEncrypted(const std::string& extension);
	void RunEncryptedFile(const std::string& filename, const std::string& extension);

	void RunPrecompiledFile(const std::string& filename);
	void RunFile(const std::string& filename);
    // Actions handled in client

    //! Change current working directory in client
    //! \param[in] dir Fully qualified path of the new directory
    void OnChangeDir( const std::string& dir);
    //! Refreshes directories in client
    void OnRefreshDirs();
	//! Prompts to save before exiting
    //! \param[in] returnCode Code to exit the application with
	void OnSaveOnExit(int returnCode);
	//! Gets user input
	//! \param[in]  prompt Prompt to display to the user
    //! \param[in]  type   Type, if specified
    //! \param[out] input  Input from user
	void OnUserInput( const std::string& prompt,
                      const std::string& type,
		              std::string&       input);
    //! Start pause
    //! \param[in] msg  User message to display
    //! \param[in] wait True if waiting for a keystroke input from user
    void OnPauseStart( const std::string& msg, 
                       bool               wait);
    //! End pause
    void OnPauseEnd();

    //! True if evaluator is quitting
    bool IsQuit();
    //! Sets the quit flag
    //! \param[in] val Sets to true if the evaluator is quitting
    void SetQuit( bool val);

    //! Returns true if evaluator is in experimental mode
    bool GetExperimental() const;
    //! Sets the experimental mode flag (-ex)
    //! \param[in] val Sets to true if the evaluator is in experimental mode
    void SetExperimental( bool val);

	int  GetVerbose() const;
	void SetVerbose(int val);

    //! Gets signal handler
    SignalHandlerBase* GetSignalHandler() const;

	int      GetBaseEnvHandle();
	int		 GetCurrentEnvHandle();
	int		 GetNewEnvHandle();
	Currency GetEnvValue(int handle, std::string varname);
	void     SetEnvValue(int handle, std::string varname, const Currency& new_val);
	void     RemoveEnvValue(int handle, std::string varname);	
	void     ImportEnv(int handle1, int handle2);
	void     DeleteEnv(int handle1);

	void     RegisterFunc(OMLTree* def);

    //! Updates function list in language
    void OnUpdateFuncList();
    //! Suspend function list updates in language
    void SuspendFuncListUpdate();
    //! Unsuspends function list updates, call OnUpdateFuncList to refresh
    void UnsuspendFuncListUpdate();

    //! True if an interrupt has been requested
    bool IsInterrupt() const;
    //! Sets the interrupt flag
    //! \parm True if evaluator is interrupted
    void SetInterrupt(bool);

	void RegisterOMLDecryptor(const std::string& extension, ENCRPTR ptr);

	void RegisterDLL(void*);

    //! Gets the application directory
    std::string GetApplicationDir() const;
    //! Sets the application directory
    void SetApplicationDir( const std::string& dir);

    //! Returns true if the given string is a keyword
    //! \param[in] in Given input
    bool IsKeyword( const std::string& in) const;
    //! Returns true if the given string is an operator
    //! \param[in] in Given input
    bool IsOperator( const std::string& in) const;

	void RegisterChildEvaluator(ExprTreeEvaluator*);
	void RemoveChildEvaluator();

	void SetNumberOfThreads(int threads);

	Currency GetProfileData() const;
	void     ClearProfileData();
	void     Profile(bool on);

    std::string GetFunctionArgumentName(int index);

    void CacheBCIPointer(OMLImplBase* ptr);
    void BCIGarbageCollect();

    //!
    //! Returns diary filestream
    //!
    std::ofstream& GetDiary();

private:
    ExprTreeEvaluator* eval;
    bool delete_eval;
};

#endif
