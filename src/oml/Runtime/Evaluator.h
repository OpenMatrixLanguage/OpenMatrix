/**
* @file Evaluator.h
* @date August 2013
* Copyright (C) 2013-2018 Altair Engineering, Inc.  
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

#ifndef __ExprTreeEvaluator_h
#define __ExprTreeEvaluator_h
#define FIRST_USER_FILE 3
#define NUM_MANDATORY_PATHS 0
#if OS_WIN
#define DIRECTORY_DELIM '\\'
#else
#define DIRECTORY_DELIM '/'
#endif

#include "Currency.h"
#include "MemoryScope.h"
#include "OutputFormat.h"

#include <map>
#include <iostream>
#include <fstream>
using namespace std;

#include "EvaluatorInt.h"

class BoundClassInfo;         // Interface for (swig) bound classes
class ClassInfo;              // Interface for external classes in OML language
class ExprTreeEvaluator;
class EvaluatorDebugInterface;
class FunctionInfo;
class SignalHandlerBase;      // Base implementation for handling client signals
class OMLTree;

// forward declarations for ANTLR stuff
struct ANTLR3_BASE_TREE_struct;
typedef ANTLR3_BASE_TREE_struct* pANTLR3_BASE_TREE;
struct ANTLR3_COMMON_TOKEN_STREAM_struct;
typedef	ANTLR3_COMMON_TOKEN_STREAM_struct* pANTLR3_COMMON_TOKEN_STREAM;
struct ANTLR3_VECTOR_struct;
typedef ANTLR3_VECTOR_struct* pANTLR3_VECTOR;

struct UserFile
{
    bool        isopen;
    std::FILE*  file;
    std::string name;
    std::string mode;
    UserFile    (std::FILE* filepointer, std::string filename, std::string filemode) :
                isopen(true), file(filepointer), name(filename), mode(filemode) {}
};

struct BuiltinFunc
{
    FUNCPTR     fptr;
	ALT_FUNCPTR alt_fptr;

	FunctionMetaData md;
    BuiltinFunc(FUNCPTR func, FunctionMetaData in_md) : fptr(func), alt_fptr(nullptr), md(in_md) {}
    BuiltinFunc(ALT_FUNCPTR func) : fptr(nullptr), alt_fptr(func) {}
    BuiltinFunc() : fptr(nullptr), alt_fptr(nullptr) {}
};

struct UserFunc
{
    FunctionInfo* fi;
    bool locked;
    UserFunc(FunctionInfo* func) : fi(func), locked(false) {}
    UserFunc() : fi(nullptr), locked(false) {}
    ~UserFunc();
};

class DebugInfo
{
public:
	DebugInfo(int line_num, const std::string* file_name) : lineNum(line_num), fileName(file_name) {}

	bool               Matches(OMLTree* tree) const;
	int                LineNum() const { return lineNum; }
	std::string        Filename() const { return *fileName; }
	const std::string* FilenamePtr() const {return fileName;}

	static DebugInfo DebugInfoFromTree(OMLTree* tree);

private:
	int                lineNum;
	const std::string* fileName;
};

typedef Currency (ExprTreeEvaluator::*TREE_FPTR)(OMLTree*);

class ExprTreeEvaluator
{
public:
	ExprTreeEvaluator();
	ExprTreeEvaluator(const ExprTreeEvaluator* source);

	~ExprTreeEvaluator();

	void     PreprocessAST(pANTLR3_BASE_TREE, pANTLR3_COMMON_TOKEN_STREAM);
	Currency RunTree(OMLTree*);

	static TREE_FPTR GetFuncPointerFromType(int type);

	const Currency& GetValue(std::string varname, int offset=0) const;
	const Currency& GetGlobalValue(std::string varname) const;
	bool            SetGlobalValue(std::string varname, const Currency& value);
	bool            SetValue(std::string varname, const Currency& value);
    bool            IsGlobal(const std::string& varname) const { return msm->IsGlobal(varname); }
	void            Clear(const std::string& varname);
    bool            ClearFromFunctions(const std::string& name);
    bool            ClearFromGlobals  (const std::string& name);
    bool            ClearFromVariables(const std::string& name);
	void            Clear(const std::regex& varname);
    bool            ClearFromFunctions(const std::regex& name);
    bool            ClearFromGlobals  (const std::regex& name);
    bool            ClearFromVariables(const std::regex& name);
    void            ClearFunctions();
    void            ClearGlobals();
    void            ClearVariables();
	int             RenameVariable(std::string old_name, std::string new_name);
	bool            Contains(const std::string& varname) const;

    void            LockCurrent();
    void            UnlockCurrent();
    void            Unlock(const std::string& fname);
    bool            IsCurrentLocked() const;
    bool            IsLocked(const std::string& fname) const;

    std::string     GetCurrentFilename() const;

	std::vector<std::string> GetVariableNames(int offset) const;
	std::vector<std::string> GetBuiltinFunctionNames() const;
	std::vector<std::string> GetFunctionNames() const;
	std::vector<std::string> GetKeywords() const;
	std::vector<std::string> GetOperators() const;

	int   GetContextEndValue() const;
	int   GetNargoutValue() const;
	int   GetNarginValue() const;
	void  PushNargValues(int argsin, int argsout);
	void  PopNargValues();

	void  Mark();
	void  Unmark();
	void  Restore();

	bool  RegisterBuiltInFunction(const std::string& func_name, FUNCPTR fp, int nargin, int nargout);
	bool  RegisterBuiltInFunction(const std::string& func_name, FUNCPTR fp, FunctionMetaData fmd);
	bool  RegisterBuiltInFunction(const std::string& func_name, ALT_FUNCPTR fp);
	void  RegisterExternalVariableFunction(EXTPTR ep) { ext_var_ptr = ep; }
	void  ClearTemporaryExternalVariables() { _externals.clear(); }

    //! Returns true if successful in registering (swig) bound class
    //! \param[in] name Bound class name
    //! \param[in] info Bound class info
    bool RegisterBoundClass( const std::string& name,
                             BoundClassInfo*    info);

	int                         GetStackDepth() const { return msm->GetStackDepth(); }
	std::vector<DebugStateInfo> GetCallStack() const { return msm->GetCallStack(); }
 	void                        GetStackInfo(int level, std::string& filename, int& line_number);

    std::vector<Currency> DoMultiReturnFunctionCall(FunctionInfo* fi, const std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<std::string>* out_vars = nullptr);
    std::vector<Currency> DoMultiReturnFunctionCall(FUNCPTR fptr, const std::string& func_name, const std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<std::string>* out_vars);

	void DoMultiReturnFunctionCall(FUNCPTR fptr,     const std::string& func_name, const std::vector<Currency>& param_values, int num_ins, bool suppress_output, OMLTree* out_tree);
	void DoMultiReturnFunctionCall(ALT_FUNCPTR fptr,     const std::string& func_name, const std::vector<Currency>& param_values, int num_ins, bool suppress_output, OMLTree* out_tree);
	void DoMultiReturnFunctionCall(FunctionInfo* fi, const std::vector<Currency>& param_values, int num_ins, bool suppress_output, OMLTree* out_tree);

	std::vector<Currency> DoAnonymousMultiReturnFunctionCall(FunctionInfo* fi, const std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<std::string>* out_vars = nullptr);

    bool     FindFunctionByName(const std::string&, FunctionInfo**, FUNCPTR*);

    Currency CallFunction(const std::string&, const std::vector<Currency>&);
	Currency CallFunction(const std::string*, const std::vector<Currency>&);
    void     CallFunction(const std::string&, const std::vector<Currency>&, std::vector<Currency>&, int nargout = -1);
    Currency CallBuiltinFunction(FUNCPTR, const std::string&, const std::vector<Currency>&);
    Currency CallBuiltinFunction(ALT_FUNCPTR, const std::string&, const std::vector<Currency>&);
	Currency CallInternalFunction(FunctionInfo*, const std::vector<Currency>&);
	Currency CallInternalFunction(FunctionInfo*, const std::vector<Currency>&, const std::string&, const Currency&);

    HML_CELLARRAY* CreateVararginCell(const std::vector<Currency>& params, int start_index);

    inline bool IsUserFunction(const std::string& func_name) { return functions->find(func_name) != functions->end(); }
    inline bool IsStdFunction(const std::string& func_name)  { return std_functions->find(func_name) != std_functions->end(); }
    FUNCPTR GetStdFunction(const std::string& func_name) const;
	std::string GetHelpModule(const std::string& func_name);

	bool IsA(const Currency& target, const std::string& classname) const;
	FunctionInfo* GetBaseClassFunctionInfo(const std::string* func_name, std::string* base_name, Currency* base_val);

	bool IsKeyword(const std::string& func_name) const;
	bool IsOperator(const std::string& func_name) const;

	std::vector<Currency>	GetOutputResults() const;
	void					ClearOutputResults();
	Currency				GetLastResult() const;
	void				    PushResult(const Currency& res, bool to_output=true);
	void				    PrintResult(Currency res);

	bool					CloseAllFiles();
	int						AddFile(std::FILE *newfile, const std::string &fname, const std::string &fmode);
	int						GetNumFiles() { return (int) userFileStreams->size(); }
    void                    AddPath(std::string pathname, bool end);
	void                    AddPath2(const std::string& pathname, const std::vector<std::string>& func_names);
    bool                    RemovePath(std::string &pathname);
    inline void             ClearPath() { paths->erase(paths->begin() + NUM_MANDATORY_PATHS, paths->end()); }
    bool                    FindFileInPath(const std::string& file_plus_ext, std::string& filepath) const;
    void                    RestorePath();

    bool                    IsInPaths   (const std::string &str) const;
    inline const std::vector<std::string>& GetPaths() const { return *paths; }

	// these assume the index is valid
	bool					CloseFile   (int i);
    inline std::FILE*		GetFile     (int i) { return userFileStreams->at(i).file; }
	inline std::string		GetFileName (int i) { return userFileStreams->at(i).name; }
	inline std::string		GetFileMode (int i) { return userFileStreams->at(i).mode; }
	std::vector<int>		GetFileIndices(int start);

    inline const OutputFormat* GetOutputFormat() const { return format; }
    inline void  SetOutputFormat(const OutputFormat& fmt) { *format = fmt; } 
    
    inline void ResetFuncSearchCache() { not_found_functions->clear(); } 

    //! True if using the debugger
    bool IsDebugging() const { return (debug_listener ? true : false); }
    //! Sets the debug listener
	void RegisterDebugListener(EvaluatorDebugInterface* edi) { debug_listener = edi; }

    Currency EqualityOperator   (const Currency& lhs, const Currency& rhs, int op);
    Currency LessThanOperator   (const Currency& lhs, const Currency& rhs);
	Currency GreaterThanOperator(const Currency& lhs, const Currency& rhs);
	Currency LessEqualOperator   (const Currency& lhs, const Currency& rhs);
	Currency GreaterEqualOperator(const Currency& lhs, const Currency& rhs);
    Currency LogicalOperator    (const Currency& lhs, const Currency& rhs, int op);
    Currency BinaryOperator     (const Currency& lhs, const Currency& rhs, int op);
    Currency UnaryOperator      (const Currency& operand, int op);

	std::string FormatErrorMessage(const std::string& base_message);
    /* set/get the last error message */
    static std::string GetLastErrorMessage();
    static void SetLastErrorMessage(const std::string& lasterr);
    static std::string GetLastWarning();
    static void SetLastWarning(const std::string& lastwarn);
	void ErrorCleanup();

    inline bool IsUsedForEvalin() const { return is_for_evalin; }
    ExprTreeEvaluator* MakeContextCopy(bool base, bool pop_nargs) const;

    static hwMatrix*       allocateMatrix();
    static hwMatrix*       allocateMatrix(const hwMatrix*);
    static hwMatrix*       allocateMatrix(int m, int n, void* data, hwMatrix::DataType type);
    static hwMatrix*       allocateMatrix(int m, int n, hwMatrix::DataType type);
    static hwMatrix*       allocateMatrix(int m, int n, double value);
    static hwMatrix*       allocateMatrix(int m, int n, hwComplex& value);
    static hwMatrixI*      allocateMatrix(int m, int n, hwMatrixI::DataType type);
    static hwMatrixI*      allocateMatrix(int m, int n, int val);

    static const hwMatrix* allocateColumn(const hwMatrix* mtx, int col);

    static hwMatrixN*      allocateMatrixN();
	static hwMatrixN*      allocateMatrixN(const std::vector<int>& dims, const hwMatrixN::DataType& dataType);
	static hwMatrixN*      allocateMatrixN(const hwMatrixN*);

    static HML_CELLARRAY*  allocateCellArray();
    static HML_CELLARRAY*  allocateCellArray(int m, int n);
    static HML_CELLARRAY*  allocateCellArray(const HML_CELLARRAY*);

	static HML_ND_CELLARRAY*  allocateNDCellArray();
	static HML_ND_CELLARRAY*  allocateNDCellArray(std::vector<int> dims);
	static HML_ND_CELLARRAY*  allocateNDCellArray(const HML_ND_CELLARRAY*);

    static StructData*     allocateStruct(const StructData*);
    static StructData*     allocateStruct();

    bool HasBuiltin(const std::string& func_name) const;
    int NargoutFor(const std::string& func_name) const;
    int NarginFor(const std::string& func_name) const;

	void     SetInterrupt(bool);
	bool	 IsInterrupt();

    //!
    //! Returns true if there is a pause request pending
    //!
    bool IsPauseRequestPending() const { return _pauseRequestPending; }
    //!
    //! Requests/resets pause. Pause request/reset will be processed after 
    //! execution of current statement
    //! \param val True if a pause request is intiated, false otherwise
    //!
    void SetPauseRequestPending(bool val) { _pauseRequestPending = val; }
    //!
    //! Returns true if interpreter is paused
    //!
    bool IsPaused() const { return _paused; }

	void     SetDiary(bool);
	void     SetDiary(std::string filename);
	bool     IsDiaryOpen();
	void     WriteToDiaryFile(std::string filename);

	void StoreSuppressedResults(bool store) { _store_suppressed = store; }

    //! Returns true if evaluator is quitting
    bool IsQuit() const { return _quit; }
    //! Sets the quit flag
    //! \param[in] val Sets to true if the evaluator is quitting
    void SetQuit( bool val) { _quit = val; }

    //! Returns true if evaluator is in experimental mode
    bool GetExperimental() const;
    //! Sets the experimental mode flag (-ex)
    //! \param[in] val Sets to true if the evaluator is in experimental mode
    void SetExperimental( bool val);

	Currency  VariableIndex(const Currency&, const std::vector<Currency>&);
    Currency  CellValueHelper(const Currency&, const std::vector<Currency>&);
	Currency  NDCellValueHelper(const Currency&, const std::vector<Currency>&);

	void AssignHelper(Currency& target, const std::vector<Currency>& indices, const Currency& value, int refcnt_target=1);
    void CellAssignmentHelper(Currency& target, const std::vector<Currency>& params, const Currency& value);
    void NDAssignmetHelper(Currency& target, const std::vector<Currency>& params, const Currency& value, int refcnt_target=1);
	void NDCellAssignmetHelper(Currency& target, const std::vector<Currency>& params, const Currency& value, int refcnt_target = 1);

	void SetScriptName(std::string script) { _script_name = script; }

	FunctionInfo* FunctionInfoFromString(const std::string&);

    //! Gets signal handler
    SignalHandlerBase* GetSignalHandler() const { return _signalHandler; }
    //! Sets signal handler, if not null or to null if already set
    //! \param[in] handler Signal handler
    void SetSignalHandler( SignalHandlerBase* handler);

    // Actions handled in the client

    //! Clear results stored in client
    void OnClearResults();
    //! Print the given currency
    //! \param[in] cur Currency to print
    void OnPrintResult( const Currency& cur);

    //! Start pause
    //! \param[in] msg  User message to display
    //! \param[in] wait True if waiting for a keystroke input from user
    void OnPauseStart( const std::string& msg, 
                       bool               wait);
    //! End pause
    void OnPauseEnd();

    //! Change current working directory in client
    //! \param[in] dir Fully qualified path of the new directory
    void OnChangeDir( const std::string& dir);
    //! Refreshes directories in client
    void OnRefreshDirs();

    //! Prompts to save before exiting
	void OnSaveOnExit();
    //! Update function list in language
    void OnUpdateFuncList();
    //! Gets user input
    //! \param[in]  prompt Prompt to display to user
    //! \param[in]  type   Type, if specified
    //! \param[out] input  Input from user
    void OnUserInput( const std::string& prompt,
                      const std::string& type,
                      std::string&       input);

	std::string GetCurrentDebugFile() const;
	int         GetCurrentDebugLine() const;
	void        SetDebugInfo(const std::string* file, int lin);
	void        CacheLineInfomation();
	void        UncacheLineInfomation();

	int      GetBaseEnvHandle();
	int      GetCurrentEnvHandle();
	int      GetNewEnvHandle();
	Currency GetEnvValue(int, std::string);
	void     RemoveEnvValue(int, std::string);
	void     SetEnvValue(int, std::string, const Currency&);
	void     ImportEnv(int source, int dest);

    //! Suspend function list updates in language
    void SuspendFuncListUpdate() { _suspendFunclistUpdate = true; }
    //! Unsuspends function list updates, call OnUpdateFuncList to refresh
    void UnsuspendFuncListUpdate() { _suspendFunclistUpdate = false; }

	void EncryptOMLFile(const std::string& in_file, const std::string& out_file);
	void RegisterOMLDecryptor(const std::string& extenstion, ENCRPTR ptr);
	void RunEncryptedFile(const std::string& extension, const std::string& file_name);
	bool IsExtensionEncrypted(const std::string& extension);

	void RunPrecompiledFile(const std::string& extension);
	void RunFile(const std::string& filename);

	std::string IsValidString(const std::string& in) const;

    //! Gets the application directory
    std::string GetApplicationDir() const { return _appdir; }
    //! Sets the application directory
    void SetApplicationDir( const std::string& dir);

	void WritePFile(const std::string& infile, const std::string& outfile);
	
	Currency Analyze(const std::string& infile);
	Currency GetMetadata(const std::string& infile);

private:
	Currency AddOperator(const Currency&, const Currency&);
	Currency SubtractOperator(const Currency&, const Currency&);
	Currency MultiplyOperator(const Currency&, const Currency&);
	Currency EntrywiseMultiplyOperator(const Currency&, const Currency&);
	Currency DivideOperator(const Currency&, const Currency&);
	Currency EntrywiseDivideOperator(const Currency&, const Currency&);
	Currency LeftDivideOperator(const Currency&, const Currency&);
	Currency EntrywiseLeftDivideOperator(const Currency&, const Currency&);

	Currency PowOperator(const Currency&, const Currency&);
	Currency EntrywisePowOperator(const Currency&, const Currency&);
	Currency NegateOperator(const Currency&);
	Currency NotOperator(const Currency&);

	Currency Number(OMLTree* tree);
	Currency HexNumber(OMLTree* tree);
	Currency TextString(OMLTree* tree);
	Currency Identifier(OMLTree* tree);
	Currency Nothing(OMLTree* tree);
	Currency Return(OMLTree* tree);
	Currency AssignOperator(OMLTree* tree);
	Currency CellAssignOperator(OMLTree* tree);
	Currency UnaryOperator(OMLTree* tree);
	Currency BinaryOperator(OMLTree* tree);
    Currency EqualityOperator(OMLTree* tree);
	Currency EqualityOperatorEx(const Currency& lhs, const Currency& rhs);
	bool     EqualityHelper(const Currency& lhs, const Currency& rhs);
	Currency LessThanOperator(OMLTree* tree);
	Currency GreaterThanOperator(OMLTree* tree);
	Currency LessEqualOperator(OMLTree* tree);
	Currency GreaterEqualOperator(OMLTree* tree);
    Currency LogicalOperator(OMLTree* tree);
    Currency LogicalOperatorEx(const Currency& lhs, const Currency& rhs, int op);
	bool     LogicalHelper(double lhs, double rhs, int op);
    Currency ShortCircuitAndOperator(OMLTree* tree);
    Currency ShortCircuitOrOperator(OMLTree* tree);
	bool     ShortCircuitHelper(OMLTree* tree);
	Currency FunctionCall(OMLTree* tree);	
	Currency MultiReturnFunctionCall(OMLTree* tree);
	Currency Conditional(OMLTree* tree);
	Currency SwitchCase(OMLTree* tree);
	Currency WhileLoop(OMLTree* tree);
	Currency ForLoop(OMLTree* tree);
	Currency FunctionDefinition(OMLTree* tree);
	Currency MatrixCreation(OMLTree* tree);
	Currency StatementList(OMLTree* tree);
	Currency Statement(OMLTree* tree);
	Currency RangeOperator(OMLTree* tree);
	Currency GlobalReference(OMLTree* tree);
	Currency PersistentReference(OMLTree* tree);
	Currency TransposeOperator(OMLTree* tree);
	Currency ConjTransposeOperator(OMLTree* tree);
	Currency AnonymousFunctionDefinition(OMLTree* tree);
	Currency CellArrayCreation(OMLTree* tree);
	Currency CellValue(OMLTree* tree);
	void     CellAssignmentHelper(Currency& target, OMLTree* indices, const Currency& value);
	Currency InlineIndex(OMLTree* tree);
	Currency InlineIndexCell(OMLTree* tree);
	Currency StructValue(OMLTree* tree);
    Currency StructValueHelper(const Currency* parent, OMLTree* indices, OMLTree* field_tree);
	Currency ObjectMethodCall(Currency* parent, OMLTree* indices, OMLTree* field_tree);
	Currency MRObjectMethodCall(OMLTree* tree);
	Currency Clear(OMLTree* tree);
	Currency TryCatch(OMLTree* tree);
	Currency CellExtraction(OMLTree* tree);
	Currency InPlaceExpansion(OMLTree* tree);
	Currency ClassDefinition(OMLTree* tree);

	Currency AssignmentUtility(OMLTree* tree, const Currency& value);

	hwMatrix* SubmatrixSingleIndexHelper(Currency& target, const Currency& indices, const Currency& value, int target_refcnt=1);
	hwMatrix* SubmatrixDoubleIndexHelper(Currency& target, const Currency& index1, const Currency& index2, const Currency& value, int target_refcnt=1);
	void      ReplaceMatrixElementHelper(hwMatrix*& target, const int index1, const int index2, const double value);
	void      ReplaceMatrixElementHelper(hwMatrix*& target, int index1, int index2, const hwComplex& value);
	void      ReplaceMatrixElementHelper(hwMatrix*& target, int index1, double value);
	void      ReplaceMatrixElementHelper(hwMatrix*& target, int index1, const hwComplex& value);
	Currency  VariableIndex(const std::string*, const std::vector<Currency>&);

	Currency  CallOverloadedOperator(const std::string&, const Currency& op);	
	Currency  CallOverloadedOperator(const std::string&, const Currency& op1, const Currency& op2);
	bool      HasOverloadedFunction(const Currency&, const std::string&);
	Currency  CallOverloadedFunction(const std::string&, const std::vector<Currency>& params);

	Currency  CheckForExternalVariable(const std::string&);

	Currency  PostFunctionIndexHelper(const Currency& target, OMLTree* tree);
	
	MemoryScope* GetCurrentScope() const;		

	void OpenScope(FunctionInfo* info = NULL);
	void CloseScope();

    void ImportPathNames   (const ExprTreeEvaluator*);
	void ImportMemoryScope (const ExprTreeEvaluator*);
	void ImportFunctionList(const ExprTreeEvaluator*);
	void ImportUserFileList(const ExprTreeEvaluator*);
	
	bool FindFunction(const std::string& func_name, std::string& file_name);
	bool FindPrecompiledFunction(const std::string& func_name, std::string& file_name);
	bool FindEncryptedFunction(const std::string& func_name, std::string& file_name, std::string& extension);
	bool ParseAndRunFile(const std::string& file_name, bool allow_script);
	bool ParseAndRunString(const std::string& str, const std::string& use_filename);
    bool RemoveFuncs(std::map<std::string, UserFunc*>& funcs, const std::regex& name);
    void RemoveFunctionsInPath(const std::string& path);
	bool CheckForFunctionInAST(const std::string& func_name);

    void Lock(const std::string& fname);
    void SetLock(const std::string& fname, bool state);
    UserFunc* GetUserFunc(const std::string& fname) const;

	FunctionInfo* FunctionInfoFromTree(OMLTree* tree);

	std::vector<Currency> DoMultiReturnFunctionCall(FunctionInfo* fi, const std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<const std::string*>* out_vars = nullptr);
	std::vector<Currency> DoMultiReturnFunctionCall(FUNCPTR fp, const std::string& func_name, const std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<const std::string*>* out_vars = nullptr);
	std::vector<Currency> DoAnonymousMultiReturnFunctionCall(FunctionInfo* fi, const std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<const std::string*>* out_vars = nullptr);

	bool FindFunctionByName(const std::string*, FunctionInfo**, FUNCPTR*);

	bool CheckForPreviousTrailingDot(pANTLR3_BASE_TREE tree, pANTLR3_COMMON_TOKEN_STREAM tokens);
	bool ValidateMatrix(const std::vector<std::vector<Currency>>& currencies, unsigned int*, unsigned int*);
	bool PadStringsIfNecessary(std::vector<std::vector<Currency>>& currencies);
	bool ValidateFunction(OMLTree* tree);
    void ValidateExtractionIndices(const hwMatrix& target, int index1, int index2);
	void ValidateRowIndex(const hwMatrix& target, int index);
	void ValidateColumnIndex(const hwMatrix& target, int index);

    //! Called when result needs to be processed for print - sets format/prints
    //! \param[in] cur      Currency to print
    //! \param[in] toOutput True if this has an output
    void ProcessResultForPrint( const Currency& cur, 
                                bool            toOutput);

    //! Gets result after assignment to a bound class/bound class property
    //! \param[in] in     Input currency
    //! \param[in] field  Field tree
    //! \param[in] value  Value to assign
    Currency BoundClassAssign( const Currency*   in, 
                               OMLTree* field,
                               const Currency&   value);
    //! Gets result after calling a method in a bound class
    //! \param[in] in    Input currencty
    //! \param[in] field Field info
    Currency BoundClassMethod( const Currency*   in, 
                               OMLTree* field);

private:

	//data members
	MemoryScopeManager* msm;
    SignalHandlerBase*  _signalHandler;    //! Handles signals

	std::map<std::string, UserFunc*> *functions;
    std::map<std::string, BuiltinFunc> *std_functions;
	std::vector<std::string> *preregistered_functions;

    std::vector<std::string> *not_found_functions;

	std::string _script_name;

	std::map<std::string, ClassInfo*>      *class_info_map;      //classdef classes
    std::map<std::string, BoundClassInfo*> _boundclassinfo; //! Bound classes

	std::map<std::string, ENCRPTR> _decryptors;

    void OnEvaluationEnd( const Currency& currency);
    bool PathMatches(const std::string& s1, const std::string& s2) const;

	std::string* end_context_varname;
	int          end_context_index;
	Currency*    end_context_currency;

	std::vector<int> nargin_values;
	std::vector<int> nargout_values;
	std::vector<std::string> user_func_calls; // used for keeping track of which funcs to lock/unlock
	int              assignment_nargout;

	int         current_statement_index;
	OMLTree*    current_tree;

	int         nested_function_marker;
	bool        _owns_msm;
    bool        _owns_functions;
    bool        _owns_userfiles;
    bool        _owns_pathnames;
    bool        _owns_format;
	bool        suppress_multi_ret_output;
	bool        _diary_on;
	bool        _store_suppressed;
    bool        _suspendFunclistUpdate; //! Suspends function list update

	std::string builtin_error_scope;

	ofstream       diaryFile;

	std::vector<Currency> results;
	Currency              last_suppressed_result;

    bool is_for_evalin;
	bool _lhs_eval;
	bool _break_on_continue;

	bool _interrupt;
    bool _quit;                 //! True if evaluator is quitting

	bool _pauseRequestPending;  //!< True if pause request is pending
    bool _paused;               //!< True if interpreter is paused


    OutputFormat* format;

	std::vector<MemoryScope*> marks;
	int          mark_narg_size;
	
    std::vector<UserFile>* userFileStreams;
    std::vector<std::string>* paths;

	EXTPTR ext_var_ptr;
	std::vector<Currency> _externals;

	EvaluatorDebugInterface* debug_listener;

    std::string _appdir;     //! Default directory for the application
	std::string cached_filename;
	int         cached_line;

    static std::string lasterrormsg;
    static std::string lastwarning;

	std::vector<hwSliceArg> slices;
	std::vector<int> indices;
};

#endif
