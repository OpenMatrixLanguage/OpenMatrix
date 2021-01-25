/**
* @file MemoryScope.h
* @date August 2013
* Copyright (C) 2013-2019 Altair Engineering, Inc.  
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

#ifndef __MemoryScope_h
#define __MemoryScope_h

#include "Currency.h"
#include "EvaluatorDebug.h"
#include <map>
#include <unordered_map>
#include <set>
#include <regex>

class FunctionInfo;

class MemoryScope
{
	friend class MemoryScopeManager;

public:
    //!
    //! Constructor
    //! \param fi Function info pointer
    //!
	MemoryScope();
    MemoryScope(FunctionInfo* info);
	~MemoryScope();
	MemoryScope(const MemoryScope&);
	MemoryScope(const MemoryScope&, FunctionInfo* fi);

	void Initialize(FunctionInfo* info);

	// I would like to hide these as well, but I need to provide access to them
	// in one particular case - namely anonymous functions
	const Currency& GetValue(const std::string& varname) const;
	const Currency& GetValue(const std::string* var_ptr) const;
	void            SetValue(const std::string& varname, const Currency& new_val);
	void            SetValue(const std::string* var_ptr, const Currency& new_val);
	bool            IsGlobal(const std::string& varname) const;
	bool            IsPersistent(const std::string* varname) const;
	bool            Contains(const std::string* var_ptr) const;
	bool            IsEmpty() const;

	bool            IsAnonymous() const;
	bool            IsNested() const;
	std::string     FunctionName() const;

	std::vector<std::string> GetVariableNames() const;
	std::vector<const std::string*> GetVariableNamePtrs() const;

	FunctionInfo*   GetNestedFunction(const std::string* func_name);
	FunctionInfo*   GetLocalFunction(const std::string* func_name);

	FunctionInfo*   GetFunctionInfo() const { return fi; }
	void            SetFunctionInfo(FunctionInfo* new_fi) { fi = new_fi; }

	void               SetDebugInfo(const std::string* fname, int line) { debug_filename = fname; debug_line = line; }
	const std::string& GetFilename() const;
	const std::string* GetFilenamePtr() const { return debug_filename; }
	int                GetLineNumber() const { return debug_line; }
	void               RegisterNestedFunction(FunctionInfo* fi);

	void			   Remove(const std::string& varname);

	void               AddGlobalReference(const std::string& varname);
	void			   Reset();

protected:
	void            RemoveGlobalReference(const std::string& varname);
	void            AddPersistentReference(const std::string& varname);
	bool            Remove(const std::regex& varname);
	void            HideGlobal(const std::string& varname);
	void            ClearLocals();
	void            ClearGlobals();
	void            ClearFromGlobals(const std::string& varname);
	bool            ClearFromGlobals(const std::regex& varname);
	Currency&       GetMutableValue(const std::string& varname);
	Currency&       GetMutableValue(const std::string* var_ptr);
	Currency*       GetMutablePointer(const std::string* var_ptr);

private:
	std::map<const std::string*, Currency> scope;
	std::set<std::string> global_names;

	static std::map<std::string, Currency> globals;

	FunctionInfo*       fi;
    const std::string*  debug_filename;
	int                 debug_line;
};

class MemoryScopeManager
{
public:
	MemoryScopeManager();
	~MemoryScopeManager();

	MemoryScope* GetCurrentScope() const;	
	MemoryScope* GetScope(int offset) const;
	MemoryScope* GetParentScope() const;	

	void CopyFullScope(MemoryScopeManager* ref);

	void OpenScope(FunctionInfo* fi);
	void CloseScope();
	void CloseTo(MemoryScope* target);

	const Currency& GetValue(const std::string& varname, int offset=0) const;
	const Currency& GetValue(const std::string* var_ptr, int offset=0) const;
	void            SetValue(const std::string& varname, const Currency& new_val);
	void            SetValue(const std::string* var_ptr, const Currency& new_val);
	void            SetParentValue(const std::string& varname, const Currency& new_val);
	void            SetParentValue(const std::string* var_ptr, const Currency& new_val);
	void            AddGlobalReference(const std::string& varname);
	void            RemoveGlobalReference(const std::string& varname);
	void            AddPersistentReference(const std::string& varname);
	void            Remove(const std::string& varname);
	bool            Remove(const std::regex& varname);
	void            ClearLocals();
	void            ClearGlobals();
	void            ClearFromGlobals(const std::string& varname);
	bool            ClearFromGlobals(const std::regex& varname);
	void            HideGlobal(const std::string& varname);
	Currency&       GetMutableValue(const std::string& varname);
	Currency&       GetMutableValue(const std::string* var_ptr);
	Currency*       GetMutablePointer(const std::string* var_ptr);
	Currency&       GetMutableParentValue(const std::string& varname);
	bool            IsGlobal(const std::string& varname) const;
	bool            IsPersistent(const std::string* var_ptr) const;
	bool            Contains(const std::string* var_ptr) const;

	int                         GetStackDepth() const { return (int) memory_stack.size(); }
	std::vector<DebugStateInfo> GetCallStack() const;
	MemoryScopeManager*         MakeContextCopy(bool base) const;
	bool                        IsFunctionInCallStack(const std::string&);

	FunctionInfo*   GetNestedFunction(const std::string* func_name) const;
	FunctionInfo*   GetLocalFunction(const std::string* func_name) const;

	void     RegisterNestedFunction(FunctionInfo* fi);

	int      GetBaseEnvHandle();
	int      GetCurrentEnvHandle();
	int      GetNewEnvHandle();
	Currency GetEnvValue(int, std::string);
	void     RemoveEnvValue(int, std::string);
	void     SetEnvValue(int, std::string, const Currency&);
	void     ImportEnv(int source, int dest);
	void     ClearEnv(MemoryScope*);
	void     ClearEnv(int source);

private:
	std::vector<MemoryScope*> memory_stack;
	bool delete_scopes;

	std::map<int, MemoryScope*> _envs;
	std::map<MemoryScope*, int> _rev_envs;
	std::vector<MemoryScope*> _allocated_envs;
	static int _env_counter;

	int first_scope_with_nested;

	MemoryScope* stack_pool;
	std::vector<MemoryScope*> unused_scopes;
	int pool_in_use;
};
#endif
