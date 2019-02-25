/**
* @file MemoryScope.cpp
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

#include "MemoryScope.h"
#include "FunctionInfo.h"
#include "OML_Error.h"
#include <algorithm>

std::map<std::string, Currency> MemoryScope::globals;

int MemoryScopeManager::_env_counter = 0;

MemoryScope::MemoryScope(const MemoryScope& in)
{
	scope            = in.scope;
	global_names     = in.global_names;
	fi               = in.fi;
	debug_filename   = in.debug_filename;
	debug_line       = in.debug_line;
	nested_functions = in.nested_functions;
}

MemoryScope::MemoryScope(const MemoryScope& in, FunctionInfo* finfo)
{
	scope            = in.scope;
	global_names     = in.global_names;
	fi               = finfo;
	debug_filename   = in.debug_filename;
	debug_line       = in.debug_line;

	if (finfo)
		nested_functions = in.nested_functions;
}

MemoryScope::~MemoryScope()
{
	std::unordered_map<const std::string*, FunctionInfo*>::const_iterator iter;

	for (iter = nested_functions.begin(); iter != nested_functions.end(); iter++)
		delete iter->second;
}

const Currency& MemoryScope::GetValue(const std::string& varname) const
{ 
	static Currency _not_found(-1.0, Currency::TYPE_NOTHING);

	if (std::find(global_names.begin(), global_names.end(), varname) != global_names.end())
		return globals[varname];

	if (fi && fi->Persistent())
	{
		const Currency& temp_cur  = fi->Persistent()->GetValue(varname);

		if (!temp_cur.IsNothing())
			return temp_cur;
	}

	std::map<const std::string*, Currency>::const_iterator temp;
	const std::string* var_ptr = Currency::vm.GetStringPointer(varname);
	temp = scope.find(var_ptr);

	if (temp == scope.end())
		return _not_found;

	return temp->second;
}

const Currency& MemoryScope::GetValue(const std::string* var_ptr) const
{ 
	static Currency _not_found(-1.0, Currency::TYPE_NOTHING);

	MemoryScope* anon_scope = NULL;
	
	if (fi)
		anon_scope = fi->AnonScope();

	if (anon_scope && anon_scope->global_names.size())
	{
		if (std::find(anon_scope->global_names.begin(), anon_scope->global_names.end(), *var_ptr) != anon_scope->global_names.end())
			return globals[*var_ptr];
	}
	else
	{
		if (std::find(global_names.begin(), global_names.end(), *var_ptr) != global_names.end())
			return globals[*var_ptr];
	}

	if (fi && fi->Persistent())
	{
		if (fi->Persistent()->Contains(var_ptr))
		{
			const Currency& temp_cur = fi->Persistent()->GetValue(var_ptr);

			if (!temp_cur.IsNothing())
				return temp_cur;
		}
	}

	std::map<const std::string*, Currency>::const_iterator temp;
	temp = scope.find(var_ptr);

	if (temp == scope.end())
	{
		if (fi && anon_scope && anon_scope->scope.size())
		{
			temp = anon_scope->scope.find(var_ptr);
			
			if (temp == anon_scope->scope.end())
				return _not_found;
			else
				return temp->second;
		}
		else
		{
			return _not_found;
		}
	}

	return temp->second;
}

Currency& MemoryScope::GetMutableValue(const std::string& varname)
{
	if (std::find(global_names.begin(), global_names.end(), varname) != global_names.end())
		return globals[varname];

	const std::string* var_ptr = Currency::vm.GetStringPointer(varname);

	if (fi && fi->Persistent())
	{
		if (fi->Persistent()->Contains(var_ptr))
			return fi->Persistent()->GetMutableValue(var_ptr);
	}

	std::map<const std::string*, Currency>::iterator temp = scope.find(var_ptr);

	if (temp == scope.end())
	{
		Currency empty;
		scope[var_ptr] = empty;
		temp = scope.find(var_ptr);
	}

	return temp->second;
}

Currency& MemoryScope::GetMutableValue(const std::string* var_ptr)
{
	if (std::find(global_names.begin(), global_names.end(), *var_ptr) != global_names.end())
		return globals[*var_ptr];

	if (fi && fi->Persistent())
	{
		if (fi->Persistent()->Contains(var_ptr))
			return fi->Persistent()->GetMutableValue(var_ptr);
	}

	std::map<const std::string*, Currency>::iterator temp = scope.find(var_ptr);

	if (temp == scope.end())
	{
		Currency empty;
		scope[var_ptr] = empty;
		temp = scope.find(var_ptr);
	}

	return temp->second;
}

Currency* MemoryScope::GetMutablePointer(const std::string* var_ptr)
{
	if (std::find(global_names.begin(), global_names.end(), *var_ptr) != global_names.end())
		return &(globals[*var_ptr]);

	if (fi && fi->Persistent())
	{
		if (fi->Persistent()->Contains(var_ptr))
			return fi->Persistent()->GetMutablePointer(var_ptr);
	}

	std::map<const std::string*, Currency>::iterator temp = scope.find(var_ptr);

	if (temp == scope.end())
	{
		Currency empty;
		scope[var_ptr] = empty;
		temp = scope.find(var_ptr);
	}

	temp->second.SetOutputName(var_ptr);
	return &(temp->second);
}

bool MemoryScope::Contains(const std::string* var_ptr) const
{
	std::map<const std::string*, Currency>::const_iterator temp = scope.find(var_ptr);

	if (temp == scope.end())
		return false;

	return true;
}

void MemoryScope::SetValue(const std::string& varname, const Currency& new_val)
{ 
	new_val.SetOutputName(varname);

	if (std::find(global_names.begin(), global_names.end(), varname) != global_names.end())
	{
		globals[varname] = new_val;
		return;
	}
	
	if (fi && fi->Persistent())
	{
		const Currency& temp_cur = fi->Persistent()->GetValue(varname);

		if (!temp_cur.IsNothing())
		{
			fi->Persistent()->SetValue(varname, new_val); 
			return;
		}
	}

	const std::string* var_ptr = Currency::vm.GetStringPointer(varname);
	std::map<const std::string*, Currency>::iterator iter = scope.find(var_ptr);
	if (iter == scope.end())
	{
		std::pair<const std::string*, Currency> new_pair(var_ptr, new_val);
		scope.insert(new_pair);
	}
	else
	{
		iter->second = new_val;
	}
}

void MemoryScope::SetValue(const std::string* var_ptr, const Currency& new_val)
{ 
	new_val.SetOutputName(var_ptr);

	if (std::find(global_names.begin(), global_names.end(), *var_ptr) != global_names.end())
	{
		globals[*var_ptr] = new_val;
		return;
	}
	
	if (fi && fi->Persistent())
	{
		if (fi->Persistent()->Contains(var_ptr))
		{
			fi->Persistent()->SetValue(*var_ptr, new_val); 
			return;
		}
	}

	std::map<const std::string*, Currency>::iterator iter = scope.find(var_ptr);
	if (iter == scope.end())
	{
		std::pair<const std::string*, Currency> new_pair(var_ptr, new_val);
		scope.insert(new_pair);
	}
	else
	{
		iter->second = new_val;
	}
}

void MemoryScope::AddGlobalReference(const std::string& varname)
{
	if (std::find(global_names.begin(), global_names.end(), varname) == global_names.end())
	{
		global_names.insert(varname);

		if (globals.find(varname) == globals.end()) // someone else might already have declared it
		{
			hwMatrix* empty = ExprTreeEvaluator::allocateMatrix();
			globals[varname] = empty;
		}
	}
}

void MemoryScope::RemoveGlobalReference(const std::string& varname)
{
	global_names.erase(varname);
}

void MemoryScope::AddPersistentReference(const std::string& varname)
{
	if (!fi || !fi->Persistent())
		throw OML_Error("Invalid use of persistent keyword"); // can't have persistent at the outermost level

	const Currency& temp_cur = fi->Persistent()->GetValue(varname);

	if (!temp_cur.IsNothing())
		return;

	fi->Persistent()->SetValue(varname, ExprTreeEvaluator::allocateMatrix());
}

void MemoryScope::Remove(const std::string& varname)
{
	const std::string* var_ptr = Currency::vm.GetStringPointer(varname);
	std::map<const std::string*, Currency>::iterator iter = scope.find(var_ptr);

	if (iter != scope.end())
		scope.erase(iter);
}

bool MemoryScope::Remove(const std::regex& varname)
{
    bool rv = false;

    std::smatch results;
    for (auto iter = scope.begin(); iter != scope.end();)
    {
        if (std::regex_match(*(iter->first), results, varname))
        {
            rv = true;
            scope.erase(iter++);
        }
        else
        {
            ++iter;
        }
    }
    return rv;
}

void MemoryScope::ClearLocals()
{
	scope.clear();
}

void MemoryScope::ClearGlobals()
{
	globals.clear();
    global_names.clear();
}

void MemoryScope::ClearFromGlobals(const std::string& str)
{
    globals.erase(str);
    global_names.erase(str);
}

bool MemoryScope::ClearFromGlobals(const std::regex& varname)
{
    bool rv = false;
    std::smatch results;
    for (auto iter = globals.begin(); iter != globals.end();)
    {
        if (std::regex_match(iter->first, results, varname))
        {
            rv = true;
            global_names.erase(iter->first);
            globals.erase(iter++);
        }
        else
        {
            ++iter;
        }
    }
    return rv;
}

void MemoryScope::HideGlobal(const std::string& str)
{
    global_names.erase(str);
}

bool MemoryScope::IsNested() const
{
	if (fi && fi->IsNested())
		return true;
	else
		return false;
}

bool MemoryScope::IsAnonymous() const
{
	if (fi && fi->IsAnonymous())
		return true;
	else
		return false;
}

bool MemoryScope::IsGlobal(const std::string& varname) const
{
	return global_names.find(varname) != global_names.end();
}

bool MemoryScope::IsPersistent(const std::string* varname) const
{
	if (fi && fi->Persistent())
	{
		MemoryScope* pers = fi->Persistent();

		if (pers->Contains(varname))
			return true;
	}

	return false;
}

std::string MemoryScope::FunctionName() const
{
	if (fi)
		return fi->FunctionName();

	return "";
}

std::vector<std::string> MemoryScope::GetVariableNames() const
{
	std::map<const std::string*, Currency>::const_iterator iter;
	std::set<std::string>::const_iterator iter2;
	std::vector<std::string> varnames;

	for (iter = scope.begin(); iter != scope.end(); iter++)
		varnames.push_back(*(iter->first));

	for (iter2 = global_names.begin(); iter2 != global_names.end(); iter2++)
	{
		if (std::find(varnames.begin(), varnames.end(), *iter2) == varnames.end())
			varnames.push_back(*iter2);
	}

	std::sort(varnames.begin(), varnames.end());

	return varnames;
}

std::vector<const std::string*> MemoryScope::GetVariableNamePtrs() const
{
	std::map<const std::string*, Currency>::const_iterator iter;
	std::vector<const std::string*> varnames;

	for (iter = scope.begin(); iter != scope.end(); iter++)
		varnames.push_back(iter->first);

	return varnames;
}

FunctionInfo* MemoryScope::GetNestedFunction(const std::string* func_name)
{
	if (nested_functions.find(func_name) != nested_functions.end())
		return nested_functions[func_name];
	else
		return NULL;
}

FunctionInfo* MemoryScope::GetLocalFunction(const std::string* func_name)
{
	if (!fi)
		return NULL;

	if (fi->local_functions->find(func_name) != fi->local_functions->end())
		return (*fi->local_functions)[func_name];
	else
		return NULL;
}

const std::string& MemoryScope::GetFilename() const
{
	if (debug_filename)
		return *debug_filename;

	return *(Currency::pm.GetStringPointer("Unknown"));
}

void MemoryScope::RegisterNestedFunction(FunctionInfo* fi)
{
	const std::string* fi_name = fi->FunctionNamePtr();

	std::unordered_map<const std::string*, FunctionInfo*>::iterator iter = nested_functions.find(fi_name);
	if (iter != nested_functions.end())
		delete iter->second;

	nested_functions[fi_name] = fi; 
}

MemoryScopeManager::~MemoryScopeManager()
{
    if (delete_scopes)
    {
        for (auto iter = memory_stack.begin(); iter != memory_stack.end(); iter++)
            delete *iter;
    }
}

MemoryScope* MemoryScopeManager::GetScope(int offset) const
{
	size_t stack_size = memory_stack.size();

	if (stack_size <= offset)
		return NULL; 
	else
		return memory_stack[stack_size-offset-1];
}

MemoryScope* MemoryScopeManager::GetCurrentScope() const
{
	return GetScope(0);
}

MemoryScope* MemoryScopeManager::GetParentScope() const
{
	return GetScope(1);
}

void MemoryScopeManager::OpenScope(FunctionInfo*fi)
{
#if _DEBUG
	#define MAX_VAL 50
#else
	#define MAX_VAL 100
#endif

	if (memory_stack.size() == MAX_VAL) 
		throw OML_Error("Maximum function depth reached");

	MemoryScope* temp = new MemoryScope(fi);
	memory_stack.push_back(temp);
}

void MemoryScopeManager::CloseScope()
{
	if (memory_stack.empty())
		return;

	// need to delete the last one
	MemoryScope* temp = memory_stack.back();

	if (temp->IsNested())
	{
		// copy local variables into parent scope
		// Evaluator will handle removal of the function from the table
		size_t       stack_size = memory_stack.size();
		MemoryScope* parent = memory_stack[stack_size-2];

		for (auto iter=temp->scope.begin(); iter != temp->scope.end(); iter++)
		{
			std::string varname = (*iter->first);

			const std::string* var_ptr = Currency::vm.GetStringPointer(varname);

			if (temp->fi->IsInputParameter(var_ptr))
				continue;

			if (temp->fi->IsReturnValue(var_ptr))
				continue;

			if (parent->fi->IsReferenced(var_ptr))
				parent->SetValue(varname, iter->second);
		}
	}

	ClearEnv(temp);

	delete temp;
	memory_stack.pop_back();

	if (first_scope_with_nested == GetStackDepth())
		first_scope_with_nested = -1;
}

void MemoryScopeManager::CopyFullScope(MemoryScopeManager* ref)
{
	std::vector<MemoryScope*>::const_iterator iter;

	for (iter = ref->memory_stack.begin(); iter != ref->memory_stack.end(); iter++)
	{
		MemoryScope* temp_scope = new MemoryScope(*(*iter));
		memory_stack.push_back(temp_scope);
	}
}

void MemoryScopeManager::CloseTo(MemoryScope* target)
{
	while (!memory_stack.empty())
	{
		MemoryScope* current = GetCurrentScope();

		if (current == target)
			break;

		CloseScope();
	}
}

const Currency& MemoryScopeManager::GetValue(const std::string& varname, int offset) const
{
	static Currency _not_found(-1.0, Currency::TYPE_NOTHING);

	MemoryScope* scope = GetScope(offset);

	const Currency& temp_cur = scope->GetValue(varname);

	if (!temp_cur.IsNothing())
		return temp_cur;
	
	if (scope->IsNested())
	{
		// if we didn't find it, and the scope is nested, look
		// in the parent's scope as well (and repeat until
		// we hit a non-nested scope.  The -2/-1 are because we
		// already checked the "current" scope and didn't find it.
		size_t stack_size = memory_stack.size();

		for (int j=0; j<stack_size-1; j++)
		{
			int index = (int)stack_size-2-j;

			MemoryScope* local = memory_stack[index];

			const Currency& loc_cur = local->GetValue(varname);

			if (!loc_cur.IsNothing())
				return loc_cur;

			if (local->IsNested())
				break;
		}
	}

	return _not_found;
}


const Currency& MemoryScopeManager::GetValue(const std::string* var_ptr, int offset) const
{
	static Currency _not_found(-1.0, Currency::TYPE_NOTHING);

	MemoryScope* scope = GetScope(offset);

	const Currency& temp_cur = scope->GetValue(var_ptr);

	if (!temp_cur.IsNothing())
		return temp_cur;
	
	if (scope->IsNested())
	{
		// if we didn't find it, and the scope is nested, look
		// in the parent's scope as well (and repeat until
		// we hit a non-nested scope.  The -2/-1 are because we
		// already checked the "current" scope and didn't find it.
		size_t stack_size = memory_stack.size();

		for (int j=0; j<stack_size-1; j++)
		{
			int index = (int)stack_size-2-j;

			MemoryScope* local = memory_stack[index];

			const Currency& loc_cur = local->GetValue(var_ptr);

			if (!loc_cur.IsNothing())
				return loc_cur;

			if (!local->IsNested())
				break;
		}
	}

	return _not_found;
}

void MemoryScopeManager::SetValue(const std::string& varname, const Currency& new_val)
{
	GetCurrentScope()->SetValue(varname, new_val);
}

void MemoryScopeManager::SetValue(const std::string* var_ptr, const Currency& new_val)
{
	GetCurrentScope()->SetValue(var_ptr, new_val);
}

void MemoryScopeManager::SetParentValue(const std::string& varname, const Currency& new_val)
{
	GetParentScope()->SetValue(varname, new_val);
}

void MemoryScopeManager::SetParentValue(const std::string* var_ptr, const Currency& new_val)
{
	GetParentScope()->SetValue(var_ptr, new_val);
}

void MemoryScopeManager::AddGlobalReference(const std::string& varname)
{
	GetCurrentScope()->AddGlobalReference(varname);
}

void MemoryScopeManager::RemoveGlobalReference(const std::string& varname)
{
	GetCurrentScope()->RemoveGlobalReference(varname);
}

void MemoryScopeManager::AddPersistentReference(const std::string& varname)
{
	GetCurrentScope()->AddPersistentReference(varname);
}

void MemoryScopeManager::Remove(const std::string& varname)
{
	GetCurrentScope()->Remove(varname);
}

bool MemoryScopeManager::Remove(const std::regex& varname)
{
	return GetCurrentScope()->Remove(varname);
}

void MemoryScopeManager::ClearLocals()
{
	GetCurrentScope()->ClearLocals();
}

void MemoryScopeManager::ClearGlobals()
{
    GetCurrentScope()->ClearGlobals();
}

void MemoryScopeManager::ClearFromGlobals(const std::string& varname)
{
    GetCurrentScope()->ClearFromGlobals(varname);
}

bool MemoryScopeManager::ClearFromGlobals(const std::regex& varname)
{
    return GetCurrentScope()->ClearFromGlobals(varname);
}

void MemoryScopeManager::HideGlobal(const std::string& varname)
{
    GetCurrentScope()->HideGlobal(varname);
}

Currency& MemoryScopeManager::GetMutableValue(const std::string& varname)
{
	return GetCurrentScope()->GetMutableValue(varname);
}

Currency& MemoryScopeManager::GetMutableValue(const std::string* var_ptr)
{
	return GetCurrentScope()->GetMutableValue(var_ptr);
}

Currency* MemoryScopeManager::GetMutablePointer(const std::string* var_ptr)
{
	return GetCurrentScope()->GetMutablePointer(var_ptr);
}

Currency& MemoryScopeManager::GetMutableParentValue(const std::string& varname)
{
	return GetParentScope()->GetMutableValue(varname);
}

bool MemoryScopeManager::Contains(const std::string* var_ptr) const
{
	return GetCurrentScope()->Contains(var_ptr);
}

bool MemoryScopeManager::IsGlobal(const std::string& varname) const
{
	return GetCurrentScope()->IsGlobal(varname);
}

bool MemoryScopeManager::IsPersistent(const std::string* var_ptr) const
{
	return GetCurrentScope()->IsPersistent(var_ptr);
}

std::vector<DebugStateInfo> MemoryScopeManager::GetCallStack() const
{
	std::vector<DebugStateInfo> call_stack;
	std::vector<MemoryScope*>::const_iterator iter;

	for (iter=memory_stack.begin(); iter!=memory_stack.end(); iter++)
	{
		std::string  fancy_call;
		MemoryScope* ms = (*iter);

		if (ms->FunctionName().size())
			fancy_call += ms->FunctionName();
		else
			fancy_call += "<main>";
		fancy_call += "(";

		int           num_params = 0;
		FunctionInfo* my_fi = ms->fi;

		std::vector<const std::string*> parameters; 

		if (my_fi)
		{
			parameters = my_fi->Parameters();
			num_params = (int)parameters.size();
		}

		for (int j=0; j<num_params; j++)
		{
			const std::string* param_name = parameters[j];
			fancy_call += ms->GetValue(param_name).GetOutputString(NULL);

			if (j !=(num_params-1))
				fancy_call += ", ";
		}

		fancy_call += ")";

		call_stack.push_back(DebugStateInfo(fancy_call, *ms->debug_filename, ms->debug_line));
	}

	return call_stack;
}

MemoryScopeManager* MemoryScopeManager::MakeContextCopy(bool base) const
{
    MemoryScopeManager* msm = new MemoryScopeManager();
    msm->delete_scopes = false;
    if (memory_stack.size())
    {
        if (base || memory_stack.size() == 1)
        {
            // grab the base scope
            msm->memory_stack.push_back(GetScope((int)memory_stack.size() - 1));
        }
        else
        {
            int end_scope = 1;
            while (end_scope < memory_stack.size() - 1 && GetScope(end_scope)->IsNested())
                ++end_scope;

            for (int i = end_scope; i >= 1; --i)
                msm->memory_stack.push_back(GetScope(i));
        }
    }

    return msm;
}

FunctionInfo* MemoryScopeManager::GetNestedFunction(const std::string* func_name) const
{
	if (first_scope_with_nested == -1)
		return NULL;

	int stack_depth = GetStackDepth();

	for (int j=0; j<stack_depth-first_scope_with_nested; j++)
	{
		MemoryScope* scope = GetScope(j);

		FunctionInfo* fi = scope->GetNestedFunction(func_name);

		if (fi)
			return fi;
	}

	return NULL;
}

FunctionInfo* MemoryScopeManager::GetLocalFunction(const std::string* func_name) const
{
	int stack_depth = GetStackDepth();

	for (int j=0; j<stack_depth; j++)
	{
		MemoryScope* scope = GetScope(j);

		FunctionInfo* fi = scope->GetLocalFunction(func_name);

		if (fi)
			return fi;
	}

	return NULL;
}

void MemoryScopeManager::RegisterNestedFunction(FunctionInfo* fi)
{
	GetCurrentScope()->RegisterNestedFunction(fi);

	if (first_scope_with_nested == -1)
		first_scope_with_nested = GetStackDepth()-1;
}

// fenv functions
int MemoryScopeManager::GetBaseEnvHandle()
{
	int depth = GetStackDepth();
	MemoryScope* temp = GetScope(depth-1);	
	_envs[0] = temp;
	_rev_envs[temp] = 0;
	
	return 0;
}

int MemoryScopeManager::GetCurrentEnvHandle()
{
	_env_counter++;
	MemoryScope* temp   = GetCurrentScope();

	if (_rev_envs.find(temp) != _rev_envs.end())
		return _rev_envs[temp];

	_envs[_env_counter] = temp;
	_rev_envs[temp]     = _env_counter;

	return _env_counter;
}

int MemoryScopeManager::GetNewEnvHandle()
{
	_env_counter++;

	MemoryScope* temp      = new MemoryScope(NULL);
	_envs[-1*_env_counter] = temp;
	_rev_envs[temp]        = -1*_env_counter;

	return -1*_env_counter;
}

Currency MemoryScopeManager::GetEnvValue(int handle, std::string varname)
{
	std::map<int, MemoryScope*>::const_iterator temp;
	temp = _envs.find(handle);

	if (temp == _envs.end())
		throw OML_Error("Invalid environment handle");

	MemoryScope* ms = temp->second;
	return ms->GetValue(varname);
}

void MemoryScopeManager::RemoveEnvValue(int handle, std::string varname)
{
	std::map<int, MemoryScope*>::const_iterator temp;
	temp = _envs.find(handle);

	if (temp == _envs.end())
		throw OML_Error("Invalid environment handle");

	MemoryScope* ms = temp->second;
	ms->Remove(varname);
}

void MemoryScopeManager::SetEnvValue(int handle, std::string varname, const Currency& new_val)
{
	std::map<int, MemoryScope*>::const_iterator temp;
	temp = _envs.find(handle);

	if (temp == _envs.end())
		throw OML_Error("Invalid environment handle");

	MemoryScope* ms = temp->second;
	return ms->SetValue(varname, new_val);
}

void MemoryScopeManager::ImportEnv(int source, int dest)
{
	std::map<int, MemoryScope*>::const_iterator temp;
	temp = _envs.find(source);

	if (temp == _envs.end())
		throw OML_Error("Invalid environment handle");

	MemoryScope* source_ms = temp->second;

	temp = _envs.find(dest);

	if (temp == _envs.end())
		throw OML_Error("Invalid environment handle");

	MemoryScope* dest_ms = temp->second;

	std::vector<const std::string*> varname_ptrs;
	varname_ptrs = source_ms->GetVariableNamePtrs();

	for (int j=0; j<varname_ptrs.size(); j++)
	{
		const std::string* var_ptr = varname_ptrs[j];

		//if (!dest_ms->Contains(var_ptr))
			dest_ms->SetValue(var_ptr, source_ms->GetValue(var_ptr));
	}
}

void MemoryScopeManager::ClearEnv(MemoryScope* ms)
{
	std::map<MemoryScope*, int>::const_iterator temp;	

	temp = _rev_envs.find(ms);

	if (temp != _rev_envs.end())
	{
		int idx = temp->second;
		_rev_envs.erase(temp);

		std::map<int, MemoryScope*>::const_iterator temp2;
		temp2 = _envs.find(idx);
		_envs.erase(temp2);
	}
}

bool MemoryScopeManager::IsFunctionInCallStack(const std::string& target)
{
	for (int j=0; j<memory_stack.size(); j++)
	{
		MemoryScope*  scope = memory_stack[j];
		FunctionInfo* fi    = scope->GetFunctionInfo();
		if (fi && fi->FunctionName() == target)
			return true;
	}

	return false;
}