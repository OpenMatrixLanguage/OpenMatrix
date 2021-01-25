/**
* @file FunctionInfo.cpp
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

#include "FunctionInfo.h"
#include "MemoryScope.h"

#include "ExprCppTreeLexer.h"

#include "OMLTree.h"

#include <cassert>

FunctionInfo::FunctionInfo(std::string func_name, std::vector<const std::string*> ret_vals, std::vector<const std::string*> params, 
	                       std::map<const std::string*, Currency> default_vals, OMLTree* stmt_list, std::string file_name, std::string help_str)
{
	_function_name    = Currency::vm.GetStringPointer(func_name);
	_file_name        = Currency::pm.GetStringPointer(file_name);
	_return_values    = ret_vals;
	_parameters       = params;
	_default_values   = default_vals;
	_statements       = stmt_list;
	_builtin          = NULL;
	_alt_fptr         = NULL;
	_anon_scope       = NULL;
	_help_string      = help_str;
	_refcnt           = 1;
	_parent_fi        = NULL;

	// the persistent scope doesn't have a FunctionInfo* (and therefore its own persistent), otherwise
	// you'll be in an infinite loop looking for persistents forever
	_persistent_scope = new MemoryScope(NULL);
	_is_nested        = false;
	_is_constructor   = false;
	_is_encrypted     = false;

	local_functions  = new std::map<const std::string*, FunctionInfo*>;
	nested_functions = new std::map<const std::string*, FunctionInfo*>;
}

FunctionInfo::FunctionInfo(std::string func_name, std::vector<const std::string*> ret_vals, std::vector<const std::string*> params, OMLTree* stmt_list, std::string file_name)
{
	_function_name    = Currency::vm.GetStringPointer(func_name);
	_file_name        = Currency::pm.GetStringPointer(file_name);
	_return_values    = ret_vals;
	_parameters       = params;
	_statements       = new OMLTree(*stmt_list);
	_builtin          = NULL;
	_alt_fptr         = NULL;
	_anon_scope       = NULL;
	_refcnt           = 1;
	_parent_fi        = NULL;

	// the persistent scope doesn't have a FunctionInfo* (and therefore its own persistent), otherwise
	// you'll be in an infinite loop looking for persistents forever
	_persistent_scope = new MemoryScope(NULL);
	_is_nested        = false;
	_is_constructor   = false;
	_is_encrypted     = false;

	local_functions  = new std::map<const std::string*, FunctionInfo*>;
	nested_functions = new std::map<const std::string*, FunctionInfo*>;
}

FunctionInfo::FunctionInfo(std::string func_name, FUNCPTR builtin_func)
{
	_function_name    = Currency::vm.GetStringPointer(func_name);
	_file_name        = Currency::pm.GetStringPointer("");
	_builtin          = builtin_func;
	_alt_fptr         = NULL;
	_statements       = NULL;
	_anon_scope       = NULL;
	_persistent_scope = NULL;
	_is_nested        = false;
	_is_constructor   = false; 
	_is_encrypted     = false;
	local_functions   = NULL;
	nested_functions  = NULL;
	_refcnt           = 1;
	_parent_fi        = NULL;
}

FunctionInfo::FunctionInfo(std::string func_name, ALT_FUNCPTR alt_func)
{
	_function_name = Currency::vm.GetStringPointer(func_name);
	_file_name = Currency::pm.GetStringPointer("");
	_builtin = NULL;
	_alt_fptr = alt_func;
	_statements = NULL;
	_anon_scope = NULL;
	_persistent_scope = NULL;
	_is_nested = false;
	_is_constructor = false;
	_is_encrypted = false;
	local_functions  = NULL;
	nested_functions = NULL;
	_refcnt = 1;
	_parent_fi = NULL;
}

FunctionInfo::FunctionInfo()
{
	_function_name    = Currency::vm.GetStringPointer("");
	_file_name        = Currency::pm.GetStringPointer("");
	_statements       = NULL;
	_anon_scope       = NULL;
	_persistent_scope = NULL;
	_builtin          = NULL;
	_alt_fptr         = NULL;
	_is_nested        = false;
	_is_constructor   = false;
	_is_encrypted     = false;
	local_functions   = NULL;
	nested_functions  = NULL;
	_refcnt           = 1;
	_parent_fi        = NULL;
}

FunctionInfo::~FunctionInfo()
{
	if (GetRefCount() == 0)
	{
		delete _statements;
		delete _persistent_scope;
		delete _anon_scope;

		std::map<const std::string*, FunctionInfo*>::const_iterator iter;

		if (local_functions)
		{
			for (iter = local_functions->begin(); iter != local_functions->end(); iter++)
			{
				FunctionInfo* fi = iter->second;

				fi->DecrRefCount();

				if (fi->GetRefCount() == 0)
					delete fi;
			}

			delete local_functions;
		}

		if (nested_functions)
		{
			for (iter = nested_functions->begin(); iter != nested_functions->end(); iter++)
			{
				FunctionInfo* fi = iter->second;

				fi->DecrRefCount();

				if (fi->GetRefCount() == 0)
					delete fi;
			}

			delete nested_functions;
		}
	}
}

FunctionInfo::FunctionInfo(const FunctionInfo& in)
    : _statements (nullptr)
{
	_function_name = in._function_name;
	_file_name = in._file_name;
	_builtin = in._builtin;
	_return_values = in._return_values;
	_parameters = in._parameters;
	_refcnt = 1;
	_parent_fi = in._parent_fi;

    if (in._statements)
    {
        _statements = new OMLTree(*in._statements);
    }

	if (in._persistent_scope)
		_persistent_scope = new MemoryScope(*in._persistent_scope);
	else
		_persistent_scope = NULL;

	if (in._anon_scope)
		_anon_scope = new MemoryScope(*in._anon_scope);
	else
		_anon_scope = NULL;

	_is_constructor = in._is_constructor;
	_is_nested      = in._is_nested;
	_is_encrypted   = in._is_encrypted;

	_help_string = in._help_string;

	local_functions = in.local_functions;
	nested_functions = in.nested_functions;

    _alt_fptr = in._alt_fptr;
}

bool FunctionInfo::IsReturnValue(const std::string* varname) const
{
	std::vector<const std::string*>::const_iterator iter;
	
	for (iter = _return_values.begin(); iter != _return_values.end(); iter++)
	{
		if (*iter == varname)
			return true;
	}

	return false;
}

bool FunctionInfo::IsInputParameter(const std::string* varname) const
{
	std::vector<const std::string*>::const_iterator iter;
	
	for (iter = _parameters.begin(); iter != _parameters.end(); iter++)
	{
		if (*iter == varname)
			return true;
	}

	return false;
}

bool CheckForReference(const std::string* varname, OMLTree* tree)
{
	int child_count = tree->ChildCount();

	for (int j=0; j<child_count; j++)
	{
		OMLTree*    inner_tree = tree->GetChild(j);

		if (inner_tree->GetType() == IDENT)
		{
			std::string ident_name = inner_tree->GetText();

			if (ident_name == *varname)
				return true;
		}
		else if (inner_tree->GetType() == FUNC_DEF)
		{
			continue; // ignore references in nested functions
		}

		if (inner_tree->ChildCount())
		{
			bool temp = CheckForReference(varname, inner_tree);

			if (temp)
				return true;
		}
	}

	return false;
}

bool FunctionInfo::IsReferenced(const std::string* varname) const
{
	return CheckForReference(varname, _statements);
}

int FunctionInfo::Nargin() const
{
	int nargin = (int)_parameters.size();
	if (nargin && (*_parameters.back() == "varargin"))
		nargin *= -1;
	return nargin;
}

int FunctionInfo::Nargout() const
{
	int nargout = (int)_return_values.size();
	if (nargout && (*_return_values.back() == "varargout"))
		nargout *= -1;
	return nargout;
}

bool FunctionInfo::IsLocalFunction(std::string script_name) const
{
	if (*_file_name == script_name)
		return false;

	size_t slash_index = _file_name->find_last_of("/\\");
	std::string root_file = _file_name->substr(slash_index+1);
	std::string just_file = root_file.substr(0, root_file.find('.'));

	if (just_file == *_function_name)
		return false;
	
	return true;
}

std::string FunctionInfo::RedirectedFunction() const
{
	OMLTree*  statements = _statements;
	int num_children = statements->ChildCount();

	if (num_children == 1)
	{
		OMLTree* child     = statements->GetChild(0);
		
		if (child->GetType() == FUNC)
		{
			OMLTree*    func     = child->GetChild(0);
			return func->GetText();
		}
	}

	return "";
}

int FunctionInfo::NumRedirectedInputs() const
{
	OMLTree*  statements = _statements;
	int num_children = statements->ChildCount();

	if (num_children == 1)
	{
		OMLTree*    child     = statements->GetChild(0);
		
		if (child->GetType() == FUNC)
		{
			OMLTree*    param_list     = child->GetChild(1);

			if (param_list->GetType() == PARAM_LIST)
				return param_list->ChildCount();
		}
	}


	return -1;
}

int FunctionInfo::MinimumInputs() const
{
	int ret = (int)_parameters.size();

	if (ret)
	{
		if (*_parameters[ret-1] == "varargin")
			ret--;
	}

	return ret;
}

OMLTree* FunctionInfo::RedirectedInput(int index) const
{
	OMLTree*  statements = _statements;
	int num_children = statements->ChildCount();

	if (num_children == 1)
	{
		OMLTree*    child     = statements->GetChild(0);
		
		if (child->GetType() == FUNC)
		{
			OMLTree*    param_list     = child->GetChild(1);

			if (param_list->GetType() == PARAM_LIST)
			{
				OMLTree* input = param_list->GetChild(index);
				return input;
			}
		}
	}

	return NULL;
}

std::string LocalGetAST(OMLTree* tree)
{
	if (!tree)
		return "";

	int child_count = tree->ChildCount();
	std::string result;

	if (tree->GetType() == DUMMY)
		return "";
	
	if (child_count)
		result += "(";

	if (tree->GetText().size())
	{
		if (tree->GetType() != QUOTE)
		{
			result += tree->GetText();
		}
		else
		{
			std::string temp = tree->GetText();
			
			size_t pos = temp.find('\'', 0);
			std::string replace = "''";

			while (pos != string::npos)
			{
				temp.replace(pos, 1, replace);
				pos = temp.find('\'', pos+2);
			}

			result += "'";
			result += temp;
			result += "'";
		}
	}
	else
	{
		result += "STATEMENT_LIST";
	}

	if (child_count)
		result += " ";

	for (int j=0; j<child_count; j++)
	{
		int index = j;

		if (tree->GetType() == FUNC_DEF)
		{
			if (j==1)
				index = 2;
			else if (j==2)
				index =1;
		}

		result += LocalGetAST(tree->GetChild(index));

		if (j != (child_count-1))
			result += " ";
	}

	if (child_count)
		result += ")";

	return result;
}

std::string FunctionInfo::GetAST() const
{
	std::string result = "(FUNC_DEF ";

	result += FunctionName();
	result += " (ID_LIST";

	for (int j=0; j<_parameters.size(); j++)
	{
		result += " ";
		result += *_parameters[j];
	}

	result += ") (ID_LIST";

	for (int j=0; j<_return_values.size(); j++)
	{
		result += " ";
		result += *_return_values[j];
	}

	result += ") ";

	if (_statements)
		result += LocalGetAST(_statements);

	result += ")";
	
	return result;
}

OMLTree* FunctionInfo::Statements() const 
{ 
	return _statements;
}

void FunctionInfo::SetAnonymous(MemoryScope* dummy)
{ 
	std::vector<std::string> idents;

	_statements->GetListOfIdents(idents);

	_anon_scope = new MemoryScope(NULL);

	// may want to ignore parameters here
	for (int j=0; j<idents.size(); j++)
	{
		const std::string* str_ptr = Currency::vm.GetStringPointer(idents[j]);

		if (dummy->Contains(str_ptr))
			_anon_scope->SetValue(str_ptr, dummy->GetValue(str_ptr));
		else if (dummy->IsGlobal(*str_ptr))
			_anon_scope->AddGlobalReference(*str_ptr);
	}
}

bool FunctionInfo::HasDefaultValue(const std::string* param_name) const
{
	return (_default_values.find(param_name) != _default_values.end());
}

Currency FunctionInfo::GetDefaultValue(const std::string* param_name) const
{
	static Currency _not_found(-1.0, Currency::TYPE_NOTHING);

	std::map<const std::string*, Currency>::const_iterator temp;
	temp = _default_values.find(param_name);

	if (temp != _default_values.end())
		return temp->second;
	else
		return _not_found;
}

void FunctionInfo::SetHelpString(const std::string& str)
{
	_help_string = str;
}