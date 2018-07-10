/**
* @file Evaluator.cpp
* @date August 2013
* Copyright (C) 2013-2018 Altair Engineering, Inc.  
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
#include "BoundClassInfo.h"
#include "BuiltInFuncsUtils.h"
#include "MemoryScope.h"
#include "Interpreter.h"
#include "FunctionInfo.h"
#include "SignalHandlerBase.h"
#include "ClassInfo.h"
#include "OML_Error.h"
#include "hwComplex.h"
#include "StructData.h"
#include "BuiltInFuncs.h"
#include "MatrixNUtils.h"
#include "ExprCppTreeLexer.h"
#include "ExprCppTreeParser.h"
#include "EvaluatorDebug.h"
#include "ANTLRData.h"
#include "OMLInterface.h"
#include "OMLTree.h"
#include <sys/stat.h>

#include <cassert>
#include <sstream>

#undef GetMessage // need to undefine this because Microsoft defines it somewhere in windows headers to GetMessageA

std::string ExprTreeEvaluator::lasterrormsg;
std::string ExprTreeEvaluator::lastwarning;

UserFunc::~UserFunc()  { if (fi) delete fi; }

#define RUN(tree) (this->*(tree->func_ptr))(tree)

inline const char* getText(pANTLR3_BASE_TREE tree)
{
	pANTLR3_COMMON_TOKEN tok = tree->getToken(tree);

	pANTLR3_STRING str = tok->getText(tok);

	if (!tok->tokText.text)
	{
		tok->textState = ANTLR3_TEXT_STRING;
		tok->tokText.text = str;
	}

	return (const char*) str->chars;
} 

//------------------------------------------------------------------------------
//! Constructor
//------------------------------------------------------------------------------
ExprTreeEvaluator::ExprTreeEvaluator() : nested_function_marker(0), assignment_nargout(0), is_for_evalin(false)
	, format(new OutputFormat())
    , _quit  (false), ext_var_ptr(NULL)
	, current_tree (NULL)
    , _signalHandler (NULL)
	, end_context_index (-1)  // Init so that end_context_index evaluates to correct value
	, end_context_varname(NULL)
    , _suspendFunclistUpdate (false)
{
	msm                     = new MemoryScopeManager;
	userFileStreams         = new std::vector<UserFile>;
    functions               = new std::map<std::string, UserFunc*>;
    std_functions           = new std::map<std::string, BuiltinFunc>;
	paths                   = new std::vector<std::string>;
	class_info_map          = new std::map<std::string, ClassInfo*>;
	not_found_functions     = new std::vector<std::string>;
	preregistered_functions = new std::vector<std::string>;

	mapBuiltInFuncs(std_functions);

	_owns_msm                 = true;
	_owns_userfiles           = true;
	_owns_functions           = true;
	_owns_pathnames           = true;
	_owns_format              = true;
	suppress_multi_ret_output = false;
	_lhs_eval                 = false;
	_interrupt                = false;
	_pause                    = false;
	_store_suppressed         = false;

	debug_listener       = NULL;
	end_context_currency = NULL;

	msm->OpenScope(NULL);
	
	AddFile(stdin,  "stdin",  "r");
	AddFile(stdout, "stdout", "w");
	AddFile(stderr, "stderr", "w");

	RestorePath();
}
//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
ExprTreeEvaluator::~ExprTreeEvaluator()
{
    // Don't delete signal handler as it is created in the client and will be
    // deleted there

	if (_owns_functions)
	{
		std::map<std::string, UserFunc*>::iterator iter;
		for (iter = functions->begin(); iter != functions->end(); iter++)
			delete (*iter).second;

		delete functions;
		delete not_found_functions;
		delete preregistered_functions;

		delete std_functions;
		delete class_info_map;
	}

	if (_owns_msm)
	{
		if (!is_for_evalin)
		msm->CloseScope();
		delete msm;
	}

	if (_owns_userfiles)
	{
		CloseAllFiles();
		delete userFileStreams;
	}

	if (_owns_pathnames)
	{
		delete paths;
	}

	if (_owns_format)
	{
		delete format;
	}
}
//------------------------------------------------------------------------------
//! Copy constructor
ExprTreeEvaluator::ExprTreeEvaluator(const ExprTreeEvaluator* source) : format(source->format), _owns_format(false),
	assignment_nargout(0), is_for_evalin(source->is_for_evalin)
    , _quit (false), ext_var_ptr(NULL)
    , current_tree (NULL)
    , _signalHandler (NULL)
	, end_context_varname(NULL)
    , end_context_index (-1) // Init so that end_context_index evaluates to correct value
    , _suspendFunclistUpdate (source->_suspendFunclistUpdate)
    , _appdir (source->_appdir)
{
	ImportUserFileList(source);
	ImportMemoryScope(source);
	ImportFunctionList(source);
	ImportPathNames(source);

	// We don't want debug listener chaining
	debug_listener       = NULL; 
	_interrupt           = false;
	_pause               = false;
    _paused              = false;
	_lhs_eval            = false;
	_store_suppressed    = false;

	nested_function_marker = source->nested_function_marker;

	suppress_multi_ret_output = source->suppress_multi_ret_output;

    if (source->_signalHandler)
        _signalHandler = source->_signalHandler->CreateClone();

	_decryptors = source->_decryptors;
}

Currency ExprTreeEvaluator::RunTree(OMLTree* tree)
{
	current_statement_index = 0;
    current_tree = tree;
	Currency ret = RUN(tree);
	current_tree = NULL;
	return ret;
}

Currency ExprTreeEvaluator::Nothing(OMLTree* tree)
{
	return Currency(-1.0, Currency::TYPE_NOTHING);
}

Currency ExprTreeEvaluator::Return(OMLTree* tree)
{
	return Currency(0.0, Currency::TYPE_RETURN);
}

Currency ExprTreeEvaluator::Number(OMLTree* tree)
{
	if (!tree->u)
	{
		std::string str = tree->GetText();
				
        // Linux and VS2015 will not process numbers with double scientific 
        // notation which use 'D' instead of 'E'
        std::replace(str.begin(), str.end(), 'D', 'E');
        std::replace(str.begin(), str.end(), 'd', 'e');
        const char* s = str.c_str();

        char* dummy = NULL;
		double val = strtod(s, &dummy);

		if ((*dummy == 'i') || (*dummy == 'j') || (*dummy == 'I') || (*dummy == 'J'))
			tree->u = new Currency(hwComplex(0.0, val));
		else
			tree->u = new Currency(val);
	}

	Currency* pCur = (Currency*)tree->u;
	return *pCur;
}

Currency ExprTreeEvaluator::HexNumber(OMLTree* tree)
{
	if (!tree->u)
	{
		std::string str = tree->GetText();
		const char* s = str.c_str();
		int val = 0;

		std::stringstream ss;
		ss << std::hex << s;
		ss >> val;

		tree->u = new Currency(val);
	}

	Currency* pCur = (Currency*)tree->u;
	return *pCur;
}

Currency ExprTreeEvaluator::TextString(OMLTree* tree)
{
	if (!tree->u)
	{
		OMLTree* child = NULL;
				
		if (tree->ChildCount())
			child = tree->GetChild(0);

		if (child)
			tree->u = new Currency(std::string(child->GetText()));
		else
			tree->u = new Currency("");
	}

	Currency* pCur = (Currency*)tree->u;
	return *pCur;
}

Currency ExprTreeEvaluator::Identifier(OMLTree* tree)
{
	if (!tree->u)
		tree->u = (void*)Currency::vm.GetStringPointer(tree->GetText());

	std::string* pString = (std::string*)tree->u;

	const Currency& temp_cur = msm->GetValue(pString);

	if (!temp_cur.IsNothing())
	{
		return temp_cur; 
	}
	else
	{
		std::vector<Currency> dummy;
		return CallFunction(pString, dummy);
	}
}

Currency ExprTreeEvaluator::CallFunction(const std::string& func_name, const std::vector<Currency>& params)
{
	const std::string* ptr = Currency::vm.GetStringPointer(func_name);
	return CallFunction(ptr, params);
}

Currency ExprTreeEvaluator::CallFunction(const std::string* func_name, const std::vector<Currency>& params)
{
	int ret_val;

	MemoryScope* scope = msm->GetCurrentScope();

	const std::string* current_filename = scope->GetFilenamePtr();
	std::string my_file;
	std::string my_extension;

	int stack_depth = msm->GetStackDepth();

	FunctionInfo* nested_fi = msm->GetNestedFunction(func_name);
	if (nested_fi)
		return CallInternalFunction(nested_fi, params);

	FunctionInfo* local_fi = msm->GetLocalFunction(func_name);
	if (local_fi)
		return CallInternalFunction(local_fi, params);

	if (functions->find(*func_name) != functions->end())
	{
		FunctionInfo* fi = (*functions)[*func_name]->fi;
		return CallInternalFunction(fi, params);
	}
	else if (CheckForFunctionInAST(*func_name))
	{
		FunctionInfo* fi = (*functions)[*func_name]->fi;
		return CallInternalFunction(fi, params);
	}
	else if (FindPrecompiledFunction(*func_name, my_file))
	{
		OMLTree* tree = OMLTree::ReadTreeFromFile(my_file);

		RUN(tree);

		delete tree;

		UserFunc*     uf = (*functions)[*func_name];
		FunctionInfo* fi = NULL;
		
		if (uf)
			fi = uf->fi;

		if (!fi)
		{
			if (func_name->length() < 2048)
			{
				char buffer[2048];
				sprintf(buffer, "Unknown function: %s", func_name->c_str());
				throw OML_Error(buffer);
			}
		}

		return CallInternalFunction(fi, params);
	}
	else if (FindFunction(*func_name, my_file) && (my_file != *current_filename))
	{
		bool ret_val = ParseAndRunFile(my_file, params.size() == 0);

		if (ret_val)
		{
			if (functions->find(*func_name) == functions->end())
			{
				if (func_name->length() < 2048)
				{
					char buffer[2048];
					sprintf(buffer, "Invalid function: %s", func_name->c_str());
					throw OML_Error(buffer);
				}		
			}
			
			UserFunc* uf = (*functions)[*func_name];
			FunctionInfo* fi = uf->fi;

			if (!fi)
			{
				if (func_name->length() < 2048)
				{
					char buffer[2048];
					sprintf(buffer, "Unknown function: %s", func_name->c_str());
					throw OML_Error(buffer);
				}
			}
			return CallInternalFunction(fi, params);
		}
		else // FindFunction already ran the script
		{
			return Currency(-1.0, Currency::TYPE_NOTHING);
		}
	}
	else if (FindEncryptedFunction(*func_name, my_file, my_extension))
	{
		RunEncryptedFile(my_extension, my_file);

		UserFunc* uf = NULL;

		if (functions->find(*func_name) != functions->end())
			uf = (*functions)[*func_name];

		if (uf)
		{
			FunctionInfo* fi = uf->fi;
			return CallInternalFunction(fi, params);
		}
		else
		{
			return Currency(-1.0, Currency::TYPE_NOTHING);
		}
	}
	else if (std_functions->find(*func_name) != std_functions->end())
	{
		FUNCPTR fptr = (*std_functions)[*func_name].fptr;

		if (fptr)
			return CallBuiltinFunction(fptr, *func_name, params);
		else
		{
			ALT_FUNCPTR alt_fptr = (*std_functions)[*func_name].alt_fptr;
			return CallBuiltinFunction(alt_fptr, *func_name, params);
		}
	}
	else
	{
		if (ext_var_ptr)
		{
			Currency temp = CheckForExternalVariable(*func_name);

			if (!temp.IsNothing())
				return temp;
		}

		if (func_name->length() < 2048)
		{
			if (my_file != *current_filename)
			{
				char buffer[2048];
				sprintf(buffer, "Unknown function: %s", func_name->c_str());
				throw OML_Error(buffer);
			}
			else
			{
				char buffer[2048];
				sprintf(buffer, "Illegal recursive call detected with script: %s", func_name->c_str());
				throw OML_Error(buffer);
			}
		}
		else
		{
			throw OML_Error(HW_ERROR_VARIABLENAMETOOLONG);
		}
	}
	return 0.0;
}

Currency ExprTreeEvaluator::CallBuiltinFunction(FUNCPTR fptr, const std::string& func_name, const std::vector<Currency>& params)
{
	// Check here for overloaded functions if the first param is an object
	// Don't call the overload if there isn't one.  Most functions will fail either way,
	// but a few built-ins (e.g. class) will work with objects and in fact we need them to.
	if (params.size() && params[0].IsObject())
	{
		if (HasOverloadedFunction(params[0], func_name))
			return CallOverloadedFunction(func_name, params);
	}

	if ((func_name != "nargout") && (func_name != "nargin"))
		PushNargValues((int)params.size(), assignment_nargout);

	std::string old_builtin_scope = builtin_error_scope;

	if ((func_name != "warning") && (func_name != "error"))
		builtin_error_scope = func_name;

	std::vector<Currency> ret;
	fptr(EvaluatorInterface(this), params, ret);
	builtin_error_scope = old_builtin_scope;

	if ((func_name != "nargout") && (func_name != "nargin"))
		PopNargValues();

	if (ret.size() > 0)
		return ret[0];
	else
		return Currency(-1.0, Currency::TYPE_NOTHING);
}

Currency ExprTreeEvaluator::CallBuiltinFunction(ALT_FUNCPTR alt_fptr, const std::string& func_name, const std::vector<Currency>& params)
{
	if ((func_name != "nargout") && (func_name != "nargin"))
		PushNargValues((int)params.size(), assignment_nargout);

	std::string old_builtin_scope = builtin_error_scope;

	OMLCurrencyListImpl in_list;
	OMLCurrencyListImpl out_list;

	for (int j=0; j<params.size(); j++)
	{
		if (params[j].IsScalar())
			in_list.AddScalar(params[j].Scalar());
		else if (params[j].IsString())
			in_list.AddString(params[j].StringVal().c_str());
		else if (params[j].IsCellArray())
			in_list.AddCellArray(params[j].CellArray());
		else if (params[j].IsMatrix())
			in_list.AddMatrix(params[j].Matrix());
		else if (params[j].IsNDMatrix())
			in_list.AddNDMatrix(params[j].MatrixN());
		else if (params[j].IsStruct())
			in_list.AddStruct(params[j].Struct());
		else if (params[j].IsFunctionHandle())
			in_list.AddFunctionHandle(params[j].FunctionHandle());
	}

	EvaluatorInterface ei(this);
	OMLInterfaceImpl impl(&ei);

	alt_fptr(&impl, &in_list, &out_list);

	builtin_error_scope = old_builtin_scope;

	if ((func_name != "nargout") && (func_name != "nargin"))
		PopNargValues();

	if (out_list.Size() > 0)
	{
		const OMLCurrencyImpl* temp = (OMLCurrencyImpl*)out_list.Get(0);
		Currency cur = temp->GetCurrency();

		OMLCurrencyImpl::GarbageCollect();
		OMLComplexImpl::GarbageCollect();
		OMLMatrixImpl::GarbageCollect();
		OMLNDMatrixImpl::GarbageCollect();
		OMLCellArrayImpl::GarbageCollect();
		OMLStructImpl::GarbageCollect();
		OMLFunctionHandleImpl::GarbageCollect();

		return cur;
	}
	else
	{
		return Currency(-1.0, Currency::TYPE_NOTHING);
	}
}

void ExprTreeEvaluator::CallFunction(const std::string& func_name, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int nargout)
{
	FunctionInfo* fi = NULL;

	const Currency& temp_cur = msm->GetValue(func_name);

	if (!temp_cur.IsNothing())
	{
		if (temp_cur.IsFunctionHandle())
			fi = temp_cur.FunctionHandle();
	}

	if (!fi)
		fi = (*functions)[func_name]->fi;

	std::vector<Currency> loc_inputs = inputs;

	std::vector<const std::string*> return_values = fi->ReturnValues();

	int my_nargout = (int)return_values.size();

	if (nargout != -1)
		my_nargout = nargout;

	if (fi)
		outputs = DoMultiReturnFunctionCall(fi, loc_inputs, (int)loc_inputs.size(), my_nargout, true, &return_values);
}

Currency ExprTreeEvaluator::VariableIndex(const std::string* var_name, const std::vector<Currency>& params)
{
	Currency target = msm->GetValue(var_name);

	return VariableIndex(target, params);
}

Currency ExprTreeEvaluator::VariableIndex(const Currency& target, const std::vector<Currency>& params)
{
	if (target.IsMatrix() && (params.size() == 1))
	{
		int index;

		const hwMatrix* data = target.Matrix();

		if (params[0].IsPositiveInteger())
		{
			index = static_cast<int>(params[0].Scalar());
		}
		else if (params[0].IsPositiveVector())
		{
			if (params[0].GetMask() == Currency::MASK_DOUBLE)
			{
				const hwMatrix* mtx = params[0].Matrix(); // this must be real
				hwMatrix* ret = allocateMatrix(mtx->M(), mtx->N(), data->Type());

				if (data->IsVector() && (data->N() == 1) && (mtx->M() == 1))
					ret->Transpose();

				int max_index = data->Size();

				for (int j=0; j<(int)mtx->Size(); j++)
				{
					int mtx_index = (int)((*mtx)(j)-1);

					if ((mtx_index < 0) || (mtx_index >= max_index))
					{
						delete ret;
						throw OML_Error(HW_ERROR_INDEXRANGE);
					}

					if (data->IsReal())
						(*ret)(j) = (*data)((int)((*mtx)(j)-1));
					else
						ret->z(j) = data->z((int)((*mtx)(j)-1));
				}
				
				return ret;
			}		
		}
		else if (params[0].IsVector() && params[0].IsLogical())
		{
			const hwMatrix* mtx = params[0].Matrix(); // this must be real
			int mtx_size = mtx->Size();

			if (data->IsReal())
			{
				std::vector<double> ret;
				ret.reserve(mtx_size);

				for (int j=0; j<mtx_size; j++)
				{
					if ((*mtx)(j))
						ret.push_back((*data)(j));
				}

				int new_rows = 1;
				int new_cols = 1;

				if ((data->M() != 1) && (data->N() == 1))
					new_rows = (int)ret.size();
				else
					new_cols = (int)ret.size();

				hwMatrix* ret_mat = allocateMatrix(new_rows, new_cols, hwMatrix::REAL);

				for (int j=0; j<ret.size(); j++)
					(*ret_mat)(j) = ret[j];

				return ret_mat;
			}
			else
			{
				std::vector<hwComplex> ret;
				ret.reserve(mtx_size);

				for (int j=0; j<mtx_size; j++)
				{
					if ((*mtx)(j))
						ret.push_back(data->z(j));
				}

				hwMatrix* ret_mat = allocateMatrix((int)ret.size(), 1, hwMatrix::REAL);

				for (int j=0; j<ret.size(); j++)
				{
					if (ret[j].IsReal() && ret_mat->IsReal())
					{
						(*ret_mat)(j) = ret[j].Real();
					}
					else
					{
						ret_mat->MakeComplex();
						ret_mat->z(j) = ret[j];
					}
				}

				return ret_mat;
			}
		}
		else if (params[0].IsMatrix())
		{
			const hwMatrix* val = params[0].Matrix();
			if (val->IsEmpty())
			{
				return params[0];
			}
			else if (params[0].IsLogical())
			{
				if ((data->M() == val->M()) && (data->N() == val->N()))
				{
					std::vector<double> temp_vec;
					temp_vec.reserve(data->Size());

					for (int j=0; j<data->Size(); j++)
					{
						if ((*val)(j) == 1.0)
							temp_vec.push_back((*data)(j));
					}

					return temp_vec;
				}
				else
				{
					throw OML_Error(HW_ERROR_MATRIXDIM);
				}
			}
			else
			{
				if (params[0].IsVector())
					throw OML_Error(HW_ERROR_INVIND);

				throw OML_Error(HW_ERROR_INVINDTYPE);
			}
		}
		else if (params[0].IsColon())
		{
			hwMatrix* ret = NULL;
			
			if (data)
			{
				ret = allocateMatrix(data->Size(), 1, data->Type());

				for (int j=0; j<data->Size(); j++)
				{
					if (data->IsReal())
						(*ret)(j) = (*data)(j);
					else
						ret->z(j) = data->z(j);
				}
			}

			return ret;
		}
        else if (params[0].IsNDMatrix())
        {
			hwMatrixN* ret = allocateMatrixN();
			const hwMatrixN* idxm = params[0].MatrixN();

            if (!idxm->IsReal())
    			throw OML_Error(HW_ERROR_INDEXPOSINT);

            if (data->IsReal())
                ret->Dimension(idxm->Dimensions(), hwMatrixN::REAL);
            else
                ret->Dimension(idxm->Dimensions(), hwMatrixN::COMPLEX);

            int size_idx = idxm->Size();
            int size_tgt = data->Size();

            for (int i = 0; i < size_idx; ++i)
            {
                double idxd = (*idxm)(i);
                int idx = static_cast<int>(idxd) - 1;

                if (idx < 0 || idx >= size_tgt)
                    throw OML_Error(OML_ERR_INVALID_RANGE);

                if (ret->IsReal())
                    (*ret)(i) = (*data)(idx);
                else
                    ret->z(i) = data->z(idx);
            }

			return ret;
        }
        else
		{
			throw OML_Error(HW_ERROR_INDEXPOSINT);
		}

		if (!data || (index > data->Size()) || (index <= 0))
			throw OML_Error(HW_ERROR_INDEXRANGE);

		if (data->IsReal())
		{
			return (*data)(index-1);
		}
		else
		{
			const hwComplex temp = data->z(index-1);
			return temp;
		}
	}
	else if ((target.IsMatrixOrString() && (params.size() == 2)))
	{
		const hwMatrix* data = target.Matrix();

		if (params[0].IsScalar() && params[1].IsScalar())
		{
			Currency c1 = params[0];
			Currency c2 = params[1];

			int index1;
			int index2;

			if (c1.IsPositiveInteger())
				index1 = static_cast<int>(c1.Scalar());
			else
				throw OML_Error(HW_ERROR_INDEXPOSINT);

			if (c2.IsPositiveInteger())
				index2 = static_cast<int>(c2.Scalar());
			else
				throw OML_Error(HW_ERROR_INDEXPOSINT);

			if ((index1 > data->M()) || (index1 <= 0))
				throw OML_Error(HW_ERROR_INDEXRANGE);

			if ((index2 > data->N()) || (index2 <= 0))
				throw OML_Error(HW_ERROR_INDEXRANGE);

			if (data->IsReal())
			{
				Currency ret_val = (*data)(index1-1, index2-1);

				if (target.IsLogical())
					ret_val.SetMask(Currency::MASK_LOGICAL);

				if (target.IsString() && ret_val.IsScalar())
				{
					double   val = ret_val.Scalar();

					Currency new_ret(EvaluatorInterface::allocateMatrix(1, 1, val));
					new_ret.SetMask(Currency::MASK_STRING);

					return new_ret;
				}

				return ret_val;
			}
			else
			{
				return data->z(index1-1, index2-1);
			}
		}
		else if (params[0].IsScalar() && params[1].IsColon())  // get the ith row
		{
			int       idx;

			if (params[0].IsPositiveInteger())
				idx = static_cast<int>(params[0].Scalar());
			else
				throw OML_Error(HW_ERROR_INDEXPOSINT);


			hwMatrix* ret = allocateMatrix(1, data->N(), data->Type());
			ValidateRowIndex(*data, idx-1);
			data->ReadRow(idx-1, *ret);
			Currency ret_val(ret);

			if (target.IsString())
				ret_val.SetMask(Currency::MASK_STRING);

			return ret_val;
		}
		else if (params[0].IsColon() && params[1].IsScalar())  // get the ith column
		{
			int idx;

			if (params[1].IsPositiveInteger())
				idx = static_cast<int>(params[1].Scalar());
			else
				throw OML_Error(HW_ERROR_INDEXPOSINT);


			hwMatrix* ret = allocateMatrix(data->M(),1,data->Type());
			ValidateColumnIndex(*data, idx-1);
			data->ReadColumn(idx-1, *ret);

			Currency ret_val(ret);

			if (target.IsString())
				ret_val.SetMask(Currency::MASK_STRING);

			return ret_val;
		}
		else if (params[0].IsVector() || params[1].IsVector())
		{
			std::vector<double> first_range;
			std::vector<double> second_range;

			if (params[0].IsPositiveVector())
			{
				// remove this copy when practical
				first_range = params[0].Vector();
			}
			else if (params[0].IsPositiveInteger())
			{
				first_range.push_back(params[0].Scalar());
			}
			else if (params[0].IsVector() && params[0].IsLogical())
			{
				std::vector<double> logicals = params[0].Vector();

				for (size_t j=0; j<logicals.size(); j++)
				{
					if (logicals[j])
						first_range.push_back((double)j+1);
				}
			}
			else if (params[0].IsColon())
			{
				for (size_t j=1; j<=data->M(); j++)
					first_range.push_back((double)j);
			}
			else if (params[0].IsEmpty()) 
			{
			}
			else
			{
				throw OML_Error(HW_ERROR_INV1STIND);
			}

			if (params[1].IsPositiveVector())
			{
				// remove this copy when practical
				second_range = params[1].Vector();
			}
			else if (params[1].IsVector() && params[1].IsLogical())
			{
				std::vector<double> logicals = params[1].Vector();

				for (size_t j=0; j<logicals.size(); j++)
				{
					if (logicals[j])
						second_range.push_back((double)j+1);
				}
			}
			else if (params[1].IsPositiveInteger())
			{
				second_range.push_back(params[1].Scalar());
			}
			else if (params[1].IsColon())
			{
				for (size_t j=1; j<=data->N(); j++)
					second_range.push_back((double)j);
			}
			else if (params[1].IsEmpty())
			{
			}
			else
			{
				throw OML_Error(HW_ERROR_INV2NDIND);
			}

			hwMatrix* ret = allocateMatrix((int)first_range.size(), (int)second_range.size(), data->Type());

			int r1 = (int)first_range.size();
			int r2 = (int)second_range.size();

			int data_M = data->M();
			int data_N = data->N();

			for (int k=0; k<r2; k++)
			{
				int index_2 = (int)second_range[k]-1;

				if ((index_2 < 0) || (index_2 >= data_N))
					throw OML_Error(HW_ERROR_INDEXRANGE);

			    for (int j=0; j<r1; j++)
			    {
				    int index_1 = (int)first_range[j]-1;

			    	if ((index_1 < 0) || (index_1 >= data_M))
				    	throw OML_Error(HW_ERROR_INDEXRANGE);

					if (ret->IsReal())
						(*ret)(j,k) = (*data)(index_1, index_2);
					else
						ret->z(j,k) = data->z(index_1, index_2);
				}
			}
            Currency out (ret);
            out.SetMask(target.GetMask());
			return out;
		}
		else if (params[0].IsEmpty() || params[1].IsEmpty())
		{
			int new_rows = 0;
			int new_cols = 0;

			if (params[0].IsEmpty())
			{
				if (params[1].IsPositiveInteger())
					new_cols = 1;
				else if (params[1].IsColon())
					new_cols = data->N();
			}
			else
			{
				if (params[0].IsPositiveInteger())
					new_rows = 1;
				else if (params[0].IsColon())
					new_rows = data->M();
			}

			hwMatrix* ret = allocateMatrix(new_rows, new_cols, data->Type());
			Currency ret_cur(ret);

			if (target.IsString())
				ret_cur.SetMask(Currency::MASK_STRING);

			return ret_cur;
		}
		else if (params[0].IsColon() && params[1].IsColon())
		{
			return target;
		}
		else
		{
			throw OML_Error(HW_ERROR_UNSUPINDOP);
		}
	}
	else if (target.IsNDMatrix() || (target.IsMatrix() && params.size() > 2))
	{
		if (params.size() == 1)
		{
			int index;

			if (params[0].IsPositiveInteger())
			{
				index = (int)params[0].Scalar() - 1;
				const hwMatrixN* mat_n = target.MatrixN();

				if ((index >= mat_n->Size()) || (index < 0))
					throw OML_Error(HW_ERROR_INDEXRANGE);

				if (mat_n->IsReal())
					return (*mat_n)(index);
				else
					return mat_n->z(index);
			}
			else if (params[0].IsColon())
			{
				const hwMatrixN* mat_n = target.MatrixN();
				std::vector<hwSliceArg> slices;
				slices.push_back(hwSliceArg());
				hwMatrixN* temp_slice = allocateMatrixN();
				mat_n->SliceRHS(slices, *temp_slice);
                Currency ret = ConvertNDto2D(temp_slice);
                ret.SetMask(target.GetMask());
                return ret;
			}
			else if (params[0].IsMatrix())
			{
				const hwMatrixN* mat_n = target.MatrixN();  // values
				const hwMatrix* idxm = params[0].Matrix();  // indices

                if (!idxm->IsReal())
        			throw OML_Error(HW_ERROR_INDEXPOSINT);

    			hwMatrix* retm = nullptr;

                if (mat_n->IsReal())
    			    retm = allocateMatrix(idxm->M(), idxm->N(), hwMatrix::REAL);
                else
    			    retm = allocateMatrix(idxm->M(), idxm->N(), hwMatrix::COMPLEX);

                int size_idx = idxm->Size();
                int size_tgt = mat_n->Size();

                for (int i = 0; i < size_idx; ++i)
                {
                    double idxd = (*idxm)(i);
                    int idx = static_cast<int>(idxd) - 1;

                    if (idx < 0 || idx >= size_tgt)
                        throw OML_Error(OML_ERR_INVALID_RANGE);

                    if (mat_n->IsReal())
                        (*retm)(i) = (*mat_n)(idx);
                    else
                        retm->z(i) = mat_n->z(idx);
                }

                Currency ret(retm);
                ret.SetMask(target.GetMask());
                return ret;
			}
			else if (params[0].IsNDMatrix())
			{
				const hwMatrixN* mat_n = target.MatrixN();    // values
				const hwMatrixN* idxm = params[0].MatrixN();  // indices

                if (!idxm->IsReal())
        			throw OML_Error(HW_ERROR_INDEXPOSINT);

    			hwMatrixN* retm = allocateMatrixN();

   			    retm->Dimension(idxm->Dimensions(), mat_n->Type());

                int size_idx = idxm->Size();
                int size_tgt = mat_n->Size();

                for (int i = 0; i < size_idx; ++i)
                {
                    double idxd = (*idxm)(i);
                    int idx = static_cast<int>(idxd) - 1;

                    if (idx < 0 || idx >= size_tgt)
                        throw OML_Error(OML_ERR_INVALID_RANGE);

                    if (mat_n->IsReal())
                        (*retm)(i) = (*mat_n)(idx);
                    else
                        retm->z(i) = mat_n->z(idx);
                }

                Currency ret(retm);
                ret.SetMask(target.GetMask());
                return ret;
			}
		}
		else
		{
			indices.clear();
			bool all_scalars = true;

			for (int j=0; j<params.size(); j++)
			{
				if (params[j].IsPositiveInteger())
				{
					int index = (int)params[j].Scalar() - 1;
					indices.push_back(index);
				}
				else
				{
					all_scalars = false;
					break;
				}
			}
				
			const hwMatrixN* mat_n = NULL;
			
			if (target.IsNDMatrix())
			{
				mat_n = target.MatrixN();
			}
			else if (target.IsMatrix())
			{
				const hwMatrix* mat = target.Matrix();
				mat_n = Convert2DtoND(mat);
			}

			if (all_scalars)
			{
				try
				{
					mat_n->BoundsCheckRHS(indices);
				}
				catch (hwMathException&)
				{
					throw OML_Error(HW_ERROR_INDEXRANGE);
				}

				return (*mat_n)(indices);
			}
			else
			{
				slices.clear();

				for (int j=0; j<params.size(); j++)
				{
					if (params[j].IsScalar())
					{
						slices.push_back((int)params[j].Scalar()-1);
					}
					else if (params[j].IsColon())
					{
						slices.push_back(hwSliceArg());
					}
					else if (params[j].IsPositiveVector())
					{
						std::vector<double> temp = params[j].Vector();
						std::vector<int>    temp_int;

						for (int k=0; k<temp.size(); k++)
							temp_int.push_back((int)(temp[k]-1));
						slices.push_back(hwSliceArg(temp_int));
					}
					else
					{
						throw OML_Error("Unsupported slice type");
					}
				}

				hwMatrixN* temp_slice = allocateMatrixN();
				mat_n->SliceRHS(slices, *temp_slice);

				std::vector<int> dims = temp_slice->Dimensions();

				if (dims.size() == 2)
				{
					hwMatrix* degenerate = ConvertNDto2D(temp_slice);
					delete temp_slice;
                    Currency ret(degenerate);
                    ret.SetMask(target.GetMask());
                    return ret;
				}
				else
				{
                    Currency ret(temp_slice);
                    ret.SetMask(target.GetMask());
                    return ret;
				}
			}
		}
	}
	else if ((target.IsString() && (params.size() == 1)))
	{
		const hwMatrix* mtx = target.Matrix();
		std::string     result;

		if (params[0].IsPositiveInteger())
		{
			int index = (int)(params[0].Scalar()-1);

			if (index >= mtx->Size())
				throw OML_Error(HW_ERROR_INDEXRANGE);

			char character = char((*mtx)(index));

			result += character;
		}
		else if (params[0].IsPositiveVector())
		{
			std::vector<double> indices = params[0].Vector();

			for (unsigned int j=0; j<indices.size(); j++)
			{
				if (indices[j] > mtx->Size())
					throw OML_Error(HW_ERROR_INDEXRANGE);

				char character = char((*mtx)((int)indices[j]-1));

				result += character;
			}
		}
		else if (params[0].IsLogical() && params[0].IsVector())
		{
			std::vector<double> indices = params[0].Vector();

			for (unsigned int j=0; j<indices.size(); j++)
			{
				if (indices[j])
				{
					char character = char((*mtx)(j));

					result += character;
				}
			}
		}
		else if (params[0].IsScalar() || (params[0].IsMatrix() && !params[0].IsEmpty()))
		{
			// i.e.  it's got a negative
			throw OML_Error(HW_ERROR_INDEXRANGE);
		}

		return result;
	}
	else if (target.IsFunctionHandle())
	{
		FunctionInfo* safe_temp = new FunctionInfo(*target.FunctionHandle());

		std::vector<Currency> expanded_params;

		for (int j=0; j<params.size(); j++)
		{
			Currency param = params[j];

			if (param.IsCellList())
			{
				HML_CELLARRAY* cells = param.CellArray();

				for (int k=0; k<cells->Size(); k++)
					expanded_params.push_back((*cells)(k));
			}
			else
			{
				expanded_params.push_back(param);

				if (expanded_params.back().IsNothing())
					expanded_params.pop_back();
			}
		}

		Currency ret = CallInternalFunction(safe_temp, expanded_params);
		delete safe_temp;
		return ret;
	}
	else if (target.IsScalar())
	{
		if (params.size() == 1)
		{
			if (params[0].IsScalar())
			{
				if (params[0].Scalar() == 1)
				{
					return target.Scalar();
				}
				else if (params[0].IsLogical())
				{
					if (params[0].Scalar() == 0)
						return Currency(); // i.e. empty matrix
				}
				else
				{
					throw OML_Error(HW_ERROR_INDEXRANGE);
				}
			}
			else if (params[0].IsColon())
			{
				return target.Scalar();
			}
			else if (params[0].IsEmpty())
			{
				return params[0];
			}
			else
			{
				throw OML_Error(HW_ERROR_NOTIMP);
			}
		}
		else if (params.size() == 2)
		{
			if (params[0].IsScalar() && params[1].IsColon())
			{
				if (params[0].Scalar() == 1)
					return target.Scalar();
			}
			else if (params[1].IsScalar() && params[0].IsColon())
			{
				if (params[1].Scalar() == 1)
					return target.Scalar();			
			}
			else if (params[0].IsScalar() && params[1].IsScalar())
			{
				if ((params[0].Scalar() == 1) && (params[1].Scalar() == 1))
					return target.Scalar();	
			}
			else if (params[0].IsEmpty() || params[1].IsEmpty())
			{
				int rows = 0;
				int cols = 0;

				if (params[0].IsScalar())
				{
					if (params[0].Scalar() == 1.0)
						rows = 1;
					else
						throw OML_Error(HW_ERROR_INDEXRANGE);
				}
				else if (params[0].IsColon())
				{
					rows = 1;
				}

				if (params[1].IsScalar())
				{
					if (params[1].Scalar() == 1.0)
						cols = 1;
					else
						throw OML_Error(HW_ERROR_INDEXRANGE);
				}
				else if (params[1].IsColon())
				{
					cols = 1;
				}

				return allocateMatrix(rows, cols, 0.0);
			}

			throw OML_Error(HW_ERROR_INDEXRANGE);
		}
		else
		{
			size_t num_params = params.size();

			for (size_t j=0; j<num_params; j++)
			{
				if (!params[j].IsScalar())
					throw OML_Error(HW_ERROR_INDEXRANGE);

				if (params[j].Scalar() != 1)
					throw OML_Error(HW_ERROR_INDEXRANGE);
			}

			return target.Scalar();
		}
	}
	else if (target.IsComplex())
	{
		if (params.size() == 1)
		{
			if (params[0].IsScalar())
			{
				if (params[0].Scalar() == 1)
				{
					return target.Complex();
				}
				else if (params[0].IsLogical())
				{
					if (params[0].Scalar() == 0)
						return Currency(); // i.e. empty matrix
				}
				else
				{
					throw OML_Error(HW_ERROR_INDEXRANGE);
				}
			}
			else if (params[0].IsColon())
			{
				return target.Complex();
			}	
		}
		else
		{
			throw OML_Error(HW_ERROR_NOTIMP);
		}
	}
	else if (target.IsStruct())
	{
		if (params.size() == 1)
		{
			if (params[0].IsScalar())
			{
				StructData* sd = target.Struct();
				return sd->GetElement((int)params[0].Scalar(), -1);
			}
			else if (params[0].IsColon())
			{
				StructData* sd = target.Struct();
				StructData* result = new StructData(*sd);
				result->Reshape(sd->M()*sd->N(), 1);
				return result;
			}
			else if (params[0].IsPositiveVector())
			{
				StructData* sd     = target.Struct();	
				StructData* new_sd = new StructData;

				std::vector<double> indices = params[0].Vector();
				new_sd->DimensionNew((int)indices.size(), 1);

				for (int j=0; j<indices.size(); j++)
				{
					int index = (int)indices[j];
					StructData* temp = sd->GetElement(index, -1);
					new_sd->SetElement(j+1, -1, temp);
				}

				return new_sd;
			}
			else
			{
				throw OML_Error(OML_ERR_POS_INTEGER_VEC_MTX);
			}
		}	
		else if (params.size() == 2)
		{
			if (params[0].IsScalar() && params[1].IsScalar())
			{
				StructData* sd = target.Struct();
				return sd->GetElement((int)params[0].Scalar(), (int)params[1].Scalar());
			}
			else if (params[0].IsPositiveInteger() && params[1].IsColon())
			{
				StructData* sd     = target.Struct();	
				StructData* new_sd = new StructData;

				int index = (int)params[0].Scalar();
				
				new_sd->DimensionNew(1, sd->N());
				
				for (int j=0; j<sd->N(); j++)
				{
					StructData* temp = sd->GetElement(index, j+1);
					new_sd->SetElement(j+1, -1, temp);
				}

				return new_sd;
			}
			else if (params[0].IsColon() && params[1].IsPositiveInteger())
			{
				StructData* sd     = target.Struct();	
				StructData* new_sd = new StructData;

				int index = (int)params[1].Scalar();
				
				new_sd->DimensionNew(sd->M(), 1);
				
				for (int j=0; j<sd->M(); j++)
				{
					StructData* temp = sd->GetElement(j+1, index);
					new_sd->SetElement(j+1, -1, temp);
				}

				return new_sd;
			}
			else
			{
				throw OML_Error(OML_MSG_INVALID_INDEX);
			}
		}
		else if (params.size() == 0)
		{
			return target;
		}
	}
	else if (target.IsObject())
	{
		if (params.size() == 1)
		{
			if (params[0].IsScalar())
			{
				StructData* sd = target.Struct();
				Currency cur = sd->GetElement((int)params[0].Scalar(), -1);
				cur.SetClass(target.GetClassname());
				return cur;
			}
			else if (params[0].IsColon())
			{
				throw OML_Error("Unsupported operation");
			}
		}	
		else if ((params.size() == 2)  && (params[0].IsScalar()) && (params[1].IsScalar()))
		{
			StructData* sd = target.Struct();
			Currency cur = sd->GetElement((int)params[0].Scalar(), (int)params[1].Scalar());
			cur.SetClass(target.GetClassname());
			return cur;
		}	
	}
	else if (target.IsCellArray())
	{
		if (params.size() == 0)
		{
			return target;
		}

		if (target.IsCellList())
			throw OML_Error(HW_ERROR_INVCELLIND);

		if (params.size() == 1)
		{
			if (params[0].IsPositiveInteger())
			{
				HML_CELLARRAY* cells = target.CellArray();

				int index = (int)params[0].Scalar()-1;

				if ((index < cells->Size()) && (index >= 0))
				{
					HML_CELLARRAY* ret = allocateCellArray(1, 1);
					(*ret)(0) = (*cells)(index);
					return ret;
				}
				else
				{
					throw OML_Error(HW_ERROR_INVIND);
				}
			}
			else if (params[0].IsPositiveVector())
			{
				std::vector<double> indices = params[0].Vector();

				HML_CELLARRAY* cells = target.CellArray();
				HML_CELLARRAY* ret = NULL;
				
				if (cells->N() == 1)
					ret = allocateCellArray((int)indices.size(), 1);
				else
					ret = allocateCellArray(1, (int)indices.size());

				for (int j=0; j<indices.size(); j++)
				{
					int index = (int)(indices[j]-1);

					if (index < cells->Size())
						(*ret)(j) = (*cells)(index);
					else
						throw OML_Error(HW_ERROR_INDEXRANGE);
				}

				return ret;
			}
			else if (params[0].IsEmpty())
			{
				const hwMatrix* mtx = params[0].Matrix();
				return allocateCellArray(mtx->M(), mtx->N());
			}
			else if (params[0].IsColon())
			{
				HML_CELLARRAY* cells = target.CellArray();
				HML_CELLARRAY* ret = allocateCellArray((int)cells->Size(), 1);
				
				for (int j=0; j<cells->Size(); j++)
					(*ret)(j) = (*cells)(j);

				return ret;
			}
			else
			{
				throw OML_Error(HW_ERROR_INVIND);
			}
		}
		else if (params.size() == 2)
		{
			if ((params[0].IsPositiveInteger()) && (params[1].IsPositiveInteger()))
			{
				HML_CELLARRAY* cells = target.CellArray();

				int index1 = (int)params[0].Scalar()-1;
				int index2 = (int)params[1].Scalar()-1;

				if ((index1 < cells->M()) && (index2 < cells->N()))
				{
					HML_CELLARRAY* ret = allocateCellArray(1, 1);
					(*ret)(0) = (*cells)(index1, index2);
					return ret;
				}
				else
				{
					throw OML_Error(HW_ERROR_INVIND);
				}
			}
			else if ((params[0].IsPositiveInteger()) && (params[1].IsColon()))
			{
				HML_CELLARRAY* cells = target.CellArray();

				int size = cells->N();

				HML_CELLARRAY* ret = allocateCellArray(1, size);

				int index1 = (int)params[0].Scalar()-1;

				for (int j=0; j<size; j++)
					(*ret)(j) = (*cells)(index1, j);

				return ret;
			}
			else if ((params[0].IsColon()) && (params[1].IsPositiveInteger()))
			{
				HML_CELLARRAY* cells = target.CellArray();

				int size = cells->M();

				HML_CELLARRAY* ret = allocateCellArray(size, 1);

				int index1 = (int)params[1].Scalar()-1;

				for (int j=0; j<size; j++)
					(*ret)(j) = (*cells)(j, index1);

				return ret;
			}
			else if (params[0].IsVector() || params[1].IsVector())
			{
				HML_CELLARRAY* cells = target.CellArray();

				std::vector<double> first_range;
				std::vector<double> second_range;

				if (params[0].IsPositiveVector())
				{
					// remove this copy when practical
					first_range = params[0].Vector();
				}
				else if (params[0].IsPositiveInteger())
				{
					first_range.push_back(params[0].Scalar());
				}
				else if (params[0].IsVector() && params[0].IsLogical())
				{
					std::vector<double> logicals = params[0].Vector();

					for (size_t j=0; j<logicals.size(); j++)
					{
						if (logicals[j])
							first_range.push_back((double)j+1);
					}
				}
				else if (params[0].IsColon())
				{
					for (size_t j=1; j<=cells->M(); j++)
						first_range.push_back((double)j);
				}
				else if (params[0].IsEmpty()) 
				{
				}
				else
				{
					throw OML_Error(HW_ERROR_INV1STIND);
				}

				if (params[1].IsPositiveVector())
				{
					// remove this copy when practical
					second_range = params[1].Vector();
				}
				else if (params[1].IsVector() && params[1].IsLogical())
				{
					std::vector<double> logicals = params[1].Vector();

					for (size_t j=0; j<logicals.size(); j++)
					{
						if (logicals[j])
							second_range.push_back((double)j+1);
					}
				}
				else if (params[1].IsPositiveInteger())
				{
					second_range.push_back(params[1].Scalar());
				}
				else if (params[1].IsColon())
				{
					for (size_t j=1; j<=cells->N(); j++)
						second_range.push_back((double)j);
				}
				else if (params[1].IsEmpty())
				{
				}
				else
				{
					throw OML_Error(HW_ERROR_INV2NDIND);
				}

				HML_CELLARRAY* ret = allocateCellArray((int)first_range.size(), (int)second_range.size());

				int r1 = (int)first_range.size();
				int r2 = (int)second_range.size();

				int data_M = cells->M();
				int data_N = cells->N();

				for (int k=0; k<r2; k++)
				{
					int index_2 = (int)second_range[k]-1;

					if ((index_2 < 0) || (index_2 >= data_N))
						throw OML_Error(HW_ERROR_INDEXRANGE);

					for (int j=0; j<r1; j++)
					{
						int index_1 = (int)first_range[j]-1;

			    		if ((index_1 < 0) || (index_1 >= data_M))
				    		throw OML_Error(HW_ERROR_INDEXRANGE);

						(*ret)(j,k) = (*cells)(index_1, index_2);
					}
				}
				Currency out (ret);
				out.SetMask(target.GetMask());
				return out;
			}
			else
			{
				throw OML_Error(HW_ERROR_INVIND);
			}
		}
	}
	else
	{
		throw OML_Error(HW_ERROR_INVINDOP);
	}

	return 0.0;
}

HML_CELLARRAY* ExprTreeEvaluator::CreateVararginCell(const std::vector<Currency>& params, int start_index)
{
    HML_CELLARRAY* param_cells = new HML_CELLARRAY;

	if (params.size()-start_index)
	{
		param_cells->Dimension(1, (int)params.size()-start_index, HML_CELLARRAY::REAL);

		for (int k = start_index; k < params.size(); k++)
			(*param_cells)(0,k-start_index) = params[k];	
	}
	else
	{
		param_cells->Dimension(0, 0, HML_CELLARRAY::REAL);
	}

    return param_cells;
}

Currency ExprTreeEvaluator::CallInternalFunction(FunctionInfo*fi, const std::vector<Currency>& param_values)
{
	class NestedFunctionHelper
	{
	public:
		NestedFunctionHelper(ExprTreeEvaluator* eval) : _eval(eval) { _eval->nested_function_marker++; }
		~NestedFunctionHelper() { _eval->nested_function_marker--; }
	private:
		ExprTreeEvaluator* _eval;
	};

	Currency r;

	if (debug_listener)
		debug_listener->PreFunctionCall(fi->FunctionName());

	if (fi->IsBuiltIn())
	{
		return CallBuiltinFunction(fi->Builtin(), fi->FunctionName(), param_values);
	}
	else
	{
		NestedFunctionHelper nfh(this);

		MemoryScope* parent = GetCurrentScope();

		FunctionInfo* temp = fi;

		OpenScope(temp);
		MemoryScope* memory = GetCurrentScope();
		
		int args_used = 0;

		std::vector<const std::string*> parameters = fi->Parameters();

		for (int j=0; j<parameters.size(); j++)
		{
			if (*parameters[j] == "varargin")
			{
                HML_CELLARRAY *param_cells = CreateVararginCell(param_values, j);
                args_used += param_cells->Size();
				memory->SetValue("varargin", param_cells);
				break;
			}
			else
			{
				if (j < param_values.size())
				{
					memory->SetValue(parameters[j], param_values[j]);
					args_used++;
				}
			}
		}

		if (GetExperimental())
		{
			// handle default parameters if necessary
			for (int j=0; j<parameters.size(); j++)
			{
				if (!memory->Contains(parameters[j]))
				{
					if (fi->HasDefaultValue(parameters[j]))
						memory->SetValue(parameters[j], fi->GetDefaultValue(parameters[j]));
				}
			}
		}
			
		if (args_used < param_values.size())
		{
			CloseScope();
			throw OML_Error(OML_ERR_NUMARGIN);
		}

		PushNargValues((int)param_values.size(), assignment_nargout);
        user_func_calls.push_back(fi->FunctionName());

		// before running the function, make sure we register all the nested functions into the proper scope
		OMLTree* statements = fi->Statements();

		if (statements)
		{
			int num_stmts = statements->ChildCount();

			for (int j=0; j<num_stmts; j++)
			{
				OMLTree*   tree = statements->GetChild(j);

				if (tree->GetType() == FUNC_DEF)
					RUN(tree);
			}

			// now we can run the statements.  Yes we're duplicating the nested function registration, but there could be cross-dependency
			r = RUN(statements);
		}

        user_func_calls.pop_back();

		PopNargValues();

		std::vector<const std::string*> return_values = fi->ReturnValues();
			
		if (!fi->IsAnonymous())
		{
			r = Currency(-1.0, Currency::TYPE_NOTHING);

			if (return_values.size())
			{
				const Currency& temp_cur = memory->GetValue(return_values[0]);

				if (!temp_cur.IsNothing())
					r = temp_cur;
			}
		}

		if (return_values.size() && (*return_values[0] == "varargout") && (r.IsCellArray()))
		{
			HML_CELLARRAY* cells = r.CellArray();

			if (cells->Size())
				r = (*cells)(0);
		}

		// need to call this before we close the scope
		if (debug_listener)
			debug_listener->PreFunctionReturn(fi->FunctionName());

		CloseScope();
	}

	r.ClearOutputName();

	if (fi->IsConstructor())
	{
		if (r.IsStruct())
			r.SetClass(fi->FunctionName());
	}
	return r;
}

FUNCPTR ExprTreeEvaluator::GetStdFunction(const std::string& func_name) const
{
    std::map<std::string, BuiltinFunc>::const_iterator iter = std_functions->find(func_name);
    return iter == std_functions->cend() ? nullptr : iter->second.fptr;
}

std::string ExprTreeEvaluator::GetHelpModule(const std::string& func_name)
{
	if (std_functions->find(func_name) != std_functions->end())
	{
		BuiltinFunc bif = (*std_functions)[func_name];
		return bif.md.module;
	}
	else 
	{
		FunctionInfo* fi = NULL;

		FindFunctionByName(func_name, &fi, NULL);

		if (fi)
		{
			std::string filename   = fi->FileName();
			size_t lastslash       = filename.find_last_of("/\\");
			std::string base_dir   = filename.substr(0, lastslash);
			size_t      start_pos  = base_dir.find_last_of("/\\");
			std::string result     = base_dir.substr(start_pos+1, base_dir.length());

			return result;
		}
	}

	return "";
}

MemoryScope* ExprTreeEvaluator::GetCurrentScope() const
{
	return msm->GetCurrentScope();
}

void ExprTreeEvaluator::OpenScope(FunctionInfo* info)
{
	return msm->OpenScope(info);
}

void ExprTreeEvaluator::CloseScope()
{
	MemoryScope* scope = msm->GetCurrentScope();
	return msm->CloseScope();
}

void ExprTreeEvaluator::ImportPathNames(const ExprTreeEvaluator *source)
{
	paths = source->paths;
	_owns_pathnames = false;
}

void ExprTreeEvaluator::ImportUserFileList(const ExprTreeEvaluator *source)
{
	userFileStreams = source->userFileStreams;
	_owns_userfiles = false;
}

void ExprTreeEvaluator::ImportMemoryScope(const ExprTreeEvaluator *source)
{
	msm = source->msm;
	_owns_msm = false;
}

void ExprTreeEvaluator::ImportFunctionList(const ExprTreeEvaluator *source)
{
	std_functions   = source->std_functions;
	functions       = source->functions;
	_owns_functions = false;
	class_info_map  = source->class_info_map;

	not_found_functions     = source->not_found_functions;
	preregistered_functions = source->preregistered_functions;
}

Currency ExprTreeEvaluator::UnaryOperator(OMLTree* tree)
{
	int	op = tree->GetType();

	OMLTree* child = tree->GetChild(0);
	Currency operand = RUN(child);

	return UnaryOperator(operand, op);
}

Currency ExprTreeEvaluator::UnaryOperator(const Currency& operand, int op)
{
	switch(op) 
	{
		case UMINUS:
			{
				return NegateOperator(operand);
			}
		case NEGATE:
			{
				Currency ret = NotOperator(operand);
				ret.SetMask(Currency::MASK_LOGICAL);
				return ret;
			}
		default:
			{
				return 0.0;
			}
	}
}

Currency ExprTreeEvaluator::BinaryOperator(OMLTree* tree)
{
	int	oper = tree->GetType();

	int goo = tree->ChildCount();

	OMLTree* child = tree->GetChild(0);

	Currency op1 = RUN(child);
	Currency op2;

	// Handle a parsing irregularity
	if (tree->ChildCount() == 3)
	{
		OMLTree* extra_node = tree->GetChild(1);
		OMLTree* child2 = tree->GetChild(2);
		op2 = RUN(child2);

		if (extra_node->GetType() == MINUS)
			op2 = NegateOperator(op2);
	}
	else
	{
		OMLTree* child1 = tree->GetChild(1);
		op2 = RUN(child1);
	}

	if (op1.GetMask() == Currency::MASK_STRING)
		op1.SetMask(Currency::MASK_NONE);

	if (op2.GetMask() == Currency::MASK_STRING)
		op2.SetMask(Currency::MASK_NONE);

	return BinaryOperator(op1, op2, oper);
}

Currency ExprTreeEvaluator::BinaryOperator(const Currency& lhs, const Currency& rhs, int oper)
{
	switch(oper) 
	{
		case PLUS:
			return AddOperator(lhs, rhs);
		case MINUS:
			return SubtractOperator(lhs, rhs);
		case TIMES:
			return MultiplyOperator(lhs, rhs);
		case ETIMES:
			return EntrywiseMultiplyOperator(lhs, rhs);
		case DIV:		
			return DivideOperator(lhs, rhs);
		case EDIV:
			return EntrywiseDivideOperator(lhs, rhs);
		case LDIV:					
			return LeftDivideOperator(lhs, rhs);
		case ELDIV:
			return EntrywiseLeftDivideOperator(lhs, rhs);
		case POW:
			return PowOperator(lhs, rhs);
		case DOTPOW:
			return EntrywisePowOperator(lhs, rhs);
		default:
			return Currency();
	}
}

Currency ExprTreeEvaluator::AddOperator(const Currency& op1,const Currency& op2)
{
	if (op1.IsScalar() && op2.IsScalar())
	{
		return op1.Scalar()+op2.Scalar();
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();

		hwMatrix* ret = allocateMatrix();
		ret->Add (*m1, op1.Scalar());
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix();
		ret->Add (*m1, op2.Scalar());
		return ret;
	}
	else if (op1.IsScalar() && op2.IsComplex())
	{
		hwComplex temp = op2.Complex();
		temp += op1.Scalar();
		return Currency(temp);
	}
	else if (op1.IsComplex() && op2.IsScalar())
	{
		hwComplex temp = op1.Complex();
		temp += op2.Scalar();
		return Currency(temp);
	}
	else if (op1.IsComplex() && op2.IsComplex())
	{
		hwComplex temp1 = op1.Complex();
		hwComplex temp2 = op2.Complex();
		temp1 += temp2;

		if (temp1.Imag() == 0.0)
			return Currency(temp1.Real());
		else
			return Currency(temp1);
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();

		hwMatrix* ret = allocateMatrix();
		ret->Add (*m1, op1.Complex());
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsComplex())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix();
		ret->Add (*m1, op2.Complex());
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		hwMathStatus stat;

		hwMatrix* ret = allocateMatrix();
		stat = ret->Add(*m1,*m2);

		if (stat.IsOk())
			return ret;
		else
			throw OML_Error(HW_ERROR_INCOMPDIM);
	}
	else if (op1.IsNDMatrix() || op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::AddOperator);
	}
	else if (op1.IsObject() || op2.IsObject())
	{
		return CallOverloadedOperator("plus", op1, op2);
	}
	else 
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
	}
}

Currency ExprTreeEvaluator::SubtractOperator(const Currency& op1, const Currency& op2)
{
	if (op1.IsScalar() && op2.IsScalar())
	{
		return op1.Scalar()-op2.Scalar();
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();

		hwMatrix* ret = allocateMatrix();
		ret->Subtr(op1.Scalar(), *m1);
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix();
		ret->Subtr(*m1, op2.Scalar());
		return ret;
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();

		hwMatrix* ret = allocateMatrix();
		ret->Subtr(op1.Complex(), *m1);
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsComplex())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix();
		ret->Subtr(*m1, op2.Complex());
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		hwMathStatus stat;

		hwMatrix* ret = allocateMatrix();
		stat = ret->Subtr(*m1,*m2);

		if (stat.IsOk())
			return ret;
		else
			throw OML_Error(HW_ERROR_INCOMPDIM);
	}
	else if (op1.IsScalar() && op2.IsComplex())
	{
		hwComplex temp(op1.Scalar()-op2.Real(), -op2.Imag());
		return temp;
	}
	else if (op1.IsComplex() && op2.IsScalar())
	{
		hwComplex temp = op1.Complex();
		temp -= op2.Scalar();
		return Currency(temp);
	}
	else if (op1.IsComplex() && op2.IsComplex())
	{
		hwComplex temp1 = op1.Complex();
		hwComplex temp2 = op2.Complex();
		temp1 -= temp2;

		if (temp1.Imag() == 0.0)
			return Currency(temp1.Real());
		else
			return Currency(temp1);
	}
	else if (op1.IsNDMatrix() || op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::SubtractOperator);
	}
	else if (op1.IsObject() || op2.IsObject())
	{
		return CallOverloadedOperator("minus", op1, op2);
	}
	else
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
	}
}

Currency ExprTreeEvaluator::MultiplyOperator(const Currency& op1, const Currency& op2)
{
	if (op1.IsScalar() && op2.IsScalar())
	{
		return op1.Scalar()*op2.Scalar();
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->MultEquals(op1.Scalar());
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->MultEquals(op2.Scalar());
		return ret;
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m2 = op2.Matrix();
		hwComplex c1 = op1.Complex();

		hwMatrix* ret = allocateMatrix(m2);
		ret->MultEquals(c1);
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsComplex())
	{
		const hwMatrix* m1 = op1.Matrix();
		hwComplex c2 = op2.Complex();

		hwMatrix* ret = allocateMatrix(m1);
		ret->MultEquals(c2);
		return ret;
	}
	else if (op1.IsScalar() && op2.IsComplex())
	{
		hwComplex temp = op2.Complex()*op1.Scalar();
		return Currency(temp);
	}
	else if (op1.IsComplex() && op2.IsScalar())
	{
		hwComplex temp = op1.Complex()*op2.Scalar();
		return Currency(temp);
	}
	else if (op1.IsComplex() && op2.IsComplex())
	{
		hwComplex temp = op1.Complex() * op2.Complex();

		if (temp.Imag() == 0.0)
			return Currency(temp.Real());
		else
			return Currency(temp);
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		int dim1 = m1->N();
		int dim2 = m2->M();

		if (m1->N() == m2->M())
		{
			hwMatrix* ret = allocateMatrix();
			ret->Mult(*m1,*m2);
			return ret;
		}
		
		throw OML_Error(HW_ERROR_INCOMPDIM);
	}
	else if (op1.IsNDMatrix() && !op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::MultiplyOperator);
	}
	else if (!op1.IsNDMatrix() && op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::MultiplyOperator);
	}
    else if (op1.IsObject() || op2.IsObject())
	{
		return CallOverloadedOperator("mtimes", op1, op2);
	}
	else
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
	}
}

Currency ExprTreeEvaluator::EntrywiseMultiplyOperator(const Currency& op1, const Currency& op2)
{
	if (op1.IsScalar() && op2.IsScalar())
	{
		return op1.Scalar()*op2.Scalar();
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->MultEquals(op1.Scalar());
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->MultEquals(op2.Scalar());
		return ret;
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->MultEquals(op1.Complex());
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsComplex())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->MultEquals(op2.Complex());
		return ret;
	}
	else if (op1.IsComplex() && op2.IsScalar())
	{
		return op1.Complex()*op2.Scalar();
	}
	else if (op1.IsScalar() && op2.IsComplex())
	{
		return op1.Scalar()*op2.Complex();
	}
	else if (op1.IsComplex() && op2.IsComplex())
	{
		return op1.Complex()*op2.Complex();
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		if ((m1->M() == m2->M()) && (m1->N() == m2->N()))
		{
			hwMatrix* ret = allocateMatrix();
			ret->MultByElems(*m1,*m2);
			return ret;
		}
		
		throw OML_Error(HW_ERROR_INCOMPDIM);
	}
	else if (op1.IsNDMatrix() || op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::EntrywiseMultiplyOperator);
	}
    else if (op1.IsObject() || op2.IsObject())
	{
		return CallOverloadedOperator("times", op1, op2);
	}
	else
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
	}
}

Currency ExprTreeEvaluator::DivideOperator(const Currency& op1, const Currency& op2)
{
	hwMathStatus stat;

	if (op1.IsScalar() && op2.IsScalar())
	{
		return op1.Scalar()/op2.Scalar();
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->DivideEquals(op2.Scalar());
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsComplex())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->DivideEquals(op2.Complex());
		return ret;
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();
		hwMatrix* ret = allocateMatrix();

		if (m1->N() == 1)
		{
			hwMatrix* ret = allocateMatrix();
			hwMatrix* temp = allocateMatrix(1, 1, hwMatrix::REAL);
			temp->SetElements(op1.Scalar());
			stat = ret->DivideRight(*temp, *m1);
			delete temp;
			return ret;
		}
		else
		{
			throw OML_Error(HW_ERROR_INCOMPDIM);
		}
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();
		hwMatrix* ret = allocateMatrix();

		if (m1->N() == 1)
		{
			hwMatrix* ret = allocateMatrix();
			hwMatrix* temp = allocateMatrix(1, 1, hwMatrix::COMPLEX);
			temp->SetElements(op1.Complex());
			stat = ret->DivideRight(*temp,*m1);
			delete temp;
			return ret;
		}
		else
		{
			throw OML_Error(HW_ERROR_INCOMPDIM);
		}
	}
	else if (op1.IsScalar() && op2.IsComplex())
	{
		hwComplex temp = op1.Scalar() / op2.Complex();
		return Currency(temp);
	}
	else if (op1.IsComplex() && op2.IsScalar())
	{
		hwComplex temp = op1.Complex() / op2.Scalar();
		return Currency(temp);
	}
	else if (op1.IsComplex() && op2.IsComplex())
	{
		hwComplex temp = op1.Complex() / op2.Complex();

		if (temp.Imag() == 0.0)
			return Currency(temp.Real());
		else
			return Currency(temp);
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		hwMatrix* ret = allocateMatrix();
		stat = ret->DivideRight(*m1,*m2);

		if (stat.IsOk())
		{
			return ret;
		}
		else if (stat.IsWarning())
		{
			stat.ResetArgs();   // suppress argument numbers because A/b has no argument list
			EvaluatorInterface eval(this);
			BuiltInFuncsUtils::SetWarning(eval, stat.GetMessage());
			return ret;
		}
		else
			throw OML_Error(HW_ERROR_INCOMPDIM);
	}
	else if (op1.IsNDMatrix() && !op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::DivideOperator);
	}
    else if (op1.IsObject() || op2.IsObject())
	{
		return CallOverloadedOperator("mrdivide", op1, op2);
	}
	else
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
	}
}

Currency ExprTreeEvaluator::EntrywiseDivideOperator(const Currency& op1, const Currency& op2)
{
	hwMathStatus stat;

	if (op1.IsScalar() && op2.IsScalar())
	{
		return op1.Scalar() / op2.Scalar();
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();
		hwMatrix* ret = allocateMatrix();

		stat = ret->Divide(op1.Scalar(), (*m1));
		return ret;
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		stat = ret->Divide(op1.Complex(), (*m1));
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsComplex())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->DivideEquals(op2.Complex());
		return ret;
	}
	else if (op1.IsComplex() && op2.IsScalar())
	{
		return op1.Complex()/op2.Scalar();
	}
	else if (op1.IsScalar() && op2.IsComplex())
	{
		return op1.Scalar()/op2.Complex();
	}
	else if (op1.IsComplex() && op2.IsComplex())
	{
		return op1.Complex()/op2.Complex();
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->DivideEquals(op2.Scalar());
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), m1->Type());
		stat = result->DivideByElems(*m1, *m2);

		if (stat.IsOk())
			return result;
		else
			throw OML_Error(HW_ERROR_INCOMPDIM);
	}
	else if (op1.IsNDMatrix() || op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::EntrywiseDivideOperator);
	}
    else if (op1.IsObject() || op2.IsObject())
	{
		return CallOverloadedOperator("rdivide", op1, op2);
	}
	else
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
	}
}

Currency ExprTreeEvaluator::LeftDivideOperator(const Currency& op1, const Currency& op2)
{
	hwMathStatus stat;

	if (op1.IsScalar() && op2.IsScalar())
	{
		return op2.Scalar()/op1.Scalar();
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->DivideEquals(op1.Scalar());

		return ret;
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->DivideEquals(op1.Complex());

		return ret;
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix();
		hwMatrix* temp = allocateMatrix(1, 1, hwMatrix::REAL);
		temp->SetElements(op2.Scalar());
		stat = ret->DivideLeft(*m1, *temp);
		delete temp;

		if (stat.IsOk())
			return ret;
		else
			throw OML_Error(HW_ERROR_INCOMPDIM);
	}
	else if (op1.IsMatrix() && op2.IsComplex())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix();
		hwMatrix* temp = allocateMatrix(1, 1, hwMatrix::COMPLEX);
		temp->SetElements(op2.Complex());
		stat = ret->DivideLeft(*m1, *temp);
		delete temp;

		if (stat.IsOk())
			return ret;
		else
			throw OML_Error(HW_ERROR_INCOMPDIM);
	}
	else if (op1.IsScalar() && op2.IsComplex())
	{
		hwComplex temp = op2.Complex() / op1.Scalar();
		return Currency(temp);
	}
	else if (op1.IsComplex() && op2.IsScalar())
	{
		hwComplex temp = op2.Scalar() / op1.Complex();
		return Currency(temp);
	}
	else if (op1.IsComplex() && op2.IsComplex())
	{
		hwComplex temp = op2.Complex() / op1.Complex();

		if (temp.Imag() == 0.0)
			return temp.Real();
		else
			return temp;
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		hwMatrix* ret = allocateMatrix();
		stat = ret->DivideLeft(*m1,*m2);

		if (stat.IsOk())
		{
			return ret;
		}
		else if (stat.IsWarning())
		{
			stat.ResetArgs();   // suppress argument numbers because A\b has no argument list
			EvaluatorInterface eval(this);
			BuiltInFuncsUtils::SetWarning(eval, stat.GetMessage());
			return ret;
		}
		else
			throw OML_Error(HW_ERROR_INCOMPDIM);
	}
	else if (!op1.IsNDMatrix() && op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::LeftDivideOperator);
	}
    else if (op1.IsObject() || op2.IsObject())
	{
		return CallOverloadedOperator("mldivide", op1, op2);
	}
	else
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
	}
}

Currency ExprTreeEvaluator::EntrywiseLeftDivideOperator(const Currency& op1, const Currency& op2)
{
	hwMathStatus stat;

	if (op1.IsScalar() && op2.IsScalar())
	{
		return op2.Scalar() / op1.Scalar();
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), m1->Type());
		stat = result->DivideByElems(*m2, *m1);

		if (stat.IsOk())
			return result;
		else
			throw OML_Error(HW_ERROR_INCOMPDIM);
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();
		hwMatrix* ret = allocateMatrix(m1);

		ret->DivideEquals(op1.Scalar());
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix();
		ret->Divide(op2.Scalar(), *m1);
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->MultEquals(op2.Scalar());
		return ret;
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->DivideEquals(op1.Complex());
		return ret;
	}
	else if (op1.IsMatrix() && op2.IsComplex())
	{
		const hwMatrix* m1 = op1.Matrix();

		hwMatrix* ret = allocateMatrix(m1);
		ret->Divide(op2.Complex(), *m1);
		return ret;
	}
	else if (op1.IsComplex() && op2.IsScalar())
	{
		return op2.Scalar()/op1.Complex();
	}
	else if (op1.IsScalar() && op2.IsComplex())
	{
		return op2.Complex()/op1.Scalar();
	}
	else if (op1.IsComplex() && op2.IsComplex())
	{
		return op2.Complex()/op1.Complex();
	}
	else if (op1.IsNDMatrix() || op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::EntrywiseLeftDivideOperator);
	}
	else if (op1.IsObject() || op2.IsObject())
	{
		return CallOverloadedOperator("ldivide", op1, op2);
	}
	else
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
	}
}

bool ExprTreeEvaluator::EqualityHelper(const Currency& lhs, const Currency& rhs)
{
	bool result = false;

	if (lhs.IsScalar() && rhs.IsScalar())
	{
		if (lhs.Scalar() == rhs.Scalar())
			result = true;
	}
	else if (lhs.IsScalar() && rhs.IsComplex())
	{
		if (rhs.Complex().IsReal() && (lhs.Scalar() == rhs.Complex().Real()))
			result = true;
	}
	else if (lhs.IsComplex() && rhs.IsScalar())
	{
		if (lhs.Complex().IsReal() && (rhs.Scalar() == lhs.Complex().Real()))
			result = true;
	}
	else if (lhs.IsComplex() && rhs.IsComplex())
	{
		hwComplex c1 = lhs.Complex();
		hwComplex c2 = rhs.Complex();
		return (c1 == c2);
	}
	else
	{
		throw OML_Error(HW_ERROR_UNSUPCOMP);
	}

	return result;
}

Currency ExprTreeEvaluator::EqualityOperatorEx(const Currency& op1, const Currency& op2)
{
	if (op1.IsScalar() && op2.IsScalar())
	{
		return op1.Scalar() == op2.Scalar();
	}
	if (op1.IsScalar() && op2.IsComplex())
	{
		return EqualityHelper(op1, op2);
	}
	if (op1.IsComplex() && op2.IsScalar())
	{
		return EqualityHelper(op1, op2);
	}
	else if (op1.IsScalar() && op2.IsMatrixOrString())
	{
		const hwMatrix* m2 = op2.Matrix();
		hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

		if (m2->IsReal())
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = (op1.Scalar() == (*m2)(j));
		}
		else
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = EqualityHelper(op1, m2->z(j));
		}

		return result;
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m2 = op2.Matrix();
		hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

		for (int j=0; j<m2->Size(); j++)
		{
			if (m2->IsReal())
				(*result)(j) = EqualityHelper(op1, (*m2)(j));
			else
				(*result)(j) = EqualityHelper(op1, m2->z(j));
		}

		return result;
	}
	else if (op1.IsMatrixOrString() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();

		if (!m1)
			return allocateMatrix();

		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

		int max_idx = m1->Size();
		double oper2  = op2.Scalar();

		if (m1->IsReal())
		{
			for (int j=0; j<max_idx; j++)
				(*result)(j) = ((*m1)(j) == oper2);
		}
		else
		{
			for (int j=0; j<max_idx; j++)
				(*result)(j) = EqualityHelper(m1->z(j), op2);
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsComplex())
	{
		const hwMatrix* m1 = op1.Matrix();
		
		if (!m1)
			return allocateMatrix();

		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

		for (int j=0; j<m1->Size(); j++)
		{
			if (m1->IsReal())
				(*result)(j) = EqualityHelper((*m1)(j), op2);
			else
				(*result)(j) = EqualityHelper(m1->z(j), op2);
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		if ((m1->M() == m2->M()) && (m1->N() == m2->N()))
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			if (m1->IsReal() && m2->IsReal())
			{
				for (int j=0; j<m1->Size(); j++)
					(*result)(j) = ((*m1)(j) == (*m2)(j));
			}
			else if (m1->IsReal() && !m2->IsReal())
			{
				for (int j=0; j<m1->Size(); j++)			
					(*result)(j) = EqualityHelper((*m1)(j), m2->z(j));
			}
			else if (!m1->IsReal() && m2->IsReal())
			{
				for (int j=0; j<m1->Size(); j++)			
					(*result)(j) = EqualityHelper((*m2)(j), m1->z(j));
			}
			else if (!m1->IsReal() && !m2->IsReal())
			{
				for (int j=0; j<m1->Size(); j++)			
					(*result)(j) = EqualityHelper(m1->z(j), m2->z(j));
			}

			return result;
		}
		else
		{
			throw OML_Error(HW_ERROR_INCOMPDIM);
		}
	}
	else if (op1.IsComplex() && op2.IsComplex())
	{
		hwComplex c1 = op1.Complex();
		hwComplex c2 = op2.Complex();

		return EqualityHelper(c1, c2);
	}
	else if (op1.IsNDMatrix() || op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::EqualityOperatorEx);
	}
	else if (op1.IsMatrixOrString() && op2.IsMatrixOrString())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		if ((m1->M() == m2->M()) && (m1->N() == m2->N()))
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			for (int j=0; j<m1->Size(); j++)
			{
				Currency lhs;
				Currency rhs;

				if (m1->IsReal())
					lhs = (*m1)(j);
				else
					lhs = m1->z(j).Real();

				if (m2->IsReal())
					rhs = (*m2)(j);
				else
					rhs = m2->z(j).Real();

				(*result)(j) = EqualityHelper(lhs, rhs);
			}

			return result;
		}
		else if (m1->Size() == 1)
		{
			hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = EqualityHelper((*m1)(0), (*m2)(j));

			return result;
		}
		else if (m2->Size() == 1)
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), m1->Type());

			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = EqualityHelper((*m1)(j), (*m2)(0));

			return result;
		}
		else
		{
			throw OML_Error(HW_ERROR_INCOMPDIM);
		}
	}
	else
	{
		throw OML_Error(HW_ERROR_INCOMPDIM);
	}

	throw OML_Error(HW_ERROR_UNSUPOP);
}

Currency ExprTreeEvaluator::EqualityOperator(OMLTree* tree)
{
	int	op = tree->GetType();

	OMLTree* child_0 = tree->GetChild(0);
	OMLTree* child_1 = tree->GetChild(1);

	Currency lhs = RUN(child_0);
	Currency rhs = RUN(child_1);

	if (lhs.IsObject() || rhs.IsObject())
	{
		if (op == EQUAL) 
			return CallOverloadedOperator("eq", lhs, rhs);
		else if (op == NEQUAL)
			return CallOverloadedOperator("ne", lhs, rhs);
		else
			throw OML_Error(HW_ERROR_UNSUPOP);
	}

	return EqualityOperator(lhs, rhs, op);
}

Currency ExprTreeEvaluator::EqualityOperator(const Currency& lhs, const Currency& rhs, int op)
{
	bool opposite = false;

	if (op == NEQUAL)
		opposite = true;

	Currency result = EqualityOperatorEx(lhs, rhs);

	if (opposite)
		result = NotOperator(result);

	result.SetMask(Currency::MASK_LOGICAL);

	return result;
}

Currency ExprTreeEvaluator::LessThanOperator(OMLTree* tree)
{
	int old_num_args = assignment_nargout;
	assignment_nargout = 1;

	OMLTree* child_0 = tree->GetChild(0);
	OMLTree* child_1 = tree->GetChild(1);

	Currency lhs = RUN(child_0);
	Currency rhs = RUN(child_1);

	assignment_nargout = old_num_args;

	Currency ret = LessThanOperator(lhs, rhs);
	ret.SetMask(Currency::MASK_LOGICAL);
	return ret;
}

Currency ExprTreeEvaluator::LessThanOperator(const Currency& op1, const Currency& op2)
{
	if (op1.IsScalar() && op2.IsScalar())
	{
		return op1.Scalar() < op2.Scalar();
	}
	if (op1.IsScalar() && op2.IsComplex())
	{
		return op1.Scalar() < op2.Complex().Mag();
	}
	if (op1.IsComplex() && op2.IsScalar())
	{
		return op1.Complex().Mag() < op2.Scalar();
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m2 = op2.Matrix();
		hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

		double lhs = op1.Scalar();

		if (m2->IsReal())
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs < (*m2)(j);
		}
		else
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs < m2->z(j).Mag();
		}

		return result;
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m2 = op2.Matrix();
		hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

		double lhs = op1.Complex().Mag();

		if (m2->IsReal())		
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs < (*m2)(j);
		}
		else
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs < m2->z(j).Mag();
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();
		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

		double rhs = op2.Scalar();

		if (m1->IsReal())
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) < rhs;
		}
		else
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = m1->z(j).Mag() < rhs;
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsComplex())
	{
		const hwMatrix* m1 = op1.Matrix();
		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

		double rhs = op2.Complex().Mag();

		if (m1->IsReal())
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) < rhs;
		}
		else
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = m1->z(j).Mag() < rhs;
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		if ((m1->M() == m2->M()) && (m1->N() == m2->N()))
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			if (m1->IsReal())
			{
				if (m2->IsReal())
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = (*m1)(j) < (*m2)(j);
				}
				else
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = (*m1)(j) < m2->z(j).Mag();
				}
			}
			else
			{
				if (m2->IsReal())
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = m1->z(j).Mag() < (*m2)(j);
				}
				else
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = m1->z(j).Mag() < m2->z(j).Mag();
				}
			}

			return result;
		}
		else
		{
			throw OML_Error(HW_ERROR_INCOMPDIM);
		}
	}
	else if (op1.IsComplex() && op2.IsComplex())
	{
		hwComplex c1 = op1.Complex();
		hwComplex c2 = op2.Complex();

		return c1.Mag() < c2.Mag();
	}
	else if (op1.IsNDMatrix() || op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::LessThanOperator);
	}
	else if (op1.IsString() && op2.IsString())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		if ((m1->M() == m2->M()) && (m1->N() == m2->N()))
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) < (*m2)(j);

			return result;
		}
		else if (m1->Size() == 1)
		{
			hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

			double lhs = (*m1)(0);

			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs < (*m2)(j);

			return result;
		}
		else if (m2->Size() == 1)
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			double rhs = (*m2)(0);

			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) < rhs;

			return result;
		}
		else
		{
			throw OML_Error(HW_ERROR_INCOMPDIM);
		}
	}
	else
	{
		throw OML_Error(HW_ERROR_INCOMPDIM);
	}

	throw OML_Error(HW_ERROR_UNSUPOP);
}

Currency ExprTreeEvaluator::GreaterThanOperator(OMLTree* tree)
{
	int old_num_args = assignment_nargout;
	assignment_nargout = 1;
	
	OMLTree* child_0 = tree->GetChild(0);
	OMLTree* child_1 = tree->GetChild(1);

	Currency lhs = RUN(child_0);
	Currency rhs = RUN(child_1);
	assignment_nargout = old_num_args;

	Currency ret = GreaterThanOperator(lhs, rhs);
	ret.SetMask(Currency::MASK_LOGICAL);
	return ret;
}

Currency ExprTreeEvaluator::GreaterThanOperator(const Currency& op1, const Currency& op2)
{
	if (op1.IsScalar() && op2.IsScalar())
	{
		return op1.Scalar() > op2.Scalar();
	}
	if (op1.IsScalar() && op2.IsComplex())
	{
		return op1.Scalar() > op2.Complex().Mag();
	}
	if (op1.IsComplex() && op2.IsScalar())
	{
		return op1.Complex().Mag() > op2.Scalar();
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m2 = op2.Matrix();
		hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

		double lhs = op1.Scalar();

		if (m2->IsReal())
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs > (*m2)(j);
		}
		else
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs > m2->z(j).Mag();
		}

		return result;
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m2 = op2.Matrix();
		hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

		double lhs = op1.Complex().Mag();

		if (m2->IsReal())		
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs > (*m2)(j);
		}
		else
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs > m2->z(j).Mag();
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();
		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

		double rhs = op2.Scalar();

		if (m1->IsReal())
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) > rhs;
		}
		else
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = m1->z(j).Mag() > rhs;
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsComplex())
	{
		const hwMatrix* m1 = op1.Matrix();
		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

		double rhs = op2.Complex().Mag();

		if (m1->IsReal())
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) > rhs;
		}
		else
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = m1->z(j).Mag() > rhs;
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		if ((m1->M() == m2->M()) && (m1->N() == m2->N()))
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			if (m1->IsReal())
			{
				if (m2->IsReal())
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = (*m1)(j) > (*m2)(j);
				}
				else
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = (*m1)(j) > m2->z(j).Mag();
				}
			}
			else
			{
				if (m2->IsReal())
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = m1->z(j).Mag() > (*m2)(j);
				}
				else
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = m1->z(j).Mag() > m2->z(j).Mag();
				}
			}

			return result;
		}
		else
		{
			throw OML_Error(HW_ERROR_INCOMPDIM);
		}
	}
	else if (op1.IsComplex() && op2.IsComplex())
	{
		hwComplex c1 = op1.Complex();
		hwComplex c2 = op2.Complex();

		return c1.Mag() > c2.Mag();
	}
	else if (op1.IsNDMatrix() || op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::GreaterThanOperator);
	}
	else if (op1.IsString() && op2.IsString())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		if ((m1->M() == m2->M()) && (m1->N() == m2->N()))
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) > (*m2)(j);

			return result;
		}
		else if (m1->Size() == 1)
		{
			hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

			double lhs = (*m1)(0);

			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs > (*m2)(j);

			return result;
		}
		else if (m2->Size() == 1)
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			double rhs = (*m2)(0);

			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) > rhs;

			return result;
		}
		else
		{
			throw OML_Error(HW_ERROR_INCOMPDIM);
		}
	}
	else
	{
		throw OML_Error(HW_ERROR_INCOMPDIM);
	}

	throw OML_Error(HW_ERROR_UNSUPOP);
}

Currency ExprTreeEvaluator::LessEqualOperator(OMLTree* tree)
{
	int old_num_args = assignment_nargout;
	assignment_nargout = 1;

	OMLTree* child_0 = tree->GetChild(0);
	OMLTree* child_1 = tree->GetChild(1);

	Currency lhs = RUN(child_0);
	Currency rhs = RUN(child_1);

	assignment_nargout = old_num_args;

	Currency ret = LessEqualOperator(lhs, rhs);
	ret.SetMask(Currency::MASK_LOGICAL);
	return ret;
}

Currency ExprTreeEvaluator::LessEqualOperator(const Currency& op1, const Currency& op2)
{
	if (op1.IsScalar() && op2.IsScalar())
	{
		return op1.Scalar() <= op2.Scalar();
	}
	if (op1.IsScalar() && op2.IsComplex())
	{
		return op1.Scalar() <= op2.Complex().Mag();
	}
	if (op1.IsComplex() && op2.IsScalar())
	{
		return op1.Complex().Mag() <= op2.Scalar();
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m2 = op2.Matrix();
		hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

		double lhs = op1.Scalar();

		if (m2->IsReal())
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs <= (*m2)(j);
		}
		else
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs <= m2->z(j).Mag();
		}

		return result;
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m2 = op2.Matrix();
		hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

		double lhs = op1.Complex().Mag();

		if (m2->IsReal())		
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs <= (*m2)(j);
		}
		else
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs <= m2->z(j).Mag();
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();
		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

		double rhs = op2.Scalar();

		if (m1->IsReal())
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) <= rhs;
		}
		else
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = m1->z(j).Mag() <= rhs;
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsComplex())
	{
		const hwMatrix* m1 = op1.Matrix();
		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

		double rhs = op2.Complex().Mag();

		if (m1->IsReal())
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) <= rhs;
		}
		else
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = m1->z(j).Mag() <= rhs;
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		if ((m1->M() == m2->M()) && (m1->N() == m2->N()))
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			if (m1->IsReal())
			{
				if (m2->IsReal())
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = (*m1)(j) <= (*m2)(j);
				}
				else
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = (*m1)(j) <= m2->z(j).Mag();
				}
			}
			else
			{
				if (m2->IsReal())
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = m1->z(j).Mag() <= (*m2)(j);
				}
				else
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = m1->z(j).Mag() <= m2->z(j).Mag();
				}
			}

			return result;
		}
		else
		{
			throw OML_Error(HW_ERROR_INCOMPDIM);
		}
	}
	else if (op1.IsComplex() && op2.IsComplex())
	{
		hwComplex c1 = op1.Complex();
		hwComplex c2 = op2.Complex();

		return c1.Mag() <= c2.Mag();
	}
	else if (op1.IsNDMatrix() || op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::LessEqualOperator);
	}
	else if (op1.IsString() && op2.IsString())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		if ((m1->M() == m2->M()) && (m1->N() == m2->N()))
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) <= (*m2)(j);

			return result;
		}
		else if (m1->Size() == 1)
		{
			hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

			double lhs = (*m1)(0);

			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs <= (*m2)(j);

			return result;
		}
		else if (m2->Size() == 1)
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			double rhs = (*m2)(0);

			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) <= rhs;

			return result;
		}
		else
		{
			throw OML_Error(HW_ERROR_INCOMPDIM);
		}
	}
	else
	{
		throw OML_Error(HW_ERROR_INCOMPDIM);
	}

	throw OML_Error(HW_ERROR_UNSUPOP);
}

Currency ExprTreeEvaluator::GreaterEqualOperator(OMLTree* tree)
{
	int old_num_args = assignment_nargout;
	assignment_nargout = 1;

	OMLTree* child_0 = tree->GetChild(0);
	OMLTree* child_1 = tree->GetChild(1);

	Currency lhs = RUN(child_0);
	Currency rhs = RUN(child_1);
	assignment_nargout = old_num_args;

	Currency ret = GreaterEqualOperator(lhs, rhs);
	ret.SetMask(Currency::MASK_LOGICAL);
	return ret;
}

Currency ExprTreeEvaluator::GreaterEqualOperator(const Currency& op1, const Currency& op2)
{
	if (op1.IsScalar() && op2.IsScalar())
	{
		return op1.Scalar() >= op2.Scalar();
	}
	if (op1.IsScalar() && op2.IsComplex())
	{
		return op1.Scalar() >= op2.Complex().Mag();
	}
	if (op1.IsComplex() && op2.IsScalar())
	{
		return op1.Complex().Mag() >= op2.Scalar();
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m2 = op2.Matrix();
		hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

		double lhs = op1.Scalar();

		if (m2->IsReal())
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs >= (*m2)(j);
		}
		else
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs >= m2->z(j).Mag();
		}

		return result;
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m2 = op2.Matrix();
		hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

		double lhs = op1.Complex().Mag();

		if (m2->IsReal())		
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs >= (*m2)(j);
		}
		else
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs >= m2->z(j).Mag();
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();
		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

		double rhs = op2.Scalar();

		if (m1->IsReal())
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) >= rhs;
		}
		else
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = m1->z(j).Mag() >= rhs;
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsComplex())
	{
		const hwMatrix* m1 = op1.Matrix();
		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

		double rhs = op2.Complex().Mag();

		if (m1->IsReal())
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) >= rhs;
		}
		else
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = m1->z(j).Mag() >= rhs;
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		if ((m1->M() == m2->M()) && (m1->N() == m2->N()))
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			if (m1->IsReal())
			{
				if (m2->IsReal())
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = (*m1)(j) >= (*m2)(j);
				}
				else
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = (*m1)(j) >= m2->z(j).Mag();
				}
			}
			else
			{
				if (m2->IsReal())
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = m1->z(j).Mag() >= (*m2)(j);
				}
				else
				{
					for (int j=0; j<m1->Size(); j++)
						(*result)(j) = m1->z(j).Mag() >= m2->z(j).Mag();
				}
			}

			return result;
		}
		else
		{
			throw OML_Error(HW_ERROR_INCOMPDIM);
		}
	}
	else if (op1.IsComplex() && op2.IsComplex())
	{
		hwComplex c1 = op1.Complex();
		hwComplex c2 = op2.Complex();

		return c1.Mag() >= c2.Mag();
	}
	else if (op1.IsNDMatrix() || op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::GreaterEqualOperator);
	}
	else if (op1.IsString() && op2.IsString())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		if ((m1->M() == m2->M()) && (m1->N() == m2->N()))
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) >= (*m2)(j);

			return result;
		}
		else if (m1->Size() == 1)
		{
			hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

			double lhs = (*m1)(0);

			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = lhs >= (*m2)(j);

			return result;
		}
		else if (m2->Size() == 1)
		{
			hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

			double rhs = (*m2)(0);

			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = (*m1)(j) >= rhs;

			return result;
		}
		else
		{
			throw OML_Error(HW_ERROR_INCOMPDIM);
		}
	}
	else
	{
		throw OML_Error(HW_ERROR_INCOMPDIM);
	}

	throw OML_Error(HW_ERROR_UNSUPOP);
}

bool ExprTreeEvaluator::LogicalHelper(double lhs, double rhs, int op)
{
	switch (op)
	{
		case AND:
			return (lhs && rhs);
		case OR:
			return (lhs || rhs);
		case LAND:
			return (lhs && rhs);
		case LOR:
			return (lhs || rhs);
	}

	return false;
}

Currency ExprTreeEvaluator::LogicalOperator(OMLTree* tree)
{
	int	op = tree->GetType();

	OMLTree* child_0 = tree->GetChild(0);
	OMLTree* child_1 = tree->GetChild(1);

	Currency op1 = RUN(child_0);
	Currency op2 = RUN(child_1);

	Currency ret = LogicalOperator(op1, op2, op);
	return ret;
}

Currency ExprTreeEvaluator::LogicalOperator(const Currency& lhs, const Currency& rhs, int op)
{
	Currency ret = LogicalOperatorEx(lhs, rhs, op);
	ret.SetMask(Currency::MASK_LOGICAL);
	return ret;
}

Currency ExprTreeEvaluator::LogicalOperatorEx(const Currency& lhs, const Currency& rhs, int op)
{
	if (lhs.IsScalar() && rhs.IsScalar())
	{
		return LogicalHelper(lhs.Scalar(), rhs.Scalar(), op);
	}
	else if (lhs.IsScalar() && rhs.IsMatrixOrString())
	{
		const hwMatrix* m2 = rhs.Matrix();
		hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

		for (int j=0; j<m2->Size(); j++)
		{
			if (m2->IsReal())
				(*result)(j) = LogicalHelper(lhs.Scalar(), (*m2)(j), op);
			else
				(*result)(j) = LogicalHelper(lhs.Scalar(), m2->z(j).Mag(), op);
		}

		return result;
	}
	else if (lhs.IsMatrixOrString() && rhs.IsScalar())
	{
		const hwMatrix* m1 = lhs.Matrix();
		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

		for (int j=0; j<m1->Size(); j++)
		{
			if (m1->IsReal())
				(*result)(j) = LogicalHelper((*m1)(j), rhs.Scalar(), op);
			else
				(*result)(j) = LogicalHelper(m1->z(j).Mag(), rhs.Scalar(), op);
		}

		return result;
	}
	else if (lhs.IsNDMatrix() || rhs.IsNDMatrix())
	{
		return oml_MatrixNUtil8(lhs, rhs, op, &ExprTreeEvaluator::LogicalOperatorEx);
	}
	else if (lhs.IsMatrixOrString() && rhs.IsMatrixOrString())
	{
		const hwMatrix* m1 = lhs.Matrix();
		const hwMatrix* m2 = rhs.Matrix();

		if (m1->IsReal() && m2->IsReal())
		{
			if ((m1->M() == m2->M()) && (m1->N() == m2->N()))
			{
				hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);

				for (int j=0; j<m1->Size(); j++)
				{
					double lhs;
					double rhs;

					if (m1->IsReal())
						lhs = (*m1)(j);
					else
						lhs = m1->z(j).Mag();

					if (m2->IsReal())
						rhs = (*m2)(j);
					else
						rhs = m2->z(j).Mag();

					(*result)(j) = LogicalHelper(lhs, rhs, op);
				}

				return result;
			}
		}
		else
		{
			throw OML_Error(HW_ERROR_COMPMATUNSUP);
		}
	}
	else if (lhs.IsComplex() && rhs.IsComplex())
	{
		hwComplex c1 = lhs.Complex();
		hwComplex c2 = rhs.Complex();

		return (LogicalHelper(c1.Mag(), c2.Mag(), op));
	}
	else if (lhs.IsScalar() && rhs.IsComplex())
	{
		hwComplex c2 = rhs.Complex();

		return (LogicalHelper(lhs.Scalar(), c2.Mag(), op));
	}
	else if (lhs.IsComplex() && rhs.IsScalar())
	{
		hwComplex c1 = lhs.Complex();

		return (LogicalHelper(c1.Mag(), rhs.Scalar(), op));
	}
    else if (lhs.IsComplex() && rhs.IsMatrixOrString())
	{
		hwComplex c1 = lhs.Complex();
		const hwMatrix* m2 = rhs.Matrix();
		hwMatrix* result = allocateMatrix(m2->M(), m2->N(), hwMatrix::REAL);

		if (m2->IsReal())
		{
			for (int j=0; j<m2->Size(); j++)
				(*result)(j) = LogicalHelper(c1.Mag(), (*m2)(j), op);
		}
		else
		{
			throw OML_Error(HW_ERROR_COMPMATUNSUP);
		}

		return result;
	}
	else if (lhs.IsMatrixOrString() && rhs.IsComplex())
	{
		const hwMatrix* m1 = lhs.Matrix();
		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);
		hwComplex c2 = rhs.Complex();

		if (m1->IsReal())
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = LogicalHelper((*m1)(j), c2.Mag(), op);
		}
		else
		{
			throw OML_Error(HW_ERROR_COMPMATUNSUP);
		}

		return result;
	}
	else
	{
		throw OML_Error(HW_ERROR_INCOMPDIM);
	}

	throw OML_Error(HW_ERROR_UNSUPOP);
}

bool ExprTreeEvaluator::ShortCircuitHelper(OMLTree* tree)
{
	Currency op1 = RUN(tree);;

	if (op1.IsScalar())
	{
		if (op1.Scalar() == 0.0)
			return false;
	}
	else if (op1.IsComplex())
	{
		hwComplex cplx = op1.Complex();

		if (cplx.Mag() == 0.0)
			return false;
	}
	else if (op1.IsMatrix())
	{
		const hwMatrix* mtx = op1.Matrix();

		if (mtx->IsEmpty())
		{
			return false;
		}
		else if (mtx->IsReal())
		{
			for (int j=0; j<mtx->Size(); j++)
			{
				if ((*mtx)(j) == 0.0)
					return false;
			}
		}
		else
		{
			for (int j=0; j<mtx->Size(); j++)
			{
				if (mtx->z(0).Mag() == 0.0)
					return false;
			}
		}
	}
	else
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
	}

	return true;
}

Currency ExprTreeEvaluator::ShortCircuitAndOperator(OMLTree* tree)
{
	if (!ShortCircuitHelper(tree->GetChild(0)))
		return false;

	return ShortCircuitHelper(tree->GetChild(1));
}

Currency ExprTreeEvaluator::ShortCircuitOrOperator(OMLTree* tree)
{
	if (ShortCircuitHelper(tree->GetChild(0)))
		return true;

	return ShortCircuitHelper(tree->GetChild(1));
}

Currency ExprTreeEvaluator::PowOperator(const Currency& op1, const Currency& op2)
{
	if (op1.IsScalar() && op2.IsScalar())
	{
		hwComplex result;
		result = hwComplex::pow_c(op1.Scalar(), op2.Scalar());

		if (result.IsReal())
			return result.Real();
		else
			return result;
	}
	else if (op1.IsComplex() && op2.IsScalar())
	{
		hwComplex result;
		result = hwComplex::pow(op1.Complex(), op2.Scalar());
		return result;
	}
	else if (op1.IsScalar() && op2.IsComplex())
	{
		hwComplex result;
		result = hwComplex::pow(op1.Scalar(), op2.Complex());
		return result;
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();
		hwMathStatus stat;

		hwMatrix* result = allocateMatrix();
		stat = result->Power(*m1, op2.Scalar());

		if (stat.IsOk())
			return result;
		else if (stat.IsWarning())
		{
			stat.ResetArgs();   // suppress argument numbers because A^b has no argument list
			EvaluatorInterface eval(this);
			BuiltInFuncsUtils::SetWarning(eval, stat.GetMessage());
			return result;
		}
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m2 = op2.Matrix();
		hwMathStatus stat;

		hwMatrix* result = allocateMatrix();
		stat = result->MatExp(*(m2) * log(op1.Scalar()));

		if (stat.IsOk())
			return result;
	}
	else if (op1.IsObject() || op2.IsObject())
	{
		return CallOverloadedOperator("mpower", op1, op2);
	}

	throw OML_Error(HW_ERROR_UNSUPOP);
}

Currency ExprTreeEvaluator::EntrywisePowOperator(const Currency& op1, const Currency& op2)
{
	if (op1.IsScalar() && op2.IsScalar())
	{
		hwComplex result = hwComplex::pow_c(op1.Scalar(), op2.Scalar());

		if (result.IsReal())
			return result.Real();
		else
			return result;
	}
	else if (op1.IsScalar() && op2.IsComplex())
	{
		hwComplex result;
		result = hwComplex::pow(op1.Scalar(), op2.Complex());
		return result;
	}
	else if (op1.IsComplex() && op2.IsScalar())
	{
		hwComplex result;
		result = hwComplex::pow(op1.Complex(), op2.Scalar());
		return result;
	}
	else if (op1.IsMatrix() && op2.IsScalar())
	{
		const hwMatrix* m1 = op1.Matrix();
		hwMathStatus stat;

		double exponent = op2.Scalar();

		if (exponent == 1.0)
		{
			hwMatrix* m_ret = (hwMatrix*)m1;
			m_ret->IncrRefCount();
			return m_ret;
		}

		hwMatrix* result = allocateMatrix();

		if (exponent == 2.0)
			stat = result->MultByElems(*m1, *m1);
		else
			stat = result->PowerByElems(*m1, exponent);

		if (stat.IsOk())
			return result;
	}
	else if (op1.IsScalar() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();
		hwMatrix* result = allocateMatrix (m1->M(), m1->N(), m1->Type());

		if (m1->IsReal())
		{
			for (int j=0; j<m1->Size(); j++)
				(*result)(j) = pow(op1.Scalar(), (*m1)(j));
		}
		else
		{
			for (int j=0; j<m1->Size(); j++)
				result->z(j) = hwComplex::pow(op1.Scalar(), m1->z(j));
		}

		return result;
	}
	else if (op1.IsComplex() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op2.Matrix();
		hwMatrix* result = allocateMatrix (m1->M(), m1->N(), m1->Type());
		result->MakeComplex();

		for (int j=0; j<m1->Size(); j++)
		{
			if (m1->IsReal())
				result->z(j) = hwComplex::pow(op1.Complex(), (*m1)(j));
			else
				result->z(j) = hwComplex::pow(op1.Complex(), m1->z(j));
		}

		return result;
	}
	else if (op1.IsMatrix() && op2.IsMatrix())
	{
		const hwMatrix* m1 = op1.Matrix();
		const hwMatrix* m2 = op2.Matrix();

		hwMathStatus stat;

		hwMatrix* result = allocateMatrix(m1->M(), m1->N(), m1->Type());
		stat = result->PowerByElems(*m1, *m2);

		if (stat.IsOk())
			return result;
		else if (stat.GetMsgCode() == HW_MATH_ERR_ARRAYSIZE)
			throw OML_Error(HW_ERROR_INCOMPDIM);
		else 
			throw OML_Error(HW_ERROR_COMPMATUNSUP);
	}
	else if (op1.IsNDMatrix() || op2.IsNDMatrix())
	{
		return oml_MatrixNUtil7(op1, op2, &ExprTreeEvaluator::EntrywisePowOperator);
	}
	else if (op1.IsObject() || op2.IsObject())
	{
		return CallOverloadedOperator("power", op1, op2);
	}

	throw OML_Error(HW_ERROR_UNSUPOP);
}

Currency ExprTreeEvaluator::NegateOperator(const Currency& op)
{
	if (op.IsScalar())
	{
		return -op.Scalar();
	}
	else if (op.IsComplex())
	{
		return -op.Complex();
	}
	else if (op.IsMatrix())
	{
		hwMatrix* ret = allocateMatrix((op.Matrix()));
		ret->MultEquals(-1.0);
		return ret;
	}
	else if (op.IsNDMatrix())
	{
		return oml_MatrixNUtil6(op, &ExprTreeEvaluator::NegateOperator);
	}
	else if (op.IsObject())
	{
		return CallOverloadedOperator("uminus", op);
	}

	return 0.0;
}	

Currency ExprTreeEvaluator::NotOperator(const Currency& op)
{
	if (op.IsScalar())
	{
		if (op.Scalar())
			return 0.0;
		else
			return 1.0;
	}
	else if (op.IsComplex())
	{
		// the IM=0 case is handled above
		throw OML_Error(HW_ERROR_LOGOPCOMP);
	}
	else if (op.IsMatrix() || op.IsString())
	{
		const hwMatrix* mtx = op.Matrix();
		hwMatrix* ret = allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);

		if (mtx->IsReal())
		{
			for (int j=0; j<ret->Size(); j++)
			{
				if ((*mtx)(j))
					(*ret)(j) = 0.0;
				else
					(*ret)(j) = 1.0;
			}
		}
		else
		{
			for (int j=0; j<ret->Size(); j++)
			{
				if (mtx->z(j).Imag() == 0.0)
				{
					if (mtx->z(j).Real())
						(*ret)(j) = 0.0;
					else
						(*ret)(j) = 1.0;
				}
				else
				{
					throw OML_Error(HW_ERROR_LOGOPCOMP);
    			}
	    	}
		}

		return ret;
	}
	else if (op.IsNDMatrix())
	{
		return oml_MatrixNUtil6(op, &ExprTreeEvaluator::NotOperator);
	}

	return 0.0;
}

Currency ExprTreeEvaluator::PostFunctionIndexHelper(const Currency& target, OMLTree* index_tree)
{
	std::vector<Currency> param_vals;

	if (index_tree->GetType() == PARAM_LIST)
	{
		int num_params = index_tree->ChildCount();
			
		param_vals.clear();
		param_vals.reserve(num_params);

		end_context_varname  = NULL;
		end_context_currency = (Currency*)&target;

		for (int j=0; j<num_params; j++)
		{
			end_context_index = j;
		
			OMLTree* child_j = index_tree->GetChild(j);
			param_vals.push_back(RUN(child_j));
						
			// Don't copy the Currency we get back from run
			// If this happens, it's an empty Currency only so it's cheap
			if (param_vals.back().IsNothing())
				param_vals.pop_back();	
		}

		return VariableIndex(target, param_vals);
	}

	throw OML_Error(HW_ERROR_INVFUNCALL);
}

Currency ExprTreeEvaluator::FunctionCall(OMLTree* tree)
{
	OMLTree* func = tree->GetChild(0);

	Currency val;
	const std::string* var_ptr = NULL;

	if (func->GetType() == IDENT)
	{
		if (!tree->u)
			tree->u = (void*)Currency::vm.GetStringPointer(func->GetText());

		 var_ptr = (const std::string*)tree->u;

		val = msm->GetValue(var_ptr);
	}
	else if (func->GetType() == STRUCT)
	{
		// it might be an object method call
		OMLTree* struct_child = func->GetChild(0);

		Currency temp = RUN(struct_child);

		if (temp.IsObject())
		{
			OMLTree* method = func->GetChild(1);
			OMLTree* params = NULL;
			
			if (tree->ChildCount() == 2)
				params = tree->GetChild(1);

			return ObjectMethodCall(&temp, params, method);
		}
		else if (temp.IsPointer())
		{
			Currency* new_temp = temp.Pointer();

			if (new_temp->IsObject())
			{
				OMLTree* method = func->GetChild(1);
				OMLTree* params = NULL;
			
				if (tree->ChildCount() == 2)
					params = tree->GetChild(1);

				return ObjectMethodCall(new_temp, params, method);
			}
			else
			{
				val = RUN(func);
			}
		}
		else
		{
			val = RUN(func);
		}
	}
	else
	{
		val = RUN(func);
	}

	std::vector<Currency> param_vals;

	int num_func_children = tree->ChildCount();
	if (num_func_children >= 2) // func + args
	{
		OMLTree* func_args = tree->GetChild(1);

		if (func_args->GetType() == PARAM_LIST)
		{
			int num_params = func_args->ChildCount();
			
			param_vals.reserve(num_params);

			if (!val.IsNothing())
			{
				// I hate having to repeat the loop code, but I can't think of a
				// better way to do it.  I need to set the end_context BEFORE evaluating
				// the children AND I need to UN-set it when I'm done - which would be
				// difficult to do otherwise because of all the places where we return. - JDS
				if (val.IsFunctionHandle())
				{
					for (int j=0; j<num_params; j++)
					{
						OMLTree* child_j = func_args->GetChild(j);
						Currency param = RUN(child_j);

						if (param.IsCellList())
						{
							HML_CELLARRAY* cells = param.CellArray();

							for (int k=0; k<cells->Size(); k++)
								param_vals.push_back((*cells)(k));
						}
						else
						{
							param_vals.push_back(param);

							if (param_vals.back().IsNothing())
								param_vals.pop_back();
						}
					}
					
					FunctionInfo* fi = new FunctionInfo(*val.FunctionHandle());
					Currency ret = CallInternalFunction(fi, param_vals);
					delete fi;
					return ret;
				}
				else
				{
					std::string* old_context  = end_context_varname;
					Currency*    old_currency = end_context_currency;

					end_context_varname = NULL;
					end_context_currency = &val;

					for (int j=0; j<num_params; j++)
					{
						if (num_params == 1)
							end_context_index = -1;
						else
							end_context_index = j;

						OMLTree* child_j = func_args->GetChild(j);

						param_vals.push_back(RUN(child_j));
						
						if (param_vals.back().IsNothing())
							param_vals.pop_back();
					}

					end_context_varname  = old_context;
					end_context_currency = old_currency;

					// it's a variable that we want to index
					Currency ret = VariableIndex(val, param_vals);

					end_context_varname = old_context;
					return ret;
				}
			}
			else
			{
                // Cache assignment_nargout during function calls
                int oldAssignmentArgs = assignment_nargout;
                assignment_nargout = 1;

				for (int j=0; j<num_params; j++)
				{
					OMLTree* child_j = func_args->GetChild(j);
					Currency param = RUN(child_j);

					if (param.IsCellList())
					{
						HML_CELLARRAY* cells = param.CellArray();

						for (int k=0; k<cells->Size(); k++)
							param_vals.push_back((*cells)(k));
					}
					else
					{
						param_vals.push_back(param);

						if (param_vals.back().IsNothing())
							param_vals.pop_back();
					}
				}
                // Restore assignment_nargout value
                assignment_nargout = oldAssignmentArgs;

				if (num_func_children == 2)
				{
					return CallFunction(var_ptr, param_vals);
				}
				else
				{
					Currency result = CallFunction(var_ptr, param_vals);

					// now we need to run the third child and use that to index into the result
					OMLTree* idx_args = tree->GetChild(2);
					return PostFunctionIndexHelper(result, idx_args);
				}
			}
		}
	}
	else if (num_func_children == 1) // func w/o args
	{
		const Currency& val = msm->GetValue(var_ptr);

		if (!val.IsNothing())
		{
			if (val.IsFunctionHandle())
			{
				FunctionInfo* fi = val.FunctionHandle();

				if (fi->Statements() || fi->IsBuiltIn())
					return CallInternalFunction(val.FunctionHandle(), param_vals);
				else if (!fi->Builtin())
					throw OML_Error("Invalid function handle");
			}

			return val;
		}

		return CallFunction(var_ptr, param_vals);
	}
	else
	{
		throw OML_Error(HW_ERROR_INVFUNCALL);
	}

	return -1;
}

bool ExprTreeEvaluator::FindFunctionByName(const std::string& func_name, FunctionInfo** fi, FUNCPTR* fptr)
{
	const std::string* var_ptr = Currency::vm.GetStringPointer(func_name);
	return FindFunctionByName(var_ptr, fi, fptr);
}

bool ExprTreeEvaluator::FindFunctionByName(const std::string* var_ptr, FunctionInfo** fi, FUNCPTR* fptr)
{
	std::string current_filename = msm->GetCurrentScope()->GetFilename();
	std::string found_file;
	bool        found = false;

	FunctionInfo* cur_fi = msm->GetCurrentScope()->GetFunctionInfo();

	const Currency& cur = msm->GetValue(var_ptr);

	if (!cur.IsNothing())
	{
		if (cur.IsFunctionHandle())
			*fi = cur.FunctionHandle();
	}
	else if (msm->GetNestedFunction(var_ptr))
	{
		*fi = msm->GetNestedFunction(var_ptr);
	}
	else if (msm->GetLocalFunction(var_ptr))
	{
		*fi = msm->GetLocalFunction(var_ptr);
	}
	else if ((functions->find(*var_ptr) != functions->end()) || (found = FindFunction(*var_ptr, found_file)))
	{
		if (found)
        {
			ParseAndRunFile(found_file, false);
            if (functions->find(*var_ptr) == functions->end())
                return false;
        }

		*fi = (*functions)[*var_ptr]->fi;
	}
	else if (std_functions->find(*var_ptr) != std_functions->end())
	{
		*fptr = (*std_functions)[*var_ptr].fptr;
	}
	else if (CheckForFunctionInAST(*var_ptr))
	{
		if (cur_fi && cur_fi->local_functions->find(var_ptr) != cur_fi->local_functions->end())
			*fi = (*cur_fi->local_functions)[var_ptr];
		else if ((functions->find(*var_ptr) != functions->end()) || (found = FindFunction(*var_ptr, found_file)))
			*fi = (*functions)[*var_ptr]->fi;
	}
	else
	{
		return false;
	}

	return true;
}

bool ExprTreeEvaluator::CheckForFunctionInAST(const std::string& func_name)
{
	int temp = nested_function_marker;
	nested_function_marker = false;

	if (current_tree)
	{
		int tree_children = current_tree->ChildCount();

		for (int j=current_statement_index+1; j<tree_children; j++)
		{
			OMLTree*    stmt = current_tree->GetChild(j);

			if (stmt->GetType() == FUNC_DEF)
			{
				OMLTree* fname = stmt->GetChild(0);

				std::string new_func = fname->GetText(); 

				if (func_name == new_func)
				{
					RUN(stmt);
					nested_function_marker = temp;
					return true;
				}
			}
		}
	}
	nested_function_marker = temp;
	return false;
}

std::vector<Currency> ExprTreeEvaluator::DoMultiReturnFunctionCall(FunctionInfo* fi_in, std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<std::string>* out_vars)
{
	std::vector<const std::string*> out_var_ptrs;

	if (out_vars)
	{
		for (int j=0; j<out_vars->size(); j++)
			out_var_ptrs.push_back(Currency::vm.GetStringPointer((*out_vars)[j]));
	}

	return DoMultiReturnFunctionCall(fi_in, param_values, num_ins, num_rets, suppress_output, &out_var_ptrs);
}

std::vector<Currency> ExprTreeEvaluator::DoMultiReturnFunctionCall(FunctionInfo* fi_in, std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<const std::string*>* out_vars)
{
	FunctionInfo* fi = fi_in;

	if (fi->IsBuiltIn())
		return DoMultiReturnFunctionCall(fi->Builtin(), fi->FunctionName(), param_values, num_ins, num_rets, suppress_output, out_vars);

	if (fi->IsAnonymous())
		return DoAnonymousMultiReturnFunctionCall(fi, param_values, num_ins, num_rets, suppress_output, out_vars);

	if (debug_listener)
		debug_listener->PreFunctionCall(fi->FunctionName());

	std::vector<const std::string*> parameters = fi->Parameters();
	
	if (parameters.size() < param_values.size())
		throw OML_Error(OML_ERR_NUMARGIN);
	else if (parameters.size() > param_values.size())
		throw OML_Error(OML_ERR_NUMARGIN);

	OpenScope(fi);

	for (size_t j=0; j<parameters.size(); j++)
		msm->SetValue(parameters[j], param_values[j]);

	PushNargValues(num_ins, num_rets);
    user_func_calls.push_back(fi->FunctionName());

	OMLTree* stmts = fi->Statements();
	RUN(stmts);

    user_func_calls.pop_back();
	PopNargValues();

	std::vector<const std::string*> return_values = fi->ReturnValues();

	int last = (int)return_values.size();

	if (last && (*return_values[last-1] == "varargout"))
	{
		Currency varargout = msm->GetValue("varargout");
		HML_CELLARRAY* varg_cells = varargout.CellArray();
		last = last-1+varg_cells->Size();
	}

	if (last < num_rets)
	{
		CloseScope();
		throw OML_Error(HW_ERROR_MISSRETURNS);
	}

	std::vector<Currency> ret;
	for (int j=0; j<num_rets; j++)
	{
		int inner_loop_size = 1;

		if (j >= return_values.size())
			break; // must have finished via varargout
		
		Currency cur = msm->GetValue(return_values[j]);
		
		if ((*return_values[j] == "varargout") && (j == (return_values.size()-1)))
		{
			HML_CELLARRAY* cells = cur.CellArray();
			inner_loop_size = cells->Size();
		}		

		for (int k=0; k<inner_loop_size; k++)
		{ 
			Currency loop_cur = cur;

			if (inner_loop_size > 1)
			{
				HML_CELLARRAY* cells = cur.CellArray();
				loop_cur = (*cells)(k);
			}
				
			ret.push_back(loop_cur);

			if (out_vars && out_vars->size())
			{
				// it's possible that not enough output variables were assigned
				if ((j+k) >= out_vars->size())
					break;

				msm->SetParentValue((*out_vars)[j+k], loop_cur);

				// this is a special case since we can't ordinarily push multiple outputs from a single statement
				if (!suppress_output)
				{
					cur.SetOutputName((*out_vars)[j+k]);
					PushResult(loop_cur); 
				}			
			}
		}
	}

	if (debug_listener)
		debug_listener->PreFunctionReturn(fi->FunctionName());

	CloseScope();
	return ret;
}

std::vector<Currency> ExprTreeEvaluator::DoMultiReturnFunctionCall(FUNCPTR fptr, const std::string& func_name, std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<std::string>* out_vars)
{
	std::vector<const std::string*> out_var_ptrs;

	if (out_vars)
	{
		for (int j=0; j<out_vars->size(); j++)
			out_var_ptrs.push_back(Currency::vm.GetStringPointer((*out_vars)[j]));
	}

	return DoMultiReturnFunctionCall(fptr, func_name, param_values, num_ins, num_rets, suppress_output, &out_var_ptrs);
}

std::vector<Currency> ExprTreeEvaluator::DoMultiReturnFunctionCall(FUNCPTR fptr, const std::string& func_name, std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<const std::string*>* out_vars)
{
	PushNargValues(num_ins, num_rets);

	std::vector<Currency> ret;

	if (func_name != "warning")
	builtin_error_scope = func_name;

	fptr(EvaluatorInterface(this), param_values, ret);
	builtin_error_scope = "";

	if (ret.size() < num_rets)
		throw OML_Error(HW_ERROR_MISSRETURNS);

	for (int j=0; j<num_rets; j++)
	{
		Currency loop_cur = ret[j];

		if (j >= out_vars->size())
			break;

		const std::string* output_name = (*out_vars)[j];

		msm->SetValue(output_name, loop_cur);

		// this is a special case since we can't ordinarily push multiple outputs from a single statement
		if (!suppress_output)
		{
			ret[j].SetOutputName(output_name);
			PushResult(ret[j]); 
		}
	}

	PopNargValues();

	return ret;
}

void ExprTreeEvaluator::DoMultiReturnFunctionCall(FUNCPTR fptr, const std::string& func_name, std::vector<Currency>& param_values, int num_ins, bool suppress_output, OMLTree* out_tree)
{
	int num_rets = out_tree->ChildCount();

	PushNargValues(num_ins, num_rets);

	std::vector<Currency> ret;

	if (func_name != "warning")
	builtin_error_scope = func_name;

	fptr(EvaluatorInterface(this), param_values, ret);
	builtin_error_scope = "";

	if (ret.size() < num_rets)
		throw OML_Error(HW_ERROR_MISSRETURNS);

	for (int j=0; j<num_rets; j++)
	{
		Currency loop_cur = ret[j];

		std::string* output_name = NULL;

		OMLTree*    out_var  = out_tree->GetChild(j);
		int         out_type = out_var->GetType();

		if ((out_type == CELL_VAL) || (out_type == FUNC))
		{
			OMLTree* indexed_var = out_var->GetChild(0);
			OMLTree* indices     = out_var->GetChild(1);
		
			if (!indexed_var->u)
				indexed_var->u = (void*)Currency::vm.GetStringPointer(indexed_var->GetText());

			output_name = (std::string*)indexed_var->u;

			Currency& target  = msm->GetMutableValue(*output_name);

			std::vector<Currency> output_indices;

			int num_indices = indices->ChildCount();

			for (int k=0; k<num_indices; k++)
			{
				OMLTree* child_k = indices->GetChild(k);
				Currency index = RUN(child_k);
				output_indices.push_back(index);
			}

			AssignHelper(target, output_indices, loop_cur);
		}
		else if (out_type == STRUCT)
		{
			OMLTree* ident_tree = out_var->GetChild(0);
			OMLTree* index_tree = NULL;
			OMLTree* field_tree = NULL;

			if (!ident_tree->u)
				ident_tree->u = (void*)Currency::vm.GetStringPointer(ident_tree->GetText());

			output_name = (std::string*)ident_tree->u;
				
			int num_children = out_var->ChildCount();

			if (num_children == 2)
			{
				field_tree = out_var->GetChild(1);
			}
			else if (num_children == 3)
			{
				index_tree = out_var->GetChild(1);					
				field_tree = out_var->GetChild(2);
			}

			Currency& target  = msm->GetMutableValue(*output_name);
			target.MakeStruct();

			StructAssignmentHelper(&target, index_tree, field_tree, loop_cur);
			msm->SetValue(output_name, target);
		}
		else if (out_type == IDENT)
		{
			if (!out_var->u)
				out_var->u = (void*)Currency::vm.GetStringPointer(out_var->GetText());

			output_name = (std::string*)out_var->u;

			msm->SetValue(output_name, loop_cur);
			}
			else
			{
			throw OML_Error("Unknown output type");
			}

			// this is a special case since we can't ordinarily push multiple outputs from a single statement
			if (!suppress_output)
			{
			ret[j].SetOutputName(output_name);
				PushResult(ret[j]); 
			}
		}

	PopNargValues();
}

std::vector<Currency> ExprTreeEvaluator::DoAnonymousMultiReturnFunctionCall(FunctionInfo* fi, std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<std::string>* out_vars)
{
	std::vector<const std::string*> out_var_ptrs;

	if (out_vars)
	{
		for (int j=0; j<out_vars->size(); j++)
			out_var_ptrs.push_back(Currency::vm.GetStringPointer((*out_vars)[j]));
	}

	return DoMultiReturnFunctionCall(fi, param_values, num_ins, num_rets, suppress_output, &out_var_ptrs);
}

std::vector<Currency> ExprTreeEvaluator::DoAnonymousMultiReturnFunctionCall(FunctionInfo* fi, std::vector<Currency>& param_values, int num_ins, int num_rets, bool suppress_output, std::vector<const std::string*>* out_vars)
{
	std::string my_func   = fi->RedirectedFunction();
	int         num_redir = fi->NumRedirectedInputs(); 

	std::vector<Currency> redirected_inputs;
	OpenScope(fi);

	std::vector<const std::string*> parameters = fi->Parameters();

	for (size_t j=0; j<parameters.size(); j++)
		msm->SetValue(parameters[j], param_values[j]);

	for (int j=0; j<num_redir; j++)
	{
		OMLTree* ri = fi->RedirectedInput(j);
		redirected_inputs.push_back(RUN(ri));
	}

	CloseScope();

	FunctionInfo* inner_fi = NULL;
	FUNCPTR       fptr     = NULL;

	if (!FindFunctionByName(my_func, &inner_fi, &fptr))
		throw OML_Error(HW_ERROR_UNKOWNFUN);

	std::vector<std::vector<Currency>> dummy;

	if (inner_fi)
		return DoMultiReturnFunctionCall(inner_fi, redirected_inputs, (int)redirected_inputs.size(), num_rets, suppress_output, out_vars);
	else if (fptr)
		return DoMultiReturnFunctionCall(fptr, my_func, redirected_inputs, (int)redirected_inputs.size(), num_rets, suppress_output, out_vars);
	else
		throw OML_Error(HW_ERROR_INVFUNCALL);
}

void ExprTreeEvaluator::DoMultiReturnFunctionCall(FunctionInfo* fi, std::vector<Currency>& param_values, int num_ins, bool suppress_output, OMLTree* out_tree)
{
	int num_rets = out_tree->ChildCount();

	if (fi->IsAnonymous())
	{
		std::vector<const std::string*> out_vars;

		for (size_t j=0; j<num_rets; j++)
		{
			OMLTree*    out_var  = out_tree->GetChild((int)j);
			int         out_type = out_var->GetType();

			if (out_type != IDENT)
				throw OML_Error("Invalid output");

			if (!out_var->u)
				out_var->u = (void*)Currency::vm.GetStringPointer(out_var->GetText());	

			const std::string* out_str = (std::string*)out_var->u;

			out_vars.push_back(out_str);
		}

		DoAnonymousMultiReturnFunctionCall(fi, param_values, num_ins, num_rets, suppress_output, &out_vars);
		return;
	}

	if (fi->IsBuiltIn())
		return DoMultiReturnFunctionCall(fi->Builtin(), fi->FunctionName(), param_values, num_ins, suppress_output, out_tree);

	if (debug_listener)
		debug_listener->PreFunctionCall(fi->FunctionName());
	
	//if (param_values.size() < fi->MinimumInputs())
	//	throw OML_Error(OML_ERR_NUMARGIN);

	OpenScope(fi);

	std::vector<const std::string*> parameters = fi->Parameters();

	for (size_t j=0; j<parameters.size(); j++)
	{
		if (j < param_values.size())
			msm->SetValue(parameters[j], param_values[j]);
		else if (*parameters[j] == "varargin")
			msm->SetValue(parameters[j], allocateCellArray());
	}

	PushNargValues(num_ins, num_rets);
    user_func_calls.push_back(fi->FunctionName());

	OMLTree* stmts = fi->Statements();
	RUN(stmts);

    user_func_calls.pop_back();
	PopNargValues();

	MemoryScope* stored_vals = new MemoryScope(fi);

	std::vector<const std::string*> return_values = fi->ReturnValues();

	for (size_t j=0; j<return_values.size(); j++)
	{
		const std::string* ret_name = return_values[j];
		stored_vals->SetValue(ret_name, msm->GetValue(ret_name));
	}

	CloseScope();

	int last = (int)return_values.size();

	if (last && (*return_values[last-1] == "varargout"))
	{
		Currency varargout = stored_vals->GetValue("varargout");

		if (varargout.IsCellArray())
		{
			HML_CELLARRAY* varg_cells = varargout.CellArray();

			if (varg_cells)
			last = last-1+varg_cells->Size();
		}
	}

	if (last < num_rets)
		throw OML_Error(HW_ERROR_MISSRETURNS);

	HML_CELLARRAY* cells;

	std::vector<Currency> ret;
	for (int j=0; j<num_rets; j++)
	{
		cells               = NULL;
		int inner_loop_size = 1;

		if (j >= return_values.size())
			break; // must have finished via varargout
		
		Currency cur = stored_vals->GetValue(return_values[j]);

		if (cur.IsNothing())
			throw OML_Error(HW_ERROR_UNASSIGNEMPTRIGHT);
		
		if ((*return_values[j] == "varargout") && (j == (return_values.size()-1)) && cur.IsCellArray())
		{
			cells = cur.CellArray();
			inner_loop_size = cells->Size();
		}		

		bool skip_var = false;

		for (int k=0; k<inner_loop_size; k++)
		{ 
			Currency loop_cur = cur;

			std::string* output_name = NULL;

			if (cells)
				loop_cur = (*cells)(k);
				
			ret.push_back(loop_cur);

				// it's possible that not enough output variables were assigned
			if ((j+k) >= num_rets)
					break;

			OMLTree*    out_var  = out_tree->GetChild(j+k);
			int         out_type = out_var->GetType();

			if ((out_type == CELL_VAL) || (out_type == FUNC))
			{
				OMLTree* indexed_var = out_var->GetChild(0);
				OMLTree* indices     = out_var->GetChild(1);

				if (!indexed_var->u)
					indexed_var->u = (void*)Currency::vm.GetStringPointer(indexed_var->GetText());

				output_name = (std::string*)indexed_var->u;

				Currency& target  = msm->GetMutableValue(*output_name);

				std::vector<Currency> output_indices;

				int num_indices = indices->ChildCount();
				
				// the indices have to be run before we call OpenScope
				for (int k=0; k<num_indices; k++)
				{
					OMLTree* child_k = indices->GetChild(k);
					Currency index = RUN(child_k);
					output_indices.push_back(index);
				}

				AssignHelper(target, output_indices, loop_cur);
			}
			else if (out_type == STRUCT)
			{
				OMLTree* ident_tree = out_var->GetChild(0);
				OMLTree* index_tree = NULL;
				OMLTree* field_tree = NULL;

				if (!ident_tree->u)
					ident_tree->u = (void*)Currency::vm.GetStringPointer(ident_tree->GetText());

				output_name = (std::string*)ident_tree->u;
				
				int num_children = out_var->ChildCount();

				if (num_children == 2)
				{
					field_tree = out_var->GetChild(1);
				}
				else if (num_children == 3)
				{
					index_tree = out_var->GetChild(1);					
					field_tree = out_var->GetChild(2);
				}

				Currency& target  = msm->GetMutableValue(*output_name);
				target.MakeStruct();

				StructAssignmentHelper(&target, index_tree, field_tree, loop_cur);
				msm->SetValue(output_name, target);
			}
			else if (out_type == IDENT)
			{
				if (!out_var->u)
					out_var->u = (void*)Currency::vm.GetStringPointer(out_var->GetText());

				output_name = (std::string*)out_var->u;

				msm->SetValue(output_name, loop_cur);
			}
			else if (out_type == NEGATE)
			{
				msm->SetValue("~", loop_cur);
				skip_var = true;
			}
			else
			{
				throw OML_Error("Unknown output type");
			}

			// this is a special case since we can't ordinarily push multiple outputs from a single statement
			if (!suppress_output && !skip_var)
			{
				cur.SetOutputName(output_name);
				PushResult(loop_cur); 
			}
		}
	}

	if (debug_listener)
		debug_listener->PreFunctionReturn(fi->FunctionName());

	delete stored_vals;

	return;

}

Currency ExprTreeEvaluator::MultiReturnFunctionCall(OMLTree* tree)
{
	OMLTree* func = tree->GetChild(1);
	
	if (func->GetType() == FUNC)
	{
		OMLTree* func_tree = func->GetChild(0);	

		if (func_tree->GetType() == IDENT)
		{
			if (!func->u)
				func->u = (void*)Currency::vm.GetStringPointer(func_tree->GetText());
		}
		else if (func_tree->GetType() == STRUCT)
		{
			return MRObjectMethodCall(tree);
		}
	}
	else if (func->GetType() == IDENT)
	{
		if (!func->u)
			func->u = (void*)Currency::vm.GetStringPointer(func->GetText());
	}
	else if (func->GetType() == CELL_VAL)
	{
		return CellExtraction(tree);
	}
	else if (func->GetType() == STRUCT)
	{
		return CellExtraction(tree); // it looks funny but it works
	}

	const std::string* my_func = (const std::string*)func->u;

	FunctionInfo*   fi   = NULL;
	FUNCPTR         fptr = NULL;

	bool            suppress_output = suppress_multi_ret_output; // running the statements in the function will change this
	
	OMLTree* ret_args = tree->GetChild(0);

	std::vector<Currency> param_values;
	int num_params = 0;

	OMLTree* func_args = NULL;
	
	if (func->ChildCount() == 2)
	{
		func_args  = func->GetChild(1);
		num_params = func_args->ChildCount();
	}
	param_values.reserve(num_params);

	std::vector<const std::string*> parameters;

	OMLTree* child_0 = NULL;
	Currency temp;

	if (func_args && func_args->ChildCount())
	{
		child_0 = func_args->GetChild(0);
		temp = RUN(child_0);
	}

	if (temp.IsObject())
	{
		if (HasOverloadedFunction(temp, *my_func))
		{ 
			std::string class_name = temp.GetClassname();
			ClassInfo*  ci         = (*class_info_map)[class_name];

			fi = ci->GetFunctionInfo(*my_func);
		}
	}

	if (!fi && my_func)
		FindFunctionByName(my_func, &fi, &fptr);

	if (!fi)
	{			
		std::string my_file;

		if (FindPrecompiledFunction(*my_func, my_file))
		{
			OMLTree* tree = OMLTree::ReadTreeFromFile(my_file);

			RUN(tree);

			delete tree;

			UserFunc* uf = (*functions)[*my_func];
			fi = uf->fi;
		}
	}

	if (!fi)
	{
		std::string my_file;
		std::string my_extension;

		if (FindEncryptedFunction(*my_func, my_file, my_extension))
		{
			RunEncryptedFile(my_extension, my_file);

			UserFunc* uf = NULL;

			if (functions->find(*my_func) != functions->end())
				uf = (*functions)[*my_func];

			if (uf)
				fi = uf->fi;
			else
				throw OML_Error("Error: unknown function: " + *my_func);
		}
	}

	if (fi)
		parameters = fi->Parameters();

	for (int j=0; j<num_params; j++)
	{
		if (j < parameters.size() && (*parameters[j] == "varargin"))
		{
			std::vector<Currency> vararg_params;

			for (int k=j; k<num_params; k++)
			{
				OMLTree* child_k = func_args->GetChild(k);
				Currency param = RUN(child_k);

				if (param.IsCellList())
				{
					HML_CELLARRAY* cells = param.CellArray();

					for (int k=0; k<cells->Size(); k++)
						vararg_params.push_back((*cells)(k));
				}
				else
				{
					vararg_params.push_back(param);

					if (vararg_params.back().IsNothing())
						vararg_params.pop_back();
				}
			}

			HML_CELLARRAY *param_cells = CreateVararginCell(vararg_params, 0);
			param_values.push_back(param_cells);
			break;
		}
		else
		{
			OMLTree* child_j = func_args->GetChild(j);
			param_values.push_back(RUN(child_j));
		}
	}

	if (fi)
	{
		DoMultiReturnFunctionCall(fi, param_values, num_params, suppress_output, ret_args);
	}
	else if (fptr)
	{
		DoMultiReturnFunctionCall(fptr, *my_func, param_values, num_params, suppress_output, ret_args);
	}
	else
	{
		std::string err_str = "Unknown function: ";
		err_str += *my_func;
		throw OML_Error(err_str);
	}
	
	return Currency(-1.0, Currency::TYPE_NOTHING);
}

Currency ExprTreeEvaluator::MRObjectMethodCall(OMLTree* tree)
{
	OMLTree* id_list   = tree->GetChild(0);
	OMLTree* func_tree = tree->GetChild(1);

	if (func_tree->ChildCount() == 2)
	{
		OMLTree*    struct_tree = func_tree->GetChild(0);
		OMLTree*    params_tree = func_tree->GetChild(1);

		if (params_tree->GetType() == PARAM_LIST)
		{
			OMLTree* struct_child_0 = struct_tree->GetChild(0);
			OMLTree* struct_child_1 = struct_tree->GetChild(1);

			Currency cur = RUN(struct_child_0);

			if (cur.IsStruct())
			{
				StructData* sd = cur.Struct();
				std::string field_name = struct_child_1->GetText();
				Currency fptr = sd->GetValue(0, -1, field_name);

				if (fptr.IsFunctionHandle())
				{
					std::vector<Currency> param_values;

					for (int j=0; j< (int)params_tree->ChildCount(); j++)
					{
						OMLTree* child_j = params_tree->GetChild(j);
						param_values.push_back(RUN(child_j));
					}

					FunctionInfo* fi = fptr.FunctionHandle();

					if (fi->IsBuiltIn())
						DoMultiReturnFunctionCall(fi->Builtin(), fi->FunctionName(), param_values, (int)param_values.size(), suppress_multi_ret_output, id_list);
					else
						DoMultiReturnFunctionCall(fi, param_values, (int)param_values.size(), suppress_multi_ret_output, id_list);
				}
			}
			else if (cur.IsObject())
			{
				ClassInfo* ci = (*class_info_map)[cur.GetClassname()];
				std::string field_name = struct_child_1->GetText();
				FunctionInfo* fi = ci->GetFunctionInfo(field_name);

				std::vector<Currency> param_values;

				param_values.push_back(cur);

				for (int j=0; j< (int)params_tree->ChildCount(); j++)
				{
					OMLTree* child_j = params_tree->GetChild(j);
					param_values.push_back(RUN(child_j));
				}

				DoMultiReturnFunctionCall(fi, param_values, (int)param_values.size(), suppress_multi_ret_output, id_list);
			}
		}
	}

	return Currency(-1, Currency::TYPE_NOTHING);
}

double GetTestVal(Currency cond_value)
{
	double test_val = 0.0;

	if (cond_value.IsScalar())
	{
		test_val = cond_value.Scalar();
	}
	else if (cond_value.IsComplex())
	{
		test_val = cond_value.Complex().Mag();
	}
	else if (cond_value.IsMatrix())
	{
		const hwMatrix* mtx = cond_value.Matrix();
		test_val = 1.0;

		if (mtx->IsReal() && !mtx->IsEmpty())
		{
			for (int j=0; j<mtx->Size(); j++)
			{
				if ((*mtx)(j) == 0.0)
				{
					test_val = 0.0;
					break;
				}
			}
		}
		else
		{
			test_val = 0.0;
		}
	}
	else
	{
		throw OML_Error("Invalid conditional value");
	}
	
	return test_val;
}

Currency ExprTreeEvaluator::Conditional(OMLTree* tree)
{
	int num_children = tree->ChildCount(); // should be at least 2 (one IF and one ELSE)
	bool found_one = false;
	Currency ret;

	const std::string* dbg_filename;
	int                dbg_linenum;

	for (int j=0; j<num_children-1; j++)
	{
		OMLTree* if_tree = tree->GetChild(j);
		
		DebugInfo stored_dbg = DebugInfo::DebugInfoFromTree(if_tree);
		dbg_filename = stored_dbg.FilenamePtr();
		dbg_linenum  = stored_dbg.LineNum();
		
		msm->GetCurrentScope()->SetDebugInfo(dbg_filename, dbg_linenum);

		OMLTree* child_0 = if_tree->GetChild(0);
		Currency cond_value = RUN(child_0);
		double test_val = GetTestVal(cond_value);
		
		if (test_val != 0.0)
		{
			if (if_tree->ChildCount() == 2)
			{
				OMLTree* child_1 = if_tree->GetChild(1);
				ret = RUN(child_1);
			}
			found_one = true;
			break;
		}
	}

	if (!found_one)
	{
		OMLTree* else_tree = tree->GetChild(num_children-1);

		if (else_tree->ChildCount() == 1)
		{
			OMLTree* child_0 = else_tree->GetChild(0);
			ret = RUN(child_0);
		}
	}

	if (ret.IsReturn() || ret.IsBreak() || ret.IsContinue())
		return ret;
	else
		return Currency(-1.0, Currency::TYPE_NOTHING);
}

Currency ExprTreeEvaluator::SwitchCase(OMLTree* tree)
{
	int      num_children = tree->ChildCount();

	OMLTree* child_0 = tree->GetChild(0);
	Currency switch_val       = RUN(child_0);
	Currency ret;

	if (!(switch_val.IsString() || switch_val.IsScalar()))
		throw OML_Error(HW_ERROR_INVSWITCH);

	for (int j=1; j<num_children; j++)
	{
		OMLTree*    case_tree = tree->GetChild(j);

		if (case_tree->GetType() == CASE)
		{
			OMLTree* case_0 = case_tree->GetChild(0);
			Currency case_val = RUN(case_0);
						
			if (case_val.IsScalar() && switch_val.IsScalar())
			{
				if (switch_val.Scalar() == case_val.Scalar())
				{
					if (case_tree->ChildCount() == 2)
					{
						OMLTree* case_1 = case_tree->GetChild(1);
						ret = RUN(case_1);
					}
					break;
				}
			}
			else if (case_val.IsString() && switch_val.IsString())
			{
				if (switch_val.StringVal() == case_val.StringVal())
				{
					if (case_tree->ChildCount() == 2)
					{
						OMLTree* case_1 = case_tree->GetChild(1);
						ret = RUN(case_1);
					}
					break;
				}
			}
			else if (case_val.IsCellArray())
			{
				// cycle through the cell array
				HML_CELLARRAY* cell              = case_val.CellArray();
				int            num_cell_children = cell->Size();
				bool           found             = false;

				for (int i=0; i<num_cell_children; i++)
				{
					Currency cell_val = (*cell)(i); // current cell val

					if (cell_val.IsScalar() && switch_val.IsScalar())
					{
						if (switch_val.Scalar() == cell_val.Scalar())
						{
							OMLTree* case_1 = case_tree->GetChild(1);
							ret = RUN(case_1);
							found = true;
							break;
						}
					}
					else if (cell_val.IsString() && switch_val.IsString())
					{
						if (switch_val.StringVal() == cell_val.StringVal())
						{
							OMLTree* case_1 = case_tree->GetChild(1);
							ret = RUN(case_1);
							found = true;
							break;
						}
					}
				}

				if (found == true)
					break;
			}
		}
		else if (case_tree->GetType() == OTHERWISE)
		{
			OMLTree* case_0 = case_tree->GetChild(0);
			ret = RUN(case_0);
			break;
		}
	}

	if (ret.IsReturn())
		return ret;

	return Currency(-1.0, Currency::TYPE_NOTHING);
}

Currency ExprTreeEvaluator::WhileLoop(OMLTree* tree)
{
	Currency ret;

	OMLTree* condition  = tree->GetChild(0);
	OMLTree* statements = tree->GetChild(1);

	bool old_val      = _store_suppressed;
	_store_suppressed = false;

	while (1)
	{
		Currency conditional = RUN(condition);

		double test_val = GetTestVal(conditional);

		if (test_val != 0)
			ret = RUN(statements);
		else
			break;

		if (ret.IsReturn())
			return ret;

		if (ret.IsBreak())
			break;

		if (ret.IsContinue())
			continue;
	}

	_store_suppressed = old_val;

	return Currency(-1, Currency::TYPE_NOTHING);
}

class LoopHelper
{
public:
	LoopHelper() { _mtx = NULL; _num_steps = 0; _is_range = false; }

	void      SetMatrix(const hwMatrix* mtx);
	void      SetRange(double start, double end, double incr);
	bool      Done();
	bool      IsRange() { return _is_range; }
	bool      IsVector();
	bool      IsMatrix();
	bool      IsReal();
	double    NextRealValue();
	double    FirstRealValue();
	double    LastRealValue();
	hwComplex NextComplexValue();
	hwComplex FirstComplexValue();
	hwMatrix* FirstColumn();
	hwMatrix* NextColumn();

private:
	long long        _counter;
	long long        _num_steps;
	double           _start;
	double           _incr;
	const hwMatrix*  _mtx;
	bool             _is_range;
};

void LoopHelper::SetMatrix(const hwMatrix* mtx)
{
	_counter   = 0;
	_mtx       = mtx;
	_is_range  = false;

	if (IsVector())
		_num_steps = mtx->Size()-1;
	else
		_num_steps = mtx->N()-1;
}

void LoopHelper::SetRange(double start, double end, double incr=1.0)
{
	_counter   = 0;
	_mtx       = NULL;
	_start     = start;
	_num_steps = (long long)((end - start) / incr);
	_incr      = incr;
	_is_range  = true;

	double num_steps_prime = (end - start)/incr + 1e-12;

	if (num_steps_prime >= 0.0)
	_num_steps = (long long)num_steps_prime;

	if ((num_steps_prime - _num_steps) > 1.0)
		throw OML_Error(HW_ERROR_UNSUPRANGEOP);
}

bool LoopHelper::Done()
{
	if (_num_steps < 0)
		return true;

	if (_counter <= _num_steps)
		return false;

	return true;
}

bool LoopHelper::IsVector()
{
	if (_mtx)
	{
		if ((_mtx->M() == 1) && (_mtx->N() != 1))
			return true;
		else if ((_mtx->M() != 1) && (_mtx->N() == 1)) 
			return true;
	}
	
	return false;
}

bool LoopHelper::IsMatrix()
{
	if (_mtx)
	{
		if ((_mtx->M() != 1) && (_mtx->N() != 1))
			return true;
	}
	
	return false;
}

bool LoopHelper::IsReal()
{
	if (_is_range)
		return true;
	else if (_mtx)
		return _mtx->IsReal();
	
	return true;
}

double LoopHelper::FirstRealValue()
{
	if (!_is_range)
		return (*_mtx)(0);
	else
		return _start;
}

double LoopHelper::LastRealValue()
{
	_counter = _num_steps+1;

	if (!_is_range)
		return (*_mtx)(_mtx->Size()-1);
	else
		return _start + _num_steps*_incr;
}

double LoopHelper::NextRealValue()
{
	if (!_is_range)
	{
		double ret = (*_mtx)((int)_counter);
		_counter++;
		return ret;
	}
	else
	{
		double ret = _start + _incr*_counter;
		_counter++;
		return ret;
	}
}

hwComplex LoopHelper::FirstComplexValue()
{
	return _mtx->z(0);
}

hwComplex LoopHelper::NextComplexValue()
{
	hwComplex ret = _mtx->z((int)_counter);
	_counter++;
	return ret;
}

hwMatrix* LoopHelper::FirstColumn()
{
	hwMatrix* new_mat = new hwMatrix;
	_mtx->ReadColumn(0, *new_mat);
	return new_mat;
}

hwMatrix* LoopHelper::NextColumn()
{
	hwMatrix* new_mat = new hwMatrix;
	_mtx->ReadColumn((int)_counter, *new_mat);
	_counter++;
	return new_mat;
}

Currency ExprTreeEvaluator::ForLoop(OMLTree* tree)
{
	OMLTree* loop_tree = tree->GetChild(0);

	if (!loop_tree->u)
		loop_tree->u = (void*)Currency::vm.GetStringPointer(loop_tree->GetText());
	std::string* loop_var = (std::string*)loop_tree->u;

	OMLTree* test = tree->GetChild(1);

	LoopHelper lh;

	bool old_val      = _store_suppressed;
	_store_suppressed = false;

	if (test->GetType() == COLON)
	{
		OMLTree* left = test->GetChild(0);
		OMLTree* right = test->GetChild(1);

		if (left->GetType() != COLON)
		{
			Currency start = RUN(left);
			Currency end   = RUN(right);

			if (start.IsScalar() && end.IsScalar())
			{
				double start_val = start.Scalar();
				double end_val   = end.Scalar();

				lh.SetRange(start_val, end_val);
			}
		}
		else
		{
			OMLTree* left2  = left->GetChild(0);
			OMLTree* right2 = left->GetChild(1);

			Currency start = RUN(left2);
			Currency incr  = RUN(right2);
			Currency end   = RUN(right);

			if (start.IsScalar() && end.IsScalar() && incr.IsScalar())
			{
				double start_val = start.Scalar();
				double end_val   = end.Scalar();
				double incr_val  = incr.Scalar();

				lh.SetRange(start_val, end_val, incr_val);
			}
		}
	}

	Currency loop_vals;

	if (!lh.IsRange())
	{
		loop_vals = RUN(test);

		if (loop_vals.IsScalar())
			lh.SetRange(loop_vals.Scalar(), loop_vals.Scalar());
		else if (loop_vals.IsMatrix())
			lh.SetMatrix(loop_vals.Matrix());
	}

	Currency ret;

	OMLTree* run_tree = NULL;
					
	if (tree->ChildCount() == 3)
		run_tree = tree->GetChild(2);

	if (lh.IsVector() || lh.IsRange())
	{
		if (!lh.Done())
		{
			if (lh.IsReal())
				msm->SetValue(loop_var, lh.FirstRealValue());
			else
				msm->SetValue(loop_var, lh.FirstComplexValue());

			if (!run_tree)
			{
				if (lh.IsReal())
					msm->SetValue(loop_var, lh.LastRealValue());
			}

			while (!lh.Done())
			{
				if (!msm->Contains(loop_var))
					msm->SetValue(loop_var, Currency());

				Currency& loop_cur = msm->GetMutableValue(loop_var);

				if (lh.IsReal())
					loop_cur.ReplaceScalar(lh.NextRealValue());
				else
					loop_cur.ReplaceComplex(lh.NextComplexValue());

				if (run_tree)
					ret = RUN(run_tree);

				if (ret.IsBreak())
				{
					break;
				}
				else if (ret.IsReturn())
				{
					return ret;
				}
				else if (ret.IsContinue())
				{
					continue;
				}
			}
		}
	}
	else if (lh.IsMatrix())
	{
		if (!lh.Done())
		{
			msm->SetValue(loop_var, lh.FirstColumn());

			Currency& loop_cur = msm->GetMutableValue(loop_var);

			while (!lh.Done())
			{
				loop_cur = msm->GetMutableValue(loop_var);

				loop_cur.ReplaceMatrix(lh.NextColumn());

				if (run_tree)
					ret = RUN(run_tree);

				if (ret.IsBreak())
				{
					break;
				}
				else if (ret.IsReturn())
				{
					return ret;
				}
				else if (ret.IsContinue())
				{
					continue;
				}
			}
		}
	}
	else if (loop_vals.IsCellArray())
	{
		HML_CELLARRAY* loop_values = loop_vals.CellArray();

		size_t loop_size = loop_values->N();
		
		if (loop_size)
		{
			HML_CELLARRAY* cell_col = new HML_CELLARRAY;
			loop_values->ReadColumn(0, *cell_col);
			msm->SetValue(loop_var, cell_col);
		}

		Currency& loop_cur = msm->GetMutableValue(loop_var);

		for (int j=0; j<loop_size; j++)
		{
			HML_CELLARRAY* loop_val = new HML_CELLARRAY;
			loop_values->ReadColumn(j, *loop_val);
			msm->SetValue(loop_var, loop_val);

			if (run_tree)
				ret = RUN(run_tree);

			if (ret.IsBreak())
			{
				break;
			}
			else if (ret.IsReturn())
			{
				return ret;
			}
			else if (ret.IsContinue())
			{
				continue;
			}
		}
	}
	else if (loop_vals.IsStruct())
	{
		StructData* loop_values = loop_vals.Struct();

		size_t loop_size = loop_values->Size();
		
		if (loop_size)
		{
			StructData* loop_val = loop_values->GetElement(1, -1);
			msm->SetValue(loop_var, loop_val);
		}

		Currency& loop_cur = msm->GetMutableValue(loop_var);

		for (int j=0; j<loop_size; j++)
		{
			StructData* loop_val = loop_values->GetElement(j+1, -1);
			msm->SetValue(loop_var, loop_val);

			if (run_tree)
				ret = RUN(run_tree);

			if (ret.IsBreak())
			{
				break;
			}
			else if (ret.IsReturn())
			{
				return ret;
			}
			else if (ret.IsContinue())
			{
				continue;
			}
		}
	}
	else
	{
		throw OML_Error("Invalid loop variable assignment");
	}

	_store_suppressed = old_val;

	return Currency(-1, Currency::TYPE_NOTHING);
}

FunctionInfo* ExprTreeEvaluator::FunctionInfoFromTree(OMLTree* tree)
{
	OMLTree* name = tree->GetChild(0);
	std::string func_name(name->GetText());

	OMLTree* ret_vals = tree->GetChild(1);

	int token_type = ret_vals->GetType(); 
	int num_rets   = ret_vals->ChildCount();

	std::vector<const std::string*> returns;

	for (int j=0; j<num_rets; j++)
	{
		OMLTree* child = ret_vals->GetChild(j);

		if (!child->u)
			child->u = (void*)Currency::vm.GetStringPointer(child->GetText());

		returns.push_back((const std::string*)child->u);
	}

	OMLTree* func_args = tree->GetChild(2);

	token_type = func_args->GetType(); // should be ID_LIST

	int num_params = func_args->ChildCount();

	std::vector<const std::string*>        params;
	std::map<const std::string*, Currency> default_values;

	for (int j=0; j<num_params; j++)
	{
		OMLTree* child = func_args->GetChild(j);

		if (!child->u)
			child->u = (void*)Currency::vm.GetStringPointer(child->GetText());

		const std::string* param_ptr = (const std::string*)child->u;

		params.push_back(param_ptr);

		if (GetExperimental())
		{
			if (child->ChildCount() == 1)
			{
				OMLTree* default_val = child->GetChild(0);
				Currency default_cur = RunTree(default_val);
				default_values[param_ptr] = default_cur;
			}
		}
	}

	OMLTree* copy_func_stmts = NULL;
	std::string       help_string;

	OMLTree*    func_stmts       = tree->GetChild(3);

	if (!func_stmts)
		return NULL; // already been taken care of

	int func_stmt_type   = func_stmts->GetType();

	if ((func_stmt_type == STATEMENT_LIST) || (func_stmt_type == STMT))
	{
		if ((!nested_function_marker) && (func_name != "anonymous"))
			copy_func_stmts = tree->DetachChild(3);
		else
			copy_func_stmts = new OMLTree(*tree->GetChild(3));

		int statement_count = func_stmts->ChildCount();

		for (int j=0; j<statement_count; j++)
		{
			OMLTree* test_stmt   = func_stmts->GetChild(j);

			if (test_stmt->GetType() == DUMMY)
			{
				int help_count = test_stmt->ChildCount();

				if (help_count)
				{
					OMLTree* help_tree = test_stmt->GetChild(0);

					if (help_string.size())
					{
						help_string += '\n';
						help_string += '\r';
					}

					help_string += help_tree->GetText();
				}
			}
			else
			{
				break;
			}
		}
	}

	DebugInfo dbg_info = DebugInfo::DebugInfoFromTree(tree);

	return new FunctionInfo(func_name, returns, params, default_values, copy_func_stmts, dbg_info.Filename(), help_string);
}

Currency ExprTreeEvaluator::FunctionDefinition(OMLTree* tree)
{
	FunctionInfo* fi      = FunctionInfoFromTree(tree);
	bool          defined = false;

	if (fi)
	{
		fi->IsNested(nested_function_marker);

		if (fi->IsNested())
		{
			msm->RegisterNestedFunction(fi);
			defined = true;
		}
		else if (fi->IsLocalFunction(_script_name))
		{
			std::string parent_filename = fi->FileName();
			size_t slash_index          = parent_filename.find_last_of("/\\");
			std::string root_file       = parent_filename.substr(slash_index+1);
			std::string parent_function = root_file.substr(0, root_file.find('.'));

			std::map<std::string, UserFunc*>::iterator iter = functions->find(parent_function);

			if (iter != functions->end())
			{
				UserFunc* uf            = iter->second;
				FunctionInfo* parent_fi = uf->fi;
				const std::string* func_name = Currency::vm.GetStringPointer(fi->FunctionName());
				(*parent_fi->local_functions)[func_name] = new FunctionInfo(*fi);
				defined = true;
			}
			else
			{
				// there is no parent so it must be a global function
			}
		}
		
		if (!defined)
		{
			int depth = msm->GetStackDepth();

			for (int j=0; j<depth; j++)
			{
				if (msm->GetScope(j)->FunctionName() == fi->FunctionName())
					throw OML_Error("Invalid function redefinition");
			}

			std::map<std::string, UserFunc*>::iterator iter = functions->find(fi->FunctionName());
			if (iter != functions->end())
				delete iter->second;
		
			(*functions)[fi->FunctionName()] = new UserFunc(fi);
		
			// don't do this for nested or local functions
			OnUpdateFuncList();
		}
	}

	return Currency(-1.0, Currency::TYPE_NOTHING);
}

Currency ExprTreeEvaluator::AnonymousFunctionDefinition(OMLTree* tree)
{
	int num_children = tree->ChildCount();
	
	FunctionInfo* fi = NULL;

	if (num_children == 1)
	{
		 //it's a simple remapping
		std::string func_name = tree->GetChild(0)->GetText();
		std::string my_file;
		std::string my_extension;

		FUNCPTR fptr = NULL;
		if (FindFunctionByName(func_name, &fi, &fptr))
		{
			if (fptr)
				fi = new FunctionInfo(func_name, fptr);
			else if (fi)
				fi = new FunctionInfo(*fi);
			else
				throw OML_Error(func_name + " is not a function handle");
		}
		else if (FindEncryptedFunction(func_name, my_file, my_extension))
		{
			RunEncryptedFile(my_extension, my_file);

			UserFunc* uf = NULL;

			if (functions->find(func_name) != functions->end())
				uf = (*functions)[func_name];

			if (uf)
				fi = new FunctionInfo(*uf->fi);
			else
				throw OML_Error("Error: unknown function: " + func_name);
		}
		else
		{
			throw OML_Error("Error: unknown function: " + func_name);
		}
	}
	else
	{
		fi = FunctionInfoFromTree(tree);
		fi->SetAnonymous(GetCurrentScope());
	}

	return fi;
}

bool ExprTreeEvaluator::ValidateMatrix(const std::vector<std::vector<Currency>>& currencies, unsigned int* out_rows, unsigned int* out_cols )
{
	unsigned int num_rows = (unsigned int)currencies.size();

	unsigned int      current_column_count = 0;
	unsigned int      current_row_height   = 0;
	unsigned int      total_row_height     = 0;
	unsigned int      target_column_count  = 0;

	bool has_cells = false;

	for (unsigned int j=0; j<num_rows; j++)
	{
		unsigned int row_cols = (unsigned int)currencies[j].size();

		current_column_count = 0;

		for (unsigned int k=0; k<row_cols; k++)
		{
			Currency temp_cur = currencies[j][k];

			if (k == 0)
			{
				if (temp_cur.IsMatrix() || temp_cur.IsString())
				{
					const hwMatrix* mtx = temp_cur.Matrix();
					
					current_row_height = mtx->M();
				}
				else if (temp_cur.IsCellArray())
				{
					HML_CELLARRAY* cells = temp_cur.CellArray();

					current_row_height = cells->M();
				}
				else if (temp_cur.IsScalar() || temp_cur.IsComplex())
				{
					current_row_height = 1;
				}
				else if (temp_cur.IsStruct())
				{
					StructData* sd = temp_cur.Struct();
					current_row_height = sd->M();
				}
				else
				{
					throw OML_Error(HW_ERROR_INVMATRIXINP);
				}
			}

			if (temp_cur.IsMatrix() || temp_cur.IsString())
			{
				const hwMatrix* mtx = temp_cur.Matrix();

				if (mtx)
				{
					if (mtx->M() != current_row_height)
						throw OML_Error(HW_ERROR_NOTCREATEMATRIX);

					if (has_cells)
						current_column_count++;
					else
						current_column_count += mtx->N();
				}
			}
			else if (temp_cur.IsCellArray())
			{
				HML_CELLARRAY* cells = temp_cur.CellArray();

				has_cells = true;

				if (cells)
				{
					if (cells->M() != current_row_height)
						throw OML_Error(HW_ERROR_NOTCREATEMATRIX);

					current_column_count += cells->N();
				}
			}
			else if (temp_cur.IsStruct())
			{
				StructData* sd = temp_cur.Struct();

				if (sd)
				{
					if (sd->M() != current_row_height)
						throw OML_Error(HW_ERROR_NOTCREATEMATRIX);

					current_column_count += sd->N();
				}
			}
			else if (temp_cur.IsScalar() || temp_cur.IsComplex())
			{
				// should really make sure the struct is a single and not an array
				if (current_row_height != 1)
					throw OML_Error(HW_ERROR_NOTCREATEMATRIX);

				current_column_count += 1;
			}
			else
			{
				throw OML_Error(HW_ERROR_NOTCREATEMATRIX);
			}
		}

		if (j==0)
			target_column_count = current_column_count;
		
		if (target_column_count != current_column_count)
		{
			if (current_column_count == 0)
				continue;

			throw OML_Error(HW_ERROR_NOTCREATEMATRIX);
		}

		total_row_height += current_row_height;
	}

	*out_rows = total_row_height;
	*out_cols = target_column_count;

	return true;
}

bool ExprTreeEvaluator::PadStringsIfNecessary(std::vector<std::vector<Currency>>& currencies)
{
	int string_length = 0;

	for (int j=0; j<currencies.size(); j++)
	{
		std::vector<Currency> cur_vec = currencies[j];

		if (cur_vec.size() != 1)
			return false;

		Currency cur = cur_vec[0];

		if (!cur.IsString())
			return false;

		const hwMatrix* mtx = cur.Matrix();

		if (string_length < mtx->N())
			string_length = (int)mtx->N();
	}

	for (int j=0; j<currencies.size(); j++)
	{
		std::vector<Currency> cur_vec = currencies[j];

		Currency cur = cur_vec[0];

		const hwMatrix* orig_mtx = cur.Matrix();

		if (orig_mtx->N() < string_length)
		{
			// need to pad the string
			hwMatrix* mtx = cur.GetWritableMatrix();

			int old_length = mtx->N();

			mtx->Resize(1, string_length, false);

			for (int k=old_length; k<string_length; k++)
				(*mtx)(k) = 32; // 32 is space

			cur_vec[0] = cur;
			currencies[j] = cur_vec;
		}
	}

	return true;
}

Currency ExprTreeEvaluator::MatrixCreation(OMLTree* tree)
{
// Non-trivial matrix creation (involving vectors and matrices as well as scalars) can
// be done two ways.  I've chosen to first validate the matrix and compute the proper 
// output size, and then fill the matrix.  The trade-off is that each tree node has to
// be run twice - once during validation and again during filling.  This operation can
// be done in a single pass, but that requires that the matrix be resized (reallocated)
// over and over.  It also makes the code quite a bit more complex b/c of the real/complex
// nature of the matrix.  Also, with the two pass approach, I don't have to worry about
// cleaning up the allocation if an error is thrown.  Otherwise, I'd have to catch the 
// error, free the matrix, and then re-throw it.

// addendum to above - and then I got smart and just stored the output currencies from
// running the tree the first time :-)

	unsigned int target_rows;
	unsigned int target_cols;

	bool         is_cell              = false;
	bool         is_string            = false;
	bool         is_struct            = false;
	bool         all_elements_logical = false;

	std::vector<std::vector<Currency>> currencies;

	unsigned int num_tree_rows = tree->ChildCount();
	for (unsigned int j=0; j<num_tree_rows; j++)
	{
		std::vector<Currency> temp;

		OMLTree* row_tree = tree->GetChild(j);
		unsigned int num_row_cols = row_tree->ChildCount();

		for (unsigned int k=0; k<num_row_cols; k++)
		{
			OMLTree* child_k = row_tree->GetChild(k);
			Currency temp_cur = RUN(child_k);

			if (temp_cur.IsMatrixOrString())
			{
				const hwMatrix* mtx = temp_cur.Matrix();

				if (temp_cur.IsString())
					is_string = true;

				if (!mtx)
					continue;

				if ((mtx->M() == 0) && (mtx->N() == 0))
					continue;
			}
			else if (temp_cur.IsCellArray())
			{
				const hwMatrix* mtx = NULL;

				if (temp_cur.IsCellList())
					mtx = temp_cur.ConvertToMatrix();

				if (!mtx)
				{
					HML_CELLARRAY* loc_cells = temp_cur.CellArray();			

					if (!loc_cells)
						continue;

					if ((loc_cells->M() == 0) && (loc_cells->N() == 0))
					{
						is_cell = true;
						continue;
					}
				}
			}
			else if (temp_cur.IsStruct())
			{
				StructData* sd = temp_cur.Struct();			

				if (!sd)
					continue;

				if ((sd->M() == 0) && (sd->N() == 0))
				{
					is_struct = true;
					continue;
				}
			}

			temp.push_back(temp_cur);
		}

		if (temp.size())
			currencies.push_back(temp);
	}

	if (currencies.size() == 1)
	{
		if (currencies[0].size() == 1)
		{
			if (currencies[0][0].IsScalar()) // this check may not be necessary; need to think about it more -- JDS
				return currencies[0][0];
		}
	}

	PadStringsIfNecessary(currencies);

	ValidateMatrix(currencies, &target_rows, &target_cols);

	hwMatrix* output_mtx = allocateMatrix();
	hwMathStatus status  = output_mtx->Dimension(target_rows, target_cols, hwMatrix::REAL);
	StructData* sd       = NULL;
	HML_CELLARRAY* cells = NULL;

	unsigned int row_offset = 0;
	unsigned int col_offset;

	for (unsigned int j=0; j<currencies.size(); j++)
	{
		if (j==0)
			all_elements_logical = true;

		unsigned int row_height = 0;

		col_offset = 0;

		bool init_flag = true;

		for (unsigned int k=0; k<currencies[j].size(); k++)
		{
			Currency temp_cur = currencies[j][k];

			if (!temp_cur.IsLogical())
				all_elements_logical = false;

			if (temp_cur.IsMatrix() || temp_cur.IsString())
			{
				const hwMatrix* mtx = temp_cur.Matrix();

				if (!mtx || mtx->IsEmpty())
				{
					if (init_flag && temp_cur.IsString())
						is_string = true;

					continue;
				}

				if (init_flag)
				{
					init_flag = false;

					row_height = mtx->M();

					if (temp_cur.IsString())
						is_string = true;
				}

				// If any row is not a string, then the output will not be a string
				if (temp_cur.IsString() && is_string)
					is_string = true;
				
				if (temp_cur.IsString() && !output_mtx->IsReal())
					throw OML_Error(HW_ERROR_NOTMIXCOMPSTR);
				
				if (cells)
				{
					(*cells)(row_offset, col_offset) = temp_cur;
					col_offset++;
				}
				else
				{
					output_mtx->WriteSubmatrix(row_offset, col_offset, *mtx);
					col_offset += mtx->N();
				}
			}
			else if (temp_cur.IsStruct())
			{
				// we're making an array of structs!
				is_struct = true;

				if (!sd)
				{
					sd = new StructData;
					sd->DimensionNew(target_rows, target_cols);
				}

				// need to validate the field names between sd & temp_cur here
				std::map<std::string, int> field_names = sd->GetFieldNames();
				std::map<std::string, int>::iterator iter;

				for (iter = field_names.begin(); iter != field_names.end(); iter++)
				{
					if (!temp_cur.Struct()->Contains(iter->first))
						throw OML_Error(HW_ERROR_STRUCTMUSTSAMEFIELDNAME);
				}

				StructData* loc_sd = temp_cur.Struct();

				if (init_flag)
				{
					init_flag = false;

					row_height = loc_sd->M();
				}

				for (int y=0; y<loc_sd->M(); y++)
				{
					for (int z=0; z<loc_sd->N(); z++)					
						sd->SetElement(row_offset+y+1, col_offset+z+1, loc_sd->GetElement(y+1, z+1));
				}

				col_offset += loc_sd->N();
			}
			else if (temp_cur.IsCellArray())
			{
				if (!is_cell)
				{
					if ((j != 0) || (k != 0))
						throw OML_Error(HW_ERROR_NOTCREATEMATRIX);

					// we're making an array of cells!
					is_cell = true;
				}

				if (!cells)
					cells = allocateCellArray(target_rows, target_cols);

				HML_CELLARRAY* loc_cells = temp_cur.CellArray();

				if (init_flag)
				{
					init_flag = false;

					row_height = loc_cells->M();
				}

				for (int y=0; y<loc_cells->M(); y++)
				{
					for (int z=0; z<loc_cells->N(); z++)					
						(*cells)(y+row_offset, z+col_offset) = (*loc_cells)(y,z);
				}

				col_offset += loc_cells->N();
			}
			else // must be scalar or complex
			{
				if (k==0)
					row_height = 1;

				if (temp_cur.IsScalar())
				{
					if (output_mtx->IsReal())
						(*output_mtx)(row_offset, col_offset) = temp_cur.Scalar();
					else
						output_mtx->z(row_offset, col_offset) = temp_cur.Scalar();
				}
				else if (temp_cur.IsComplex())
				{
					if (is_string)
						throw OML_Error(HW_ERROR_NOTMIXCOMPSTR);

					if (output_mtx->IsReal())
						output_mtx->MakeComplex();

					output_mtx->z(row_offset, col_offset).Set(temp_cur.Real(), temp_cur.Imag());
				}

				col_offset += 1;
			}
		}

		row_offset += row_height;
	}

	Currency ret(output_mtx);

	if (all_elements_logical)
		ret.SetMask(Currency::MASK_LOGICAL);

	if (is_string)
		ret.SetMask(Currency::MASK_STRING);

	if (is_struct)
		return sd;

	if (is_cell)
	{
		if (!cells)
			cells = allocateCellArray();

		return cells;
	}

	return ret;
}
#if OS_WIN
BOOL FileExists(LPCTSTR szPath)
{
  DWORD dwAttrib = GetFileAttributes(szPath);

  return (dwAttrib != INVALID_FILE_ATTRIBUTES && 
		 !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY));
}
#else
bool FileExists(const char *szPath)
{
	struct stat st;
	return (!stat(szPath, &st));
}
#endif

void ExprTreeEvaluator::SetInterrupt(bool val)
{
	_interrupt = val;
}

bool ExprTreeEvaluator::IsInterrupt()
{
	return _interrupt;
}

void ExprTreeEvaluator::SetPause(bool val)
{
	_pause = val;
}

bool ExprTreeEvaluator::IsPause()
{
	return _pause;
}

void ExprTreeEvaluator::SetDiary(bool val)
{
	if(val)
		diaryFile.open("diary.txt", ios::trunc);
	else
		diaryFile.close();
}
void ExprTreeEvaluator::SetDiary(std::string newfilename)
{
	diaryFile.open(newfilename, ios::trunc);
}

bool ExprTreeEvaluator::IsDiaryOpen()
{
	if(diaryFile.is_open())
		return true;
	else
		return false;
}

void ExprTreeEvaluator::WriteToDiaryFile(std::string input)
{
	diaryFile << input;
}

Currency ExprTreeEvaluator::StatementList(OMLTree* tree)
{
	int      num_statements = tree->ChildCount();
	Currency r;

	int starting_statement = 0;

	for (int i = starting_statement; i < num_statements; ++i)
	{
		if (tree == current_tree)
			current_statement_index = i;

		OMLTree* stmt = tree->GetChild(i);

		if (IsInterrupt())
			return Currency(-1,Currency::TYPE_BREAK);

		while (IsPause()) {
            _paused = true;
			sleep(0.5); // burns less cpu
        }
		
        _paused = false;

		int tok_type = stmt->GetType();
		if ((tok_type != FUNC_DEF) && (tok_type != DUMMY))
		{
			DebugInfo stored_dbg = DebugInfo::DebugInfoFromTree(stmt);
			msm->GetCurrentScope()->SetDebugInfo(stored_dbg.FilenamePtr(), stored_dbg.LineNum());

			if (debug_listener)
				debug_listener->PreStatement(DebugStateInfo(stored_dbg.Filename(), stored_dbg.LineNum()));
		}

		r = RUN(stmt);

		if (r.IsBreak() || r.IsReturn() || r.IsError() || r.IsContinue())
			break;
	}

	return r;
}

Currency ExprTreeEvaluator::Statement(OMLTree* tree)
{
	int      num_statements = tree->ChildCount();

	// this loop is to account for multiple statements on one 
	// line (e.g. a=3,b=4;)
	for (int i = 0; i < num_statements; i++)
	{
		OMLTree*    stmt      = tree->GetChild(i);
		OMLTree*    next_stmt = NULL;

		if (IsInterrupt())
			return Currency(-1,Currency::TYPE_BREAK);

		while (IsPause()) {
            _paused = true;
			sleep(0.5); // burns less cpu
        }
        _paused = false;
		
        bool suppress = false;

		if (i+1 < num_statements)
			next_stmt = tree->GetChild(i+1);

		if (next_stmt)
		{			
			if (next_stmt->GetType() == SEMIC)
			{
				suppress = true;
				suppress_multi_ret_output = true; // special case for output of assignments of a multi-return function call

			}
		}

		const Currency r = RUN(stmt);

		suppress_multi_ret_output = false; 

		if (r.IsBreak() || r.IsReturn() || r.IsError() || r.IsContinue())
			return r;

		PushResult(r, !suppress);
	}

	return Currency();
}

Currency ExprTreeEvaluator::AssignOperator(OMLTree* tree)
{
	OMLTree* lhs = tree->GetChild(0);

	if (!lhs->u)
	{
		lhs->u = (void*)Currency::vm.GetStringPointer(lhs->GetText());

		std::string* pString = (std::string*)lhs->u;

		if ((*pString == "spmd") || (*pString == "break") || (*pString == "continue"))
			throw OML_Error(HW_ERROR_NOTUSEKEY);
	}

	std::string* pString = (std::string*)lhs->u;

	if (lhs->GetType() == STRUCT)
	{
		Currency ret = StructAssignment(tree);

		if (ret.GetOutputName().size() == 0)
		{
			OMLTree* child = lhs->GetChild(0);
			std::string output_name;

			if (child->GetType() == CELL_VAL)
				output_name = child->GetChild(0)->GetText();
			ret.SetOutputName(output_name);
		}
		return ret;
	}

	int num_children = tree->ChildCount();
	
	int old_assignment_nargs = assignment_nargout;
	assignment_nargout = 1;

	OMLTree* child_1 = tree->GetChild(1);
	Currency value = RUN(child_1);

	if (value.IsNothing())
		throw OML_Error(HW_ERROR_UNASSIGNEMPTRIGHT);

	assignment_nargout = old_assignment_nargs;

	if (num_children == 3) // assignment with indices
	{
		std::vector<Currency> indices;
		
		OMLTree* index_tree = tree->GetChild(2);
		int num_indices = index_tree->ChildCount();

		indices.reserve(num_indices);

		end_context_varname = (std::string*)pString;

		for (int j=0; j<num_indices; j++)
		{
			if (num_indices == 1)
				end_context_index = -1;
			else
				end_context_index = j;

			OMLTree* child_j = index_tree->GetChild(j);
			const Currency& idx = RUN(child_j);
			indices.push_back(idx);
		}

		end_context_varname = NULL;

		if (lhs->GetType() == CELL_VAL)
		{
			_lhs_eval = true;
			Currency lhs_cur = RUN(lhs);
			_lhs_eval = false;

			AssignHelper(lhs_cur, indices, value);
			std::string parent_name = lhs->GetChild(0)->GetText();
			return msm->GetValue(parent_name);
		}

		bool      check  = msm->Contains(pString);
		Currency& temp   = msm->GetMutableValue(pString);

		try
		{
            if (value.IsNDMatrix())
			    NDAssignmetHelper(temp, indices, value);
            else
			    AssignHelper(temp, indices, value);

			temp.SetOutputName(pString);
			return temp;
		}
		catch (OML_Error& e)
		{
			if (!check)
				msm->Remove(*pString);
			throw e;
		}
	}
	else // assignment without indices
	{
		if (value.IsCellList()) // special case for cell lists
		{
			HML_CELLARRAY* cells = value.CellArray();

			if (cells->Size())
				value = (*cells)(0);
			else
				value = allocateCellArray();
		}

		if (value.IsFunctionHandle())
		{
			FunctionInfo* fi = value.FunctionHandle();
			fi->ClearAnonymousVariable(pString);
		}

		msm->SetValue(pString, value);

		return value;
	}
	return 0.0;
}

hwMatrix* AdjustMatrixSize(hwMatrix* temp_mtx, int index1)
{
	if (temp_mtx->GetRefCount() != 1)
		assert(0); // detach the matrix before calling this function

	if (index1 >= temp_mtx->Size())
	{
		if (temp_mtx->IsEmpty())
		{
			temp_mtx->Dimension(1, index1+1, hwMatrix::REAL);
			temp_mtx->SetElements(0.0);
		}
		else if (temp_mtx->M() == 1)
		{
			temp_mtx->Resize(1, index1+1, true); // make all the newly-created elements 0
		}
		else if (temp_mtx->N() == 1)
		{
			temp_mtx->Resize(index1+1, 1, true); // make all the newly-created elements 0
		}
		else
		{
			throw OML_Error(HW_ERROR_INVMATRESIZE);
		}
	}

	return temp_mtx;
}

hwMatrix* AdjustMatrixSize(hwMatrix* temp_mtx, const Currency& idx1, const Currency& idx2, int index1, int index2)
{
	if (!idx1.IsColon() && !idx2.IsColon())
	{
		if (((index1 >= temp_mtx->M()) || (index2 >= temp_mtx->N())))
		{   
			int size1 = temp_mtx->M();
			if (index1 >= size1)
				size1 = index1+1;

			int size2 = temp_mtx->N();
			if (index2 >= size2)
				size2 = index2+1;

			// This is unpleasant because of CoW.  We can't just resize the matrix, because we may
			// be sharing it with someone.  But we allocated it here, so we need to clean it up
			// later (since it never goes into a Currency).

			if (temp_mtx->IsEmpty())
			{
				hwMatrix* new_matrix = EvaluatorInterface::allocateMatrix(temp_mtx);
				hwMathStatus stat = new_matrix->Dimension(size1, size2, hwMatrix::REAL);

				if (!stat.IsOk())
					throw OML_Error(stat);

				new_matrix->SetElements(0.0);
				temp_mtx = new_matrix;
			}
			else
			{
				hwMatrix* new_matrix = temp_mtx;
				if (new_matrix->GetRefCount() != 1)
					new_matrix = EvaluatorInterface::allocateMatrix(temp_mtx);

				hwMathStatus stat = new_matrix->Resize(size1, size2, true); // make all the newly-created elements 0

				if (!stat.IsOk())
					throw OML_Error(stat);

				temp_mtx = new_matrix;
			}
		}
	}

	return temp_mtx;
}

void ExprTreeEvaluator::AssignHelper(Currency& target, const std::vector<Currency>& indices, const Currency& value)
{
	int  num_indices      = (int)indices.size();

	if (num_indices) // assignment with indices
	{
		if (num_indices > 2) // creation of a 3-D matrix where the target is not 3-D
			return NDAssignmetHelper(target, indices, value);

		if (target.IsNDMatrix())
			return NDAssignmetHelper(target, indices, value);

		if (target.IsEmpty() && value.IsStruct())
			target.MakeStruct();

		if (target.IsMatrix() || target.IsString() || target.IsPointer())
		{
			hwMatrix* temp_mtx = NULL;
			
			if (target.IsPointer())
				temp_mtx = target.Pointer()->GetWritableMatrix();
			else
				temp_mtx = target.GetWritableMatrix();

			Currency idx1 = indices[0];
			Currency idx2;

			int index1 = 0;
			int index2 = 0;

			if (idx1.IsLogical())
				idx1.ConvertToMatrix();
			
			if (idx1.IsPositiveInteger())
			{
				index1 = static_cast<int>(idx1.Scalar() - 1);
			}
			else if (idx1.IsMatrix())
			{
				if (num_indices == 1)
				{
					hwMatrix* new_mtx = SubmatrixSingleIndexHelper(target, idx1, value);

					target.ReplaceMatrix(new_mtx);

					return;
				}
				else if (num_indices == 2)
				{
					idx2 = indices[1];
					hwMatrix* new_mtx = SubmatrixDoubleIndexHelper(target, idx1, idx2, value);

					target.ReplaceMatrix(new_mtx);

					return;
				}
			}
			else if (!idx1.IsColon())
			{
				throw OML_Error(HW_ERROR_INDEXPOSINT);
			}

			if (num_indices == 2)
			{
				idx2 = indices[1];

				if (idx2.IsPositiveInteger())
				{
					index2 = static_cast<int>(idx2.Scalar() - 1); // remember we want to be 1-based
				}
				else if (idx2.IsMatrix())
				{
					hwMatrix* new_mtx = SubmatrixDoubleIndexHelper(target, idx1, idx2, value);
						
					target.ReplaceMatrix(new_mtx);

					return;
				}
				else if (!idx2.IsColon())
				{
					throw OML_Error(HW_ERROR_INDEXPOSINT);
				}
			}

			if (idx1.IsColon() || idx2.IsColon())
			{
				if (value.IsMatrix())
				{
					const hwMatrix* rhs_mtx = value.Matrix();

					if (rhs_mtx->IsEmpty()) // removing a row or column
					{		
						hwMatrix* new_matrix = allocateMatrix();

						if (idx1.IsColon()) 
						{
							if (idx2.IsScalar())
							{
								int col = (int)idx2.Scalar();				

								if ((col < 1) || (col > temp_mtx->N()))
									throw OML_Error(HW_ERROR_INDEXRANGE);

								hwMathStatus stat = new_matrix->DeleteColumns(*temp_mtx, col-1, 1);
							}
							else if (idx2.IsVector())
							{
								new_matrix = SubmatrixDoubleIndexHelper(target, idx1, idx2, value);
							}
						}
						else
						{
							if (idx1.IsScalar())
							{
								int row = (int)idx1.Scalar();

								if ((row < 1) || (row > temp_mtx->M()))
									throw OML_Error(HW_ERROR_INDEXRANGE);

								hwMathStatus stat = new_matrix->DeleteRows(*temp_mtx, row-1, 1);
							}
							else if (idx1.IsVector())
							{
								if (temp_mtx)
									target.ReplaceMatrix(temp_mtx);

								new_matrix = SubmatrixDoubleIndexHelper(target, idx1, idx2, value);
							}
						}

						target.ReplaceMatrix(new_matrix);
						return;
					}
					else
					{
						if (temp_mtx)
							target.ReplaceMatrix(temp_mtx);

						hwMatrix* new_mtx = NULL;
						
						if (num_indices == 2)
							new_mtx = SubmatrixDoubleIndexHelper(target, idx1, idx2, value);
						else if (num_indices == 1)
							new_mtx = SubmatrixSingleIndexHelper(target, idx1, value);

						target.ReplaceMatrix(new_mtx);
						return;
					}
				}
				else if (value.IsScalar() || value.IsComplex() || value.IsEmpty())
				{
					if (temp_mtx)
						target.ReplaceMatrix(temp_mtx);

					hwMatrix* new_mtx = NULL;
						
					if (num_indices == 2)
						new_mtx = SubmatrixDoubleIndexHelper(target, idx1, idx2, value);
					else if (num_indices == 1)
						new_mtx = SubmatrixSingleIndexHelper(target, idx1, value);

					target.ReplaceMatrix(new_mtx);
					return;
				}
			}
			else
			{
				hwMatrix* new_matrix = temp_mtx;

				if (temp_mtx->GetRefCount() != 1)
					new_matrix = allocateMatrix(temp_mtx);

				if (num_indices == 2)
				{
					if (value.IsScalar())
					{
						new_matrix = AdjustMatrixSize(new_matrix, idx1, idx2, index1, index2);

						if (new_matrix->IsReal())
							(*new_matrix)(index1, index2) = value.Scalar();
						else
							new_matrix->z(index1, index2) = value.Scalar();
					}
					else if (value.IsComplex())
					{
						new_matrix = AdjustMatrixSize(new_matrix, idx1, idx2, index1, index2);
						new_matrix->MakeComplex();
						new_matrix->z(index1, index2) = value.Complex();
					}
					else if (value.IsCharacter())
					{
						std::string str = value.StringVal();
						new_matrix = AdjustMatrixSize(new_matrix, idx1, idx2, index1, index2);
						(*new_matrix)(index1, index2) = (unsigned int)str[0];
					}
					else
					{
						throw OML_Error(HW_ERROR_INVRIGHT);
					}
				}
				else if (num_indices == 1)
				{
					if (value.IsScalar())
					{
						new_matrix = AdjustMatrixSize(new_matrix, index1);

						if (new_matrix->IsReal())
							(*new_matrix)(index1) = value.Scalar();
						else
							new_matrix->z(index1) = value.Scalar();
					}
					else if (value.IsComplex())
					{
						new_matrix = AdjustMatrixSize(new_matrix, index1);
						new_matrix->MakeComplex();
						new_matrix->z(index1) = value.Complex();
					}
					else if (value.IsEmpty())
					{
						if (new_matrix->M() == 1)
						{
							new_matrix->DeleteColumns(index1, 1);
						}
						else if (new_matrix->N() == 1)
						{
							new_matrix->DeleteRows(index1, 1);
						}
						else
						{
							new_matrix->Reshape(1, new_matrix->Size());
							new_matrix->DeleteColumns(index1, 1);
						}
					}
					else if (value.IsString())
					{
                        const hwMatrix* mtx = value.Matrix();
                        if (!mtx || mtx->IsEmpty())
                        {
                            new_matrix->DeleteElements(index1);
                        }
                        else if (mtx->Size() == 1)
                        {
						    new_matrix = AdjustMatrixSize(new_matrix, index1);
                            (*new_matrix)(index1) = static_cast<unsigned int>((*mtx)(0));
                        }
						else
						{
							throw OML_Error(HW_ERROR_INVRIGHT);
						}
					}
					else
					{
						throw OML_Error(HW_ERROR_INVRIGHT);
					}
				}

				target.ReplaceMatrix(new_matrix);  // make sure this doesn't affect the mask
				return;
			}
		}
		else if (target.IsCellArray())
		{
			if (value.IsCellArray())
			{
				HML_CELLARRAY* cells = value.CellArray();

				if (indices[0].IsScalar())
				{
					if ((cells->M() ==1) && (cells->N()==1))
					{
						HML_CELLARRAY* dest_cells = target.CellArray();

						if (dest_cells->GetRefCount() != 1)
						{
							dest_cells = allocateCellArray(dest_cells);
							target.ReplaceCellArray(dest_cells);
						}

						int index1 = static_cast<int>(indices[0].Scalar() - 1); // remember we want to be 1-based

						if (indices.size() > 1)
						{
							int index2 = static_cast<int>(indices[1].Scalar() - 1); // remember we want to be 1-based

							int old_m = dest_cells->M();
							int old_n = dest_cells->N();
							int new_m = old_m;
							int new_n = old_n;

							if (index1 >= old_m)
								new_m = index1+1;

							if (index2 >= old_n)
								new_n = index2+1;

							if ((old_m != new_m) || (old_n != new_n))
							{
								if (dest_cells->IsEmpty())
									dest_cells->Dimension(new_m, new_n, HML_CELLARRAY::REAL);
								else
									dest_cells->Resize(new_m, new_n);
							}

							(*dest_cells)(index1,index2) = (*cells)(0,0);
						}
						else
						{
							if (dest_cells->IsEmpty())
								dest_cells->Dimension(index1+1, HML_CELLARRAY::REAL);
							else if (index1 >= dest_cells->Size())
								dest_cells->Resize(index1+1);

							(*dest_cells)(index1) = (*cells)(0,0);
						}
					}
					else if ((indices.size() > 1) && (indices[1].IsColon()))
					{
						if (value.IsCellArray())
						{
							HML_CELLARRAY* dest_cells = target.CellArray();
							HML_CELLARRAY* cells = value.CellArray();

							if (dest_cells->GetRefCount() != 1)
							{
								dest_cells = allocateCellArray(dest_cells);
								target.ReplaceCellArray(dest_cells);
							}

							if ((dest_cells->N() == cells->N()) && (cells->M() == 1))
							{
								int row_index = static_cast<int>(indices[0].Scalar()-1);

								if (row_index >= dest_cells->M())
								{
									if (dest_cells->IsEmpty())
										dest_cells->Dimension(row_index+1, dest_cells->N(), HML_CELLARRAY::REAL);
									else
										dest_cells->Resize(row_index+1, dest_cells->N());
								}

								for (int j=0; j<cells->N(); j++)
									(*dest_cells)(row_index, j) = (*cells)(0, j);
							}
							else
							{
								throw OML_Error(HW_ERROR_INVRHS);
							}
						}
					}
					else if (cells->IsEmpty())
					{
						HML_CELLARRAY* dest_cells = target.CellArray();

						if (dest_cells->GetRefCount() != 1)
						{
							dest_cells = allocateCellArray(dest_cells);
							target.ReplaceCellArray(dest_cells);
						}
						
						if (indices.size() == 1)
						{
							int index1 = static_cast<int>(indices[0].Scalar() - 1); // remember we want to be 1-based

							if (index1 >= dest_cells->Size())
								throw OML_Error(HW_ERROR_INDEXRANGE);
							else
								(*dest_cells)(index1) = Currency(allocateMatrix());
						}
						else if (indices.size() == 2)
						{
							int index1 = static_cast<int>(indices[0].Scalar() - 1); // remember we want to be 1-based
							int index2 = static_cast<int>(indices[1].Scalar() - 1); // remember we want to be 1-based
							
							if (index1 >= dest_cells->M())
								throw OML_Error(HW_ERROR_INDEXRANGE);
							else if (index2 >= dest_cells->N())
								throw OML_Error(HW_ERROR_INDEXRANGE);
							else
								(*dest_cells)(index1,index2) = Currency(allocateMatrix());
						}
					}
				}
				else if (indices[0].IsColon())
				{
					if (value.IsCellArray())
					{
						HML_CELLARRAY* dest_cells = target.CellArray();

						if (dest_cells->GetRefCount() != 1)
						{
							dest_cells = allocateCellArray(dest_cells);
							target.ReplaceCellArray(dest_cells);
						}

						HML_CELLARRAY* cells = value.CellArray();

						if ((dest_cells->M() == cells->M()) && (cells->N() == 1))
						{
							int col_index = static_cast<int>(indices[1].Scalar()-1);

							if (col_index >= dest_cells->N())
								dest_cells->Resize(dest_cells->M(), col_index+1);

							for (int j=0; j<cells->M(); j++)
								(*dest_cells)(j, col_index) = (*cells)(j, 0);
						}
						else
						{
							throw OML_Error(HW_ERROR_INVRHSSIZE);
						}
					}
				}
				else
				{
					CellAssignmentHelper(target, indices, value);
					return;
				}
			}
			else if (value.IsEmpty())
			{
				if (indices.size() == 1)
				{
					Currency idx1 = indices[0];
					if (idx1.IsPositiveInteger())
					{
						int elem = (int)idx1.Scalar();
						HML_CELLARRAY* temp_mtx   = target.CellArray();
						HML_CELLARRAY* new_matrix = allocateCellArray();

						if (temp_mtx->M() == 1)
						{
							if ((elem > temp_mtx->N()))
								throw OML_Error(HW_ERROR_INDEXRANGE);

							hwMathStatus stat = new_matrix->DeleteColumns(*temp_mtx, elem-1, 1);
						}
						else if (temp_mtx->N() == 1)
						{
							if ((elem > temp_mtx->M()))
								throw OML_Error(HW_ERROR_INDEXRANGE);

							hwMathStatus stat = new_matrix->DeleteRows(*temp_mtx, elem-1, 1);
						}

						target.ReplaceCellArray(new_matrix);
						return;
					}
					else
					{
						throw OML_Error(HW_ERROR_INVIND);
					}
				}
				else if (indices.size() == 2)
				{
					if (indices[0].IsColon())
					{
						if (indices[1].IsPositiveInteger())
						{
							int col = (int)indices[1].Scalar();
							HML_CELLARRAY* temp_mtx   = target.CellArray();
							HML_CELLARRAY* new_matrix = allocateCellArray();

							if ((col > temp_mtx->N()))
								throw OML_Error(HW_ERROR_INDEXRANGE);

							hwMathStatus stat = new_matrix->DeleteColumns(*temp_mtx, col-1, 1);
							target.ReplaceCellArray(new_matrix);
							return;
						}
						else
						{
							throw OML_Error(HW_ERROR_INVIND);
						}
					}
					else if (indices[1].IsColon())
					{
						if (indices[0].IsPositiveInteger())
						{
							int row = (int)indices[0].Scalar();
							HML_CELLARRAY* temp_mtx   = target.CellArray();
							HML_CELLARRAY* new_matrix = allocateCellArray();

							if ((row > temp_mtx->M()))
								throw OML_Error(HW_ERROR_INDEXRANGE);

							hwMathStatus stat = new_matrix->DeleteRows(*temp_mtx, row-1, 1);
							target.ReplaceCellArray(new_matrix);
							return;
						}
						else
						{
							throw OML_Error(HW_ERROR_INVIND);
						}
					}
					else
					{
						throw OML_Error(HW_ERROR_INVCELLDEL);
					}
				}
			}
			else
			{
				CellAssignmentHelper(target, indices, value);
				return;
			}
		}
		else if (target.IsStruct())
		{
			if (value.IsStruct())
			{
				StructData* sd  = target.Struct();
				StructData* sd2 = value.Struct();

				if (indices.size() == 1)
				{
					int index1 = static_cast<int>(indices[0].Scalar());

					if (index1 <= 0)
						throw OML_Error("Invalid index");

					if (index1 >= sd->M()*sd->N())
						sd->DimensionNew(index1);
					sd->SetElement(index1, -1, sd2);
				}
				else if (indices.size() == 2)
				{
					int index1 = static_cast<int>(indices[0].Scalar());
					int index2 = static_cast<int>(indices[1].Scalar());

					if ((index1 >= sd->M()) || (index2 >= sd->N()))
					{
						int num_rows = sd->M();
						int num_cols = sd->N();

						if (index1 >= sd->M())
							num_rows = index1;

						if (index2 >= sd->N())
							num_cols = index2;

						sd->DimensionNew(num_rows, num_cols);
					}
					sd->SetElement(index1, index2, sd2);
				}
			}
			else
			{
				throw OML_Error("Invalid assignment");
			}
		}
		else if (target.IsObject())
		{
			if (value.IsObject())
			{
				StructData* sd  = target.Struct();
				StructData* sd2 = value.Struct();

				int old_size = sd->Size();

				if (indices.size() == 1)
				{
					int index1 = static_cast<int>(indices[0].Scalar());

					if (index1 >= sd->M()*sd->N())
						sd->DimensionNew(index1);
					sd->SetElement(index1, -1, sd2);
				}
				else if (indices.size() == 2)
				{
					int index1 = static_cast<int>(indices[0].Scalar());
					int index2 = static_cast<int>(indices[1].Scalar());

					if ((index1 >= sd->M()) || (index2 >= sd->N()))
					{
						int num_rows = sd->M();
						int num_cols = sd->N();

						if (index1 >= sd->M())
							num_rows = index1;

						if (index2 >= sd->N())
							num_cols = index2;

						sd->DimensionNew(num_rows, num_cols);
					}
					sd->SetElement(index1, index2, sd2);
				}

				// start with old_size
				for (int j=1; j<=sd->Size(); j++)
				{
					StructData* temp = sd->GetElement(j, -1);
					
					if (temp->IsEmpty())
					{
						std::string class_name = target.GetClassname();
						std::vector<Currency> empty_params;
						Currency default_inst = CallFunction(class_name, empty_params);

						if (default_inst.IsObject())
							sd->SetElement(j, -1, default_inst.Struct());
					}
				}
			}
		}
		else if (target.IsScalar() || target.IsComplex())
		{
			hwMatrix* mtx = allocateMatrix();

			mtx->Dimension(1, hwMatrix::REAL);

			if (target.IsComplex())
			{
				mtx->MakeComplex();
				mtx->z(0) = target.Complex();
			}
			else
			{
				(*mtx)(0) = target.Scalar();
			}

			if (num_indices == 1)
			{
				 if (indices[0].IsEmpty())
				 {
					 target.ReplaceMatrix(mtx);
					 return;
				 }
				 else if (indices[0].IsPositiveInteger())
				 {
					 int idx = (int)indices[0].Scalar();

					 if (value.IsComplex())
					 {
						 ReplaceMatrixElementHelper(mtx, idx-1, value.Complex());
					 }
					 else if (value.IsScalar())
					 {
						 ReplaceMatrixElementHelper(mtx, idx-1, value.Scalar());
					 }
					 else if (value.IsEmpty())
					 {
						 mtx->DeleteColumns(0);
					 }
					 else
					 {
						 throw OML_Error(HW_ERROR_INVRHS);
					 }

					 target.ReplaceMatrix(mtx);
					 return;
				 }
				 else if (indices[0].IsLogical())
				 {
					 int idx = (int)indices[0].Scalar();

					 if (!idx)
						 return;

					 if (value.IsComplex())
					 {
						 ReplaceMatrixElementHelper(mtx, idx-1, value.Complex());
					 }
					 else if (value.IsScalar())
					 {
						 ReplaceMatrixElementHelper(mtx, idx-1, value.Scalar());
					 }
					 else
					 {
						 throw OML_Error(HW_ERROR_INVRHS);
					 }

					 target.ReplaceMatrix(mtx);
					 return;
				 }
				 else if (indices[0].IsPositiveVector())
				 {
					 if (value.IsVector())
					 {
						 const hwMatrix* vals = value.Matrix();
						 std::vector<double> ind = indices[0].Vector();

						 if (vals->Size() == ind.size())
						 {
							 if (vals->IsReal())
							 {
								 for (int j=0; j<vals->Size(); j++)
									 ReplaceMatrixElementHelper(mtx, (int)(ind[j]-1), (*vals)(j));
							 }
							 else
							 {
								 for (int j=0; j<vals->Size(); j++)
									 ReplaceMatrixElementHelper(mtx, (int)(ind[j]-1), vals->z(j));
							 }

							 target.ReplaceMatrix(mtx);
							 return;
						 }
						 else
						 {
							 throw OML_Error(HM_ERROR_INDSAMESIZE);
						 }
					 }
				 }
				 else
				 {
					 throw OML_Error(HW_ERROR_INVIND);
				 }
			 }
			 else if (num_indices == 2)
			 {
				 if (indices[0].IsPositiveInteger() && indices[1].IsPositiveInteger())
				 {
					 int idx1 = (int)indices[0].Scalar();
					 int idx2 = (int)indices[1].Scalar();

					 if (value.IsComplex())
					 {
						 ReplaceMatrixElementHelper(mtx, idx1-1, idx2-1, value.Complex());
					 }
					 else if (value.IsScalar())
					 {
						 ReplaceMatrixElementHelper(mtx, idx1-1, idx2-1, value.Scalar());
					 }
					 else
					 {
						 throw OML_Error(HW_ERROR_INVRHS);
					 }

					 target = Currency(mtx);
					 return;
				 }
				 else if (indices[0].IsColon() && indices[1].IsScalar())
				 {
					 int idx2 = (int)indices[1].Scalar();

					 if (value.IsScalar())
					 {
						 ReplaceMatrixElementHelper(mtx, 0, idx2-1, value.Scalar());
					 }
					 else if (value.IsComplex())
					 {
						 ReplaceMatrixElementHelper(mtx, 0, idx2-1, value.Complex());
					 }
					 else if (value.IsMatrix())
					 {
						 if (value.IsEmpty())
						 {
							 if (idx2 == 1)
							 {
								 if (target.IsScalar())
							         mtx = allocateMatrix(0, 1, hwMatrix::REAL);
								 else
									 mtx = allocateMatrix(0, 1, hwMatrix::COMPLEX);
							 }
							 else
							 {
								 throw OML_Error(HW_ERROR_INDEXRANGE);						 
							 }
						 }
						 else
						 {
							throw OML_Error(HW_ERROR_NOTIMP);
						 }
					 }

					 target = Currency(mtx);
					 return;
				 }
				 else if (indices[1].IsColon() && indices[0].IsScalar())
				 {
					 int idx1 = (int)indices[0].Scalar();

					 if (value.IsScalar())
					 {
						 ReplaceMatrixElementHelper(mtx, idx1-1, 0, value.Scalar());
					 }
					 else if (value.IsComplex())
					 {
						 ReplaceMatrixElementHelper(mtx, idx1-1, 0, value.Complex());
					 }
					 else if (value.IsMatrix())
					 {
						 if (value.IsEmpty())
						 {
							 if (idx1 == 1)
							 {
								 if (target.IsScalar())
							         mtx = allocateMatrix(1, 0, hwMatrix::REAL);
								 else
									 mtx = allocateMatrix(1, 0, hwMatrix::COMPLEX);
							 }
							 else
							 {
								 throw OML_Error(HW_ERROR_INDEXRANGE);						 
							 }
						 }
						 else
						 {
							throw OML_Error(HW_ERROR_NOTIMP);
						 }
					 }

					 target = Currency(mtx);
					 return;
				 }				 
				 else
				 {
					throw OML_Error(HW_ERROR_NOTIMP);
				 }
			 }
			 else if (num_indices > 2)
			 {
				 throw OML_Error(HW_ERROR_NOTIMP);
			 }
		}
		else
		{
			throw OML_Error(HW_ERROR_NOTIMP);
		}
	}

	return;
}
//------------------------------------------------------------------------------
//! Utility to convert 2D matrix to ND matrix
//! \param[in] mtx_in Given input matrix
//------------------------------------------------------------------------------
hwMatrixN* ExprTreeEvaluator::Convert2DtoND(const hwMatrix* mtx_in)
{
	if (mtx_in)
	{
		std::vector<int> dims;

		dims.push_back(mtx_in->M());
		dims.push_back(mtx_in->N());

		int size = mtx_in->Size();

		if (mtx_in->IsReal())
		{
			hwMatrixN* real =  ExprTreeEvaluator::allocateMatrixN(dims, hwMatrixN::REAL);

			for (int j=0; j<size; j++)
				(*real)(j) = (*mtx_in)(j);

			return real;
		}
		else
		{
			hwMatrixN* imag =  ExprTreeEvaluator::allocateMatrixN(dims, hwMatrixN::COMPLEX);

			for (int j=0; j<size; j++)
				imag->z(j) = mtx_in->z(j);

			return imag;
		}
	}
	else
	{
		return ExprTreeEvaluator::allocateMatrixN();
	}
}

void ExprTreeEvaluator::NDAssignmetHelper(Currency& target, const std::vector<Currency>& params, const Currency& value)
{
	int num_indices = (int)params.size();
	int num_colons = 0;
		
	hwMatrixN* LHS = NULL;

	bool       need_assign = false;

	if (target.IsMatrix())
	{
		LHS = Convert2DtoND(target.Matrix());
		need_assign = true;
	}
	else if (target.IsNDMatrix())
	{
		LHS = target.GetWritableMatrixN();
		if (LHS->GetRefCount() != 1)
		{
			LHS = allocateMatrixN(LHS);
			need_assign = true;
		}
	}
	else if (target.IsScalar())
	{
		std::vector<int> dims;
		dims.push_back(1);
		LHS = allocateMatrixN(dims, hwMatrixN::REAL);
		(*LHS)(0) = target.Scalar();
		need_assign = true;
	}
	else
	{
		throw OML_Error(OML_MSG_ARRAYSIZE);
	}

    if (num_indices == 1)
    {
        if (params[0].IsMatrix())
        {
            const hwMatrix* idxm = params[0].Matrix();  // indices
            int size_idx = idxm->Size();
            int size_lhs = LHS->Size();

            if (!idxm->IsReal())
                throw OML_Error(HW_ERROR_INDEXPOSINT);

            if (value.IsMatrix())
            {
                const hwMatrix* RHS = value.Matrix();

                if (idxm->M() != RHS->M() || idxm->N() != RHS->N())
                    throw OML_Error(HW_ERROR_INVINDOP);

                if (!RHS->IsReal() && LHS->IsReal())
                    LHS->MakeComplex();

                for (int i = 0; i < size_idx; ++i)
                {
                    double idxd = (*idxm)(i);
                    int idx = static_cast<int>(idxd) - 1;

                    if (idx < 0 || idx >= size_lhs)
                        throw OML_Error(OML_ERR_INVALID_RANGE);

                    if (LHS->IsReal())
                        (*LHS)(idx) = (*RHS)(i);
                    else
                        LHS->z(idx) = RHS->z(i);
                }
            }
            else if (value.IsScalar())
            {
                double rhs = value.Scalar();

                for (int i = 0; i < size_idx; ++i)
                {
                    double idxd = (*idxm)(i);
                    int idx = static_cast<int>(idxd) - 1;

                    if (idx < 0 || idx >= size_lhs)
                        throw OML_Error(OML_ERR_INVALID_RANGE);

                    if (LHS->IsReal())
                        (*LHS)(idx) = rhs;
                    else
                        LHS->z(idx) = rhs;
                }
            }
            else if (value.IsComplex())
            {
                hwComplex rhs = value.Complex();

                if (LHS->IsReal())
                    LHS->MakeComplex();

                for (int i = 0; i < size_idx; ++i)
                {
                    double idxd = (*idxm)(i);
                    int idx = static_cast<int>(idxd) - 1;

                    if (idx < 0 || idx >= size_lhs)
                        throw OML_Error(OML_ERR_INVALID_RANGE);

                    LHS->z(idx) = rhs;
                }
            }

            if (need_assign)
            {
                int old_mask = target.GetMask();
                target = Currency(LHS);
                target.SetMask(old_mask);
            }

            return;
        }
        else if (params[0].IsNDMatrix())
        {
            const hwMatrixN* idxm = params[0].MatrixN();  // indices
            int size_idx = idxm->Size();
            int size_lhs = LHS->Size();

            if (!idxm->IsReal())
                throw OML_Error(HW_ERROR_INDEXPOSINT);

            if (value.IsNDMatrix())
            {
                const hwMatrixN* RHS = value.MatrixN();

                if (idxm->Dimensions() != RHS->Dimensions())
                    throw OML_Error(HW_ERROR_INVINDOP);

                if (!RHS->IsReal() && LHS->IsReal())
                    LHS->MakeComplex();

                for (int i = 0; i < size_idx; ++i)
                {
                    double idxd = (*idxm)(i);
                    int idx = static_cast<int>(idxd) - 1;

                    if (idx < 0 || idx >= size_lhs)
                        throw OML_Error(OML_ERR_INVALID_RANGE);

                    if (LHS->IsReal())
                        (*LHS)(idx) = (*RHS)(i);
                    else
                        LHS->z(idx) = RHS->z(i);
                }
            }
            else if (value.IsScalar())
            {
                double rhs = value.Scalar();

                for (int i = 0; i < size_idx; ++i)
                {
                    double idxd = (*idxm)(i);
                    int idx = static_cast<int>(idxd) - 1;

                    if (idx < 0 || idx >= size_lhs)
                        throw OML_Error(OML_ERR_INVALID_RANGE);

                    if (LHS->IsReal())
                        (*LHS)(idx) = rhs;
                    else
                        LHS->z(idx) = rhs;
                }
            }
            else if (value.IsComplex())
            {
                hwComplex rhs = value.Complex();

                if (LHS->IsReal())
                    LHS->MakeComplex();

                for (int i = 0; i < size_idx; ++i)
                {
                    double idxd = (*idxm)(i);
                    int idx = static_cast<int>(idxd) - 1;

                    if (idx < 0 || idx >= size_lhs)
                        throw OML_Error(OML_ERR_INVALID_RANGE);

                    LHS->z(idx) = rhs;
                }
            }

            if (need_assign)
            {
                int old_mask = target.GetMask();
                target = Currency(LHS);
                target.SetMask(old_mask);
            }

            return;
        }
    }

    slices.clear();

	for (int j=0; j<num_indices; j++)
	{
		const Currency& temp = params[j];

		if (temp.IsColon())
		{
			slices.push_back(hwSliceArg());
			num_colons++;
		}
        else if (temp.IsScalar())
		{
			slices.push_back((int)temp.Scalar()-1);
		}
		else if (temp.IsPositiveVector())
		{
			std::vector<double> temp_vec = params[j].Vector();
			std::vector<int>    temp_int;

			for (int k=0; k<temp_vec.size(); k++)
				temp_int.push_back((int)(temp_vec[k]-1));
			slices.push_back(hwSliceArg(temp_int));
		}
		else
		{
			throw OML_Error(OML_MSG_INVALID_INDEX);
		}
	}

	if (value.IsScalar())
	{
		try
		{
			LHS->SliceLHS(slices, value.Scalar());
		}
		catch (hwMathException& hme)
		{
            throw OML_Error(hme.what());
		}	

		std::vector<int> dims = LHS->Dimensions();

		if (dims.size() < 3)
			target = Currency(ConvertNDto2D(LHS));
	}
	else if (value.IsComplex())
	{
		try
		{
			LHS->SliceLHS(slices, value.Complex());
		}
		catch (hwMathException& hme)
		{
            throw OML_Error(hme.what());
		}		
	}
	else if (value.IsMatrix())
	{
		hwMatrixN* RHS = Convert2DtoND(value.Matrix());

		try
		{
			LHS->SliceLHS(slices, *RHS);
		}
		catch (hwMathException& hme)
		{
            throw OML_Error(hme.what());
		}

		delete RHS;
	}
	else if (value.IsNDMatrix())
	{
		const hwMatrixN* RHS = value.MatrixN();

		try
		{
			LHS->SliceLHS(slices, *RHS);
		}
		catch (hwMathException& hme)
		{
            throw OML_Error(hme.what());
		}
	}
	else
	{
		throw OML_Error("Invalid RHS");
	}

	if (need_assign)
	{
		int old_mask = target.GetMask();
		target = Currency(LHS);
		target.SetMask(old_mask);
	}
}

Currency ExprTreeEvaluator::RangeOperator(OMLTree* tree)
{
	// It's more complicated than I'd like, but it works.  The root is a colon with 2 children.
	// If the first child is also a colon, then it's the a:b:c variant.  Otherwise just a:b.  
	// If the first child is a colon, and ITS first child is also a colon, then we have an error.  a:b:c:d is not legit.

	int num_statements = tree->ChildCount(); // should be 2
	Currency start;
	Currency end;
	double start_val = 0.0;
	double end_val = 0.0;
	double increment = 1.0;

	// this is a single colon to be used for accessing a complete row or column
	// need to communicate this to the caller somehow
	if (num_statements == 0)
		return Currency(0.0, Currency::TYPE_COLON);

	OMLTree* first_child = tree->GetChild(0);

	int token_type = first_child->GetType();

	if (token_type == COLON)
	{
		OMLTree* first_grandchild = first_child->GetChild(0);

		int gc_token = first_grandchild->GetType();

		if (gc_token == COLON)
			throw OML_Error(HW_ERROR_UNSUPRANGEOP);

		start = RUN(first_grandchild);

		OMLTree* child_1 = first_child->GetChild(1);
		Currency incr = RUN(child_1);

		if (incr.IsScalar())
			increment = incr.Scalar();
		else if (incr.IsComplex())
			increment = incr.Complex().Real();
		else if (incr.IsMatrix() && !incr.Matrix()->IsEmpty())
			increment = (*incr.Matrix())(0);

		OMLTree* child_end = tree->GetChild(1);
		end   = RUN(child_end);
	}
	else
	{
		start = RUN(first_child);

		OMLTree* child_1 = tree->GetChild(1);
		end   = RUN(child_1);
	}

	bool abort = false;

	if (start.IsScalar())
		start_val = start.Scalar();
	else if (start.IsComplex())
		start_val = start.Complex().Real();
	else if (start.IsMatrix() && !start.Matrix()->IsEmpty())
		start_val = (*start.Matrix())(0);
	else if (start.IsCharacter())
		start_val = (unsigned int)start.StringVal()[0];
	else if (start.IsCellArray())
		throw OML_Error(HW_ERROR_INVSTARTVAL);
	else if (start.IsStruct())
		throw OML_Error(HW_ERROR_INVSTARTVAL);
	else
		abort = true;

	if (end.IsScalar())
		end_val = end.Scalar();
	else if (end.IsComplex())
		end_val = end.Complex().Real();
	else if (end.IsMatrix() && !end.Matrix()->IsEmpty())
		end_val = (*end.Matrix())(0);
	else if (end.IsCharacter())
		end_val = (unsigned int)end.StringVal()[0];
	else if (end.IsCellArray())
		throw OML_Error(HW_ERROR_INVENDVAL);
	else if (end.IsStruct())
		throw OML_Error(HW_ERROR_INVENDVAL);
	else
		abort = true;

	if (abort)
		return allocateMatrix(1, 0, hwMatrix::REAL);

	if (IsInf_T(start_val) || IsNaN_T(start_val))
		throw OML_Error(HW_ERROR_INVSTARTVAL);

	if (IsInf_T(end_val) || IsNaN_T(end_val))
		throw OML_Error(HW_ERROR_INVENDVAL);

	std::vector<double> temp;

	// if increment is 0, num_steps will end up being negative and we
	// return an empty matrix like we should
	double num_steps_prime = (end_val - start_val)/increment + 1e-12;

	int num_steps = (int)(num_steps_prime + 1.0);

	// overflow check - if this happens, give up
	if ((num_steps_prime > 0) && (num_steps < 0) && (increment > 0.0))
	{
		throw OML_Error(HW_ERROR_OUTMEM);
	}
	else
	{
		if (IsInf_T(increment))
			num_steps = 0;

		if (num_steps > 0)
			temp.reserve(num_steps);

		if (increment != 1.0)
		{
			for (int j=0; j<num_steps; j++)
				temp.push_back(j*increment+start_val);
		}
		else
		{
			for (int j=0; j<num_steps; j++)
				temp.push_back(j+start_val);
		}
	}

	Currency cur(temp);
	if (start.IsCharacter())
		cur.SetMask(Currency::MASK_STRING);

	return cur;
}

Currency ExprTreeEvaluator::GlobalReference(OMLTree* tree)
{
	// this only handles the first variable - need to fix
	int child_count = tree->ChildCount();

	for (int j=0; j<child_count; j++)
	{
		std::string varname = tree->GetChild(j)->GetText();
		msm->AddGlobalReference(varname);
	}

	return Currency(-1.0, Currency::TYPE_NOTHING);
}
Currency ExprTreeEvaluator::PersistentReference(OMLTree* tree)
{
	// this only handles the first variable - need to fix
	std::string varname = tree->GetChild(0)->GetText();
	msm->AddPersistentReference(varname);

	return Currency(-1.0, Currency::TYPE_NOTHING);
}

Currency ExprTreeEvaluator::TransposeOperator(OMLTree* tree)
{
	OMLTree* child_0 = tree->GetChild(0);
	Currency ret = RUN(child_0);

	if (ret.IsScalar() || ret.IsComplex())
	{
		return ret;
	}
	else if (ret.IsMatrix())
	{
		const hwMatrix* old_mtx = ret.Matrix();
		hwMatrix* new_mtx = allocateMatrix(old_mtx);
		new_mtx->Transpose();
		return new_mtx;
	}
	else if (ret.IsString())
	{
		const hwMatrix* old_mtx = ret.Matrix();
		hwMatrix* new_mtx = allocateMatrix(old_mtx);
		new_mtx->Transpose();
		Currency cur(new_mtx);
		cur.SetMask(Currency::MASK_STRING);
		return cur;
	}
	else
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
	}

	return Currency();
}

Currency ExprTreeEvaluator::ConjTransposeOperator(OMLTree* tree)
{
	OMLTree* child_0 = tree->GetChild(0);
	Currency ret = RUN(child_0);

	if (ret.IsScalar())
	{
		return ret;
	}
	else if (ret.IsComplex())
	{
		hwComplex res = ret.Complex();
		hwComplex new_res = res.Conjugate();
		return new_res;
	}
	else if (ret.IsMatrix())
	{
		const hwMatrix* old_mtx = ret.Matrix();
		hwMatrix* new_mtx = allocateMatrix(old_mtx);
		new_mtx->Transpose();
		new_mtx->Conjugate();
		return new_mtx;
	}
	else if (ret.IsString())
	{
		const hwMatrix* old_mtx = ret.Matrix();
		hwMatrix* new_mtx = allocateMatrix(old_mtx);
		new_mtx->Transpose();
		Currency cur(new_mtx);
		cur.SetMask(Currency::MASK_STRING);
		return cur;
	}
	else if (ret.IsCellArray())
	{
		HML_CELLARRAY* old_mtx = ret.CellArray();
		HML_CELLARRAY* new_mtx = new HML_CELLARRAY(0, 0, HML_CELLARRAY::REAL);
		new_mtx->Transpose(*old_mtx);
		return new_mtx;
	}
	else
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
	}

	return Currency();
}

Currency ExprTreeEvaluator::CellArrayCreation(OMLTree* tree)
{
	HML_CELLARRAY* ret_val = new HML_CELLARRAY;

	unsigned int cell_rows = tree->ChildCount();

	if (!cell_rows) // empty cell
		return ret_val;

	OMLTree* row_tree = tree->GetChild(0);
	unsigned int cell_cols = row_tree->ChildCount();

	std::vector<std::vector<Currency>> currencies;

	// run all the entries and store them
	for (unsigned int j=0; j<cell_rows; j++)
	{
		std::vector<Currency> loc_vec;

		OMLTree* row_tree = tree->GetChild(j);
		size_t num_cols_in_row = row_tree->ChildCount();

		for (unsigned int k=0; k<num_cols_in_row; k++)
		{
			OMLTree* col_tree = row_tree->GetChild(k);
			loc_vec.push_back(RUN(col_tree));
		}

		currencies.push_back(loc_vec);
	}

	int dest_rows   = 0;
	int dest_cols   = 0;
	int target_cols = 0;

	// validate and figure out the final size
	for (unsigned int j=0; j<currencies.size(); j++)
	{
		std::vector<Currency> temp = currencies[j];

		dest_cols = 0;

		for (unsigned int k=0; k<temp.size(); k++)
		{
			Currency temp_cur = temp[k];

			if (temp_cur.IsCellList())
			{
				HML_CELLARRAY* cells = temp_cur.CellArray();
				dest_cols += cells->Size();
			}
			else
			{
				dest_cols++;
			}
		}

		if (j==0)
		{
			target_cols = dest_cols;
		}
		else
		{
			if (target_cols != dest_cols)
				throw OML_Error(HW_ERROR_NOTCREATECELLARRAY);
			}

		dest_rows++;
		}

	ret_val->Dimension(dest_rows, dest_cols, HML_CELLARRAY::REAL);

	int col_idx = 0;

	// put the currencies in the correct place
	for (unsigned int j=0; j<currencies.size(); j++)
	{
		std::vector<Currency> temp = currencies[j];

		col_idx = 0;

		for (unsigned int k=0; k<temp.size(); k++)
		{
			Currency temp_cur = temp[k];

			if (temp_cur.IsCellList())
			{
				HML_CELLARRAY* cells = temp_cur.CellArray();

				for (int m=0; m<cells->Size(); m++)
					(*ret_val)(j, col_idx+m) = (*cells)(m);

				col_idx += cells->Size();
		    }
		    else
		    {
				temp_cur.ClearOutputName();
				(*ret_val)(j, col_idx) = temp_cur;
				col_idx++;
	        }
		}
	}
	
	return ret_val;
}

Currency ExprTreeEvaluator::CellValue(OMLTree* tree)
{
	OMLTree* source = tree->GetChild(0);

	Currency target;

	if (_lhs_eval)
	{
		if (source->GetType() == IDENT)
		{
			const Currency& loc_target = msm->GetValue(source->GetText());
			const Currency* temp = &loc_target;
			target = Currency((Currency*)temp);
		}
		else
		{
			throw OML_Error(HW_ERROR_UNKNOWNLHS);
		}
	}
	else
	{	
		target = RUN(source);
	}

	if (!target.IsCellArray())
	{
		if (!_lhs_eval)
			throw OML_Error(HW_ERROR_NOTCELLINDNONCELL);
	}

	OMLTree* index_token = tree->GetChild(1);
		
	std::string* old_end_context_varname  = end_context_varname;
	Currency*    old_end_context_currency = end_context_currency;
	end_context_varname = NULL;
	end_context_currency = &target;

	int num_indices = index_token->ChildCount();
	std::vector<Currency> indices;

	for (int j=0; j<num_indices; j++)
	{
		if (num_indices == 1)
			end_context_index = -1;
		else
			end_context_index = j;
	
		OMLTree* child_j = index_token->GetChild(j);
		indices.push_back(RUN(child_j));
	}

	HML_CELLARRAY* temp = NULL;
	
	if (target.IsCellArray())
	{
		temp = target.CellArray();
	}
	else if (target.IsPointer())
	{
		Currency* temp_cur = target.Pointer();

		if (temp_cur->IsCellArray())
			temp = temp_cur->CellArray();
	}

	Currency ret_val;

	if (indices.size() == 1)
	{
		Currency idx  = indices[0];

		if (idx.IsScalar())
		{
			int index_1 = static_cast<int>(idx.Scalar())-1;

			if ((index_1 >= temp->Size()) || (index_1 < 0))
				throw OML_Error(HW_ERROR_CELLINDEXRANGE);

			if (_lhs_eval)
				ret_val = &((*temp)(index_1));
			else
				ret_val = (*temp)(index_1);
		}
		else if (idx.IsMatrix())
		{
			if (idx.GetMask() == Currency::MASK_DOUBLE)
			{
				const hwMatrix*      mtx = idx.Matrix();
				HML_CELLARRAY* ret = new HML_CELLARRAY;
				ret->Dimension(mtx->Size(), 1 , HML_CELLARRAY::REAL);

				int max_index = temp->Size();

				for (int j=0; j<mtx->Size(); j++)
				{
					int mtx_index = (int)((*mtx)(j)-1);

					if ((mtx_index < 0) || (mtx_index >= max_index))
					{
						delete ret;
						throw OML_Error(HW_ERROR_CELLINDEXRANGE);
					}

					(*ret)(j) = (*temp)((int)((*mtx)(j)-1));
				}

				ret_val = Currency(ret);
				ret_val.SetMask(Currency::MASK_CELL_LIST);
			}
		}
		else if (idx.IsColon())
		{
			ret_val = target;
			ret_val.SetMask(Currency::MASK_CELL_LIST);
		}
		else
		{
			throw OML_Error(HW_ERROR_INVCELLINDTYPE);
		}
	}
	else if (indices.size() == 2)
	{
		HML_CELLARRAY* temp = target.CellArray();
			
		if (indices[0].IsScalar() && indices[1].IsScalar())
		{
			end_context_index = 0;
			int index_1 = static_cast<int>(indices[0].Scalar())-1;
			end_context_index = 1;
			int index_2 = static_cast<int>(indices[1].Scalar())-1;

			if ((index_1 >= temp->M()) || (index_1 < 0))
				throw OML_Error(HW_ERROR_CELLINDEXRANGE);

			if ((index_2 >= temp->N()) || (index_2 < 0))
				throw OML_Error(HW_ERROR_CELLINDEXRANGE);

			ret_val = (*temp)(index_1, index_2);
		}
		else if (indices[0].IsScalar() && indices[1].IsColon())
		{
			int index_1 = static_cast<int>(indices[0].Scalar())-1;
			HML_CELLARRAY* ret = allocateCellArray();
			hwMathStatus stat = temp->ReadRow(index_1, *ret);

			if (stat.IsOk())
			{
				Currency cur(ret);
				cur.SetMask(Currency::MASK_CELL_LIST);
				return cur;
			}
			else
			{
				throw OML_Error(HW_ERROR_INVCELLIND);
			}
		}
		else if (indices[1].IsScalar() && indices[0].IsColon())
		{
			int index_2 = static_cast<int>(indices[1].Scalar())-1;
			HML_CELLARRAY* ret = allocateCellArray();
			hwMathStatus stat = temp->ReadColumn(index_2, *ret);

			if (stat.IsOk())
			{
				Currency cur(ret);
				cur.SetMask(Currency::MASK_CELL_LIST);
				return cur;
			}
			else
			{
				throw OML_Error(HW_ERROR_INVCELLIND);
			}
		}
		else if (indices[0].IsColon() && indices[1].IsColon())
		{
			HML_CELLARRAY* ret = allocateCellArray(temp);
			Currency cur(ret);
			cur.SetMask(Currency::MASK_CELL_LIST);
			return cur;
		}
		else
		{
			throw OML_Error(HW_ERROR_INVCELLIND);
		}
	}

	end_context_currency = old_end_context_currency;
	end_context_varname  = old_end_context_varname;

	return ret_val;
}

/*
Currency ExprTreeEvaluator::CellValue(OMLTree* tree)
{
	std::string var_name(tree->GetChild(0)->GetText());
	int num_index_groups = tree->ChildCount() - 1;

	Currency outer_target;

	std::string* old_end_context_varname = end_context_varname;

	for (int k=0; k<num_index_groups; k++)
	{	
		OMLTree* index_token = tree->GetChild(k+1);
		int num_indices = index_token->ChildCount();

		if (k == 0)
		{
			const Currency& target = msm->GetValue(var_name);

			if (_lhs_eval)
			{
				const Currency* temp = &target;
				outer_target = Currency((Currency*)temp);
			}
			else
			{
				outer_target = target;
			}
		}

		if (outer_target.IsNothing())
			throw OML_Error(HW_ERROR_UNKNOWNVAR);

		std::vector<Currency> indices;
		
		end_context_varname = (std::string*)Currency::vm.GetStringPointer(var_name);

		for (int j=0; j<num_indices; j++)
		{
			if (num_indices == 1)
				end_context_index = -1;
			else
				end_context_index = j;
	
			OMLTree* child_j = index_token->GetChild(j);
			indices.push_back(RUN(child_j));
		}
	
		if (outer_target.IsPointer())
			outer_target = CellValueHelper(*outer_target.Pointer(), indices);
		else
			outer_target = CellValueHelper(outer_target, indices);
	}

	end_context_varname = old_end_context_varname;

	return outer_target;
}
*/

Currency  ExprTreeEvaluator::CellValueHelper(const Currency& target, const std::vector<Currency>& indices)
{
	Currency ret_val;

	if (target.IsCellArray())
	{
		if (indices.size() == 1)
		{
			end_context_index = -1;

			HML_CELLARRAY* temp = target.CellArray();
			Currency       idx  = indices[0];

			if (idx.IsScalar())
			{
				int index_1 = static_cast<int>(idx.Scalar())-1;

				if ((index_1 >= temp->Size()) || (index_1 < 0))
					throw OML_Error(HW_ERROR_CELLINDEXRANGE);

				if (_lhs_eval)
					ret_val = &((*temp)(index_1));
				else
					ret_val = (*temp)(index_1);
			}
			else if (idx.IsMatrix())
			{
				if (idx.GetMask() == Currency::MASK_DOUBLE)
				{
					const hwMatrix*      mtx = idx.Matrix();
					HML_CELLARRAY* ret = new HML_CELLARRAY;
					ret->Dimension(mtx->Size(), 1 , HML_CELLARRAY::REAL);

					int max_index = temp->Size();

					for (int j=0; j<mtx->Size(); j++)
					{
						int mtx_index = (int)((*mtx)(j)-1);

						if ((mtx_index < 0) || (mtx_index >= max_index))
						{
							delete ret;
							throw OML_Error(HW_ERROR_CELLINDEXRANGE);
						}

						(*ret)(j) = (*temp)((int)((*mtx)(j)-1));
					}

					ret_val = Currency(ret);
					ret_val.SetMask(Currency::MASK_CELL_LIST);
				}
			}
			else if (idx.IsColon())
			{
				ret_val = target;
				ret_val.SetMask(Currency::MASK_CELL_LIST);
			}
			else
			{
				throw OML_Error(HW_ERROR_INVCELLINDTYPE);
			}
		}
		else if (indices.size() == 2)
		{
			HML_CELLARRAY* temp = target.CellArray();
			
			if (indices[0].IsScalar() && indices[1].IsScalar())
			{
				end_context_index = 0;
				int index_1 = static_cast<int>(indices[0].Scalar())-1;
				end_context_index = 1;
				int index_2 = static_cast<int>(indices[1].Scalar())-1;

				if ((index_1 >= temp->M()) || (index_1 < 0))
					throw OML_Error(HW_ERROR_CELLINDEXRANGE);

				if ((index_2 >= temp->N()) || (index_2 < 0))
					throw OML_Error(HW_ERROR_CELLINDEXRANGE);

				ret_val = (*temp)(index_1, index_2);
			}
			else if (indices[0].IsScalar() && indices[1].IsColon())
			{
				int index_1 = static_cast<int>(indices[0].Scalar())-1;
				HML_CELLARRAY* ret = allocateCellArray();
				hwMathStatus stat = temp->ReadRow(index_1, *ret);

				if (stat.IsOk())
				{
					Currency cur(ret);
					cur.SetMask(Currency::MASK_CELL_LIST);
					return cur;
				}
				else
				{
					throw OML_Error(HW_ERROR_INVCELLIND);
				}
			}
			else if (indices[1].IsScalar() && indices[0].IsColon())
			{
				int index_2 = static_cast<int>(indices[1].Scalar())-1;
				HML_CELLARRAY* ret = allocateCellArray();
				hwMathStatus stat = temp->ReadColumn(index_2, *ret);

				if (stat.IsOk())
				{
					Currency cur(ret);
					cur.SetMask(Currency::MASK_CELL_LIST);
					return cur;
				}
				else
				{
					throw OML_Error(HW_ERROR_INVCELLIND);
				}
			}
			else if (indices[0].IsColon() && indices[1].IsColon())
			{
				HML_CELLARRAY* ret = allocateCellArray(temp);
				Currency cur(ret);
				cur.SetMask(Currency::MASK_CELL_LIST);
				return cur;
			}
			else
			{
				throw OML_Error(HW_ERROR_INVCELLIND);
		}
	}
	}
	else
	{
		throw OML_Error(HW_ERROR_NOTCELLINDNONCELL);
	}

	ret_val.ClearOutputName();

	if (ret_val.IsMatrix())
		hwMatrix* junk = ret_val.GetWritableMatrix(); // force this to allocate something, otherwise we'll have to check for NULL matrices in a million places

	return ret_val;
}

Currency ExprTreeEvaluator::CellAssignment(OMLTree* tree)
{
	OMLTree*    first_child = tree->GetChild(0);

	if (first_child->GetType() == IDENT)
	{
		std::string var_name(first_child->GetText());

		const Currency& temp_cur = msm->GetValue(var_name);

		if (temp_cur.IsNothing())
		{
			Currency cur;
			cur.MakeCell();
			msm->SetValue(var_name, cur);
		}

		Currency& target = msm->GetMutableValue(var_name);
		
		int old_assignment_nargs = assignment_nargout;
		
		assignment_nargout = 1;
		OMLTree* child_1 = tree->GetChild(1);
		Currency value = RUN(child_1);

		if (value.IsNothing())
			throw OML_Error(HW_ERROR_MISSRETURNS);

		assignment_nargout = old_assignment_nargs;

		end_context_varname = (std::string*)Currency::vm.GetStringPointer(var_name);

		if (tree->ChildCount() == 3)
		{
			OMLTree* index_token = tree->GetChild(2);
			CellAssignmentHelper(target, index_token, value);
			Currency temp(target);
			temp.SetOutputName(var_name);
			return temp;
		}
		else if (tree->ChildCount() == 4)
		{
			OMLTree* first_index_token = tree->GetChild(2);
			OMLTree* first_index_child = first_index_token->GetChild(0);

			OMLTree* index_token = tree->GetChild(3);
			Currency first_index = RUN(first_index_child);

			if (target.IsCellArray() && first_index.IsPositiveInteger())
			{
				HML_CELLARRAY* cells = target.CellArray();
				int index = (int)(first_index.Scalar()-1);

				if (index < cells->Size())
				{
					Currency& new_target = (*cells)((int)(first_index.Scalar()-1));
					end_context_currency = &new_target;
					end_context_varname  = NULL;
					CellAssignmentHelper(new_target, index_token, value);
					Currency temp(target);
					temp.SetOutputName(var_name);
					return temp;
				}
				else
				{
					throw OML_Error(HW_ERROR_CELLINDEXRANGE);
				}
			}
			else
			{
				throw OML_Error(HW_ERROR_CELLINDEXRANGE);
			}
		}
		else
		{
			throw OML_Error(HW_ERROR_UNSUPOP);
		}
	}
	else // should only be struct
	{
		std::string ret_name;

		if (first_child->GetType() == STRUCT)
		{
			OMLTree*    first_grandchild = first_child->GetChild(0);

			if (first_grandchild->GetType() == IDENT)
				ret_name = first_grandchild->GetText();
		}

		_lhs_eval = true;
		Currency target = RUN(first_child);
		_lhs_eval = false;

		OMLTree* child_1 = tree->GetChild(1);
		Currency value = RUN(child_1);

		if (tree->ChildCount() == 3)
		{
			OMLTree* index_token = tree->GetChild(2);

			if (!target.IsPointer())
				throw OML_Error(HW_ERROR_UNKNOWNLHS);

			CellAssignmentHelper(*target.Pointer(), index_token, value);

			return msm->GetValue(ret_name);
		}
		else if (tree->ChildCount() == 4)
		{
			OMLTree* first_index_token  = tree->GetChild(2);
			OMLTree* first_index_child  = first_index_token->GetChild(0);
			OMLTree* second_index_child = first_index_token->GetChild(1);

			OMLTree* index_token = tree->GetChild(3);

			Currency first_index  = RUN(first_index_child);
			Currency second_index = RUN(second_index_child);

			int first_idx  = (int)first_index.Scalar() - 1;
			int second_idx = (int)second_index.Scalar() - 1;

			if (!target.IsPointer())
				throw OML_Error(HW_ERROR_UNKNOWNLHS);

			Currency* my_cur = target.Pointer();

			if (my_cur->IsCellArray() && first_index.IsPositiveInteger() && second_index.IsPositiveInteger())
			{
				HML_CELLARRAY* cells = my_cur->CellArray();
			
				if ((first_idx >= cells->M()) || (second_idx >= cells->N()))
				{
					int new_m = cells->M();
					int new_n = cells->N();

					if (first_idx >= new_m)
						new_m = first_idx+1;

					if (second_idx >= new_n)
						new_n = second_idx+1;

					cells->Resize(new_m, new_n);
				}

				Currency& new_target = (*cells)(first_idx, second_idx);
				CellAssignmentHelper(new_target, index_token, value);
				Currency temp(target);
				return temp;
			}
			else
			{
				throw OML_Error(HW_ERROR_CELLINDEXRANGE);
			}
		}
		else
		{
			throw OML_Error("Unsupported operation");
		}
	}

}

void ExprTreeEvaluator::CellAssignmentHelper(Currency& target, OMLTree* indices, const Currency& value)
{
	int num_indices = indices->ChildCount();
    std::vector<Currency> params;
	end_context_currency = &target;
	end_context_index = -1;

    for (int i = 0; i < num_indices; ++i)
    {
		if (num_indices > 1)
			end_context_index = i;
		OMLTree* child_i = indices->GetChild(i);
        params.push_back(RUN(child_i));
    }
	end_context_currency = NULL;
    return CellAssignmentHelper(target, params, value);
}

void ExprTreeEvaluator::CellAssignmentHelper(Currency& target, const std::vector<Currency>& params, const Currency& value)
{
    int num_indices = (int)params.size();

	if (num_indices > 2)
		throw OML_Error(HW_ERROR_UNSUPP2DIM);

	if (target.IsMatrix() && target.IsEmpty())
	{
		target.MakeCell();
	}
	else if (target.IsStruct())
	{
		StructData* sd = target.Struct();
		
		if (sd->Size() == 0)
			target.MakeStruct();
	}

	if (!target.IsCellArray())
		throw OML_Error("Invalid cell assignment");

	if (num_indices == 2)
	{
		// need to fix this later also
        if (!params[0].IsScalar() || !params[1].IsScalar())
            throw OML_Error(HW_ERROR_INVINDTYPE);

        if (static_cast<int>(params[0].Scalar()) <= 0)
            throw OML_Error(OML_ERR_INVALID_INDEX, 1, OML_VAR_INDEX);
                
        if (static_cast<int>(params[1].Scalar()) <= 0 )
            throw OML_Error(OML_ERR_INVALID_INDEX, 2, OML_VAR_INDEX);

		HML_CELLARRAY* temp = target.CellArray();

		if (temp->GetRefCount() != 1)
		{
			temp = allocateCellArray(temp);
			target.ReplaceCellArray(temp);
		}
		
		end_context_index = 0;
        int index_1 = static_cast<int>(params[0].Scalar())-1;
		end_context_index = 1;
        int index_2 = static_cast<int>(params[1].Scalar())-1;

		if ((index_1 >= temp->M()) || (index_2 >= temp->N()))
		{
			int new_m = temp->M();
			if (index_1 >= new_m)
				new_m = index_1+1;

			int new_n = temp->N();
			if (index_2 >= new_n)
				new_n = index_2+1;

			if (temp->IsEmpty())
				temp->Dimension(new_m, new_n, HML_CELLARRAY::REAL);
			else
			temp->Resize(new_m, new_n);
		}

		(*temp)(index_1, index_2) = value;
		(*temp)(index_1, index_2).ClearOutputName();
	}
	else if (num_indices == 1)
	{
        if (params[0].IsScalar())
		{
			HML_CELLARRAY* temp = target.CellArray();

			if (temp->GetRefCount() != 1)
			{
				temp = allocateCellArray(temp);
				target.ReplaceCellArray(temp);
			}

			end_context_index = -1;
			int index_1 = static_cast<int>(params[0].Scalar())-1;

			if (index_1 >= temp->Size())
			{
				if (temp->IsEmpty())
					temp->Dimension(1, index_1+1, HML_CELLARRAY::REAL); 

				if (temp->M() == 1)
					temp->Resize(1, index_1+1);
				else if (temp->N() == 1)
					temp->Resize(index_1+1, 1);
				else
					throw OML_Error(HW_ERROR_INVCELLRESIZE);
			}
			else if (index_1 < 0)
			{
				throw OML_Error(HW_ERROR_INVIND);
			}

			value.SetOutputName("");
			(*temp)(index_1) = value;
			(*temp)(index_1).ClearOutputName();
		}
		else if (params[0].IsMatrix())
		{
			HML_CELLARRAY*  temp  = target.CellArray();
			const hwMatrix* index = params[0].Matrix();

			if (value.IsCellArray())
			{
				HML_CELLARRAY*  val  = value.CellArray();

				if (temp->GetRefCount() != 1)
				{
					temp = allocateCellArray(temp);
					target.ReplaceCellArray(temp);
				}

				for (int j=0; j<index->Size(); j++)
				{
					int target_index = (int)((*index)(j)-1);

					if (target_index < 0)
						throw OML_Error(HW_ERROR_INVIND);

					if (target_index >= temp->Size())
					{
						if (temp->IsEmpty())
							temp->Dimension(1, target_index+1, HML_CELLARRAY::REAL);

						if (temp->M() == 1)
							temp->Resize(1, target_index+1, true);
						else if (temp->N() == 1)
							temp->Resize(target_index+1, 1, true);
						else
							throw OML_Error(HW_ERROR_UNSUPOP);
					}

					if (val->Size() == 1)
						(*temp)(target_index) = (*val)(0);
					else if (val->Size() == index->Size())
						(*temp)(target_index) = (*val)(j);
					else
						throw OML_Error(HW_ERROR_UNSUPOP);
				}
			}
			else
			{
				throw OML_Error(HW_ERROR_UNSUPOP);
			}
		}
		else if (params[0].IsColon())
		{
			HML_CELLARRAY* temp = target.CellArray();

			if (temp->GetRefCount() != 1)
			{
				temp = allocateCellArray(temp);
				target.ReplaceCellArray(temp);
			}

			end_context_index = -1;

			value.SetOutputName("");

			for (int j=0; j<temp->Size(); j++)
				(*temp)(j) = value;

			//(*temp)(index_1).ClearOutputName();
		}
		else
		{
			throw OML_Error(HW_ERROR_INVINDTYPE);
		}
	}

	end_context_varname = NULL;
}

Currency ExprTreeEvaluator::InlineIndex(OMLTree* tree)
{
	std::vector<Currency> param_vals;

	int num_func_children = tree->ChildCount();
	if (num_func_children == 2) // func + args
	{
		OMLTree* func_args = tree->GetChild(1);

		OMLTree* child_0 = tree->GetChild(0);
		Currency target = RUN(child_0);
		
		end_context_currency = &target;

		if (func_args->GetType() == PARAM_LIST)
		{
			int num_params = func_args->ChildCount();

			for (int j=0; j<num_params; j++)
			{
				OMLTree* child_j = func_args->GetChild(j);
				param_vals.push_back(RUN(child_j));
			}
		}
	
		Currency temp = VariableIndex(target, param_vals);
		end_context_currency = NULL;
		return temp;
	}
	else if (num_func_children == 1) // func w/o args
	{
		OMLTree* child_0 = tree->GetChild(0);
		Currency target = RUN(child_0);
		Currency temp = VariableIndex(target, param_vals);
		return temp;
	}

	return 0.0;
}

Currency ExprTreeEvaluator::InlineIndexCell(OMLTree* tree)
{
	std::vector<Currency> param_vals;

	int num_func_children = tree->ChildCount();
	if (num_func_children == 2) // func + args
	{
		OMLTree* func_args = tree->GetChild(1);

		OMLTree* child_0 = tree->GetChild(0);
		Currency target = RUN(child_0);

		end_context_currency = &target;

		if (func_args->GetType() == PARAM_LIST)
		{
			int num_params = func_args->ChildCount();

			for (int j=0; j<num_params; j++)
			{
				OMLTree* child_j = func_args->GetChild(j);
				param_vals.push_back(RUN(child_j));
			}
		}

		if (target.IsCellArray())
		{
			Currency temp = CellValueHelper(target, param_vals);
			end_context_currency = NULL;
			return temp;
		}
		else
		{
			throw OML_Error(HW_ERROR_UNSUPOP);
		}
	}

	return 0.0;
}

Currency ExprTreeEvaluator::StructValue(OMLTree* tree)
{
	OMLTree* struct_tree = tree->GetChild(0);
	OMLTree* field_tree  = tree->GetChild(1);

	Currency struct_target = RUN(struct_tree);

	StructData* sd = NULL;
	ClassInfo*  ci = NULL;

	if (struct_target.IsStruct() || struct_target.IsObject())
	{
		sd = struct_target.Struct();

		if (struct_target.IsObject())
			ci = (*class_info_map)[struct_target.GetClassname()];
	}
	else if (struct_target.IsPointer())
	{
		Currency* temp = struct_target.Pointer();

		if (temp->IsStruct() || temp->IsObject())
			sd = temp->Struct();

		if (temp->IsObject())
			ci = (*class_info_map)[temp->GetClassname()];
	}

	if (sd)
	{
		std::string struct_name = struct_tree->GetText();

		if (_lhs_eval)
		{
			if (sd && sd->GetRefCount() != 1)
			{
				StructData* new_sd  = new StructData(*sd);
				struct_target.ReplaceStruct(new_sd);
				SetValue(struct_name, struct_target);
				sd = new_sd;
			}
		}

		std::string field_name;

		if (field_tree->GetType() == IDENT)
		{
			field_name = field_tree->GetText();
		}
		else if (field_tree->GetType() == FIELD)
		{
			OMLTree* child = field_tree->GetChild(0);

			field_name = child->GetText();
		}
		else if (field_tree->GetType() == DYN_FIELD)
		{
			OMLTree* child = field_tree->GetChild(0);
			Currency temp  = RUN(child);

			if (temp.IsString())
				field_name = temp.StringVal();
			else
				throw OML_Error(HW_ERROR_INVFIELDNAME);
		}
		else // to account for field names that are keywords
		{
			field_name = field_tree->GetText();
		}

		if (ci)
		{
			if (ci->IsPropertyPrivate(field_name))
			{
				bool in_class_method = ci->IsClassMethod(msm->GetCurrentScope()->GetFunctionInfo());

				if (!in_class_method)
					throw OML_Error("Unable to access private property");
			}
		}

		if (field_name.size() && !sd->Contains(field_name))
		{
			if (_lhs_eval)
			{
				sd->addField(field_name);
				HML_CELLARRAY* empty_cell = EvaluatorInterface::allocateCellArray();
				sd->SetValue(0, -1, field_name, empty_cell);
			}
			else
			{
				throw OML_Error(HW_ERROR_INVFIELDNAME);
			}
		}

		if (_lhs_eval)
		{
			Currency* ptr = (Currency*)sd->GetPointer(0, -1, field_name);
			return Currency(ptr);
		}
		else
		{
			if (sd->Size() == 1)
			{
				Currency    result = sd->GetValue(0, -1, field_name);
				return result;
			}
			else
			{
				int            size  = sd->Size();
				HML_CELLARRAY* cells = EvaluatorInterface::allocateCellArray(1, size);

				for (int j=0; j<size; j++)
					(*cells)(j) = sd->GetValue(j, -1, field_name);

				Currency result(cells);
				result.SetMask(Currency::MASK_CELL_LIST);
				return result;
			}
		}
	}
	else
	{
		throw OML_Error(HW_ERROR_INPUTSTRUCT);
	}
}

#if JDS
Currency ExprTreeEvaluator::StructValue(OMLTree* tree)
{	
	OMLTree* child_tree = tree->GetChild(0);

	int          num_children = tree->ChildCount();

	OMLTree* field_tree = tree->GetChild(num_children - 1); // always the last child
	OMLTree* index_tree = NULL;

	if (num_children == 3)
		index_tree = tree->GetChild(1);

	if (child_tree->GetType() == IDENT)
	{
		std::string  struct_name = tree->GetChild(0)->GetText();
		const Currency& parent   = GetValue(struct_name);

		if (parent.IsNothing())
		{
			OMLTree* child_0 = tree->GetChild(0);
			Currency temp = RUN(child_0);

			if (temp.IsStruct())
				return StructValueHelper(&temp, index_tree, field_tree);
			else
				throw OML_Error(HW_ERROR_INPUTSTRUCT);
		}
		else
		{
			if (parent.IsStruct() || parent.IsObject())
			{
				StructData* sd = parent.Struct();

				if (_lhs_eval)
				{
					if (sd && sd->GetRefCount() != 1)
					{
						StructData* new_sd  = new StructData(*sd);
						Currency    new_cur = parent;
						new_cur.ReplaceStruct(new_sd);
						SetValue(struct_name, new_cur);
					}
				}
			}

			end_context_varname      = (std::string*)Currency::vm.GetStringPointer(struct_name);
			Currency temp            = StructValueHelper(&parent, index_tree, field_tree);
			end_context_varname      = NULL;
			return temp;
		}
	}
	else
	{
		OMLTree* child_0 = tree->GetChild(0);
		Currency temp = RUN(child_0);
		return StructValueHelper(&temp, index_tree, field_tree);
	}
}
#endif
Currency ExprTreeEvaluator::StructValueHelper(const Currency* parent, OMLTree* indices, OMLTree* field_tree)
{
	if (parent->IsObject() && !_lhs_eval)
	{
		return ObjectMethodCall((Currency*)parent, indices, field_tree);
	}
    else if (parent->IsBoundObject())
	{
        return BoundClassMethod(parent, field_tree);
	}
	else if (parent->IsPointer() && !_lhs_eval)
	{
		Currency* temp = parent->Pointer();

		if (temp->IsObject())
			return ObjectMethodCall(temp, indices, field_tree); 
	}

	int index1      = 0;
	int index2      = -1;
	int num_indices = 0;
	bool cell_index = false;

	Currency temp;
	Currency idx1;
	Currency idx2;

	if (indices)
	{
		int index_token_type = indices->GetType();

		if (index_token_type == CELL_PARAM_LIST)
			cell_index = true;

		num_indices = indices->ChildCount();

		if (num_indices > 0)
		{
			if (num_indices == 1)
				end_context_index = -1;
			else
				end_context_index = 0;

			OMLTree* child_0 = indices->GetChild(0);
			idx1 = RUN(child_0);

			if (idx1.IsScalar())
			{
				index1 = (int)(idx1.Scalar() - 1);
			}
			else if (idx1.IsPositiveVector())
			{
				if (field_tree)
				{
					std::vector<Currency> temp_vec;
					temp_vec.push_back(idx1);

					temp = VariableIndex(*parent, temp_vec);
					parent = &temp;
					num_indices = 0;
				}
			}
			else if (idx1.IsColon())
			{
				// it'll get sorted out below
			}
			else if (!parent->IsFunctionHandle())
			{
				throw OML_Error(OML_ERR_POS_INTEGER_VEC_MTX);
			}
		}

		if (num_indices == 2)
		{
			end_context_index = 1;
			OMLTree* child_1 = indices->GetChild(1);
			idx2 = RUN(child_1);
			index2 = (int)(idx2.Scalar() - 1);
		}
	}

	std::string field_name;
	int         num_field_children = 0;

	if (field_tree)
	{
		num_field_children = field_tree->ChildCount();

		if (field_tree->GetType() == FIELD)
		{
			field_name = field_tree->GetChild(0)->GetText();
		}
		else if (num_field_children)
		{
			if (field_tree->GetType() == DYN_FIELD)
			{
				OMLTree* child_0 = field_tree->GetChild(0);
				Currency temp = RUN(child_0);

				if (temp.IsString())
					field_name = temp.StringVal();
				else
					throw OML_Error(HW_ERROR_INVFIELDNAME);
			}
		}
	}

	StructData *sd = NULL;
	Currency new_parent;

	if (!parent->IsStruct() &&!parent->IsObject() && !parent->IsPointer())
	{
		std::vector<Currency> index_vec;

		end_context_currency = (Currency*)parent;
		end_context_varname  = NULL;

		for (int j=0; j<num_indices; j++)
		{
			OMLTree* child_j = indices->GetChild(j);
			index_vec.push_back(RUN(child_j));
		}

		if (!cell_index)
			new_parent = VariableIndex(*parent, index_vec);
		else
			new_parent = CellValueHelper(*parent, index_vec); // This may not be the best way to go.  We might prefer just to get the cell content from the above line.

		if (!field_tree)
		{
			return new_parent;
		}
		else 
		{
			if (new_parent.IsStruct())
				sd = new_parent.Struct();
			else
				throw OML_Error("Invalid struct access");
		}

		index1 = 0;
		index2 = -1;
	}
	else if (!field_tree)
	{
		std::vector<Currency> index_vec;

		end_context_currency = (Currency*)parent;
		end_context_varname  = NULL;

		for (int j=0; j<num_indices; j++)
		{
			OMLTree* child_j = indices->GetChild(j);
			index_vec.push_back(RUN(child_j));
		}

		return VariableIndex(*parent, index_vec);
	}
	else
	{
		sd = parent->Struct();
	}

	if (parent->IsPointer())
	{
		Currency* ptr = parent->Pointer();

		if (parent->Pointer()->IsObject() || parent->Pointer()->IsStruct())
			sd = parent->Pointer()->Struct();
	}

	if (field_name.size() && !sd->Contains(field_name))
	{
		if (_lhs_eval)
		{
			sd->addField(field_name);
			HML_CELLARRAY* empty_cell = EvaluatorInterface::allocateCellArray();
			sd->SetValue(0, -1, field_name, empty_cell);
		}
		else
		{
			throw OML_Error(HW_ERROR_INVFIELDNAME);
		}
	}

	// need an exception here if num_field_children == 2 && child(2) is an empty PARAM_LIST token
	// and even that's no good for a LHS
	// example a.b() should be equivalent to a.b, but only for a RHS
	if (num_field_children == 1)
	{
		if (_lhs_eval)
		{
			Currency* ptr = (Currency*)sd->GetPointer(index1, index2, field_name);
			return Currency(ptr);
		}
		else if (num_indices && !idx1.IsColon() && !idx2.IsColon())
		{
			if (field_name.size())
			{
				return sd->GetValue(index1, index2, field_name); 
			}
			else
			{
				if (index2 == -1)
					return sd->GetElement(index1+1, -1);
				else
					return sd->GetElement(index1+1, index2+1);
			}
		}
		else if (num_indices == 0)
		{
			int num_elements = sd->M()*sd->N();

			if (num_elements == 1)
			{
				return sd->GetValue(index1, index2, field_name);
			}
			else
			{
				HML_CELLARRAY* cell_list = allocateCellArray(num_elements, 1);

				for (int j=0; j<num_elements; j++)
					(*cell_list)(j) = sd->GetValue(j, -1, field_name);

				Currency ret(cell_list);
				ret.SetMask(Currency::MASK_CELL_LIST);
				return ret;
			}
		}
		else if ((num_indices == 1) && (idx1.IsColon()))
		{
			int num_elements = sd->M()*sd->N();

			if (num_elements == 1)
			{
				return sd->GetValue(index1, index2, field_name);
			}
			else
			{
				HML_CELLARRAY* cell_list = allocateCellArray(num_elements, 1);

				for (int j=0; j<num_elements; j++)
					(*cell_list)(j) = sd->GetValue(j, -1, field_name);

				Currency ret(cell_list);
				ret.SetMask(Currency::MASK_CELL_LIST);
				return ret;
			}
		}
		else if (num_indices == 2)
		{
			if (idx1.IsColon())
			{
				int num_elements = sd->N();

				HML_CELLARRAY* cell_list = allocateCellArray(1, num_elements);

				for (int j=0; j<num_elements; j++)
					(*cell_list)(j) = sd->GetValue(j, index2, field_name);

				Currency ret(cell_list);
				ret.SetMask(Currency::MASK_CELL_LIST);
				return ret;
			}
			else if (idx2.IsColon())
			{
				int num_elements = sd->M();

				HML_CELLARRAY* cell_list = allocateCellArray(num_elements, 1);

				for (int j=0; j<num_elements; j++)
					(*cell_list)(j) = sd->GetValue(index1, j, field_name);

				Currency ret(cell_list);
				ret.SetMask(Currency::MASK_CELL_LIST);
				return ret;
			}
			else
			{
				throw OML_Error("Invalid struct index");
			}
		}
	}
	else
	{
		if (num_field_children == 2)
		{
			OMLTree*    index_tree = field_tree->GetChild(1);

			if (index_tree->GetType() == PARAM_LIST)
			{
				if (index_tree->ChildCount() == 0)
				{
					int num_elements = sd->M()*sd->N();

					if (num_elements == 1)
					{
						Currency cur = sd->GetValue(index1, index2, field_name);

						if (!cur.IsFunctionHandle())
							return cur;
					}
					else
					{
						HML_CELLARRAY* cell_list = allocateCellArray(num_elements, 1);

						for (int j=0; j<num_elements; j++)
							(*cell_list)(j) = sd->GetValue(j, -1, field_name);

						Currency ret(cell_list);
						ret.SetMask(Currency::MASK_CELL_LIST);
						return ret;
					}
				}
			}
		}

		if (field_name.size() && !sd->Contains(field_name))
			throw OML_Error(OML_ERR_NUMARGIN);

		Currency new_parent = sd->GetValue(index1, index2, field_name);

		OMLTree*    index_tree = NULL;
		OMLTree*    field_arg  = NULL;

		if (field_tree)
		{
			index_tree = field_tree->GetChild(1);

			OMLTree*    last_child = field_tree->GetChild(num_field_children-1);

			if ((last_child->GetType() == FIELD) || (last_child->GetType() == DYN_FIELD))
				field_arg = last_child;
		}

		if (index_tree)
		{
			int index_type = index_tree->GetType();
			if ((index_type != PARAM_LIST) && (index_type != CELL_PARAM_LIST))
				index_tree = NULL;
		}

		return StructValueHelper(&new_parent, index_tree, field_arg);
	}

	return Currency();
}

Currency ExprTreeEvaluator::StructAssignment(OMLTree* tree)
{
	// tree is of type '='
	int num_children = tree->ChildCount();

	// lhs is of type struct
	OMLTree* lhs = tree->GetChild(0);
	OMLTree* rhs = tree->GetChild(1);

	int num_lhs_children = lhs->ChildCount();

	OMLTree* lhs_child = lhs->GetChild(0);

	Currency* alt_target = NULL;

	if (lhs_child->GetType() == CELL_VAL)
	{
		_lhs_eval = true;
		Currency lhs_cur = RUN(lhs_child);
		_lhs_eval = false;
		
		if (lhs_cur.IsPointer())
		{
			alt_target = lhs_cur.Pointer();

			if (alt_target->IsCellArray() || alt_target->IsEmpty())
				alt_target->ConvertToStruct(); // converting an empty cell into an empty struct
		}
		else
		{
			throw OML_Error("Unsupported operation");
		}
	}

	std::string struct_name = lhs_child->GetText();

	const Currency& temp_cur = msm->GetValue(struct_name);

	if (temp_cur.IsNothing())
	{
		StructData* sd = new StructData;
		msm->SetValue(struct_name, sd);
	}

	Currency target  = msm->GetValue(struct_name);

	bool is_handle = false;

	if (target.IsScalar() && HasBuiltin("ishandle")) // special case for plots and such
	{
		std::vector<Currency> dummy_input;
		dummy_input.push_back(target);
		Currency curhandle = CallFunction("ishandle", dummy_input);

		if (curhandle.IsScalar() && (curhandle.Scalar() == 1.0))
			is_handle = true;
	}

	if (target.IsMatrix() && target.IsEmpty())
	{
		target.MakeStruct();
	}
	else if (target.IsCellArray())
	{
		HML_CELLARRAY* cells = target.CellArray();
		
		if (cells->Size() == 0)
			target.MakeStruct();
	}

	if (!(target.IsStruct() || target.IsObject() || target.IsPointer()))
	{
		if (!is_handle)
			throw OML_Error("Invalid struct assignment");
	}

	OMLTree* indices = NULL;
	OMLTree* field   = NULL;

	for (int j=1; j<num_lhs_children; j++)
	{
		OMLTree* temp   = lhs->GetChild(j);

		if (temp->GetType() == PARAM_LIST)
			indices = temp;
		else if (temp->GetType() == FIELD)
			field = temp;
		else if (temp->GetType() == DYN_FIELD)
			field = temp;
	}

	Currency value = RUN(rhs);
	value.SetOutputName("");

	// do this after calling run on the rhs
	end_context_varname = (std::string*)Currency::vm.GetStringPointer(struct_name);
	end_context_index   = -1; // for now

    if (temp_cur.IsBoundObject())  // Swig bound class
    {
        return BoundClassAssign(&temp_cur, field, value);
    }

	if (alt_target)
	{
		StructAssignmentHelper(alt_target, indices, field, value);
		end_context_varname = NULL;
		return *alt_target;
	}
	else
	{
		if (target.IsStruct() || target.IsObject())
		{
			StructData* sd = target.Struct();

			if (sd->GetRefCount() != 2) // one stored in the memory scope and one in target
			{
				StructData* new_sd = new StructData(*sd);
				target.ReplaceStruct(new_sd);
			}
		}

		StructAssignmentHelper(&target, indices, field, value);
		msm->SetValue(struct_name, target);
	}

	end_context_varname = NULL;
	return target;
}
//------------------------------------------------------------------------------
//! Gets result after calling a method in a bound class
//! \param[in] in    Input currencty
//! \param[in] field Field info
//------------------------------------------------------------------------------
Currency ExprTreeEvaluator::BoundClassMethod(const Currency*   in, 
                                             OMLTree* field)
{  
    // Sanity checks
    if (!in) 
        return allocateMatrix();
	
	BoundClassInfo* info = _boundclassinfo[in->GetClassname()];
    if (!info) 
        return allocateMatrix();
    
    std::string methodname (field->GetChild(0)->GetText());
    FUNCPTR     fptr = info->GetMethod(methodname);
    if (!fptr)  
        return allocateMatrix();

	std::vector<Currency> inputs;   // Pack the inputs
    inputs.push_back(*in);          // First argument should be the bound object

	if (field->ChildCount() > 1)
	{
		OMLTree* args_tree = field->GetChild(1);
		int num_args = args_tree->ChildCount();

        // Pack rest of arguments needed for the bound method
		for (int j = 0; j < num_args; ++j)
		{
			OMLTree* child_j = args_tree->GetChild(j);
			inputs.push_back(RUN(child_j));
		}
	}

    return CallBuiltinFunction(fptr, methodname, inputs);
}
//------------------------------------------------------------------------------
//! Gets result after assignment to a bound class/bound class property
//! \param[in] in     Input currency
//! \param[in] field  Field tree
//! \param[in] value  Value to assign
//------------------------------------------------------------------------------
Currency ExprTreeEvaluator::BoundClassAssign(const Currency*   in, 
                                             OMLTree* field, 
                                             const Currency&   value)
{    	
    if (!in) 
        return allocateMatrix();

	BoundClassInfo* info = _boundclassinfo[in->GetClassname()];
    if (!info)
        return allocateMatrix();

    std::string propname (field->GetChild(0)->GetText());
    FUNCPTR fptr = info->GetPropertySetter(propname);
    if (!fptr) 
        return allocateMatrix();

	std::vector<Currency> inputs;   // Pack the inputs
	inputs.push_back(*in);          // First arg will always be bound object
    inputs.push_back(value);      

    return CallBuiltinFunction(fptr, propname, inputs);
}

void ExprTreeEvaluator::StructAssignmentHelper(Currency* parent, OMLTree* indices, OMLTree* field_tree, const Currency& rhs)
{
	StructData* sd = NULL;

	if (parent->IsStruct() || parent->IsObject())
	{
		sd = parent->Struct();
	}
	else if (parent->IsPointer() && parent->Pointer()->IsObject())
	{
		sd = parent->Pointer()->Struct();
	}
	else if (parent->IsPointer() && parent->Pointer()->IsStruct())
	{
		sd = parent->Pointer()->Struct();
	}
	else if (parent->IsCellArray())
	{
		HML_CELLARRAY* cells = parent->CellArray();

		if (cells->IsEmpty())
		{
			cells->Dimension(1, 1, HML_CELLARRAY::REAL);
			(*cells)(0) = new StructData;
			sd = (*cells)(0).Struct();
		}
	}

	std::string field_name         = "";
	int         num_field_children = 0;

	if (field_tree)
	{
		if (field_tree->GetType() == FIELD)
		{
			field_name = field_tree->GetChild(0)->GetText();
		}
		else
		{
			OMLTree* child_0 = field_tree->GetChild(0);
			Currency field_cur = RUN(child_0);

			if (field_cur.IsString())
				field_name = field_cur.StringVal();
			else
				throw OML_Error(HW_ERROR_DYNFIELD);
		}

 		num_field_children = field_tree->ChildCount();
	}

	if (parent->IsObject() || (parent->IsPointer() && parent->Pointer()->IsObject()))
	{
		if (!sd->Contains(field_name))
			throw OML_Error("Cannot add fields to class at run time");

		if (parent->IsObject())
		{
			ClassInfo* ci = (*class_info_map)[parent->GetClassname()];

			if (ci->IsPropertyPrivate(field_name))
				throw OML_Error("Unable to access private property");
		}
	}

	int num_indices = 0;

	if (indices)
		num_indices = indices->ChildCount();

	int index_1 = 0;
	int index_2 = -1;
	
	if (num_indices == 2)
	{
		end_context_index = 0;
		OMLTree* child_0 = indices->GetChild(0);
		Currency idx1 = RUN(child_0);
		end_context_index = 1;
		OMLTree* child_1 = indices->GetChild(1);
		Currency idx2 = RUN(child_1);

		index_1 = (int)(idx1.Scalar()-1);
		index_2 = (int)(idx2.Scalar()-1);
	}
	else if (num_indices == 1)
	{
		OMLTree* child_0 = indices->GetChild(0);
		Currency idx1 = RUN(child_0);

		index_1 = (int)(idx1.Scalar()-1);
	}
	else if (num_indices == 0)
	{
        if (parent->IsStruct())
        {
		    if ((parent->Struct()->M() > 1) || (parent->Struct()->N() > 1))
			    throw OML_Error(HW_ERROR_UNSUPOP);
        }
	}

	if (num_field_children == 1) // just the name
	{
		if (parent->IsScalar() && HasBuiltin("ishandle")) // special case for plots and such
		{
			std::vector<Currency> dummy_input;
			dummy_input.push_back(*parent);
			Currency is_handle = CallFunction("ishandle", dummy_input);

			if (is_handle.IsScalar() && (is_handle.Scalar() == 1.0))
			{
				std::vector<Currency> dummy_inputs;
				dummy_inputs.push_back(*parent);
				dummy_inputs.push_back(field_name);
				dummy_inputs.push_back(rhs);
				CallFunction("set", dummy_inputs);
				return;
			}
		}

		if (sd)
			sd->SetValue(index_1, index_2, field_name, rhs);
	}
	else
	{
		if (index_2 == -1)
		{
			if (index_1 >= sd->Size())
				sd->DimensionNew(index_1+1);
		}

		Currency new_parent = sd->GetValue(index_1, index_2, field_name);

		if (new_parent.IsScalar() && HasBuiltin("ishandle")) // special case for plots and such
		{
			std::vector<Currency> dummy_input;
			dummy_input.push_back(new_parent);
			Currency is_handle = CallFunction("ishandle", dummy_input);

			if (is_handle.IsScalar() && (is_handle.Scalar() == 1.0))
			{
				int field_children = field_tree->ChildCount();

				if (field_children == 2)
				{
					OMLTree* temp = field_tree->GetChild(1);
					OMLTree* temp2 = temp->GetChild(0);
					field_name = temp2->GetText();

					std::vector<Currency> dummy_inputs;
					dummy_inputs.push_back(new_parent);
					dummy_inputs.push_back(field_name);
					dummy_inputs.push_back(rhs);
					CallFunction("set", dummy_inputs);
					return;
				}
			}
		}

		if (new_parent.IsEmpty())
		{
			StructData* new_data = new StructData;
			sd->SetValue(index_1, index_2, field_name, new_data);
			new_parent = sd->GetValue(index_1, index_2, field_name);
		}
		else if (!new_parent.IsStruct())
		{
			std::vector<Currency> indices;
			int num_field_indices = field_tree->ChildCount();

			OMLTree* sub_index_tree = NULL;
			OMLTree* sub_field_tree = NULL;

			int cell_index_type = 0;

			for (int j=0; j<num_field_indices; j++)
			{
				OMLTree* temp = field_tree->GetChild(j);

				int index_type = temp->GetType();

				if ((index_type == PARAM_LIST) || (index_type == CELL_PARAM_LIST))
				{
					cell_index_type = index_type;
					sub_index_tree = temp;

					int num_target_indices = temp->ChildCount();

					end_context_varname = NULL;
					end_context_currency = &new_parent;

					for (int k=0; k<num_target_indices; k++)
					{
						end_context_index = k;
						OMLTree* idx = temp->GetChild(k);
						indices.push_back(RUN(idx));
					}

					end_context_currency = NULL;
				}
				else if ((index_type == FIELD)  || (index_type == DYN_FIELD))
				{
					if (cell_index_type == PARAM_LIST)
						throw OML_Error("Invalid field access");

					if (new_parent.IsCellArray())
					{					
						sub_field_tree = temp;

						HML_CELLARRAY* cells  = new_parent.CellArray();
						StructData*    loc_sd = NULL;

						if (indices.size() == 2)
						{
							int cell_index1 = (int)indices[0].Scalar()-1;
							int cell_index2 = (int)indices[1].Scalar()-1;

							if ((cell_index1 < 0) || (cell_index2 < 0))
								throw OML_Error("Invalid index");

							if ((cell_index1 < cells->M()) && (cell_index2 < cells->N()))
							{
								Currency cell_val = (*cells)(cell_index1, cell_index2);

								if (cell_val.IsStruct())
									loc_sd = cell_val.Struct();
								else
									throw OML_Error("Invalid assignment");
							}
							else
							{
								// we're growing the cell so make a new one
								loc_sd = new StructData;
							}
						}
						else if (indices.size() == 1)
						{
							int cell_index1 = (int)indices[0].Scalar()-1;

							if (cell_index1 < 0)
								throw OML_Error("Invalid index");

							if (cell_index1 < cells->Size())
							{
								Currency& cell_val = (*cells)(cell_index1);

								if (cell_val.IsStruct())
								{
									loc_sd = cell_val.Struct();
								}
								else if (cell_val.IsCellArray())
								{
									HML_CELLARRAY* cells = cell_val.CellArray();

									if (cells->Size() == 0)
									{
										cell_val.MakeStruct();
										loc_sd = cell_val.Struct();
									}
									else
									{
										throw OML_Error("Invalid assignment");
									}
								}
								else
								{
									throw OML_Error("Invalid assignment");
								}
							}
							else
							{
								// we're growing the cell so make a new one
								loc_sd = new StructData;
							}
						}
						else
						{
							throw OML_Error("Invalid assignment");
						}

						loc_sd->IncrRefCount();
						Currency temp_cur(loc_sd);
						StructAssignmentHelper(&temp_cur, NULL, sub_field_tree, rhs);
						AssignHelper(new_parent, indices, temp_cur);
						sd->SetValue(index_1, index_2, field_name, new_parent);
						return;
					}
				}
			}
			AssignHelper(new_parent, indices, rhs);
			sd->SetValue(index_1, index_2, field_name, new_parent);
			return;
		}

		OMLTree* indices2 = NULL;
		OMLTree* field    = NULL;

		for (int j=1; j<num_field_children; j++)
		{
			OMLTree* temp   = field_tree->GetChild(j);

			if (temp->GetType() == PARAM_LIST)
				indices2 = temp;
			else if (temp->GetType() == FIELD)
				field = temp;
			else if (temp->GetType() == DYN_FIELD)
				field = temp;
		}

		if (indices2 && !field)
		{
			std::vector<Currency> indices_vec;
			int num_field_indices = indices2->ChildCount();

			for (int j=0; j<num_field_indices; j++)
			{
				OMLTree* temp = indices2->GetChild(j);
				indices_vec.push_back(RUN(temp));
			}

			if (rhs.IsStruct())
			{
				if (!indices_vec[0].IsScalar())
					throw OML_Error(HW_ERROR_INV1STIND);

				int loc_index1 = (int)indices_vec[0].Scalar();
				int loc_index2 = -1;

				if (indices_vec.size() == 2)
				{
					if (!indices_vec[1].IsScalar())
						throw OML_Error(HW_ERROR_INV2NDIND);

					loc_index2 = (int)indices_vec[1].Scalar();
				}

				const Currency& loc_cur = sd->GetValue(index_1, index_2, field_name);
				StructData*     loc_sd  = loc_cur.Struct();

				if (loc_index2 == -1)
					loc_sd->DimensionNew(loc_index1);
				else
					loc_sd->DimensionNew(loc_index1, loc_index2);

				loc_sd->SetElement(loc_index1, loc_index2, rhs.Struct());
			}
			else
			{
				hwMatrix* loc_mtx = new hwMatrix(0, 0, hwMatrix::REAL);
				Currency  loc_temp(loc_mtx);
				AssignHelper(loc_temp, indices_vec, rhs);
				loc_temp.SetMask(rhs.GetMask());
				sd->SetValue(index_1, index_2, field_name, loc_temp);
			}
		}
		else
		{
			StructAssignmentHelper(&new_parent, indices2, field, rhs);
			parent->Struct()->SetValue(index_1, index_2, field_name, new_parent);
		}
	}
}

Currency ExprTreeEvaluator::ObjectMethodCall(Currency* parent, OMLTree* indices, OMLTree* field_tree)
{
	ClassInfo* ci = (*class_info_map)[parent->GetClassname()];

	// ignore indices for now
	std::string method_name = field_tree->GetText();

	int num_field_children = field_tree->ChildCount();

	FunctionInfo* fi = ci->GetFunctionInfo(method_name);

	std::vector<Currency> inputs;

	if (fi)
	{
		if (ci->IsSubclassOf("handle"))
			inputs.push_back(parent);
		else
			inputs.push_back(*parent);
	}

	if (indices)
	{
		int num_args = indices->ChildCount();

		for (int j=0; j<num_args; j++)
		{
			OMLTree* child_j = indices->GetChild(j);
			inputs.push_back(RUN(child_j));
		}
	}

	if (fi)
	{
		return CallInternalFunction(fi, inputs);
	}
	else
	{
        if (!ci->HasProperty(method_name))
            throw OML_Error("Unknown property");

		StructData* sd = parent->Struct();

        bool access_allowed = !ci->IsPropertyPrivate(method_name);

		if (!access_allowed)
			access_allowed = ci->IsClassMethod(msm->GetCurrentScope()->GetFunctionInfo());

		if (access_allowed)
		{
			Currency temp = sd->GetValue(0, 0, method_name);

			if (indices)
				return VariableIndex(temp, inputs);
			else
				return temp;
		}

        throw OML_Error("Unable to access private property");
	}

}

#if 0
Currency ExprTreeEvaluator::ObjectMethodCall(Currency* parent, OMLTree* indices, OMLTree* field_tree)
{
	ClassInfo* ci = (*class_info_map)[parent->GetClassname()];

	// ignore indices for now
	std::string method_name = field_tree->GetChild(0)->GetText();

	int num_field_children = field_tree->ChildCount();

	FunctionInfo* fi = ci->GetFunctionInfo(method_name);

	if (fi)
	{
		std::vector<Currency> inputs;

		if (ci->IsSubclassOf("handle"))
			inputs.push_back(parent);
		else
			inputs.push_back(*parent);

		if (field_tree->ChildCount() > 1)
		{
			OMLTree* args_tree = field_tree->GetChild(1);
			int num_args = args_tree->ChildCount();

			for (int j=0; j<num_args; j++)
			{
				OMLTree* child_j = args_tree->GetChild(j);
				inputs.push_back(RUN(child_j));
			}
		}

		return CallInternalFunction(fi, inputs);
	}
	else
	{
        if (!ci->HasProperty(method_name))
            throw OML_Error("Unknown property");

		StructData* sd = parent->Struct();

        bool isPrivate = ci->IsPropertyPrivate(method_name);
        if (!isPrivate)
		{
			Currency temp(new StructData(*sd));
			return StructValueHelper(&temp, indices, field_tree);
		}

		bool in_class_method = ci->IsClassMethod(
                               msm->GetCurrentScope()->GetFunctionInfo());

        if (in_class_method)
	        return sd->GetValue(0, 0, method_name);

        throw OML_Error("Unable to access private property");
	}

}
#endif

const Currency& ExprTreeEvaluator::GetValue(std::string varname, int offset) const 
{ 
	return msm->GetValue(varname, offset);
}

const Currency& ExprTreeEvaluator::GetGlobalValue(std::string varname) const 
{ 
	msm->AddGlobalReference(varname);
	return msm->GetValue(varname);
}

bool ExprTreeEvaluator::Contains(const std::string& varname) const
{ 
	const Currency& temp_cur = msm->GetValue(varname);
	return !temp_cur.IsNothing();
}

bool ExprTreeEvaluator::SetValue(std::string varname, const Currency& value)
{
	msm->SetValue(varname, value);
	return true;
}

bool ExprTreeEvaluator::SetGlobalValue(std::string varname, const Currency& value)
{
	msm->AddGlobalReference(varname);
	return SetValue(varname, value);
}

void ExprTreeEvaluator::Clear(const std::string& name)
{
	if (!(ClearFromVariables(name) || ClearFromGlobals(name)))
		ClearFromFunctions(name);
}

void ExprTreeEvaluator::Clear(const std::regex& name)
{
	if (!(ClearFromVariables(name) || ClearFromGlobals(name)))
		ClearFromFunctions(name);
}

bool ExprTreeEvaluator::ClearFromFunctions(const std::string& name)
{
    std::map<std::string, UserFunc*>::iterator func_iter = functions->find(name);

	if (func_iter != functions->end() && !func_iter->second->locked)
	{
		functions->erase(func_iter);

		OnUpdateFuncList();
		return true;
	}
	return false;
}

static bool matches(const std::regex& exp, const std::string& str)
{
    std::smatch sm;
    return std::regex_match(str, sm, exp);
}
//------------------------------------------------------------------------------
//! Returns true if successful in removing functions. If this function returns
//! true, caller must call OnUpdateFuncList()
//! \param[in] funcs Functions list which is updated
//! \param[in] name  Name of the function to remove
//------------------------------------------------------------------------------
bool ExprTreeEvaluator::RemoveFuncs(std::map<std::string, UserFunc*>& funcs, 
                                    const std::regex&                 name)
{
    bool rv = false;
    std::map<std::string, UserFunc*>::iterator iter;
    for (iter = funcs.begin(); iter != funcs.end();)
    {
        if (!iter->second->locked && matches(name, iter->first))
        {
            funcs.erase(iter++);
            rv = true;
        }
        else
        {
            ++iter;
        }
    }
    return rv;
}

void ExprTreeEvaluator::RemoveFunctionsInPath(const std::string& path)
{
	std::vector<std::string>                   to_clear;
    std::map<std::string, UserFunc*>::iterator iter;

    for (iter = functions->begin(); iter != functions->end(); iter++)
    {
		FunctionInfo* fi = iter->second->fi;

		std::string filename = fi->FileName();
		std::string pathname;

#ifdef OS_WIN
		size_t slash_index = filename.find_last_of("/\\");
#else
		size_t slash_index = filename.rfind('/');
#endif

		if (slash_index != std::string::npos)
			pathname = filename.substr(0, slash_index);

		if (pathname == path)
			to_clear.push_back(fi->FunctionName());
	}

	for (size_t j=0; j<to_clear.size(); j++)
		ClearFromFunctions(to_clear[j]);
}

bool ExprTreeEvaluator::ClearFromFunctions(const std::regex& name)
{
	std::string current_filename = msm->GetCurrentScope()->GetFilename();

    bool rv = RemoveFuncs(*functions, name);

    if (rv)
		OnUpdateFuncList();
	return rv;
}

bool ExprTreeEvaluator::ClearFromGlobals(const std::string& name)
{
	if (msm->IsGlobal(name))
	{
        msm->ClearFromGlobals(name);
		return true;
	}
	return false;
}

bool ExprTreeEvaluator::ClearFromGlobals(const std::regex& name)
{
    return msm->ClearFromGlobals(name);
}

bool ExprTreeEvaluator::ClearFromVariables(const std::string& name)
{
	const Currency& temp_cur = msm->GetValue(name);

	if (!temp_cur.IsNothing()) 
	{
		msm->Remove(name);
		msm->HideGlobal(name);
		return true;
	}
	return false;
}

bool ExprTreeEvaluator::ClearFromVariables(const std::regex& name)
{
    return msm->Remove(name);
}

void ExprTreeEvaluator::ClearFunctions()
{
    std::map<std::string, UserFunc*>::iterator iter;
    for (iter = functions->begin(); iter != functions->end();)
    {
		bool locked = false;

		if (iter->second && iter->second->locked)
			locked = true;

        if (!locked)
            functions->erase(iter++);
        else
            ++iter;
    }

    OnUpdateFuncList();
}

void ExprTreeEvaluator::ClearGlobals()
{
	msm->ClearGlobals();
}

void ExprTreeEvaluator::ClearVariables()
{
	msm->ClearLocals();
}

int ExprTreeEvaluator::RenameVariable(std::string old_name, std::string new_name)
{
	const Currency& temp_cur = msm->GetValue(new_name);

	if (!temp_cur.IsNothing())
		return -1;

	Currency cur = msm->GetValue(old_name);

	msm->Remove(old_name);
	msm->SetValue(new_name, cur);

	return 0;
}

void ExprTreeEvaluator::SetLock(const std::string& fname, bool state)
{
    UserFunc* ufunc = GetUserFunc(fname);
    if (ufunc)
        ufunc->locked = state;
}

void ExprTreeEvaluator::Lock(const std::string& fname)
{
    SetLock(fname, true);
}

void ExprTreeEvaluator::Unlock(const std::string& fname)
{
    SetLock(fname, false);
}

void ExprTreeEvaluator::LockCurrent()
{
    if (user_func_calls.empty())
        throw OML_Error(HW_ERROR_LOCKFAILNOTFUNC);

    Lock(user_func_calls.back());
}

void ExprTreeEvaluator::UnlockCurrent()
{
    if (user_func_calls.empty())
        throw OML_Error(HW_ERROR_UNLOCKFAILNOTFUNC);

    Unlock(user_func_calls.back());
}

bool ExprTreeEvaluator::IsCurrentLocked() const
{
    if (user_func_calls.empty())
        throw OML_Error(HW_ERROR_CANNOTCHECKLOCKNOTFUNC);

    return IsLocked(user_func_calls.back());
}

bool ExprTreeEvaluator::IsLocked(const std::string& fname) const
{
    UserFunc* ufunc = GetUserFunc(fname);
    return ufunc ? ufunc->locked : false;
}

std::string ExprTreeEvaluator::GetCurrentFilename() const
{
	std::string ret_val;

	if (msm->GetCurrentScope())
		ret_val = msm->GetCurrentScope()->GetFilename();

	return ret_val;
}

UserFunc* ExprTreeEvaluator::GetUserFunc(const std::string& fname) const
{
    std::map<std::string, UserFunc*>::const_iterator iter = functions->find(fname);
    if (iter != functions->end())
        return iter->second;
    return nullptr;
}

std::vector<std::string> ExprTreeEvaluator::GetVariableNames(int offset) const
{
	std::vector<std::string> ret_val;
	MemoryScope* memory = msm->GetScope(offset);

	if (memory)
		ret_val = memory->GetVariableNames();

	return ret_val;
}

std::vector<std::string> ExprTreeEvaluator::GetBuiltinFunctionNames() const
{
	std::vector<string> function_names;
	std::map<std::string, BuiltinFunc>::const_iterator iter;
	for (iter = std_functions->begin(); iter != std_functions->end(); ++iter)
		function_names.push_back(iter->first);
	return function_names;
}

std::vector<std::string> ExprTreeEvaluator::GetFunctionNames() const
{
	std::set<std::string> function_names;

	std::map<std::string, UserFunc*>::const_iterator iter;

	for (iter = functions->begin(); iter != functions->end(); ++iter)
		function_names.insert(iter->first);

	std::map<std::string, BuiltinFunc>::const_iterator iter2;

	for (iter2 = std_functions->begin(); iter2 != std_functions->end(); ++iter2)
		function_names.insert(iter2->first);

	std::vector<std::string>::const_iterator iter3;

	for (iter3 = preregistered_functions->begin(); iter3 != preregistered_functions->end(); ++iter3)
		function_names.insert(*iter3);

	std::vector<std::string> keywords = GetKeywords();

	std::vector<std::string>::iterator iter4;

	for (iter4 = keywords.begin(); iter4 != keywords.end(); ++iter4)
		function_names.erase(*iter4);

	std::set<std::string>::const_iterator iter5;

	std::vector<std::string> function_names_vector;

	for (iter5 = function_names.begin(); iter5 != function_names.end(); ++iter5)
		function_names_vector.push_back(*iter5);

	return function_names_vector;
}

std::vector<std::string> ExprTreeEvaluator::GetKeywords() const
{
	std::vector<string> keywords;

	keywords.push_back("if");
	keywords.push_back("while");
	keywords.push_back("for");
	keywords.push_back("switch");
	keywords.push_back("case");
	keywords.push_back("otherwise");
	keywords.push_back("global");
	keywords.push_back("persistent");
	keywords.push_back("return");
	keywords.push_back("try");
	keywords.push_back("catch");
	keywords.push_back("end");
	keywords.push_back("function");
	keywords.push_back("else");
	keywords.push_back("elseif");
	keywords.push_back("break");
	keywords.push_back("continue");
	keywords.push_back("varargin");
	keywords.push_back("varargout");
	keywords.push_back("classdef");
	return keywords;
}

std::vector<std::string> ExprTreeEvaluator::GetOperators() const
{
	std::vector<string> operators;
	operators.push_back("+");
	operators.push_back("-");
	operators.push_back("*");
	operators.push_back("/");
	operators.push_back("\\");
	operators.push_back(".*");
	operators.push_back("./");
	operators.push_back(".\\");
	operators.push_back("~");
	operators.push_back("^");
	operators.push_back(".^");
	operators.push_back("=");
	operators.push_back("==");
	operators.push_back("~=");
	operators.push_back("<");
	operators.push_back(">");
	operators.push_back("<=");
	operators.push_back(">=");
	operators.push_back("&");
	operators.push_back("&&");
	operators.push_back("|");
	operators.push_back("||");
	operators.push_back(":");
	return operators;
}

int ExprTreeEvaluator::GetNarginValue() const
{
	if (!nargin_values.size())
		throw OML_Error(HW_ERROR_INVCALLNARGIN);

	return nargin_values.back();
}

int ExprTreeEvaluator::GetNargoutValue() const
{
	if (!nargout_values.size())
		throw OML_Error(HW_ERROR_INVCALLNARGOUT);

	return nargout_values.back();
}

void ExprTreeEvaluator::PushNargValues(int args_in, int args_out)
{
	nargout_values.push_back(args_out);
	nargin_values.push_back(args_in);
}

void ExprTreeEvaluator::PopNargValues()
{
	assert(nargout_values.size());
	assert(nargin_values.size());
	
	if (nargout_values.size())
		nargout_values.pop_back();

	if (nargin_values.size())
		nargin_values.pop_back();
}

int ExprTreeEvaluator::GetContextEndValue() const
{
	Currency context;

	if (end_context_varname)
		context = msm->GetValue(end_context_varname);
	else if (end_context_currency)
		context = *end_context_currency;
	else
		throw OML_Error(HW_ERROR_NOCONTEXTENDFUNC);

	if (context.IsMatrixOrString())
	{
		const hwMatrix* mtx = context.Matrix();

		if (end_context_index == -1)
			return mtx->Size();
		else if (end_context_index == 0)
			return mtx->M();
		else if (end_context_index == 1)
			return mtx->N();
	}
	if (context.IsNDMatrix())
	{
		const hwMatrixN* mtx = context.MatrixN();

		if (end_context_index == -1)
		{
			return mtx->Size();
		}
		else
		{
			std::vector<int> dims = mtx->Dimensions();
			return dims[end_context_index];
		}
	}
	else if (context.IsScalar() || context.IsComplex())
	{
		return 1;
	}
	else if (context.IsStruct())
	{
		StructData* sd = context.Struct();

		if (end_context_index == -1)
			return (sd->M()*sd->N());
		else if (end_context_index == 0)
			return sd->M();
		else if (end_context_index == 1)
			return sd->N();
	}
	else if (context.IsCellArray())
	{
		HML_CELLARRAY* mtx = context.CellArray();

		if (end_context_index == -1)
			return mtx->Size();
		else if (end_context_index == 0)
			return mtx->M();
		else if (end_context_index == 1)
			return mtx->N();
	}

	return 0;
}

Currency ExprTreeEvaluator::TryCatch(OMLTree* tree)
{
	Mark();

	try
	{
		OMLTree* child_0 = tree->GetChild(0);
		Currency ret = RUN(child_0);
		Unmark();
		return ret;
	}
	catch (OML_Error& error)
	{
		lasterrormsg = error.GetErrorMessage();
		Restore(); // this also calls Unmark

		int num_children = tree->ChildCount();

		if (num_children > 1)
		{
			for (int j=num_children-1; j>0; j--)
			{
				OMLTree*    child = tree->GetChild(j);

				if (child->GetType() == IDENT)
				{
					std::string err_var = child->GetText();

					StructData* sd = new StructData();
					sd->DimensionNew(1, 1);
					sd->SetValue(0, 0, "message", lasterrormsg);
					msm->SetValue(err_var, sd);
				}
				else
				{
					return RUN(child);
				}
			}
		}
	}

	return Currency(-1, Currency::TYPE_NOTHING);
}

Currency ExprTreeEvaluator::CellExtraction(OMLTree* tree)
{
	bool suppress_output = suppress_multi_ret_output; 

	int  num_tree_children = tree->ChildCount();

	if (num_tree_children != 2)
		throw OML_Error(HW_ERROR_INVAST);

	OMLTree* lhs  = tree->GetChild(0);
	OMLTree* cell = tree->GetChild(1);

	int num_lhs = lhs->ChildCount();

	Currency rhs = RUN(cell);

	if (!rhs.IsCellList())
		throw OML_Error(HW_ERROR_INVCELLEXT);

	HML_CELLARRAY* cells = rhs.CellArray();

	if (num_lhs > cells->Size())
		throw OML_Error(HW_ERROR_MISSCELLVAL);

	for (int j=0; j<num_lhs; j++)
	{
		std::string varname = lhs->GetChild(j)->GetText();
		msm->SetValue(varname, (*cells)(j));
		PushResult(msm->GetValue(varname), !suppress_output);
	}

	return Currency(-1, Currency::TYPE_NOTHING);
}

Currency ExprTreeEvaluator::InPlaceExpansion(OMLTree* tree)
{
	OMLTree* var_tree = tree->GetChild(0);

	if (!var_tree->u)
		var_tree->u = (void*)Currency::vm.GetStringPointer(var_tree->GetText());

	std::string* pString = (std::string*)var_tree->u;

	if (!msm->IsGlobal(*pString) && !msm->Contains(pString))
	{
		char buffer[2048];
		sprintf(buffer, "Unknown function: %s", pString->c_str());
		throw OML_Error(buffer);
	}

	Currency& original = msm->GetMutableValue(pString);
	OMLTree* child_1 = tree->GetChild(1);
	Currency addition = RUN(child_1);

	if (original.IsScalar() || original.IsComplex())
	{
		original.ConvertToMatrix();
		hwMatrix* mat = original.GetWritableMatrix();

		const hwMatrix* mat2 = addition.ConvertToMatrix();
		if (mat->M() == mat2->M())
		{
			int start_col = mat->N();
			mat->WriteSubmatrix(0, start_col, *mat2);
			return original;
		}
	}
	else if (original.IsMatrix() || original.IsString())
	{
		hwMatrix* mat = original.GetWritableMatrix();

		if (mat->GetRefCount() != 1)
		{
			mat = allocateMatrix(mat);
			original.ReplaceMatrix(mat);
		}

		if (mat->IsEmpty())
		{
			if (addition.IsScalar() || addition.IsComplex())
				addition.ConvertToMatrix();

			if (!addition.IsMatrix() && !addition.IsString())
				throw OML_Error(HW_ERROR_ILLASSIGN);

			hwMatrix* new_mat = addition.GetWritableMatrix();
			new_mat->IncrRefCount();
			original.ReplaceMatrix(new_mat);

			if (addition.IsString())
				original.SetMask(Currency::MASK_STRING);

			return original;
		}
		else
		{
			if (addition.IsScalar())
			{
				int start_col = mat->N();
				mat->Resize(mat->M(), start_col+1, false);

				if (mat->IsReal())
					(*mat)(mat->M()-1, start_col) = addition.Scalar();
				else
					mat->z(mat->M()-1, start_col) = addition.Scalar();

				return original;
			}
			else if (addition.IsComplex())
			{
				int start_col = mat->N();
				mat->Resize(mat->M(), start_col+1, false);
				mat->MakeComplex();
				mat->z(mat->M()-1, start_col) = addition.Complex();
				return original;
			}
			else
			{
				const hwMatrix* mat2 = addition.ConvertToMatrix();
				if (mat->M() == mat2->M())
				{
					int start_col = mat->N();
					mat->WriteSubmatrix(0, start_col, *mat2);
					return original;
				}
			}
		}
	}
	else if (original.IsCellArray())
	{
		HML_CELLARRAY* cells = original.CellArray();

		if (cells->GetRefCount() != 1)
		{
			cells = allocateCellArray(cells);
			original.ReplaceCellArray(cells);
		}

		if (cells->IsEmpty())
		{
			if (!addition.IsCellArray())
				throw OML_Error(HW_ERROR_ILLASSIGN);

			HML_CELLARRAY* new_cells = addition.CellArray();
			new_cells->IncrRefCount();
			original.ReplaceCellArray(new_cells);

			return original;
		}
		else
		{
			if (addition.IsCellArray())
			{
				HML_CELLARRAY* added_cells = addition.CellArray();

				if (added_cells->M() == cells->M())
				{
					int start_col = cells->N();
					cells->Resize(cells->M(), start_col+added_cells->N(), false);

					for (int k=0; k<cells->M(); k++)
					{
						for (int j=0; j<added_cells->N(); j++)
							(*cells)(k, start_col+j) = (*added_cells)(j);
					}

					return original;
				}
			}
		}
	}

	throw OML_Error(HW_ERROR_ILLASSIGN);
}

Currency ExprTreeEvaluator::ClassDefinition(OMLTree* tree)
{
	OMLTree* var_tree = tree->GetChild(0);

	if (!var_tree->u)
		var_tree->u = (void*)Currency::vm.GetStringPointer(var_tree->GetText());

	std::string* class_name = (std::string*)var_tree->u;
	ClassInfo*   ci         = new ClassInfo(*class_name);
	
	int num_sections = tree->ChildCount();

	for (int j=1; j<num_sections; j++)
	{
		OMLTree*    my_tree = tree->GetChild(j);

		if (my_tree->GetType() == PROPERTIES)
		{
			int num_prop_children = my_tree->ChildCount();

            bool isPrivate = false;
			if (num_prop_children == 3)
			{
				OMLTree* lhs = my_tree->GetChild(1);
				OMLTree* rhs = my_tree->GetChild(2);

				std::string lhs_text = lhs->GetText();
				std::string rhs_text = rhs->GetText();

                isPrivate = ((lhs_text == "Access") && (rhs_text == "private"));
			}

			OMLTree* properties_tree = my_tree->GetChild(0);
	
			int num_properties = properties_tree->ChildCount();

			for (int j=0; j<num_properties; j++)
				ci->AddProperty(properties_tree->GetChild(j)->GetText(), isPrivate);
		}
		else if (my_tree->GetType() == FUNC_LIST)
		{
			int num_methods = my_tree->ChildCount();

			for (int j=0; j<num_methods; j++)
			{
				OMLTree* method = my_tree->GetChild(j);

				FunctionInfo* fi = FunctionInfoFromTree(method);

				if (fi->FunctionName() == *class_name) // it's the constructor
				{
					std::map<std::string, UserFunc*>::iterator iter = functions->find(fi->FunctionName());
					if (iter != functions->end())
						delete iter->second;

					fi->SetAsConstructor();
		
					(*functions)[fi->FunctionName()] = new UserFunc(fi);
				}
				else
				{
                    ci->AddClassMethod(fi->FunctionName(), fi);
				}
			}
		}
	}

	if (num_sections > 3)
	{
		OMLTree* subclass_tree = tree->GetChild(num_sections-1);

		for (int j=0; j<subclass_tree->ChildCount(); j++)
			ci->AddBaseClass(subclass_tree->GetChild(j)->GetText());
	}

	(*class_info_map)[*class_name] = ci;

	return Currency();
}
//------------------------------------------------------------------------------
// Returns true along with file name if given function name can be located
//------------------------------------------------------------------------------
bool ExprTreeEvaluator::FindFunction(const std::string& func_name,
                                     std::string&       file_name)
{
	// to do
	// locate a file based on func_name (either .m or .oml extension)
	// create a lexer/parser to get a new tree from the file
	// run the tree in this interpreter (so the functions get loaded)

	// avoid having to search the file system unnecessarily
	if (std::binary_search(not_found_functions->begin(), not_found_functions->end(), func_name))
		return false;

	std::string oml_file_name = func_name + ".oml";
	std::string m_file_name   = func_name + ".m";

	if (!(FindFileInPath(oml_file_name, file_name) || FindFileInPath(m_file_name, file_name)))
	{
		// Let FindEncryptedFunction modify the not_found_functions list
		return false;
	}

	return true;
}

//------------------------------------------------------------------------------
// Returns true along with file name if given function name can be located
//------------------------------------------------------------------------------
bool ExprTreeEvaluator::FindPrecompiledFunction(const std::string& func_name,
                                     std::string&       file_name)
{
	// avoid having to search the file system unnecessarily
	if (std::binary_search(not_found_functions->begin(), not_found_functions->end(), func_name))
		return false;

	std::string oml_file_name = func_name + ".omlp";

	if (!FindFileInPath(oml_file_name, file_name))
	{
		// Let FindEncryptedFunction modify the not_found_functions list
		return false;
	}

	return true;
}

bool ExprTreeEvaluator::FindEncryptedFunction(const std::string& func_name, std::string& file_name, std::string& extension)
{
	// avoid having to search the file system unnecessarily
	if (std::binary_search(not_found_functions->begin(), not_found_functions->end(), func_name))
		return false;

	std::map<std::string, ENCRPTR>::iterator iter;

	for (iter = _decryptors.begin(); iter != _decryptors.end(); iter++)
	{
		std::string file_to_search = func_name + "." + iter->first;

		if (FindFileInPath(file_to_search, file_name))
		{
			extension = iter->first;
			return true;
		}
	
	}

	not_found_functions->push_back(func_name);
	std::sort(not_found_functions->begin(), not_found_functions->end());
	return false;
}

bool ExprTreeEvaluator::ParseAndRunFile(const std::string& file_name, bool allow_script)
{
	pANTLR3_INPUT_STREAM input = ANTLRData::InputFromFilename(file_name);

	if (!input)
	{
		char buffer[1024];
		sprintf(buffer, "File %s is not accessible", file_name.c_str());
		throw OML_Error(buffer);
	}

	ANTLRData ad(input, true);

	pANTLR3_COMMON_TOKEN_STREAM tokens = ad.GetTokens();
	ANTLRData::PreprocessTokenStream(tokens);

	ad.CreateParser(tokens);
	pExprCppTreeParser parser = ad.GetParser();

	ExprCppTreeParser_prog_return r = parser->prog(parser);
		
	int temp = nested_function_marker;

	if (parser->pParser->rec->getNumberOfSyntaxErrors(parser->pParser->rec) == 0)
	{
		pANTLR3_BASE_TREE tree = r.tree;
		PreprocessAST(tree, tokens);
		OMLTree* oml_tree = OMLTree::ConvertTree(tree);

		nested_function_marker = 0;

		if (!ValidateFunction(oml_tree))
		{
			if (!allow_script)
			{
				throw OML_Error(HW_ERROR_FILENOTEXEC);
			}
			else
			{
				RunTree(oml_tree);
				delete oml_tree;
			}

			nested_function_marker = temp;
			return false;
		}

		RunTree(oml_tree);
		delete oml_tree;
	}
	else
	{
		char buffer[2048];
		sprintf(buffer, "Syntax error in included file %s at line number %d", file_name.c_str(), parser->pParser->rec->state->exception->line);
		throw OML_Error(buffer);
	}
	
	nested_function_marker = temp;	

	return true;
}

bool ExprTreeEvaluator::ParseAndRunString(const std::string& str, const std::string& use_filename)
{
	pANTLR3_INPUT_STREAM input = ANTLRData::InputFromExpression(str, use_filename);
	ANTLRData ad(input, true);

	pANTLR3_COMMON_TOKEN_STREAM tokens = ad.GetTokens();
	ANTLRData::PreprocessTokenStream(tokens);

	ad.CreateParser(tokens);
	pExprCppTreeParser parser = ad.GetParser();

	ExprCppTreeParser_prog_return r = parser->prog(parser);
		
	int temp = nested_function_marker;

	if (parser->pParser->rec->getNumberOfSyntaxErrors(parser->pParser->rec) == 0)
	{
		pANTLR3_BASE_TREE tree = r.tree;
		PreprocessAST(tree, tokens);
		nested_function_marker = 0;
		OMLTree* oml_tree = OMLTree::ConvertTree(tree);
		RunTree(oml_tree);
	}
	else
	{
		throw OML_Error("Syntax error");
	}
	
	nested_function_marker = temp;	

	return true;
}

bool ExprTreeEvaluator::ValidateFunction(OMLTree* tree)
{
	int                  num_children = tree->ChildCount();
	bool                 found_one    = false;

	for (int j=0; j<num_children; j++)
	{
		OMLTree*    child     = tree->GetChild(j);

		if (child->GetType() == FUNC_DEF)
			found_one = true;
		else if (child->GetType() == CLASSDEF)
			found_one = true;
		else if (child->GetType() != DUMMY)
			return false;
	}

	return found_one;
}

bool  ExprTreeEvaluator::RegisterBuiltInFunction(const std::string& func_name, FUNCPTR fp, int nargin, int nargout)
{
	(*std_functions)[func_name] = BuiltinFunc(fp, FunctionMetaData(nargin, nargout));
	OnUpdateFuncList();
	return true;
}

bool  ExprTreeEvaluator::RegisterBuiltInFunction(const std::string& func_name, FUNCPTR fp, FunctionMetaData fmd)
{
	(*std_functions)[func_name] = BuiltinFunc(fp, fmd);
	OnUpdateFuncList();
	return true;
}

bool  ExprTreeEvaluator::RegisterBuiltInFunction(const std::string& func_name, ALT_FUNCPTR fp)
{
	(*std_functions)[func_name] = BuiltinFunc(fp);
	OnUpdateFuncList();
	return true;
}

bool ExprTreeEvaluator::CheckForPreviousTrailingDot(pANTLR3_BASE_TREE tree, pANTLR3_COMMON_TOKEN_STREAM tokens)
{
	pANTLR3_COMMON_TOKEN tok        = tree->getToken(tree);
	ANTLR3_MARKER        index      = (int)tok->index;
	pANTLR3_VECTOR       vec        = tokens->getTokens(tokens);

	if (index <1)
		return false;

	pANTLR3_COMMON_TOKEN prev_tok   = (pANTLR3_COMMON_TOKEN)vec->get(vec, ANTLR3_UINT32(index-1));
	pANTLR3_STRING       str        = prev_tok->getText(prev_tok);
	const char*          c_str_prev = (char*)str->chars;
	size_t               len        = strlen(c_str_prev);

	if (c_str_prev[len-1] == '.')
		return true;

	return false;
}

void ExprTreeEvaluator::PreprocessAST(pANTLR3_BASE_TREE tree, pANTLR3_COMMON_TOKEN_STREAM tokens)
{
	// fix up POW, LDIV, RDIV nodes
	pANTLR3_COMMON_TOKEN tok = tree->getToken(tree);

	if(tok) 
	{
		if (tok->type == DIV)
		{
			if (CheckForPreviousTrailingDot(tree, tokens))
			{
				tok->type = EDIV;
				pANTLR3_STRING str = tok->getText(tok);
				str->set(str, "./");
				tok->setText(tok, str);
			}
		}
		else if (tok->type == LDIV)
		{
			if (CheckForPreviousTrailingDot(tree, tokens))
			{
				tok->type = ELDIV;
				pANTLR3_STRING str = tok->getText(tok);
				str->set(str, ".\\");
				tok->setText(tok, str);
			}
		}
		else if (tok->type == POW)
		{
			if (CheckForPreviousTrailingDot(tree, tokens))
			{
				tok->type = DOTPOW;
				pANTLR3_STRING str = tok->getText(tok);
				str->set(str, ".^");
				tok->setText(tok, str);
			}
		}
		else if (tok->type == TIMES)
		{
			if (CheckForPreviousTrailingDot(tree, tokens))
			{
				tok->type = ETIMES;
				pANTLR3_STRING str = tok->getText(tok);
				str->set(str, ".*");
				tok->setText(tok, str);
			}
		}
		else if (tok->type == UGLY1)
		{
			tok->type = POW;

			pANTLR3_BASE_TREE    first_child = (pANTLR3_BASE_TREE)tree->getChild(tree, 0);
			pANTLR3_COMMON_TOKEN first_token = first_child->getToken(first_child);
			
			if (first_token->type == TRANSP)
				first_token->type = CTRANSP;
		}
		else if (tok->type == UGLY2)
		{
			tok->type = DOTPOW;

			pANTLR3_BASE_TREE    first_child = (pANTLR3_BASE_TREE)tree->getChild(tree, 0);
			pANTLR3_COMMON_TOKEN first_token = first_child->getToken(first_child);
			
			if (first_token->type == TRANSP)
				first_token->type = CTRANSP;
		}
		else if (tok->type == UGLY3)
		{
			tok->type = POW;
		}
		else if (tok->type == UGLY4)
		{
			tok->type = DOTPOW;
		}
		else if ((tok->type == HML_STRING) || (tok->type == DUMMY))
		{
			std::string my_str;
			int num_words = tree->getChildCount(tree);

			if (tokens && num_words)
			{
				pANTLR3_VECTOR vec = tokens->getTokens(tokens);

				pANTLR3_BASE_TREE    first_child = (pANTLR3_BASE_TREE)tree->getChild(tree, 0);
				pANTLR3_COMMON_TOKEN first_token = first_child->getToken(first_child);
				ANTLR3_MARKER        first_index = first_token->getTokenIndex(first_token);

				pANTLR3_BASE_TREE    last_child  = (pANTLR3_BASE_TREE)tree->getChild(tree, num_words-1);
				pANTLR3_COMMON_TOKEN last_token  = last_child->getToken(last_child);
				ANTLR3_MARKER        last_index  = last_token->getTokenIndex(last_token);

				if (tok->type == DUMMY)
				{
					first_child = (pANTLR3_BASE_TREE)tree->getChild(tree, 0);
					first_token = first_child->getToken(first_child);
					first_index = first_token->getTokenIndex(first_token) - 1;
					last_child  = (pANTLR3_BASE_TREE)tree->getChild(tree, num_words-1);
					last_token  = last_child->getToken(last_child);
					last_index  = last_token->getTokenIndex(last_token) + 1;
				}
				else if ((first_index == 0) && (last_index == 0))
				{
					// This is to account for the alternate function syntax
					// I'm artificially inserting the 's into the tree, so there is no token 
					// for them :-|
					first_child = (pANTLR3_BASE_TREE)tree->getChild(tree, 1);
					first_token = first_child->getToken(first_child);
					first_index = first_token->getTokenIndex(first_token) - 1;
					last_child  = (pANTLR3_BASE_TREE)tree->getChild(tree, num_words-2);
					last_token  = last_child->getToken(last_child);
					last_index  = last_token->getTokenIndex(last_token) + 1;
				}

				for (int j=(int)(first_index+1); j<(int)last_index; j++)
				{
					pANTLR3_COMMON_TOKEN tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, j);

					if (tok->getType(tok) == EQUOTE)
						my_str.append("'");
					else
						my_str.append((char*)tok->getText(tok)->chars);
				}

				for (int j=num_words-1; j>0; j--)
					tree->deleteChild(tree, j);

				tree = (pANTLR3_BASE_TREE)tree->getChild(tree, 0);
				tok = tree->getToken(tree);
				pANTLR3_STRING str = tok->getText(tok);
				str->set(str, my_str.c_str());
				tok->setText(tok, str);
			}
		}
		else if (tok->type == ASSIGN) // in-place expansion
		{
			pANTLR3_BASE_TREE    first_child = (pANTLR3_BASE_TREE)tree->getChild(tree, 0);
			pANTLR3_COMMON_TOKEN first_tok   = first_child->getToken(first_child);

			if (first_tok->getType(first_tok) == IDENT)
			{
				std::string target_name = getText(first_child);

				pANTLR3_BASE_TREE    second_child = (pANTLR3_BASE_TREE)tree->getChild(tree, 1);
				pANTLR3_COMMON_TOKEN second_tok   = second_child->getToken(second_child);

				if (second_tok->getType(second_tok) == MATRIX)
				{
					if (second_child->getChildCount(second_child) == 1)
					{
						pANTLR3_BASE_TREE    vector        = (pANTLR3_BASE_TREE)second_child->getChild(second_child, 0);
						pANTLR3_BASE_TREE    first_element = (pANTLR3_BASE_TREE)vector->getChild(vector, 0);
						pANTLR3_COMMON_TOKEN first_el_tok  = first_element->getToken(first_element);

						if (vector->getChildCount(vector) > 1) // account for a = [a]
						{
							if (first_el_tok->getType(first_el_tok) == IDENT)
							{
								std::string first_element_name = getText(first_element);

								if (first_element_name == target_name)
								{
									tok->type = INPLACE;
									vector->deleteChild(vector, 0);
									pANTLR3_STRING str = tok->getText(tok);
									str->set(str, "INPLACE_EXPANSION");
									tok->setText(tok, str);
								}
							}
						}
					}
				}
			}
		}

		int num_children = tree->getChildCount(tree);

		for (int j=0; j<num_children; j++)
			PreprocessAST((pANTLR3_BASE_TREE)tree->getChild(tree, j), tokens);
	}
}

void ExprTreeEvaluator::Mark()
{
    MemoryScope* curmem = msm->GetCurrentScope();
	marks.push_back(curmem);
	mark_narg_size = (int)nargin_values.size();
}

void ExprTreeEvaluator::Unmark()
{
	marks.pop_back();
}

void ExprTreeEvaluator::Restore()
{
    if (!marks.empty())
    {
	    msm->CloseTo(marks.back());
	    Unmark();
    }

	while (nargin_values.size() > mark_narg_size)
	{
		nargin_values.pop_back();
		nargout_values.pop_back();
	}

	assignment_nargout = 0;
	end_context_varname = NULL;
	builtin_error_scope = "";
	_lhs_eval = false;
}

hwMatrix* ExprTreeEvaluator::SubmatrixSingleIndexHelper(Currency& target, const Currency& indices, const Currency& value)
{
	hwMatrix* data = target.GetWritableMatrix();
	const hwMatrix* index = NULL;
	
	if (indices.IsMatrix())
		index = indices.Matrix();

	hwMatrix* new_matrix = data;

	if (new_matrix->GetRefCount() != 1)
		new_matrix = allocateMatrix(data);

	if (index && !index->IsRealData())
		throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_INDEX);

	if (indices.IsColon()) 
	{
        // If this this a single colon as in the case below, set elements of
        // new matrix and not data
        // Matrix=ones(6,6)
        // Matrix(:)=1:36
        // The test above would have failed if data is used instead of new_matrix
		if (value.IsScalar())
		{
            assert(new_matrix);
            new_matrix->SetElements(value.Scalar()); 
		}
		else if (value.IsComplex())
		{
			new_matrix->MakeComplex();
            new_matrix->SetElements(value.Complex());
		}
		else if (value.IsEmpty())
		{
            new_matrix->Resize(0, 0, false);
		}
		else if (value.IsMatrix())
		{
			const hwMatrix* val_mtx = value.Matrix();

			if (data->Size() == val_mtx->Size())
			{
				for (int j=0; j< data->Size(); j++)
				{
					if (val_mtx->IsReal())
					{
                        (*new_matrix)(j) = (*val_mtx)(j);
					}
					else
					{
						new_matrix->MakeComplex();
						new_matrix->z(j) = val_mtx->z(j);
					}
				}
			}
			else
			{
				throw OML_Error(HW_ERROR_ILLASSIGN);
			}
		}
	}
	else if (indices.GetMask() == Currency::MASK_DOUBLE)
	{
		if (value.IsEmpty())
		{
			std::vector<int> to_delete;

			for (int j=0; j<index->Size(); j++)
				to_delete.push_back((int)((*index)(j)));

			std::sort(to_delete.begin(), to_delete.end());

			if (new_matrix->M() == 1)
				new_matrix->Reshape(1, new_matrix->Size());		
			else
				new_matrix->Reshape(new_matrix->Size(), 1);

			if (new_matrix->M() ==1 )
			{
				int range_start = -1;
				int num_cols    = 1;

				for (int j=(int)to_delete.size()-1; j>=0; j--)
				{
					if (range_start == -1)
					{
						range_start = to_delete[j]-1;
						num_cols    = 1;
					}
					else
					{
						int cur_index = to_delete[j]-1;

						if (cur_index == (range_start-1))
						{
							num_cols++;
							range_start = cur_index;
						}
						else
						{
							new_matrix->DeleteColumns(range_start, num_cols);
							range_start = cur_index;
							num_cols = 1;
						}
					}
				}

				if (range_start != -1)
					new_matrix->DeleteColumns(range_start, num_cols);
			}
			else
			{
				int range_start = -1;
				int num_rows    = 1;

				for (int j=(int)to_delete.size()-1; j>=0; j--)
				{
					if (range_start == -1)
					{
						range_start = to_delete[j]-1;
						num_rows    = 1;
					}
					else
					{
						int cur_index = to_delete[j]-1;

						if (cur_index == (range_start-1))
						{
							num_rows++;
							range_start = cur_index;
						}
						else
						{
							new_matrix->DeleteRows(range_start, num_rows);
							range_start = -1;
						}
					}
				}

				if (range_start != -1)
					new_matrix->DeleteRows(range_start, num_rows);
			}
		}
		else
		{
			bool is_scalar = value.IsScalar();
			double scalar_val;
			
			if (is_scalar)
				scalar_val = value.Scalar();

			for (int j=0; j<index->Size(); j++)
			{
				int target_index = (int)((*index)(j)-1);

				if (target_index < 0)
					throw OML_Error(HW_ERROR_INVIND);

				if (target_index >= new_matrix->Size())
				{
					if (new_matrix->IsEmpty())
						new_matrix->Dimension(1, target_index+1, hwMatrix::REAL);

					if (new_matrix->M() == 1)
						new_matrix->Resize(1, target_index+1, true);
					else if (new_matrix->N() == 1)
						new_matrix->Resize(target_index+1, 1, true);
					else
						throw OML_Error(HW_ERROR_UNSUPOP);
				}

				if (is_scalar)
				{
					(*new_matrix)(target_index) = scalar_val;
				}
				else if (value.IsComplex())
				{
					if (new_matrix->IsReal())
						new_matrix->MakeComplex();

					new_matrix->z(target_index) = value.Complex();
				}
				else if (value.IsMatrix() || value.IsString())
				{
					const hwMatrix* val_mtx = value.Matrix();
					hwMatrix*       tmp_mtx = NULL;
						
					// special case to deal with two vectors that are different orientations
					// but still valid since this is a single index
					if ((val_mtx->M() == index->N()) && (val_mtx->N() == 1) && (index->M() == 1))
					{
						tmp_mtx = allocateMatrix(val_mtx);
						tmp_mtx->Transpose();
						val_mtx = tmp_mtx;
					}

					if (val_mtx->M() != index->M())
						throw OML_Error(HW_ERROR_RHSSIZEIND);

					if (val_mtx->N() != index->N())
						throw OML_Error(HW_ERROR_RHSSIZEIND);

					if (new_matrix->IsReal())
					{
						if (val_mtx->IsReal())
						{
							(*new_matrix)(target_index) = (*val_mtx)(j);
						}
						else
						{
							new_matrix->MakeComplex();
							new_matrix->z(target_index) = val_mtx->z(j);
						}
					}
					else
					{
						if (val_mtx->IsReal())
						{
							new_matrix->z(target_index) = (*val_mtx)(j);
						}
						else
						{
							new_matrix->z(target_index) = val_mtx->z(j);
						}
					}

					delete tmp_mtx;
				}
				else
				{
					throw OML_Error(HW_ERROR_INVRIGHT);
				}
			}
		}
	}
	else if (indices.IsLogical())
	{
		if (value.IsCharacter())
		{
			const hwMatrix* values = value.Matrix();

			for (int j=0; j<index->Size(); j++)
			{
				if ((*index)(j) == 1.0)
					(*new_matrix)(j) = (*values)(0);
			}
		}
		else if (value.IsMatrix() || value.IsString())
		{
			const hwMatrix* values = value.Matrix();

			if (values->IsReal())
			{
				int next_index = 0;

				for (int j=0; j<index->Size(); j++)
				{
					if ((*index)(j) == 1.0)
					{
						ReplaceMatrixElementHelper(new_matrix, j, (*values)(next_index));
						next_index++;
					}
				}

				if (next_index != values->Size())
					throw OML_Error(HW_ERROR_SIZENOMATCH);
			}
			else
			{
				int next_index = 0;

				if (new_matrix->IsReal())
					new_matrix->MakeComplex();


				for (int j=0; j<index->Size(); j++)
				{
					if ((*index)(j) == 1.0)
					{
						ReplaceMatrixElementHelper(new_matrix, j, values->z(next_index));
						next_index++;
					}
				}

				if (next_index != values->Size())
					throw OML_Error(HW_ERROR_SIZENOMATCH);
			}
		}
		else if (value.IsScalar())
		{
			for (int j=0; j<index->Size(); j++)
			{
				if ((*index)(j) == 1.0)
					ReplaceMatrixElementHelper(new_matrix, j, value.Scalar());
			}
		}
		else if (value.IsComplex())
		{
			if (new_matrix->IsReal())
				new_matrix->MakeComplex();

			for (int j=0; j<index->Size(); j++)
			{
				if ((*index)(j) == 1.0)
					ReplaceMatrixElementHelper(new_matrix, j, value.Complex());
			}
		}
	}
	
	return new_matrix;
}

hwMatrix* ExprTreeEvaluator::SubmatrixDoubleIndexHelper(Currency& target, const Currency& index1, const Currency& index2, const Currency& value)
{
	hwMatrix* data = target.GetWritableMatrix();
	hwMatrix* new_matrix = data;

	if (new_matrix->GetRefCount() != 1)
		new_matrix = allocateMatrix(data);

	if (index1.IsColon())
	{
		if (index2.IsScalar())
		{
			int index_2 = (int)index2.Scalar()-1;
			
			if ((index_2 < 0))
				throw OML_Error(HW_ERROR_INDEXRANGE);

			if (value.IsMatrix())
			{
				const hwMatrix* rhs_mtx = value.Matrix();

				if ((new_matrix->M() != 0) && ((rhs_mtx->M() != new_matrix->M()) || (rhs_mtx->N() != 1)))
				{
					if ((new_matrix->M() == rhs_mtx->N()) && (rhs_mtx->M() == 1))
						; // allow this even though it's technically wrong
					else
					    throw OML_Error(HW_ERROR_MATSIZE);
				}

				for (int j=0; j<rhs_mtx->Size(); j++)
				{
					if (rhs_mtx->IsReal())
						ReplaceMatrixElementHelper(new_matrix, j, index_2, (*rhs_mtx)(j));
					else
						ReplaceMatrixElementHelper(new_matrix, j, index_2, rhs_mtx->z(j));
				}
			}
			else
			{
				if (data->IsEmpty())
					data->Dimension(1, 0, hwMatrix::REAL);

				for (int j=0; j<data->M(); j++)
				{
					if (value.IsScalar())
						ReplaceMatrixElementHelper(new_matrix, j, index_2, value.Scalar());
					else if (value.IsComplex())
						ReplaceMatrixElementHelper(new_matrix, j, index_2, value.Complex());
				}
			}
		}
		else if (index2.IsVector())
		{
			std::vector<double> indices = index2.Vector();

			if (value.IsScalar())
			{
				for (int j=0; j<indices.size(); j++)
				{
					if (indices[j] < 1)
						throw OML_Error(HW_ERROR_INDEXRANGE);

    				for (int k=0; k<new_matrix->M(); k++)
						ReplaceMatrixElementHelper(new_matrix, k, (int)indices[j]-1, value.Scalar());
				}
			}
			else if (value.IsComplex())
			{
				new_matrix->MakeComplex();

				for (int j=0; j<indices.size(); j++)
				{
					if (indices[j] < 1)
						throw OML_Error(HW_ERROR_INDEXRANGE);

    				for (int k=0; k<new_matrix->M(); k++)
						ReplaceMatrixElementHelper(new_matrix, k, (int)indices[j]-1, value.Complex());
				}
			}
			else if (value.IsEmpty()) // really should be Is0x0
			{
				// This only works if the values are strictly increasing -- JDS
				for (int j=(int)indices.size()-1; j>=0; j--)
				{
					if ((indices[j] < 1) || (indices[j] > new_matrix->N()))
						throw OML_Error(HW_ERROR_INDEXRANGE);

					hwMathStatus stat = new_matrix->DeleteColumns((int)indices[j]-1, 1);
				}
			}
			else if (value.IsMatrix())
			{
				const hwMatrix* val_matrix = value.Matrix();

				if (val_matrix->N() != indices.size())
					throw OML_Error(HW_ERROR_INCOMPDIM); 

				if (val_matrix->M() != new_matrix->M())
					throw OML_Error(HW_ERROR_INCOMPDIM); 

   				for (int j=0; j<indices.size(); j++)
    			{
					int index2 = (int)indices[j];

    				for (int k=0; k<val_matrix->M(); k++)
	    			{
						int index1 = k+1;

						if (!val_matrix->IsReal())
							ReplaceMatrixElementHelper(new_matrix, index1-1, index2-1, val_matrix->z(k, j));
						else
							ReplaceMatrixElementHelper(new_matrix, index1-1, index2-1, (*val_matrix)(k, j));
					}
		        }
			}
		}
		else if (index2.IsColon())
		{
			if (value.IsScalar())
			{
				new_matrix->SetElements(value.Scalar());
			}
			else if (value.IsMatrix())
			{
				const hwMatrix* rhs_mtx = value.Matrix();

				if ((new_matrix->M() == rhs_mtx->M()) && (new_matrix->N() == rhs_mtx->N()))
					*new_matrix = *rhs_mtx;
				else if (new_matrix->IsEmpty())
					*new_matrix = *rhs_mtx;
				else
					throw OML_Error(HW_ERROR_MATSIZE);
			}
		}
	}
	else if (index2.IsColon())
	{
		if (index1.IsScalar())
		{
			int index_1 = (int)index1.Scalar()-1;

			if (value.IsMatrix())
			{
				const hwMatrix* rhs_mtx = value.Matrix();
				if ((new_matrix->N() != 0) && ((rhs_mtx->N() != new_matrix->N()) || (rhs_mtx->M() != 1)))
				{
					if ((new_matrix->N() == rhs_mtx->M()) && (rhs_mtx->N() == 1))
						; // allow this even though it's technically wrong
					else
					throw OML_Error(HW_ERROR_MATSIZE);
				}

				for (int j=0; j<rhs_mtx->Size(); j++)
				{
					if ((index_1 < 0))
						throw OML_Error(HW_ERROR_INDEXRANGE);

					if (rhs_mtx->IsReal())
						ReplaceMatrixElementHelper(new_matrix, index_1, j, (*rhs_mtx)(j));
					else
						ReplaceMatrixElementHelper(new_matrix, index_1, j, rhs_mtx->z(j));
				}
			}
			else
			{
				for (int j=0; j<data->N(); j++)
				{
					if (index_1 < 0)
						throw OML_Error(HW_ERROR_INDEXRANGE);

					if (value.IsScalar())
					{
						ReplaceMatrixElementHelper(new_matrix, index_1, j, value.Scalar());
					}
					else if (value.IsComplex())
					{
						ReplaceMatrixElementHelper(new_matrix, index_1, j, value.Complex());
					}
				}
			}
		}
		else if (index1.IsVector())
		{
			std::vector<double> indices = index1.Vector();

			if (value.IsScalar())
			{
				for (int j=0; j<indices.size(); j++)
				{
					if ((indices[j]-1 < 0) || (indices[j]-1 >= new_matrix->M()))
						throw OML_Error(HW_ERROR_INDEXRANGE);

					for (int k=0; k<new_matrix->N(); k++)
						(*new_matrix)((int)indices[j]-1,k) = value.Scalar();
				}
			}
			else if (value.IsComplex())
			{
				new_matrix->MakeComplex();

				for (int j=0; j<indices.size(); j++)
				{
					if ((indices[j]-1 < 0) || (indices[j]-1 >= new_matrix->M()))
						throw OML_Error(HW_ERROR_INDEXRANGE);

					for (int k=0; k<new_matrix->N(); k++)
						new_matrix->z((int)indices[j]-1,k) = value.Complex();
				}
			}
			else if (value.IsEmpty())
			{
				// This only works if the values are strictly increasing -- JDS
				for (int j=(int)indices.size()-1; j>=0; j--)
				{
					if ((indices[j]-1 < 0) || (indices[j]-1>= new_matrix->M()))
						throw OML_Error(HW_ERROR_INDEXRANGE);

					hwMathStatus stat = new_matrix->DeleteRows((int)indices[j]-1, 1);
				}
			}
			else if (value.IsMatrix())
			{
				const hwMatrix* val_matrix = value.Matrix();

				if (val_matrix->M() != indices.size())
					throw OML_Error(HW_ERROR_INCOMPDIM); 

				if (val_matrix->N() != new_matrix->N())
					throw OML_Error(HW_ERROR_INCOMPDIM); 

   				for (int j=0; j<indices.size(); j++)
    			{
					int index1 = (int)indices[j];

    				for (int k=0; k<val_matrix->N(); k++)
	    			{
		    			int index2 = k+1;

						if (!val_matrix->IsReal())
							ReplaceMatrixElementHelper(new_matrix, index1-1, index2-1, val_matrix->z(j, k));
						else
							ReplaceMatrixElementHelper(new_matrix, index1-1, index2-1, (*val_matrix)(j, k));
					}
				}
		    }
		}
	}
	else
	{
		if (value.IsMatrix() || value.IsScalar() || value.IsComplex())
		{
			const hwMatrix* val_mtx = NULL;

			if (value.IsMatrix())
				val_mtx = value.Matrix();

			bool is_real = true;

			if (val_mtx && !val_mtx->IsReal())
				is_real = false;

			if (value.IsComplex())
				is_real = false;

			if (index1.IsPositiveVector() && index2.IsPositiveInteger())
			{
				std::vector<double> indices = index1.Vector();
			
				if (val_mtx && val_mtx->Size() != indices.size())
					throw OML_Error(HW_ERROR_INCOMPDIM);//might be dimension mismatch

				if (is_real)
				{
					double new_val;

					for (int j=0; j<(int)indices.size(); j++)
					{
						if (val_mtx)
							new_val = (*val_mtx)(j);
						else
							new_val = value.Scalar();

						ReplaceMatrixElementHelper(new_matrix, (int)indices[j]-1, (int)index2.Scalar()-1, new_val);	
					}
				}
				else
				{
					hwComplex new_val;

					for (int j=0; j<(int)indices.size(); j++)
					{
						if (val_mtx)
							new_val = val_mtx->z(j);
						else
							new_val = value.Complex();

						ReplaceMatrixElementHelper(new_matrix, (int)indices[j]-1, (int)index2.Scalar()-1, new_val);
					}
				}
			}
			else if (index1.IsPositiveInteger() && index2.IsPositiveVector())
			{
				std::vector<double> indices = index2.Vector();
			
				if (val_mtx && val_mtx->Size() != indices.size())
					throw OML_Error(HW_ERROR_INCOMPDIM);//might be mismatch
				
				if (is_real)
				{
					double new_val;

					for (int j=0; j<indices.size(); j++)
					{
						if (val_mtx)
							new_val = (*val_mtx)(j);
						else
							new_val = value.Scalar();

						ReplaceMatrixElementHelper(new_matrix, (int)index1.Scalar()-1, (int)indices[j]-1, new_val);					
					}
				}
				else
				{
					hwComplex new_val;

					for (int j=0; j<indices.size(); j++)
					{
						if (val_mtx)
							new_val = val_mtx->z(j);
						else
							new_val = value.Complex();

						ReplaceMatrixElementHelper(new_matrix, (int)index1.Scalar()-1, (int)indices[j]-1, new_val);
					}
				}
			}
			else if (index1.IsPositiveVector() && index2.IsPositiveVector())
			{
				std::vector<double> indices1 = index1.Vector();
				std::vector<double> indices2 = index2.Vector();		

				if (val_mtx && val_mtx->M() != indices1.size())
					throw OML_Error(HW_ERROR_INCOMPDIM);//might be mismatch

				if (val_mtx && val_mtx->N() != indices2.size())
					throw OML_Error(HW_ERROR_INCOMPDIM);//might be mismatch

				if (is_real)
				{
					double new_val;

					int i1_size = (int)indices1.size();
					int i2_size = (int)indices2.size();

					for (int k=0; k<i2_size; k++)
					{
						int index_2 = (int)indices2[k]-1;

						if (index_2 < 0)
							throw OML_Error(HW_ERROR_INDEXRANGE);

    					for (int j=0; j<i1_size; j++)
	    				{
		    				int index_1 = (int)indices1[j]-1;

			    			if (index_1 < 0)
				    			throw OML_Error(HW_ERROR_INDEXRANGE);

				    		if (val_mtx)
					    		new_val = (*val_mtx)(j, k);
						    else
		    					new_val = value.Scalar();

							if ((index_1 >= new_matrix->M()) || (index_2 >= new_matrix->N()))
								ReplaceMatrixElementHelper(new_matrix, index_1, index_2, new_val);	
							else if (new_matrix->IsReal())
								(*new_matrix)(index_1, index_2) = new_val;
                            else
								new_matrix->z(index_1, index_2) = new_val;
						}
					}				
				}
				else
				{
					hwComplex new_val;

					for (int k=0; k<indices2.size(); k++)
					{
    					for (int j=0; j<indices1.size(); j++)
	    				{
							if (val_mtx)
								new_val = val_mtx->z(j, k);
							else
								new_val = value.Complex();

							ReplaceMatrixElementHelper(new_matrix, (int)indices1[j]-1, (int)indices2[k]-1, new_val);		
						}
					}
				}

			}
			else if (index1.IsVector() && index1.IsLogical() && index2.IsVector() && index2.IsLogical())
			{
				const hwMatrix* idx1 = index1.Matrix();
				const hwMatrix* idx2 = index2.Matrix();

				if (value.IsScalar())
				{
					for (int j=0; j<idx1->Size(); j++)
					{
						if ((*idx1)(j))
						{
							for (int k=0; k<idx2->Size(); k++)
							{
								if ((*idx2)(k))
									ReplaceMatrixElementHelper(new_matrix, j, k, value.Scalar());
							}
						}
					}
				}
			}
		}
	}

	return new_matrix;
}

void ExprTreeEvaluator::ReplaceMatrixElementHelper(hwMatrix*& target, const int index1, const int index2, const double value)
{
	if (index1 < 0)
		throw OML_Error(HW_ERROR_INDEXRANGE);

	if (index2 < 0)
		throw OML_Error(HW_ERROR_INDEXRANGE);

	int target_M = target->M();
	int target_N = target->N();

	if ((index1 >= target_M) || (index2 >= target_N))
	{
		int size1 = target_M;
		int size2 = target_N;

		if (index1 >= target_M)
			size1 = index1+1;

		if (index2 >= target_N)
			size2 = index2+1;

		if (target->IsEmpty())
		{
			target->Dimension(size1, size2, hwMatrix::REAL);
			target->SetElements(0.0);
		}
		else
		{
			target->Resize(size1, size2, true); // make all the newly-created elements 0
		}
	}

	if (target->IsReal())
		(*target)(index1, index2) = value;
	else
		target->z(index1, index2) = value;
}

void ExprTreeEvaluator::ReplaceMatrixElementHelper(hwMatrix*& target, int index1, int index2, const hwComplex& value)
{
	if (index1 < 0)
		throw OML_Error(HW_ERROR_INDEXRANGE);

	if (index2 < 0)
		throw OML_Error(HW_ERROR_INDEXRANGE);

	if ((index1 >= target->M()) || (index2 >= target->N()))
	{
		int size1 = target->M();
		int size2 = target->N();

		if (index1 >= target->M())
			size1 = index1+1;

		if (index2 >= target->N())
			size2 = index2+1;

		if (target->IsEmpty())
		{
			target->Dimension(size1, size2, hwMatrix::COMPLEX);
			target->SetElements(0.0);
		}
		else
		{
			target->Resize(size1, size2, true); // make all the newly-created elements 0
		}
	}

	if (target->IsReal())
		target->MakeComplex();

	target->z(index1, index2) = value;
}

void ExprTreeEvaluator::ReplaceMatrixElementHelper(hwMatrix*& target, int index1, double value)
{
	if (index1 < 0)
		throw OML_Error(HW_ERROR_INDEXRANGE);

	if (target->IsEmpty())
	{
		target->Dimension(1, index1+1, hwMatrix::REAL);
		target->SetElements(0.0);
	}

	if (index1 >= target->Size())
		target->Resize(1, index1+1, true); // make all the newly-created elements 0

	if (target->IsReal())
		(*target)(index1) = value;
	else
		target->z(index1) = value;
}

void ExprTreeEvaluator::ReplaceMatrixElementHelper(hwMatrix*& target, int index1, const hwComplex& value)
{
	if (index1 < 0)
		throw OML_Error(HW_ERROR_INDEXRANGE);

	if (target->IsEmpty())
	{
		target->Dimension(1, index1+1, hwMatrix::COMPLEX);
		target->SetElements(0.0);
	}

	if (index1 >= target->Size())
		target->Resize(1, index1+1, true); // make all the newly-created elements 0

	if (target->IsReal())
		target->MakeComplex();

	target->z(index1) = value;
}

void ExprTreeEvaluator::ValidateRowIndex(const hwMatrix& target, int index)
{
	if ((index < 0) || (index >= target.M()))
		throw OML_Error(HW_ERROR_INDEXRANGE);
}

void ExprTreeEvaluator::ValidateColumnIndex(const hwMatrix& target, int index)
{
	if ((index < 0) || (index >= target.N()))
		throw OML_Error(HW_ERROR_INDEXRANGE);
}

void ExprTreeEvaluator::ValidateExtractionIndices(const hwMatrix& target, int index1, int index2)
{
	if ((index1 < 0) || (index1 >= target.M()))
		throw OML_Error(HW_ERROR_INDEXRANGE);

	if ((index2 < 0) || (index2 >= target.N()))
		throw OML_Error(HW_ERROR_INDEXRANGE);
}

std::vector<Currency> ExprTreeEvaluator::GetOutputResults() const
{
	return results;
}

void ExprTreeEvaluator::ClearOutputResults()
{
	results.clear();
	last_suppressed_result = Currency(-1, Currency::TYPE_NOTHING);
}

Currency ExprTreeEvaluator::GetLastResult() const
{
	if (results.size())
	{
		return results[results.size() - 1];
	}
	else if (!last_suppressed_result.IsNothing())
	{
		return last_suppressed_result;
	}
	return Currency(-1, Currency::TYPE_NOTHING);
}
//------------------------------------------------------------------------------
//! Pushes result
//! \param[in] cur       Given data
//! \param[in] to_output True if results need to be sent to output
void ExprTreeEvaluator::PushResult(const Currency& cur, bool to_output)
{
	if (cur.IsNothing()) return;

	MemoryScope* scope = msm->GetCurrentScope();
	if (scope && scope->IsAnonymous() && !cur.IsError()) 
		return;

	if (cur.GetOutputName().empty() && !cur.IsError())
		SetValue("ans", cur);

	if (to_output || cur.IsFormat())
	{
		if (cur.IsObject())
		{
			ClassInfo* ci = (*class_info_map)[cur.GetClassname()];

			FunctionInfo* fi = ci->GetFunctionInfo("disp");

			if (fi)
			{
				std::vector<Currency> inputs;
				std::vector<Currency> outputs;
				inputs.push_back(cur);
				CallInternalFunction(fi, inputs); 
				results.push_back(cur);
				return;
			}

		}
		else if (cur.IsCellList())
		{
			HML_CELLARRAY* cells = cur.CellArray();

			for (int k=0; k<cells->Size(); ++k)
			{
				Currency temp = (*cells)(k);
				
				if (temp.GetOutputName().empty())
					SetValue("ans", temp);

				results.push_back(temp);

                ProcessResultForPrint(temp, to_output);
			}
            return;
		}
		results.push_back(cur);
		_store_suppressed = false;
        ProcessResultForPrint(cur, to_output);
	}
	else if (!cur.IsNothing() && _store_suppressed)
	{
		last_suppressed_result = cur;
	}
}

//------------------------------------------------------------------------------
//! Called when printing result
//! \param[in] cur Currency to be printed
//------------------------------------------------------------------------------
void ExprTreeEvaluator::PrintResult(Currency cur)
{
	if (cur.IsNothing()) return;

	results.push_back(cur);
    ProcessResultForPrint(cur, true);
}

// assumes index is valid
bool ExprTreeEvaluator::CloseFile(int i)
{
	UserFile &ufi = (*userFileStreams)[i];
	if(fclose(ufi.file)) // 0 for success
		return false;

	ufi.isopen = false;
	ufi.file   = nullptr;
	ufi.name   = std::string();
	ufi.mode   = std::string();

	return true;
}

bool ExprTreeEvaluator::CloseAllFiles()
{
	bool returnval = true;
	for (int i = FIRST_USER_FILE; i < userFileStreams->size(); i++)
	{
		if ((*userFileStreams)[i].isopen && !CloseFile(i))
			returnval = false;
	}
	return returnval;
}

int ExprTreeEvaluator::AddFile(std::FILE *newfile, const std::string &fname, const std::string &fmode)
{
	UserFile ufi(newfile, fname, fmode);
	int i;
	for (i = 0; i < userFileStreams->size(); i++)
	{
		// if the index isn't currently being used
		if (!(*userFileStreams)[i].isopen)
		{
			(*userFileStreams)[i] = ufi;
			return i;
		}
	}
	userFileStreams->push_back(ufi);
	return i;
}

bool ExprTreeEvaluator::PathMatches(const std::string& s1, const std::string& s2) const
{
#ifdef OS_WIN

	std::string ls1 = s1;
	std::string ls2 = s2;
	std::transform(s1.cbegin(), s1.cend(), ls1.begin(), &tolower);
	std::transform(s2.cbegin(), s2.cend(), ls2.begin(), &tolower);
	std::replace(ls1.begin(), ls1.end(), '/', '\\');
	std::replace(ls2.begin(), ls2.end(), '/', '\\');
	return ls1 == ls2;
#else
	return s1 == s2;
#endif
}

void ExprTreeEvaluator::AddPath(std::string pathname, bool end)
{
	// avoid repeats
	if (!IsInPaths(pathname))
	{
		if (end)
			paths->push_back(pathname);
		else
			paths->insert(paths->begin() + NUM_MANDATORY_PATHS, pathname);
		ResetFuncSearchCache();
	}
}

void ExprTreeEvaluator::AddPath2(const std::string& pathname, const std::vector<std::string>& func_names)
{
	AddPath(pathname, true);

	std::vector<std::string>::const_iterator iter;

	for (iter = func_names.begin(); iter != func_names.end(); iter++)
	{
		// don't allow duplicates
		if (std::find(preregistered_functions->begin(), preregistered_functions->end(), *iter) == preregistered_functions->end())
			preregistered_functions->push_back(*iter);
	}

	OnUpdateFuncList();
}

bool ExprTreeEvaluator::RemovePath(std::string &pathname)
{
	std::vector<std::string>::iterator iter;
	for (iter = paths->begin() + NUM_MANDATORY_PATHS; iter != paths->end(); iter++)
	{
		if (PathMatches(pathname, *iter))
		{
			std::string actual_path = *iter;

			paths->erase(iter);

			RemoveFunctionsInPath(actual_path);

			return true;
		}
	}
	return false;
}

bool ExprTreeEvaluator::FindFileInPath(const std::string& file_plus_ext, std::string& filepath) const
{
	std::string cur_path = BuiltInFuncsUtils::GetCurrentWorkingDir();
	std::replace(cur_path.begin(), cur_path.end(), '\\', '/'); // slash direction is important for the debugger
	paths->insert(paths->begin(), cur_path);

	bool ret_val = false;

	std::vector<std::string>::iterator iter;
	for (iter = paths->begin(); iter != paths->end(); iter++)
	{
		std::string temp = *iter + "/" + file_plus_ext; // slash direction is important for the debugger

		if (FileExists(temp.c_str()))
		{
			filepath = temp;
			ret_val = true;
			break;
		}
	}

	paths->erase(paths->begin(), paths->begin()+1);

	return ret_val;
}

bool ExprTreeEvaluator::IsInPaths(const std::string &str) const
{
	std::vector<std::string>::const_iterator iter;
	for (iter = paths->cbegin(); iter != paths->cend(); iter++)
	{
		if (PathMatches(str, *iter))
			return true;
	}
	return false;
}

void ExprTreeEvaluator::RestorePath()
{
	paths->clear();

    if (!_appdir.empty())
		paths->push_back(_appdir + "/scripts/oml");
}

std::vector<int> ExprTreeEvaluator::GetFileIndices(int start)
{
	std::vector<int> indices;

	for (int i = start; i < userFileStreams->size(); i++)
	{
		if ((*userFileStreams)[i].isopen)
			indices.push_back(i);
	}
	return indices;
}

void ExprTreeEvaluator::GetStackInfo(int level, std::string& filename, int& line_number)
{
	MemoryScope* scope = msm->GetScope(level);
	filename           = scope->GetFilename();
	line_number        = scope->GetLineNumber();
}

bool DebugInfo::Matches(OMLTree* tree) const
{
	DebugInfo temp = DebugInfoFromTree(tree);

	if (temp.lineNum != lineNum)
		return false;

	if (temp.fileName != fileName)
		return false;

	return true;
}

DebugInfo DebugInfo::DebugInfoFromTree(OMLTree* tree)
{
	int lineNum = tree->Line();
	const std::string* fileName = NULL;

	if (lineNum != 0)
	{		
		fileName = tree->FilenamePtr();

		return DebugInfo(lineNum, fileName);
	}
	else
	{
		//if the line number is 0 then we have to go down the child tree until we find it.
		OMLTree*    stmt_child;

		int tree_child_count = tree->ChildCount();

		for (int i = 0; i < tree_child_count; i++)
		{
			stmt_child = tree->GetChild(i);
		
			if (stmt_child && stmt_child->Line() != 0)
			{
				lineNum = stmt_child->Line();
				fileName = stmt_child->FilenamePtr();
				break;
			}
		}

		if (lineNum == 0) // if it's still 0, try the grandchildren
		{
			OMLTree*    stmt = tree->GetChild(0);
			OMLTree*    stmt_child;
	
			for (int i = 0; i < stmt->ChildCount();i++)
			{
				stmt_child = stmt->GetChild(i);
		
				if (stmt_child && stmt_child->Line() != 0)
				{
					lineNum = stmt_child->Line();
					fileName = stmt_child->FilenamePtr();
					break;
				}
			}
		}

		if (lineNum == 0) // if it's still 0, try the great-grandchildren
		{
			OMLTree*    xstmt = tree->GetChild(0);

			if (xstmt->ChildCount() > 0)
			{
				OMLTree*    stmt  = xstmt->GetChild(0);
				OMLTree*    stmt_child;

				for (int i = 0; i < stmt->ChildCount();i++)
				{
					stmt_child = stmt->GetChild(i);
		
					if (stmt_child && stmt_child->Line() != 0)
					{
						lineNum = stmt_child->Line();
					    fileName = stmt_child->FilenamePtr();
						break;
					}
				}
			}
		}

		if (!fileName)
			fileName = Currency::pm.GetStringPointer("Unknown");

		return DebugInfo(lineNum, fileName);
	}
}

std::string ExprTreeEvaluator::FormatErrorMessage(const std::string& base_message)
{
	MemoryScope* scope = GetCurrentScope();
	std::string error_str = base_message;

	std::string error_file = scope->GetFilename();
	int         error_line = scope->GetLineNumber();

	if (scope)
	{
		if (builtin_error_scope.empty())
		{
			if (!scope->FunctionName().empty())
			{
				error_str += " in function ";
				error_str += scope->FunctionName();
			}
		}
		else
		{
			error_str += " in call to function ";
			error_str += builtin_error_scope;

			if (cached_filename.size())
			{
				error_file = cached_filename;
				error_line = cached_line;

				UncacheLineInfomation();
			}
		}

		std::string possible_full_path = error_file;

		error_str += " at line number ";
		std::ostringstream convertLine;
		convertLine << error_line;
		std::string lineNumberStr = convertLine.str();
		error_str += lineNumberStr;

		if (possible_full_path != "dummy")
		{
			error_str += " in file ";
			std::replace(possible_full_path.begin(), possible_full_path.end(), '/', '\\');

			size_t slash_pos = possible_full_path.rfind('\\');
			if (slash_pos == std::string::npos)
			{
				error_str += error_file;
			}
			else
			{
				std::string filename = possible_full_path.substr(slash_pos+1);
				error_str += filename;
			}
		}
	}
	return error_str;
}

std::string ExprTreeEvaluator::GetLastErrorMessage()
{
	return lasterrormsg;
}
	
void ExprTreeEvaluator::SetLastErrorMessage(const std::string& lasterr)
{
	lasterrormsg = lasterr;
}

std::string ExprTreeEvaluator::GetLastWarning()
{
	return lastwarning;
}
	
void ExprTreeEvaluator::SetLastWarning(const std::string& lastwarn)
{
	lastwarning = lastwarn;
}

void ExprTreeEvaluator::ErrorCleanup()
{
	current_tree = NULL;
}

ExprTreeEvaluator* ExprTreeEvaluator::MakeContextCopy(bool base, bool pop_nargs) const
{
	ExprTreeEvaluator* eval = new ExprTreeEvaluator(this);
	eval->nargout_values = nargout_values;
	eval->nargin_values = nargin_values;

    if (pop_nargs && nargin_values.size() && nargout_values.size())
        eval->PopNargValues();

	eval->msm = msm->MakeContextCopy(base);
	eval->_owns_msm = true;
	eval->is_for_evalin = true;
	return eval;
}

static void CheckDims(int m, int n)
{
    if (m < 0 || n < 0)
        throw OML_Error("Error: invalid dimensions for matrix, must be positive and no larger than " + std::to_string((long long) INT_MAX));
}

template <typename T1, typename T2>
static void CheckSize(const hwTMatrix<T1, T2>* mtx, int expected_size)
{
    if (mtx->Size() != expected_size)
    {
        delete mtx;
        throw OML_Error(HW_ERROR_OUTMEM);
    }
}

hwMatrix* ExprTreeEvaluator::allocateMatrix()
{
	return new hwMatrix;
}

hwMatrix* ExprTreeEvaluator::allocateMatrix(const hwMatrix* mtx)
{
    if (!mtx)
        return allocateMatrix();

	hwMatrix* newmtx = new hwMatrix(*mtx);
    CheckSize(newmtx, mtx->Size());
	return newmtx;
}

hwMatrix* ExprTreeEvaluator::allocateMatrix(int m, int n, void* data, hwMatrix::DataType type)
{
    CheckDims(m,n);
	hwMatrix* mtx = new hwMatrix(m, n, data, type);
    CheckSize(mtx, m*n);
	return mtx;
}

hwMatrix* ExprTreeEvaluator::allocateMatrix(int m, int n, hwMatrix::DataType type)
{
    CheckDims(m,n);
	hwMatrix* mtx = new hwMatrix(m, n, type);
    CheckSize(mtx, m*n);
	return mtx;
}

hwMatrix* ExprTreeEvaluator::allocateMatrix(int m, int n, double value)
{
    CheckDims(m,n);
	hwMatrix* mtx = ExprTreeEvaluator::allocateMatrix(m, n, hwMatrix::REAL);
    CheckSize(mtx, m*n);
	mtx->SetElements(value);
	return mtx;
}

hwMatrix* ExprTreeEvaluator::allocateMatrix(int m, int n, hwComplex& value)
{
    CheckDims(m,n);
	hwMatrix* mtx = ExprTreeEvaluator::allocateMatrix(m, n, hwMatrix::COMPLEX);
    CheckSize(mtx, m*n);
	mtx->SetElements(value);
	return mtx;
}

hwMatrixI* ExprTreeEvaluator::allocateMatrix(int m, int n, hwMatrixI::DataType type)
{
    CheckDims(m,n);
	hwMatrixI* mtx = new hwMatrixI(m, n, type);
    CheckSize(mtx, m*n);
	return mtx;
}

hwMatrixI* ExprTreeEvaluator::allocateMatrix(int m, int n, int val)
{
    CheckDims(m,n);
	hwMatrixI* mtx = ExprTreeEvaluator::allocateMatrix(m, n, hwMatrixI::REAL);
    CheckSize(mtx, m*n);
	mtx->SetElements(val);
	return mtx;
}

const hwMatrix* ExprTreeEvaluator::allocateColumn(const hwMatrix* mtx, int col)
{
	const hwMatrix* mtxcol;
	if (mtx->IsReal())
	{
		const double* matrixColumnPtr = &((*mtx)(0, col));
		mtxcol = new hwMatrix(mtx->M(), (void*) matrixColumnPtr, hwMatrix::REAL);
	}
	else
	{
		const hwComplex* matrixColumnPtr = &(mtx->z(0, col));
		mtxcol = new hwMatrix(mtx->M(), (void*) matrixColumnPtr, hwMatrix::COMPLEX);
	}

    CheckSize(mtxcol, mtx->M());
	return mtxcol;
}

hwMatrixN* ExprTreeEvaluator::allocateMatrixN()
{
	return new hwMatrixN;
}

hwMatrixN* ExprTreeEvaluator::allocateMatrixN(const std::vector<int>& dims, const hwMatrixN::DataType& dataType)
{
	return new hwMatrixN(dims, dataType);
}

hwMatrixN* ExprTreeEvaluator::allocateMatrixN(const hwMatrixN* in)
{
	return new hwMatrixN(*in);
}

HML_CELLARRAY* ExprTreeEvaluator::allocateCellArray()
{
	return new HML_CELLARRAY;
}

HML_CELLARRAY* ExprTreeEvaluator::allocateCellArray(int m, int n)
{
    CheckDims(m,n);
	HML_CELLARRAY* cell = new HML_CELLARRAY(m, n, HML_CELLARRAY::REAL);
    CheckSize(cell, m*n);
	return cell;
}

HML_CELLARRAY* ExprTreeEvaluator::allocateCellArray(const HML_CELLARRAY* cell)
{
	HML_CELLARRAY* newcell = new HML_CELLARRAY(*cell);
    CheckSize(cell, cell->Size());
	return newcell;
}

StructData* ExprTreeEvaluator::allocateStruct(const StructData* strct)
{
	StructData* newstruct = new StructData(*strct);
	if (newstruct->M() != strct->M() || newstruct->N() != strct->N())
	{
		delete newstruct;
		throw OML_Error(HW_ERROR_OUTMEM);
	}
	return newstruct;
}

StructData* ExprTreeEvaluator::allocateStruct()
{
	StructData* newstruct = new StructData();
	return newstruct;
}

bool ExprTreeEvaluator::HasBuiltin(const std::string& func_name) const
{
	return std_functions->count(func_name) != 0;
}

int ExprTreeEvaluator::NargoutFor(const std::string& func_name) const
{
	return std_functions->at(func_name).md.nargout;
}

int ExprTreeEvaluator::NarginFor(const std::string& func_name) const
{
	return std_functions->at(func_name).md.nargin;
}

Currency ExprTreeEvaluator::CheckForExternalVariable(const std::string& varname)
{
	if (ext_var_ptr)
	{
		Currency my_temp = (*ext_var_ptr)(varname);

		if (!my_temp.IsNothing())
		{
			_externals.push_back(my_temp);
			return _externals.back();
		}
	}

	return Currency(-1.0, Currency::TYPE_NOTHING);
}

FunctionInfo* ExprTreeEvaluator::FunctionInfoFromString(const std::string& str)
{
	pANTLR3_INPUT_STREAM input = ANTLRData::InputFromExpression(str, "dummy");
	ANTLRData ad(input, false);
	pExprCppTreeParser parser = ad.GetParser();
	ExprCppTreeParser_prog_return r = parser->prog(parser);

	OMLTree* oml_tree = OMLTree::ConvertTree(r.tree);

	OMLTree* temp1 = oml_tree->GetChild(0);
	OMLTree* temp2 = temp1->GetChild(0);
	FunctionInfo*     fi    = FunctionInfoFromTree(temp2);

	MemoryScope* temp = new MemoryScope(*GetCurrentScope(), fi);
	fi->SetAnonymous(temp);

	return fi;
}

bool ExprTreeEvaluator::IsKeyword(const std::string& func_name) const
{
	std::vector<std::string> keywords = GetKeywords();
	
	std::vector<std::string>::const_iterator iter;

	for (iter = keywords.begin(); iter != keywords.end(); iter++)
	{
		if (*iter == func_name)
			return true;
	}

	return false;
}

bool ExprTreeEvaluator::IsOperator(const std::string& func_name) const
{
	std::vector<std::string> operators = GetOperators();
	
	std::vector<std::string>::const_iterator iter;

	for (iter = operators.begin(); iter != operators.end(); iter++)
	{
		if (*iter == func_name)
			return true;
	}

	return false;
}
//------------------------------------------------------------------------------
//! Sets signal handler, if not null or to null if already set
//! \param[in] handler Signal handler
//------------------------------------------------------------------------------
void ExprTreeEvaluator::SetSignalHandler(SignalHandlerBase* handler)
{
    // Set to the new handler only if not set before or the new handler is 
    // clearing the original signal handler to null
    if (!_signalHandler || !handler)
        _signalHandler = handler;
}
//------------------------------------------------------------------------------
//! Prompts to save before exiting
//------------------------------------------------------------------------------
void ExprTreeEvaluator::OnSaveOnExit()
{
    if (_signalHandler)
        _signalHandler->OnSaveOnExitHandler();
}
//------------------------------------------------------------------------------
//! Clears display results in the client
//------------------------------------------------------------------------------
void ExprTreeEvaluator::OnClearResults()
{
    if (_signalHandler)
        _signalHandler->OnClearResultsHandler();
}
//------------------------------------------------------------------------------
//! Gets user input
//! \param[in]  prompt Prompt to display to user
//! \param[in]  type   Type, if specified
//! \param[out] input  Input from user
//------------------------------------------------------------------------------
void ExprTreeEvaluator::OnUserInput(const std::string& prompt,
                                       const std::string& type,
                                       std::string&       input)
{
    if (_signalHandler)
        _signalHandler->OnUserInputHandler(prompt, type, input);
}
//------------------------------------------------------------------------------
//! Refreshes directories in client
//------------------------------------------------------------------------------
void ExprTreeEvaluator::OnRefreshDirs()
{
    if (_signalHandler)
        _signalHandler->OnRefreshDirsHandler();
}
//------------------------------------------------------------------------------
//! Change current working directory in client
//! \param[in] dir Fully qualified path of the new directory
//------------------------------------------------------------------------------
void ExprTreeEvaluator::OnChangeDir(const std::string& dir)
{
    if (_signalHandler)
        _signalHandler->OnChangeDirHandler(dir);
}
//------------------------------------------------------------------------------
//! Start pause
//! \param[in] msg  User message to display
//! \param[in] wait True if waiting for a keystroke input from user
//------------------------------------------------------------------------------
void ExprTreeEvaluator::OnPauseStart(const std::string& msg, bool wait)
{
    if (_signalHandler)
        _signalHandler->OnPauseStartHandler(msg, wait);
}
//------------------------------------------------------------------------------
//! End pause
//------------------------------------------------------------------------------
void ExprTreeEvaluator::OnPauseEnd()
{
    if (_signalHandler)
        _signalHandler->OnPauseEndHandler();
}
//------------------------------------------------------------------------------
//! Update function list in language
//------------------------------------------------------------------------------
void ExprTreeEvaluator::OnUpdateFuncList()
{
    if (_signalHandler && !_suspendFunclistUpdate)
        _signalHandler->OnUpdateFuncListHandler();
}
//------------------------------------------------------------------------------
//! Called when result needs to be processed for print - sets format/prints
//! \param[in] cur      Currency to print
//! \param[in] toOutput True if this has an output
//------------------------------------------------------------------------------
void ExprTreeEvaluator::ProcessResultForPrint(const Currency& cur, bool toOutput)
{
    if (cur.IsFormat())
    {
        if (cur.Format())
            SetOutputFormat(*cur.Format());
        return;
    }
    
    if (toOutput)  
        OnPrintResult(cur); // Print
}

std::string ExprTreeEvaluator::GetCurrentDebugFile() const
{
	return msm->GetCurrentScope()->GetFilename();
}

int ExprTreeEvaluator::GetCurrentDebugLine() const
{
	return msm->GetCurrentScope()->GetLineNumber();
}

void ExprTreeEvaluator::SetDebugInfo(const std::string* file, int line)
{
	msm->GetCurrentScope()->SetDebugInfo(file, line);
}
//------------------------------------------------------------------------------
//! Print the given currency
//! \param[in] cur Currency to print
//------------------------------------------------------------------------------
void ExprTreeEvaluator::OnPrintResult(const Currency& cur)
{
    if (_signalHandler)
        _signalHandler->OnPrintResultHandler(cur);
}

// End of file:

// fenv functions
int ExprTreeEvaluator::GetBaseEnvHandle()
{
	return msm->GetBaseEnvHandle();
}

int ExprTreeEvaluator::GetCurrentEnvHandle()
{
	return msm->GetCurrentEnvHandle();
}

int ExprTreeEvaluator::GetNewEnvHandle()
{
	return msm->GetNewEnvHandle();
}

Currency ExprTreeEvaluator::GetEnvValue(int handle, std::string varname)
{
	return msm->GetEnvValue(handle, varname);
}

void ExprTreeEvaluator::RemoveEnvValue(int handle, std::string varname)
{
	return msm->RemoveEnvValue(handle, varname);
}

void ExprTreeEvaluator::SetEnvValue(int handle, std::string varname, const Currency& new_val)
{
	return msm->SetEnvValue(handle, varname, new_val);
}

void ExprTreeEvaluator::ImportEnv(int source, int dest)
{
	return msm->ImportEnv(source, dest);
}

Currency ExprTreeEvaluator::CallOverloadedOperator(const std::string& op_name, const Currency& op)
{
	std::string class_name = op.GetClassname();
	ClassInfo*  ci         = (*class_info_map)[class_name];

	FunctionInfo* fi = ci->GetFunctionInfo(op_name);

	if (fi)
	{
		std::vector<Currency> inputs;
		inputs.push_back(op);
		return CallInternalFunction(fi, inputs);
	}
	else
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
	}
}

Currency ExprTreeEvaluator::CallOverloadedOperator(const std::string& op_name, const Currency& op1, const Currency& op2)
{
	std::string class_name;
	
	if (op1.IsObject())
		class_name = op1.GetClassname();
	else if (op2.IsObject())
		class_name = op2.GetClassname();

	ClassInfo*  ci = (*class_info_map)[class_name];

	FunctionInfo* fi = ci->GetFunctionInfo(op_name);

	if (fi)
	{
		std::vector<Currency> inputs;
		inputs.push_back(op1);
		inputs.push_back(op2);
		return CallInternalFunction(fi, inputs);
	}
	else
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
	}
}

bool ExprTreeEvaluator::HasOverloadedFunction(const Currency& obj, const std::string& func_name)
{
	std::string class_name = obj.GetClassname();
	ClassInfo*  ci         = (*class_info_map)[class_name];

	FunctionInfo* fi = ci->GetFunctionInfo(func_name);

	if (fi)
		return true;
	else
		return false;
}

Currency ExprTreeEvaluator::CallOverloadedFunction(const std::string& func_name, const std::vector<Currency>& inputs)
{
	std::string class_name = inputs[0].GetClassname();
	ClassInfo*  ci         = (*class_info_map)[class_name];

	FunctionInfo* fi = ci->GetFunctionInfo(func_name);

	if (fi)
		return CallInternalFunction(fi, inputs);
	else
		throw OML_Error(HW_ERROR_UNSUPOP);
}
//------------------------------------------------------------------------------
//! Returns true if evaluator is in experimental mode
//------------------------------------------------------------------------------
bool ExprTreeEvaluator::GetExperimental() const
{
    return Currency::GetExperimental(); 
}
//------------------------------------------------------------------------------
//! Sets the experimental mode flag (-ex)
//! \param[in] val Sets to true if the evaluator is in experimental mode
//------------------------------------------------------------------------------
void ExprTreeEvaluator::SetExperimental(bool val)
{
    Currency::SetExperimental(val); 
}

void ExprTreeEvaluator::CacheLineInfomation()
{
	MemoryScope *ms = msm->GetCurrentScope();
	cached_filename = ms->GetFilename();
	cached_line     = ms->GetLineNumber();
}

void ExprTreeEvaluator::UncacheLineInfomation()
{
	cached_filename = "";
	cached_line     = 0;
}
//------------------------------------------------------------------------------
//! Returns true if successful in registering (swig) bound class
//! \param[in] name Bound class name
//! \param[in] info Bound class info
//------------------------------------------------------------------------------
bool ExprTreeEvaluator::RegisterBoundClass(const std::string& name,
                                           BoundClassInfo*    info)
{
    // Set warnings, if needed
    std::string warn;
    std::string warnprefix ("Warning: Cannot register class ; ");

    if (name.empty())
        warn = warnprefix + "name is empty ";

    else if (!info)
        warn = warnprefix + "information is empty for [" + name + "] ";

    if (!warn.empty())
    {
        std::string msg(FormatErrorMessage(warn));
        Currency tmp(msg);
        tmp.DispOutput();
        PrintResult(tmp);
        SetLastWarning(msg);
        return false;
    }

    _boundclassinfo[name] = info;
    OnUpdateFuncList();
    return true;
}

bool ExprTreeEvaluator::IsExtensionEncrypted(const std::string& extension)
{
	if (_decryptors.find(extension) == _decryptors.end())
		return false;

	return true;
}

void ExprTreeEvaluator::RunEncryptedFile(const std::string& extension, const std::string& file_name)
{
	if (_decryptors.find(extension) == _decryptors.end())
		throw OML_Error("Missing decryptor");

	ENCRPTR decrypt = _decryptors[extension];

	if (decrypt)
	{
		char* decrypted = (char*)decrypt(file_name);

		ParseAndRunString(decrypted, file_name);

		delete [] decrypted;
	}
}

void ExprTreeEvaluator::RegisterOMLDecryptor(const std::string& extension, ENCRPTR ptr)
{
	_decryptors[extension] = ptr;
}

void ExprTreeEvaluator::RunPrecompiledFile(const std::string& filename)
{
	OMLTree* tree = OMLTree::ReadTreeFromFile(filename);
	RUN(tree);
	delete tree;
}

void ExprTreeEvaluator::RunFile(const std::string& filename)
{
	pANTLR3_INPUT_STREAM input = ANTLRData::InputFromFilename(filename);

	ANTLRData ad(input, true);
	pANTLR3_COMMON_TOKEN_STREAM tokens = ad.GetTokens();                                         	

	ANTLRData::PreprocessTokenStream(tokens);
	pANTLR3_VECTOR vec = tokens->getTokens(tokens);

	bool keep_going = true;
	
	if (vec->size(vec) == 0)
		keep_going = false;

	OMLTree* oml_tree = NULL;

	std::vector<OMLTree*> statements;

	while (keep_going)
	{
		pExprCppTreeParser parser = ExprCppTreeParserNew(tokens);

		pANTLR3_COMMON_TOKEN tok = (pANTLR3_COMMON_TOKEN)vec->get(vec, tokens->p);

		ExprCppTreeParser_stmt_return r = parser->stmt(parser);
	
		if (parser->pParser->rec->getNumberOfSyntaxErrors(parser->pParser->rec) > 0)
		{		
			throw OML_Error("Invalid oml file");
			return;
		}
		else
		{
			pANTLR3_BASE_TREE tree = r.tree;

			PreprocessAST(tree, tokens);
			oml_tree = OMLTree::ConvertTree(tree);

			if (tok->input)
				oml_tree->SetDebugInfo((const char*)tok->input->fileName->chars, tok->line);

			if (oml_tree->GetType() == FUNC_DEF)
			{
				RunTree(oml_tree);
				delete oml_tree;
			}
			else
			{
				statements.push_back(oml_tree);
			}

			if (tokens->p == vec->size(vec))
				keep_going = false;
		}

		parser->free(parser);
	}

	std::string error;
	Currency    rr;
    bool formatmsg = true;

	try
	{
		if (input && input->fileName)
			SetScriptName((char*)input->fileName->chars);

		for (int j=0; j<statements.size(); j++)
		{
			oml_tree = statements[j];
			SetDebugInfo(oml_tree->FilenamePtr(), oml_tree->Line());
			rr = RunTree(oml_tree);
			delete oml_tree;
		}

	}
	catch (OML_Error& e)
	{
		error = e.GetErrorMessage();
        formatmsg = e.GetFormatMessage();
		ErrorCleanup();
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
		ErrorCleanup();
	}
	if (error.size())
		throw OML_Error(error, formatmsg);

	return;
}
//------------------------------------------------------------------------------
//! Sets application directory and restores paths
//------------------------------------------------------------------------------
void ExprTreeEvaluator::SetApplicationDir(const std::string& dir)
{
    _appdir = dir;
    RestorePath();
}

std::string ExprTreeEvaluator::IsValidString(const std::string& in) const
{
	pANTLR3_INPUT_STREAM input = ANTLRData::InputFromExpression(in, "dummy");
	ANTLRData ad(input, true);

	pANTLR3_COMMON_TOKEN_STREAM tokens = ad.GetTokens();
	ANTLRData::PreprocessTokenStream(tokens);

	ad.CreateParser(tokens);
	pExprCppTreeParser parser = ad.GetParser();

	ExprCppTreeParser_prog_return r = parser->prog(parser);

	if (parser->pParser->rec->getNumberOfSyntaxErrors(parser->pParser->rec) == 0)
		return "";

	char buffer[2048];
	sprintf(buffer, "Syntax error at line number %d near character position %d", parser->pParser->rec->state->exception->line, parser->pParser->rec->state->exception->charPositionInLine);
	return buffer;
}

TREE_FPTR ExprTreeEvaluator::GetFuncPointerFromType(int type)
{
	switch(type) 
	{
		case NUMBER:
			return &ExprTreeEvaluator::Number;
		case HEXVAL: 
			return &ExprTreeEvaluator::HexNumber;
		case HML_STRING:
			return &ExprTreeEvaluator::TextString;
		case IDENT: 
			return &ExprTreeEvaluator::Identifier;
		case PLUS:
		case MINUS:
		case TIMES:
		case ETIMES:
		case DIV:
		case EDIV:
		case LDIV:
		case ELDIV:
		case POW:
		case DOTPOW:
			return &ExprTreeEvaluator::BinaryOperator;
		case EQUAL: 
		case NEQUAL:
			return &ExprTreeEvaluator::EqualityOperator;
		case LTHAN:
			return &ExprTreeEvaluator::LessThanOperator;
		case GTHAN:
			return &ExprTreeEvaluator::GreaterThanOperator;
		case GEQ:
			return &ExprTreeEvaluator::GreaterEqualOperator;
		case LEQ:
			return &ExprTreeEvaluator::LessEqualOperator;
		case AND:
		case OR:
			return &ExprTreeEvaluator::LogicalOperator;
		case LAND:
			return &ExprTreeEvaluator::ShortCircuitAndOperator;
		case LOR:
			return &ExprTreeEvaluator::ShortCircuitOrOperator;
		case UMINUS:			
		case NEGATE:
			return &ExprTreeEvaluator::UnaryOperator;
		case ASSIGN: 
			return &ExprTreeEvaluator::AssignOperator;
		case FUNC:
			return &ExprTreeEvaluator::FunctionCall;
		case MR_FUNC:
			return &ExprTreeEvaluator::MultiReturnFunctionCall;
		case PARAM_LIST:
			return &ExprTreeEvaluator::Nothing;
		case CONDITIONAL:
			return &ExprTreeEvaluator::Conditional;
		case SWITCH: 
			return &ExprTreeEvaluator::SwitchCase;
		case IF: 
		case ELSE:
		case DUMMY: 
			return &ExprTreeEvaluator::Nothing;
		case VECTOR:
			return &ExprTreeEvaluator::Nothing;
		case WHILE: 
			return &ExprTreeEvaluator::WhileLoop;
		case FOR: 
			return &ExprTreeEvaluator::ForLoop;
		case STATEMENT_LIST: 
			return &ExprTreeEvaluator::StatementList;
		case STMT: 
			return &ExprTreeEvaluator::Statement;
		case FUNC_DEF: 
			return &ExprTreeEvaluator::FunctionDefinition;
		case FUNC_HANDLE:
			return &ExprTreeEvaluator::AnonymousFunctionDefinition;
		case MATRIX:
			return &ExprTreeEvaluator::MatrixCreation;
		case COLON:
			return &ExprTreeEvaluator::RangeOperator;
		case GLOBAL:
			return &ExprTreeEvaluator::GlobalReference;
		case PERSISTENT:
			return &ExprTreeEvaluator::PersistentReference;
		case CELL_ARRAY:
			return &ExprTreeEvaluator::CellArrayCreation;
		case CELL_VAL:
			return &ExprTreeEvaluator::CellValue;
		case CELL_ASSIGN:
			return &ExprTreeEvaluator::CellAssignment;
		case RETURN:
			return &ExprTreeEvaluator::Return;
		case TRANSP:
			return &ExprTreeEvaluator::TransposeOperator;
		case CTRANSP:
			return &ExprTreeEvaluator::ConjTransposeOperator;
		case INLINE_INDEX:
			return &ExprTreeEvaluator::InlineIndex;
		case INLINE_INDEX_CELL:
			return &ExprTreeEvaluator::InlineIndexCell;
		case STRUCT:
			return &ExprTreeEvaluator::StructValue;
		case TRY:
			return &ExprTreeEvaluator::TryCatch;
		case SEMIC:
		case COMMA:
			return &ExprTreeEvaluator::Nothing;
		case INPLACE:
			return &ExprTreeEvaluator::InPlaceExpansion;
		case CLASSDEF:
			return &ExprTreeEvaluator::ClassDefinition;	
	}

	return NULL;
}

void ExprTreeEvaluator::WritePFile(const std::string& infile, const std::string& outfile)
{
	pANTLR3_INPUT_STREAM input = ANTLRData::InputFromFilename(infile);
	ANTLRData ad(input, true);

	pANTLR3_COMMON_TOKEN_STREAM tokens = ad.GetTokens();
	ANTLRData::PreprocessTokenStream(tokens);

	ad.CreateParser(tokens);
	pExprCppTreeParser parser = ad.GetParser();

	ExprCppTreeParser_prog_return r = parser->prog(parser);

	if (parser->pParser->rec->getNumberOfSyntaxErrors(parser->pParser->rec) == 0)
	{
		pANTLR3_BASE_TREE tree = r.tree;
		PreprocessAST(tree, tokens);
		OMLTree* oml_tree = OMLTree::ConvertTree(tree);

		oml_tree->WriteTreeToFile(outfile);
	}
}

void AnalyzeHelper(OMLTree* tree, std::set<std::string>& required, std::set<std::string>& defined)
{
	int type = tree->GetType();

	if (type == ASSIGN)
	{
		OMLTree* assignee = tree->GetChild(0);

		if (assignee->GetType() == IDENT)
			defined.insert(assignee->GetText());

		OMLTree* rhs = tree->GetChild(1);

		AnalyzeHelper(rhs, required, defined);
	}
	else if (type == FUNC)
	{
		OMLTree* func_name = tree->GetChild(0);

		if (func_name->GetType() == IDENT)
			required.insert(func_name->GetText());
	}
	else if (type == IDENT)
	{
		std::string text = tree->GetText();

		if (defined.find(text) == defined.end())
			required.insert(tree->GetText());
	}
	else if (type == FUNC_HANDLE)
	{
		OMLTree* child   = tree->GetChild(0);
		std::string text = child->GetText();

		if (text != "anonymous")
			required.insert(text);
	}
	else if (type == GLOBAL)
	{
		for (int j=0; j<tree->ChildCount(); j++)
		{
			OMLTree* child = tree->GetChild(j);
			defined.insert(child->GetText());
		}
	}
	else if (type == FUNC_DEF)
	{
		// need to set up the input vars as defined
		// and then just run the helper for the statements

		OMLTree* child = tree->GetChild(0); // 0 is the token for the name
		defined.insert(child->GetText());

		std::set<std::string> local_required;
		std::set<std::string> local_defined;

		child = tree->GetChild(2); // 2 is the token for the inputs
		for (int j=0; j<child->ChildCount(); j++)
		{
			OMLTree* input = child->GetChild(j);
			local_defined.insert(input->GetText());
		}

		child = tree->GetChild(3); // 3 is the token for the statements
		AnalyzeHelper(child, local_required, local_defined);

		std::set<std::string>::iterator iter;
		for (iter = local_required.begin(); iter != local_required.end(); iter++)
		{
			std::string test = *iter;

			if (local_defined.find(test) == local_defined.end())
				required.insert(test);
		}
	}
	else if (type != DUMMY)
	{
		int child_count = tree->ChildCount();

		for (int j=0; j<child_count; j++)
		{
			OMLTree* child = tree->GetChild(j);
			AnalyzeHelper(child, required, defined);
		}
	}
}

Currency ExprTreeEvaluator::Analyze(const std::string& infile)
{
	pANTLR3_INPUT_STREAM input = ANTLRData::InputFromFilename(infile);

	if (!input)
	{
		char buffer[1024];
		sprintf(buffer, "File %s is not accessible", infile.c_str());
		throw OML_Error(buffer);
	}

	ANTLRData ad(input, true);

	pANTLR3_COMMON_TOKEN_STREAM tokens = ad.GetTokens();
	ANTLRData::PreprocessTokenStream(tokens);

	ad.CreateParser(tokens);
	pExprCppTreeParser parser = ad.GetParser();

	ExprCppTreeParser_prog_return r = parser->prog(parser);
		
	int temp = nested_function_marker;

	if (parser->pParser->rec->getNumberOfSyntaxErrors(parser->pParser->rec) != 0)
		throw OML_Error("Invalid file");

	pANTLR3_BASE_TREE tree = r.tree;
	PreprocessAST(tree, tokens);
	OMLTree* oml_tree = OMLTree::ConvertTree(tree);

	// finally have the tree
	std::set<std::string> required;
	std::set<std::string> defined;

	std::vector<std::string> needed;

	AnalyzeHelper(oml_tree, required, defined);

	std::set<std::string>::iterator iter;
	for (iter = required.begin(); iter != required.end(); iter++)
	{
		std::string test = *iter;

		if (defined.find(test) == defined.end())
		{
			FunctionInfo* fi   = NULL;
			FUNCPTR       fptr = NULL;

			FindFunctionByName(test, &fi, &fptr);

			// it's a built-in function so don't include it
			if (!fptr)
				needed.push_back(test);
		}
	}

	HML_CELLARRAY* ans_cell = ExprTreeEvaluator::allocateCellArray(1, (int)needed.size());

	for (int j=0; j<needed.size(); j++)
		(*ans_cell)(j) = needed[j];

	return ans_cell;

	return 0.0;
}