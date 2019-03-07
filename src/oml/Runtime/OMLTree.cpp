/**
* @file OMLTree.cpp
* @date October 2017
* Copyright (C) 2017-2018 Altair Engineering, Inc.  
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
#include "OMLTree.h"
#include "OML_Error.h"
#include "ANTLRData.h"
#include <climits>

inline pANTLR3_COMMON_TOKEN getToken(pANTLR3_BASE_TREE tree)
{
	return (pANTLR3_COMMON_TOKEN)(((pANTLR3_COMMON_TREE)(tree->super))->token);
}

inline pANTLR3_BASE_TREE getChild(pANTLR3_BASE_TREE tree, unsigned i)
{
	return (pANTLR3_BASE_TREE)tree->children->elements[i].element;
}

inline const char* getText(pANTLR3_BASE_TREE tree)
{
	pANTLR3_COMMON_TOKEN tok = getToken(tree);

	pANTLR3_STRING str = tok->getText(tok);

	if (!tok->tokText.text)
	{
		tok->textState = ANTLR3_TEXT_STRING;
		tok->tokText.text = str;
	}

	return (const char*) str->chars;
} 

OMLTree::OMLTree(int type, const std::string& text, const std::string* filename, int line, int num_children) :  _text(text), _type(type), _line(line), _filename(filename)
{
	func_ptr = ExprTreeEvaluator::GetFuncPointerFromType(type);
	u        = NULL;

	_children.reserve(num_children);
}

OMLTree::~OMLTree()
{
	if ((_type == NUMBER) || (_type == HML_STRING) || (_type == HEXVAL))
	{
		Currency* temp = (Currency*)u;
		delete temp;
	}

	for (int j=0; j<ChildCount(); j++)
		delete GetChild(j);
}

OMLTree::OMLTree(const OMLTree& tree)
{
	_text     = tree._text;
	_line     = tree._line;
	_filename = tree._filename;
	_type     = tree._type;

	func_ptr = ExprTreeEvaluator::GetFuncPointerFromType(_type);
	u        = NULL;

	for (int j=0; j<tree._children.size(); j++)
		_children.push_back(new OMLTree(*tree._children[j]));
}

void OMLTree::AddChild(OMLTree* tree)
{
	_children.push_back(tree);
}

void OMLTree::SetDebugInfo(const char* filename, int line_number)
{
	_filename = Currency::pm.GetStringPointer(filename);
	_line     = line_number;
}

OMLTree* OMLTree::ConvertTree(void* antlr_tree)
{
	pANTLR3_BASE_TREE tree = (pANTLR3_BASE_TREE)antlr_tree;

	if (!tree)
		return NULL;

	pANTLR3_COMMON_TOKEN tok = tree->getToken(tree);

	std::string tree_text;
	int         tok_type = tok->getType(tok);

	// uncommenting this line will make the conversion faster, but it will break the ast function 
	// that Activate uses for code generation.  This does not store the text for tokens that are not variable
	// (i.e. +), but we need to put that back somehow when generating the AST
    //if ((tok_type == IDENT) || (tok_type == NUMBER) || (tok_type == HEXVAL) || (tok_type == QUOTE))
		tree_text = getText(tree);

	std::string filename;

	if (tok->input)
		filename =(char *)tok->input->fileName->chars;

	const std::string* file_ptr = Currency::pm.GetStringPointer(filename);

	int num_children = tree->getChildCount(tree);

	OMLTree* loc_tree = new OMLTree(tok_type, tree_text, file_ptr, tok->line, num_children);

	if (num_children != 0)
	{
		for (int j=0; j<num_children; j++)
		{
			pANTLR3_BASE_TREE child_j = (pANTLR3_BASE_TREE)getChild(tree, j);
			pANTLR3_COMMON_TOKEN child_j_tok = child_j->getToken(child_j);
			loc_tree->AddChild(ConvertTree(child_j));
		}
	}

	return loc_tree;
}

std::string OMLTree::DumpAST()
{
	int child_count = ChildCount();
	std::string result;

	int tok = GetType();

	if (tok == DUMMY)
		return "";
	
	if (child_count)
		result += "(";

	if (tok != QUOTE)
	{
		result += _text;
	}
	else
	{
		std::string temp = _text;
			
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

	if (child_count)
		result += " ";

	for (int j=0; j<child_count; j++)
	{
		result += GetChild(j)->DumpAST();

		if (j != (child_count-1))
			result += " ";
	}

	if (child_count)
		result += ")";

	return result;
}

OMLTree* OMLTree::DetachChild(int index)
{
	OMLTree* temp = _children[index];
	_children[index] = NULL;
	return temp;
}

OMLTree* OMLTree::ReadTreeFromFile(const std::string& filename)
{
	FILE* file = fopen(filename.c_str(), "rb");
	OMLTree* tree = NULL;

	if (file)
	{
		short version;
		size_t res = fread(&version, sizeof(short), 1, file);

		if (version != 3)
			throw OML_Error("Invalid omlp file");

		const std::string* file_ptr = Currency::pm.GetStringPointer(filename);

		tree = new OMLTree(0, "", file_ptr, 0, 0);
		tree->ReadFromBinaryFile(file, file_ptr);
	}

	return tree;
}

void OMLTree::ReadFromBinaryFile(FILE* outfile, const std::string* filename)
{
	// read the type
	short dummy;
	size_t res = fread(&dummy, sizeof(short), 1, outfile);
	_type = dummy;

	func_ptr = ExprTreeEvaluator::GetFuncPointerFromType(_type);

	// read the text
	short num_chars;
	res = fread(&num_chars, sizeof(short), 1, outfile);

	if (num_chars)
	{
		char* my_string = new char[num_chars+1];
		res = fread(my_string, sizeof(char), num_chars, outfile);

		my_string[num_chars] = '\0';
		_text = my_string;

		delete my_string;
	}

	// read the line
	res = fread(&dummy, sizeof(short), 1, outfile);
	_line = dummy;

	// read the children
	int num_children = 0;

	res = fread(&num_children, sizeof(int), 1, outfile);
	_children.reserve(num_children);

	for (int j=0; j<num_children; j++)
	{
		OMLTree* child = new OMLTree(0, "", filename, 0, 0);
		child->ReadFromBinaryFile(outfile, filename);
		_children.push_back(child);
	}
}

void OMLTree::WriteTreeToFile(const std::string& filename)
{
	FILE* file = fopen(filename.c_str(), "wb");

	short version = 3;
	size_t res = fwrite(&version, sizeof(short), 1, file);

	if (file)
		WriteToBinaryFile(file);

	fclose(file);
}

void OMLTree::WriteToBinaryFile(FILE* outfile)
{
	// write the type
	size_t res = fwrite(&_type, sizeof(short), 1, outfile);

	// write the string
	int num_chars = 0;
	if ((_type == IDENT) || (_type == NUMBER) || (_type == HEXVAL) || (_type == QUOTE))	
		num_chars = (short)_text.length();

	res = fwrite(&num_chars, sizeof(short), 1, outfile);

	if (num_chars)
		res = fwrite(_text.c_str(), sizeof(char), num_chars, outfile);

	// write the line
	res = fwrite(&_line, sizeof(short), 1, outfile);

	// write the children
	int num_children = (int)_children.size();

	res = fwrite(&num_children, sizeof(int), 1, outfile);

	for (int j=0; j < num_children; j++)
		_children[j]->WriteToBinaryFile(outfile);
}

void OMLTree::GetListOfIdents(std::vector<std::string>& idents) const
{
	if (_type == IDENT)
		idents.push_back(_text);

	for (int j=0; j<_children.size(); j++)
		GetChild(j)->GetListOfIdents(idents);
}	

// may want to change this to return a string pointer and use the u member (assuming it's set)
const std::string* OMLTree::GetLeadingIdent()
{
	if (_type == IDENT)
	{
		if (!u)
			u = (void*)Currency::vm.GetStringPointer(_text);
		return (const std::string*)u;
	}

	for (int j=0; j<_children.size(); j++)
	{
		OMLTree* child = GetChild(j);

		const std::string* ans = child->GetLeadingIdent();

		if (ans != NULL)
			return ans;
	}

	return NULL;
}