/**
* @file OMLTree.h
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

#ifndef __OMLTree_h
#define __OMLTree_h

#include "OMLDll.h"
#include <vector>
#include "Currency.h"
#include "Evaluator.h"

class OMLTree;

class OMLDLL_DECLS OMLTree
{
public:
	OMLTree(int type, const std::string& text, const std::string* filename, int line, int num_children);
	~OMLTree();

	OMLTree(const OMLTree&);

	int        ChildCount() const { return (int)_children.size(); }
	OMLTree*   GetChild(int idx) const { return _children[idx]; }
	OMLTree*   DetachChild(int idx);
	int        GetType() const { return _type; }
	Currency   Run();

	std::string        Filename() const { return *_filename; }
	const std::string* FilenamePtr() const { return _filename; }
	int                Line() const { return _line; } 
	std::string        GetText() const { return _text; }

	void        AddChild(OMLTree*);

	void        SetDebugInfo(const char* filename, int line_number);

	static OMLTree* ConvertTree(void* ANTLR_tree);

	static OMLTree* ReadTreeFromFile(const std::string& filename);
	void            WriteTreeToFile(const std::string& filename);
	void            ReadFromBinaryFile(FILE* outfile, const std::string* filename);
	void            WriteToBinaryFile(FILE* outfile);

	void               GetListOfIdents(std::vector<std::string>&) const;
    const std::string* GetLeadingIdent();

	std::string DumpAST() const;

	std::string GetStringRepresentation() const;

	const OMLTree* FindParentOf(std::string& ident_name) const;

	bool IsBroadcastOutput() const;

	void*      u; // user data
	TREE_FPTR  func_ptr;

private:
	std::vector<OMLTree*> _children;
	std::string           _text;
	int                   _type;
	int                   _line;
	const std::string*    _filename;
};

#endif
