/**
* @file ANTLRData.h
* @date April 2015
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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

#ifndef __ANTLRdata_h
#define __ANTLRdata_h

#include "ExprCppTreeLexer.h"
#include "ExprCppTreeParser.h"
#include "Hml2Dll.h"
#include <string>

class HML2DLL_DECLS ANTLRData
{
public:
	ANTLRData(pANTLR3_INPUT_STREAM input, bool delay_parser = false);
	~ANTLRData();

	void CreateParser(pANTLR3_COMMON_TOKEN_STREAM);

	pExprCppTreeParser          GetParser()  { return parser; }
	pANTLR3_COMMON_TOKEN_STREAM GetTokens()  { return tokens; }

	void DumpTokenStream(bool include_hidden = false);

	static pANTLR3_INPUT_STREAM InputFromFilename(const std::string&);
	static pANTLR3_INPUT_STREAM InputFromExpression(const std::string& expr, const std::string& use_filename);
	static void PreprocessTokenStream(pANTLR3_COMMON_TOKEN_STREAM& tokens);

private:
	pANTLR3_INPUT_STREAM input;
	pExprCppTreeLexer lex;
	pANTLR3_COMMON_TOKEN_STREAM tokens;
	pExprCppTreeParser parser;
};

#endif
