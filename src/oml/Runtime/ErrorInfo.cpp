/**
* @file ErrorInfo.cpp
* @date March 2014
* Copyright (C) 2014-2018 Altair Engineering, Inc.  
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

#include "ErrorInfo.h"
#include <iostream>
#include <sstream>
#include <ostream>
#include <string>

ErrorInfo::ErrorInfo()
{
}

ErrorInfo::ErrorInfo(std::string file_name, int line_num, int ch_pos)
{
	_fileName = file_name;
	_lineNum = line_num;
	_chPos = ch_pos;
}

ErrorInfo::ErrorInfo(const ErrorInfo& other)
{ 
	_fileName = other._fileName; 
	_lineNum = other._lineNum; 
	_chPos = other._chPos;
}

ErrorInfo& ErrorInfo::operator= (const ErrorInfo& other)
{
	_fileName = other._fileName; 
	_lineNum = other._lineNum; 
	_chPos = other._chPos;
	return *this;
}

std::string ErrorInfo::GetErrStr() const
{ 			
	//gotta convert some stuff from int to std::string

	int lineNumber = _lineNum;
	std::ostringstream convertLine;
	convertLine << lineNumber;
	std::string lineNumberStr = convertLine.str();
	int chPosNumber = _chPos;
	std::ostringstream convertChPos;
	convertChPos << chPosNumber;
	std::string chPosNumberStr = convertChPos.str();
	std::string errStr = "Syntax error at line number " + lineNumberStr;
	
	if (_fileName != "dummy")
	{
		errStr += " in file ";
		errStr += _fileName;
	}

	errStr += " near character position ";
	errStr += chPosNumberStr;

	return errStr;
}

int ErrorInfo::GetLine() const
{
	return _lineNum;
}

std::string ErrorInfo::GetFilename() const
{
	return _fileName;
}