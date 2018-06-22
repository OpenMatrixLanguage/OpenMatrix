/**
* @file hwMathException.cxx
* @date June 2014
* Copyright (C) 2014-2018 Altair Engineering, Inc.  
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

#include <string.h>

// strcpy wrapper
#ifndef OS_WIN
  #ifndef strcpy_s
    static void* strcpy_s_wrapper(char* dest, size_t destsz, const char* src)
    {
        return strcpy(dest, src);
    }
    #define strcpy_s strcpy_s_wrapper
  #endif
#endif

#include <hwMathException.h>

//! Explicit constructor
hwMathException::hwMathException(hwMathMsgCode mathCode, int arg1, int arg2)
    : exception(), m_status(mathCode, arg1, arg2)
{
    m_msg = new char[80];
}

//! No throw destructor
hwMathException::~hwMathException() throw()
{
    if (m_msg)
        delete [] m_msg;
}

//! Copy constructor
hwMathException::hwMathException(const hwMathException& mathException)
    : m_status(mathException.m_status)
{
    m_msg = new char[80];
}

//! Return the exception description, overriding the std::exception version
const char* hwMathException::what() const throw()
{
    strcpy_s(m_msg, 80, m_status.GetMessage().c_str());

    return m_msg;
}
