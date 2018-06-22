/**
* @file hwMathStatus.cxx
* @date January 2009
* Copyright (C) 2009-2018 Altair Engineering, Inc.  
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

#include "hwMathStatus.h"

#if defined(_DARWIN) || defined(LINUX)
  #include <stdio.h>
  #define _itoa_s(value, str, size, base) ::sprintf((str), "%d", value)
#endif

//*******************************************************************
//                      hwMathStatus functions
//*******************************************************************

//! Empty constructor
hwMathStatus::hwMathStatus()
    : m_pContent(NULL)
{
}

//! Constructor containing information
hwMathStatus::hwMathStatus(hwMathMsgCode mathCode,
                           int arg1, int arg2)
{
    m_pContent = new hwContent(mathCode, arg1, arg2);
}

//! Copy constructor
hwMathStatus::hwMathStatus(const hwMathStatus& status)
    : m_pContent(status.m_pContent)
{
    if (this == &status)
        return;

    if (m_pContent)
        status.m_pContent = NULL;
}

//! Destructor
hwMathStatus::~hwMathStatus()
{
    if (m_pContent)
        delete m_pContent;
}

//! Set index of first argument involved in a failed function call
void hwMathStatus::SetArg1(int arg1)
{
    if (m_pContent)
        m_pContent->SetArg1(arg1);
}

//! Set index of second argument involved in a failed function call
void hwMathStatus::SetArg2(int arg2)
{
    if (m_pContent)
        m_pContent->SetArg2(arg2);
}

//! Get index of first argument involved in a failed function call
int hwMathStatus::GetArg1() const
{
    if (m_pContent)
        return m_pContent->GetArg1();
    else
        return -1;
}

//! Get index of second argument involved in a failed function call
int hwMathStatus::GetArg2() const
{
    if (m_pContent)
        return m_pContent->GetArg2();
    else
        return -1;
}

//! Reset member variables to default values
void hwMathStatus::ResetArgs()
{
    if (m_pContent)
        return m_pContent->ResetArgs();
}

//! Set the error message code for a failed function call
void hwMathStatus::SetMsgCode(hwMathMsgCode mathCode)
{
    if (m_pContent)
        m_pContent->SetMsgCode(mathCode);
    else
        m_pContent = new hwContent(mathCode);
}

//! Get the error message code for a failed function call
hwMathMsgCode hwMathStatus::GetMsgCode() const
{
    if (m_pContent)
        return m_pContent->GetMsgCode();
    else
        return HW_MATH_ERR_NONE;
}

//! Set the toolbox function name involved in a failed function call
void hwMathStatus::SetTboxFuncName(const std::string& tboxFuncName)
{
    if (!m_pContent)
        m_pContent = new hwContent;

    m_pContent->SetTboxFuncName(tboxFuncName);
}

//! Set the user function name involved in a failed function call
void hwMathStatus::SetUserFuncName(const std::string& userFuncName)
{
    if (!m_pContent)
        m_pContent = new hwContent;

    m_pContent->SetUserFuncName(userFuncName);
}

//! Set all member variables to record a failed function call
hwMathStatus& hwMathStatus::operator()(hwMathMsgCode mathCode,
                                       int arg1, int arg2)
{
    if (mathCode == HW_MATH_ERR_NONE)
    {
        if (m_pContent)
        {
            delete m_pContent;
            m_pContent = NULL;
        }
    }
    else
    {
        if (m_pContent)
            (*m_pContent)(mathCode, arg1, arg2);
        else
            m_pContent = new hwContent(mathCode, arg1, arg2);
    }

    return *this;
}

//! Assignment operator
void hwMathStatus::operator=(const hwMathStatus& status)
{
    if (this == &status)
        return;

    if (m_pContent)
        delete m_pContent;

    // content is transferred, not copied. multiple copies should not be needed
    m_pContent = status.m_pContent;

    if (m_pContent)
        status.m_pContent = NULL;
}

//! Equality operator to compare error message codes
bool hwMathStatus::operator==(hwMathMsgCode mathCode) const
{
    if (GetMsgCode() == mathCode)
        return true;
    else
        return false;
}

//! Inequality operator to compare error message codes
bool hwMathStatus::operator!=(hwMathMsgCode mathCode) const
{
    if (GetMsgCode() != mathCode)
        return true;
    else
        return false;
}

//*******************************************************************
//              private class hwContent implementation
//*******************************************************************

//! Constructor containing information
hwMathStatus::hwContent::hwContent(hwMathMsgCode mathCode,
                                   int arg1, int arg2)
    : m_mathCode(mathCode), m_arg1(arg1), m_arg2(arg2)
{
}

//! Copy constructor
hwMathStatus::hwContent::hwContent(const hwContent& content)
{
    m_mathCode = content.m_mathCode;
    m_arg1 = content.m_arg1;
    m_arg2 = content.m_arg2;
}

//! Destructor
hwMathStatus::hwContent::~hwContent()
{
}

//! Convert error code information into a string
std::string hwMathStatus::hwContent::GetMessage() const
{
    std::string message = GetHMathErrMsg(m_mathCode);

    // if (m_tboxFuncName.empty())  // not used in oml
    //     return message;

    // insert function name in message
    size_t pos = message.find(';');
    size_t len = message.length();
    std::string statusMsg = message.substr(0, pos);

    if (m_arg1 != -1)
    {
        char argChar[8];

        if (m_arg1 == -999)     // internal hopt error
        {
            // insert m_arg2
            _itoa_s(m_arg2, argChar, 8, 10);
            statusMsg += "code " + std::string(argChar);
        }
        else
        {
            // prepare to insert m_arg1
            _itoa_s(m_arg1, argChar, 8, 10);

            if (m_arg2 != -1)
            {
                // insert m_arg1
                statusMsg += " in arguments " + std::string(argChar);

                // insert m_arg2
                _itoa_s(m_arg2, argChar, 8, 10);
                statusMsg += "," + std::string(argChar);
            }
            else
            {
                // insert m_arg1
                statusMsg = message.substr(0, pos) + " in argument " + std::string(argChar);
            }
        }
    }

    // complete the string
    if (pos != std::string::npos)
        statusMsg += message.substr(pos, len-pos);

    if (!m_userFuncName.empty())
    {
        size_t pos = statusMsg.find(';');
        statusMsg.insert(pos + 2, m_userFuncName + " ");
    }
/*
    {
        pos = statusMsg.find("user function");
        statusMsg.erase(pos, 13);
        statusMsg.insert(pos, m_userFuncName);
    }
*/
    return statusMsg;
}

//! Assignment operator
void hwMathStatus::hwContent::operator=(const hwContent& content)
{
    m_mathCode = content.m_mathCode;
    m_arg1 = content.m_arg1;
    m_arg2 = content.m_arg2;
}

//! Set all member variables to record a failed function call
void hwMathStatus::hwContent::operator()(hwMathMsgCode mathCode,
                                         int arg1, int arg2)
{
    m_mathCode = mathCode;
    m_arg1 = arg1;
    m_arg2 = arg2;
}
//------------------------------------------------------------------------------
// Returns message \todo: Delete GetMessage
//------------------------------------------------------------------------------
std::string hwMathStatus::GetMessageString() const
{
    if (m_pContent)
        return m_pContent->GetMessage();

    return HW_MATH_MSG_SUCCESS;
}
//------------------------------------------------------------------------------
//! Gets generated message \todo: Redundant method, use GetMessageString
//------------------------------------------------------------------------------
std::string hwMathStatus::GetMessage() const
{
    return GetMessageString();
}
