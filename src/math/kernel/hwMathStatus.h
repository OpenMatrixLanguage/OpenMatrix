/**
* @file hwMathStatus.h
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
#ifndef _hwMathStatus_h
#define _hwMathStatus_h

#include "MathCoreExports.h"
#include "Globals.h"

#include <iostream>
#include <vector>

//------------------------------------------------------------------------------
//!
//! \brief  Math function return status class
//!
//------------------------------------------------------------------------------
class MATHCORE_DECLS hwMathStatus
{
public:
    hwMathStatus();
    hwMathStatus(hwMathMsgCode mathCode, int arg1 = -1, int arg2 = -1);
    hwMathStatus(const hwMathStatus& status);
    ~hwMathStatus();

public:
    bool IsOk() const {return (m_pContent ? false : true);}
    void SetArg1(int arg1);
    void SetArg2(int arg2);
    int GetArg1() const;
    int GetArg2() const;
    void ResetArgs();
    void SetMsgCode(hwMathMsgCode mathCode);
    hwMathMsgCode GetMsgCode() const;
    void SetTboxFuncName(const std::string& tboxFuncName);
    void SetUserFuncName(const std::string& userFuncName);
    bool IsWarning() const { return (m_pContent && GetMsgCode() < HW_MATH_WARN_EOL) ? true : false; }
    hwMathStatus& operator()(hwMathMsgCode mathCode,
                             int arg1 = -1, int arg2 = -1);
    hwMathMsgCode operator()() const { return GetMsgCode(); }   // shorthand
    void operator=(const hwMathStatus& status);
    bool operator==(hwMathMsgCode mathCode) const;
    bool operator!=(hwMathMsgCode mathCode) const;

    //!
    //! Gets generated error/status/warning message \todo: Delete GetMessage
    //!
    std::string GetMessageString() const;
    //!
    //! Gets generated message \todo: Redundant method, use GetMessageString
    //!
    std::string GetMessage() const;

private:
    class hwContent
    {
    public:
        hwContent(hwMathMsgCode mathCode = HW_MATH_ERR_NONE,
                  int arg1 = -1, int arg2 = -1);
        hwContent(const hwContent& content);
        ~hwContent();

    private:
        int m_arg1;                 // first argument related to error
        int m_arg2;                 // second argument
        hwMathMsgCode m_mathCode;   // message code enum
        std::string m_tboxFuncName; // toolbox function in which error is reported
        std::string m_userFuncName; // user function in which error is reported

    public:
        inline void SetArg1(int arg1) { m_arg1 = arg1; }
        inline void SetArg2(int arg2) { m_arg2 = arg2; }
        inline int GetArg1() const { return m_arg1; }
        inline int GetArg2() const { return m_arg2; }
        inline void ResetArgs() { m_arg1 = -1; m_arg2 = -1; }
        inline void SetMsgCode(hwMathMsgCode mathCode) { m_mathCode = mathCode; }
        inline hwMathMsgCode GetMsgCode() const { return m_mathCode; }
        inline void SetTboxFuncName(const std::string& tboxFuncName) { m_tboxFuncName = tboxFuncName; }
        inline void SetUserFuncName(const std::string& userFuncName) { m_userFuncName = userFuncName; }
        std::string GetMessage() const;

        void operator=(const hwContent& content);
        void operator()(hwMathMsgCode mathCode, int arg1, int arg2);
    };

private:
    mutable hwContent* m_pContent;
};

#endif // _hwMathStatus_h
