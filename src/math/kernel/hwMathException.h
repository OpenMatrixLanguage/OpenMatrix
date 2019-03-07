/**
* @file hwMathException.h
* @date June 2014
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

#ifndef _hwMathException_h
#define _hwMathException_h

#include <exception>
#include <hwMathStatus.h>

#pragma warning(disable : 4275)

// -------------------------------------------------------------------
//! \class hwMathException hwMathException.h math/core/hwMathException.h
//! \brief An exception class for the math core
//! 
//! The class contains an hwMathStatus object and overrides the
//! std::exception what function.
// -------------------------------------------------------------------

//! hwMathException class definition
class MATHCORE_DECLS hwMathException: public std::exception
{
public:
    //! Explicit constructor
    explicit hwMathException(hwMathMsgCode mathCode, int arg1 = -1, int arg2 = -1); 
    //! No throw destructor
    ~hwMathException() throw();
    //! Copy constructor
    hwMathException(const hwMathException& mathException);
    //! Return the exception description, overriding the std::exception version
    virtual const char* what() const throw();
    //! Return the hwMathStatus member
    hwMathStatus& Status() {return m_status;}

private:
    //! Inaccessible = operator
    void operator=(const hwMathException& mathException);

private:
    //! object that contains error coding and generates messages
    hwMathStatus m_status;
    //! string for error message
    char* m_msg;
};

#endif // _hwMathException_h
