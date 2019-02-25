/**
* @file OmlUtils.h
* @date October 2016
* Copyright (C) 2016-2018 Altair Engineering, Inc.  
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

#ifndef __OMLUTILS__
#define __OMLUTILS__

// Begin defines/includes
#include "Hml2Dll.h"

#include <string>
#include <vector>

// End defines/includes

//------------------------------------------------------------------------------
//! Class for utility functions used by evaluator/interpreter/built in func
//------------------------------------------------------------------------------
class HML2DLL_DECLS OmlUtils
{
    public:
    //! Destructor
    ~OmlUtils() {}

    //! Gets current working directory
    static std::string GetCurrentWorkingDir();

    //! Gets file names matching pattern in current directory
    //! \param[in] pattern Pattern to match
    static std::vector<std::string> GetMatchingFiles( const std::string& pattern);

private:
    //! Constructor - cannot instantiate this class
    OmlUtils() {}

};
#endif