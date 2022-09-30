/**
* @file hwSliceArg.cxx
* @date July 2014
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

//:---------------------------------------------------------------------------
//:Description
//
//  Matrix slice argument class
//  A complete slice definition is stored as std::vector<hwSliceArg>
//
//:---------------------------------------------------------------------------
#include <hwSliceArg.h>
#include <hwMathException.h>
#include <GeneralFuncs.h>

//! Constructor for colon argument
hwSliceArg::hwSliceArg(): m_type(TYPE_COLON)
{
}

//! Constructor for integer argument
hwSliceArg::hwSliceArg(int index): m_type(TYPE_SCALAR), m_index(index)
{
    if (index < 0)
        throw hwMathException(HW_MATH_ERR_SLICE_INDEX);
}

//! Constructor for vector argument
hwSliceArg::hwSliceArg(const std::vector<int>& indexVec): m_type(TYPE_VECTOR), m_indexVec(indexVec)
{
    for (int i = 0; i < indexVec.size(); ++i)
    {
        if (indexVec[i] < 0)
            throw hwMathException(HW_MATH_ERR_SLICE_INDEX);
    }
}

// utility function used when deleting slices from hwTMatrixN objects
// placed here for convenience, but needing of a better home
//! Create a vector with no repeated elements
#include <algorithm>

std::vector<int> GetUniqueVec(const std::vector<int>& input)
{
    std::vector<int> result(input);

    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());

    return result;
}

