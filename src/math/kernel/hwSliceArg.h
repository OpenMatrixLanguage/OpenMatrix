/**
* @file hwSliceArg.h
* @date July 2014
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

//:---------------------------------------------------------------------------
//:Description
//
//  Matrix slice argument class
//  A complete slice definition is stored as std::vector<hwSliceArg>
//
//:---------------------------------------------------------------------------
#ifndef __hwSliceArg_h
#define __hwSliceArg_h

#include <vector>
#include <MathCoreExports.h>

class MATHCORE_DECLS hwSliceArg
{
public:
    //! Constructor for colon argument
    hwSliceArg();
    //! Constructor for integer argument
    hwSliceArg(int index);
    //! Constructor for vector argument
    hwSliceArg(const std::vector<int>& indexVec);

    //! Query if colon type
    bool IsColon() const { return m_type == TYPE_COLON; }
    //! Query if integer type
    bool IsScalar() const { return m_type == TYPE_SCALAR; }
    //! Query if vector type
    bool IsVector() const { return m_type == TYPE_VECTOR; }

    //! Return integer value
    int Scalar() const { return m_index; }
    //! Return reference to integer value
    int& Scalar() { return m_index; }
    //! Return vector reference
    const std::vector<int>& Vector() const { return m_indexVec; }

private:
    enum { TYPE_COLON, TYPE_SCALAR, TYPE_VECTOR };
    int m_type;
    int m_index;
    std::vector<int> m_indexVec;
};

// utility function used when deleting slices from hwTMatrixN objects
// placed here for convenience, but needing of a better home
//! Create a vector with no repeated elements
MATHCORE_DECLS std::vector<int> GetUniqueVec(const std::vector<int>& input);

#endif // __hwSliceArg_h
