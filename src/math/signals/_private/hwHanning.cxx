/**
* @file  hwHanning.cxx
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
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
#include "hwHanning.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwHanning::hwHanning(bool periodic)
    : hwWindowFunc(periodic)
{
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwHanning::hwHanning(int  numPoints,
                     bool periodic)
    : hwWindowFunc(numPoints, periodic)
{
    ComputeWeights(m_weight);
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwHanning::~hwHanning()
{
}
//------------------------------------------------------------------------------
// Returns status after computing window weights
//------------------------------------------------------------------------------
hwMathStatus hwHanning::ComputeWeights(hwMatrix& weight)
{
    if (!weight.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!weight.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    int size = weight.Size();

    if (m_windowType == Periodic)
    {
        for (int j = 0; j < size; ++j)
        {
            weight(j) = 0.5 * (1.0 - cos((2*PI*j) / size));
        }
    }
    else // (m_windowType == Symmetric)
    {
        if (size == 1)
        {
            weight(0) = 1;
        }
        else
        {
            for (int j = 0; j < size; ++j)
            {
                weight(j) = 0.5 * (1.0 - cos((2*PI*j) / (size - 1.0)));
            }
        }
    }

    return hwMathStatus();
}
