/**
* @file hwChebyshev.cxx
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
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

#include <cmath>

#include "hwChebyshev.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwChebyshev::hwChebyshev(double sideLobe,
                         bool   periodic)
    : hwWindowFunc(periodic)
{
    if (sideLobe < 0.0)
        return;

    m_gamma = exp(-sideLobe * log(10.0) / 20.0);
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwChebyshev::hwChebyshev(int    num_points, 
                         double sideLobe,
                         bool   periodic)
    : hwWindowFunc(num_points, periodic)
{
    if (sideLobe < 0.0)
        return;

    m_gamma = exp(-sideLobe * log(10.0) / 20.0);

    ComputeWeights(m_weight);
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwChebyshev::~hwChebyshev()
{
}
//------------------------------------------------------------------------------
// Computes window weights
//------------------------------------------------------------------------------
hwMathStatus hwChebyshev::ComputeWeights(hwMatrix& weight)
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
    if (size == 1)
    {
        weight(0) = 1;
        return hwMathStatus();
    }

    int m;
    double beta;
    double max = 0.0;
    double value1, value2;

    if (m_windowType == Periodic)
    {
        m = (size + 1) / 2;
        beta = cosh(acosh(1.0 / m_gamma) / size);

        for (int j = 0; j < size; ++j)
        {
            value2 = 0.0;

            for (int k = 1; k <= m; ++k)
            {
                value1 = beta * cos(k * PI / (size + 1));

                value1 = (fabs(value1) > 1) ? 
                         cosh(size * acosh(value1)) :
                         cos(size * acos(value1));

                value2 += value1 * cos((2*j - size) * k * PI / (size + 1));
            }

            weight(j) = (1.0 / m_gamma + 2.0 * value2) / (size + 1);
            max = _max(max, weight(j));
        }
    }
    else // (m_windowType == Symmetric)
    {
        m = size / 2;
        beta = cosh(acosh(1.0 / m_gamma) / (size - 1.0));

        for (int j = 0; j < size; ++j)
        {
            value2 = 0.0;

            for (int k = 1; k <= m; ++k)
            {
                value1 = beta * cos(k * PI / size);

                value1 = (fabs(value1) > 1) ?
                         cosh((size - 1.0) * acosh(value1)) :
                         cos((size - 1.0) * acos(value1));

                value2 += value1 * cos((2*j - size + 1) * k * PI / size);
            }

            weight(j) = (1.0 / m_gamma + 2.0 * value2) / size;
            max = _max(max, weight(j));
        }
    }

    for (int j = 0; j < size; ++j)
    {
        weight(j) /= max;
    }
    return hwMathStatus();
}
