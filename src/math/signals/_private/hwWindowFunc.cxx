/**
* @file hwWindowFunc.h
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
#include <math.h>

#include "hwWindowFunc.h"
#include "hwMathStatus.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwWindowFunc::hwWindowFunc(bool periodic)
{
    m_windowType = (!periodic) ? Symmetric : Periodic;
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwWindowFunc::hwWindowFunc(int  numPoints,
                           bool periodic)
{
    if (numPoints < 1)
    {
        return;
    }

    hwMathStatus status = m_weight.Dimension(numPoints, hwMatrix::REAL);
    if (!status.IsOk())
    {
        return;
    }
    m_windowType = (!periodic) ? Symmetric : Periodic;
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwWindowFunc::hwWindowFunc(const hwMatrix& weight)
{
    m_weight = weight;
}
//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
hwWindowFunc::hwWindowFunc(const hwWindowFunc& src)
{
    Copy(src);
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwWindowFunc::~hwWindowFunc()
{
}
//------------------------------------------------------------------------------
// Assignment operator
//------------------------------------------------------------------------------
hwWindowFunc& hwWindowFunc::operator=(const hwWindowFunc& src)
{
    Copy(src);
    return *this;
}
//------------------------------------------------------------------------------
// Returns status after setting window length
//------------------------------------------------------------------------------
hwMathStatus hwWindowFunc::SetSize(int numPoints)
{
    if (numPoints == m_weight.Size())    // nothing to update
    {
        return hwMathStatus();
    }

    hwMathStatus status = m_weight.Dimension(numPoints, hwMatrix::REAL);
    if (!status.IsOk())
    {
        return status;
    }

    ComputeWeights(m_weight);
    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Copy function
//------------------------------------------------------------------------------
hwMathStatus hwWindowFunc::Copy(const hwWindowFunc& src)
{
    hwMathStatus status = m_weight.Dimension(m_weight.M(), m_weight.N(),
                                             hwMatrix::REAL);
    if (!status.IsOk())
    {
        return status;
    }

    m_weight     = src.m_weight;
    m_windowType = src.m_windowType;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Compute bias correction constants
//------------------------------------------------------------------------------
void hwWindowFunc::Bias(const hwMatrix& signal, 
                        int             start,
                        double&         k1, 
                        double&         k2)
{
    int    size = m_weight.Size();
    double sum3 = 0.0;

    for (int j = 0; j < size; ++j)
    {
#if 0
        //sum1 += m_weight(j) * signal(start);
        //sum2 += m_weight(j);
#endif
        sum3 += m_weight(j) * m_weight(j);
        ++start;
    }
#if 0
/*
    // k1 gives the windowed data mean. It is used in Templex, but is omitted here.
    if (IsZero(sum2, 1.0e-12))
        k1 = sum1 / sum2;
    else
        k1 = 0.0;
*/    
#endif

    k2 = 1.0;
    if (!IsZero(sum3, 1.0e-12))
    {
        k2 = sqrt(size / sum3);
    }       
}
//------------------------------------------------------------------------------
// Window the signal data
//------------------------------------------------------------------------------
hwMathStatus hwWindowFunc::ApplyWindow(const hwMatrix& data_to_window, 
                                       int             start,
                                       hwMatrix&       windowed_data, 
                                       bool            bias)
{
    if (!data_to_window.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data_to_window.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    if (start < 0)
    {
        return hwMathStatus(HW_MATH_ERR_NEGATIVE, 2);
    }

    int size = m_weight.Size();
    if (data_to_window.Size() < start + size)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1);
    }

    if (bias)
    {
        double k1;
        double k2;
        Bias(data_to_window, start, k1, k2);

        for (int j = 0; j < size; ++j)
        {
            // windowed_data(j) = m_weight(j) * (data_to_window(start) - k1) * k2;
            windowed_data(j) = m_weight(j) * data_to_window(start) * k2;
            ++start;
        }
        return hwMathStatus();
    }

    for (int j = 0; j < size; ++j)
    {
        windowed_data(j) = m_weight(j) * data_to_window(start);
        ++start;
    }
    return hwMathStatus();
}
