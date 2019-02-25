/**
* @file  hwFilterSpecs.cxx
* @date April 2009
* Copyright (C) 2009-2018 Altair Engineering, Inc.  
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
#include "hwFilterSpecs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwFilterSpecs::hwFilterSpecs(int    order, 
                             double lowCutoffFreq,
                             double highCutoffFreq, 
                             bool   digital)
    : m_order(order)
{
    // determine filter type
    double nearZero = 1.0e-6;

    if (lowCutoffFreq < 0.0)
    {
        if (digital)
        {
            m_status(HW_MATH_ERR_FILTERFREQ_D, 2);
        }
        else
        {
            m_status(HW_MATH_ERR_FILTERFREQ_A, 2);
        }
        m_filterType = Unknown;
        return;
    }

    if (digital)
    {
        if (lowCutoffFreq > 1.0)
        {
            m_status(HW_MATH_ERR_FILTERFREQ_D, 2);
            m_filterType = Unknown;
            return;
        }
    }

    if (highCutoffFreq < 0.0)
    {
        if (digital)
        {
            m_status(HW_MATH_ERR_FILTERFREQ_D, 3);
        }
        else
        {
            m_status(HW_MATH_ERR_FILTERFREQ_A, 3);
        }
        m_filterType = Unknown;
        return;
    }

    if (digital)
    {
        if (highCutoffFreq > 1.0)
        {
            m_status(HW_MATH_ERR_FILTERFREQ_D, 3);
            m_filterType = Unknown;
            return;
        }
    }

    if (lowCutoffFreq < highCutoffFreq)
    {
        if (highCutoffFreq < nearZero)
        {
            m_status(HW_MATH_ERR_FILTERSPEC_E, 3);
            m_filterType = Unknown;
            return;
        }

        if (lowCutoffFreq == 0.0)
        {
            if (digital && highCutoffFreq == 1.0)
            {
                m_status(HW_MATH_ERR_FILTERFREQ_D, 2, 3);
                m_filterType = Unknown;
                return;
            }

            m_filterType = LowPass;
            m_lowerCorner = 0.0;            // not used
            m_upperCorner = highCutoffFreq;
        }
        else if (digital && highCutoffFreq == 1.0)
        {
            m_filterType = HighPass;
            m_lowerCorner = lowCutoffFreq;
            m_upperCorner = 0.0;            // not used
        }
        else
        {
            m_filterType = BandPass;
            m_lowerCorner = lowCutoffFreq;
            m_upperCorner = highCutoffFreq;
        }
    }
    else if (highCutoffFreq < lowCutoffFreq)
    {
        if (lowCutoffFreq < nearZero)
        {
            m_status(HW_MATH_ERR_FILTERSPEC_E, 2);
            m_filterType = Unknown;
            return;
        }

        if (highCutoffFreq == 0.0)
        {
            if (digital && lowCutoffFreq == 1.0)
            {
                m_status(HW_MATH_ERR_FILTERFREQ_D, 2, 3);
                m_filterType = Unknown;
                return;
            }

            m_filterType = HighPass;
            m_lowerCorner = lowCutoffFreq;
            m_upperCorner = 0.0;            // not used
        }
        else if (digital && lowCutoffFreq == 1.0)
        {
            m_filterType = LowPass;
            m_lowerCorner = 0.0;            // not used
            m_upperCorner = highCutoffFreq;
        }
        else
        {
            m_filterType = BandStop;
            m_lowerCorner = highCutoffFreq;
            m_upperCorner = lowCutoffFreq;
        }
    }
    else
    {
        m_status(HW_MATH_ERR_FILTERFREQS_EQ, 2, 3);
        m_filterType = Unknown;
    }

    if (digital)
    {
        // The input frequencies are assumed to be normalized so that
        // F_normalized = F_cutoff / (F_sampling/2). However the algorithm
        // uses the F_normalized = F_cutoff / F_sampling convention. Therefore
        // the inputs must be multiplied by a factor of 1/2.
        m_lowerCorner *= 0.5;
        m_upperCorner *= 0.5;
    }
    else
    {
        // convert from hertz to radians/sec
        m_lowerCorner *= 2.0 * PI;
        m_upperCorner *= 2.0 * PI;
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwFilterSpecs::~hwFilterSpecs()
{
}
