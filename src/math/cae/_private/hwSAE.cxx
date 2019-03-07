/**
* @file hwSAE.cxx
* @date October 1994
* Copyright (C) 1994-2018 Altair Engineering, Inc.  
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
#include "hwSAE.h"

#include <math.h>

#include "hwMatrix.h"
#include "FourierFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwSAEFilter::hwSAEFilter(double SAE_class, 
                         double sampFreq, 
                         int    fftSize)
    : m_SAE_class(SAE_class)
    , m_sampFreq(sampFreq)
    , m_fftSize(fftSize)
{
    if (m_SAE_class <= 0.0)
    {
        m_status(HW_MATH_ERR_FILTERCLASS, 1);
        return; 
    }
    if (m_sampFreq <= 0.0)
    {
        m_status(HW_MATH_ERR_NONPOSITIVE, 2);
        return;
    }
    if (fftSize < 0)
    {
        m_status(HW_MATH_ERR_NONPOSINT, 3);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwSAEFilter::~hwSAEFilter()
{
}
//------------------------------------------------------------------------------
// Filters the signal with SAE Filter
//------------------------------------------------------------------------------
hwMathStatus hwSAEFilter::Compute(const hwMatrix& input,
                                  hwMatrix&       output)
{
    double f_low  = 0.0;
    double f_high = 1.65*m_SAE_class;

    if ((f_low < 0.0) || (f_high > (m_sampFreq/2)))
    {
        return m_status(HW_MATH_ERR_FREQCLASS, 3, 4);
    }

    // reworked to use a one sided spectrum
    hwMatrix temp;
    m_status = Fft(input, temp, m_fftSize);
    if (!m_status.IsOk())
    {
        return m_status;
    }

    hwComplex* data = temp.GetComplexData();
    int size = temp.Size();
    int start = static_cast<int>(ceil(size * f_high / m_sampFreq));
    double factor;
    double delta_f = m_sampFreq / size;
    double freq    = start * delta_f;  // one sided frequency
    double log10b2 = 1.0 / log(2.0);

    for (int j = start; j < (size+1)/2; ++j)
    {
        factor        = pow(f_high/freq, log10b2);
        data[j]      *= factor;
        data[size-j] *= factor;
        freq         += delta_f;
    }

    if (size && size%2 == 0)
    {
        factor = pow(f_high/freq, log10b2);
        data[size/2] *= factor;
    }

    m_status = Ifft(temp, output, m_fftSize);

    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
    }
    return m_status;
}
