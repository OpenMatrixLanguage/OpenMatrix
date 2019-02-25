/**
* @file hwPSD.cxx
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
#include "hwMatrix.h"

#include "hwPSD.h"
#include "hwCSD.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwPSD::hwPSD(double sampFreq, 
             int    fftSize)
    : m_sampFreq(sampFreq)
    , m_fft(fftSize)
{
    if (sampFreq <= 0.0)
    {
        m_status(HW_MATH_ERR_NONPOSITIVE, 1);
        return;
    }

    m_status = m_fft.Status();
    if (!m_status.IsOk())
    {
        m_status.SetArg1(2);
    }
}
//------------------------------------------------------------------------------
// Returns status after computing Power Spectral Density for const input
//------------------------------------------------------------------------------
hwMathStatus hwPSD::Compute(const hwMatrix& input,
                            hwMatrix&       psd)
{
    hwMatrix     fft_output;
    hwMathStatus status = m_fft.Compute(input, fft_output);
    if (!status.IsOk())
    {
        return status;
    }

    status = psd.AbsSq(fft_output);
    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    int fftSize  = m_fft.Size();
    int dataSize = input.Size();

    double norm;
    if (dataSize < fftSize)
        norm = 1.0 / ((double) dataSize * m_sampFreq);
    else
        norm = 1.0 / ((double) fftSize * m_sampFreq);

    psd *= norm;

    return status;
}
//------------------------------------------------------------------------------
// Returns status after computing Power Spectral Density for non-const input
//------------------------------------------------------------------------------
hwMathStatus hwPSD::Compute(hwMatrix& input,
                            hwMatrix& psd)
{
    hwMatrix     fft_output;
    hwMathStatus status = m_fft.Compute(input, fft_output);
    if (!status.IsOk())
    {
        return status;
    }

    status = psd.AbsSq(fft_output);
     if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    int fftSize  = m_fft.Size();
    int dataSize = input.Size();

    double norm;
    if (dataSize < fftSize)
    {
        norm = 1.0 / ((double) dataSize * m_sampFreq);
    }
    else
    {
        norm = 1.0 / ((double) fftSize * m_sampFreq);
    }

    psd *= norm;
    return status;
}
