/**
* @file hwPSD.h
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
#ifndef _Signals_PSD_h
#define _Signals_PSD_h

#include "hwFFT.h"

//------------------------------------------------------------------------------
//!
//! \class hwPSD
//! \brief Power spectral density function classes
//!
//------------------------------------------------------------------------------
class hwPSD
{
public:
    //!
    //! Constructor
    //! \param sampFreq Sampling frequency
    //! \param fftSize
    //!
    hwPSD(double sampFreq, 
          int    fftSize = 0);
    //!
    //! Destructor
    //!
    ~hwPSD() {}

    //!
    //! Returns status after computing Power Spectral Density for const input
    //! \param input    Const input
    //! \param dataSize Utilized input length
    //! \param psd      Power spectral density of input
    //!
    hwMathStatus Compute(const hwMatrix& input,
                         int             dataSize,
                         hwMatrix&       psd);
    //!
    //! Returns status after computing Power Spectral Density for non-const input
    //! \param input    Non-const input
    //! \param dataSize Utilized input length
    //! \param psd      Power spectral density of input
    //!
    hwMathStatus Compute(hwMatrix& input,
                         int       dataSize,
                         hwMatrix& psd);
    //!
    //! Returns the status
    //!
    const hwMathStatus& Status() const { return m_status; }

private:
    double       m_sampFreq;    //!< Sampling frequency
    hwFFT_f      m_fft;         //!< FFT object
    hwMathStatus m_status;      //!< Status
};

#endif // _Signals_PSD_h
