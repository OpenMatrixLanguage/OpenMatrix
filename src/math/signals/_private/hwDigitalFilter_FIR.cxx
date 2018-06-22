/**
* @file  hwDigitalFilter_FIR.cxx
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
#include "hwMatrix.h"
#include "hwDigitalFilter_FIR.h"
#include "hwFilterSpecs.h"
#include "hwDigitalFilterGen_FIR.h"

//------------------------------------------------------------------------------
// Constructor for single band FIR digital filters
//------------------------------------------------------------------------------
hwDigitalFilter_FIR::hwDigitalFilter_FIR(int             order, 
                                         double          lowCutoffFreq,
                                         double          highCutoffFreq, 
                                         const hwMatrix* window,
                                         bool            normalize, 
                                         bool            multibandStop)
{
	hwFilterSpecs filterSpecs(order, lowCutoffFreq, highCutoffFreq, true);

    hwDigitalFilterGen_FIR filterGen(filterSpecs, *this, window);

    m_status = filterGen.Status();
    if (!m_status.IsOk())
    {
        return;
    }

    switch (filterSpecs.BandType())
    {
        case LowPass:  filterGen.CreateLowPassFilter(normalize);  break;
        case HighPass: filterGen.CreateHighPassFilter(normalize); break;
        case BandPass: filterGen.CreateBandPassFilter(normalize); break;
        case BandStop: filterGen.CreateBandStopFilter(normalize, multibandStop); break;
        default: break;
    }
}
//------------------------------------------------------------------------------
// Constructor for multiband FIR digital filters
//------------------------------------------------------------------------------
hwDigitalFilter_FIR::hwDigitalFilter_FIR(int             order, 
                                         const hwMatrix& cutoffFreq,
                                         const hwMatrix* window, 
                                         const char*     type,
                                         bool            normalize)
{
    if (!cutoffFreq.IsReal())
    {
        m_status(HW_MATH_ERR_COMPLEX, 2);
        return;
    }

    int numBands = (cutoffFreq.Size() + 1) / 2;
    double lowCutoffFreq;
    double highCutoffFreq;

    m_pNumerCoef = new hwMatrix(order+1, hwMatrix::REAL);
    m_pNumerCoef->SetElements(0.0);

    for (int i = 0; i < cutoffFreq.Size()-1; ++i)
    {
        if (cutoffFreq(i) >= cutoffFreq(i+1))
        {
            m_status(HW_MATH_ERR_NONINCREASE, 2);
            return;
        }
    }

    if (!type || !strcmp(type, "DC-0"))
    {
        for (int i = 0; i < numBands; ++i)
        {
            if (i < numBands - 1)
            {
                lowCutoffFreq = cutoffFreq(2*i);
                highCutoffFreq = cutoffFreq(2*i+1);
            }
            else if (cutoffFreq.Size()%2 == 0)  // last band is band pass
            {
                lowCutoffFreq = cutoffFreq(2*i);
                highCutoffFreq = cutoffFreq(2*i+1);
            }
            else  // last band is high pass
            {
                lowCutoffFreq = cutoffFreq(cutoffFreq.Size()-1);
                highCutoffFreq = 1.0;
            }
#if 0 // Commented code \todo: Remove
            /*
            // This method adds the individual unnormalized bandpass filters
            // together, and then normalizes the multiband filter based on the first band. The
            // drawback is that remaining bands can have different gains than if they are created
            // as individuals.
            hwDigitalFilter_FIR band(order, lowCutoffFreq, highCutoffFreq, window, false);
            */
#endif
            // Normalized bandpass filters are added together, and then the
            // multiband is normalized based on the first band. This causes the 
            // remaining bands to have gains
            hwDigitalFilter_FIR band(order, lowCutoffFreq, highCutoffFreq, 
                                     window, normalize);

            m_status = band.Status();

            if (!m_status.IsOk())
            {
                return;
            }

            (*m_pNumerCoef) += (*band.GetNumerCoefs());
        }

        if (normalize)
        {
            double scale;
            RespMag(PI * 0.5 * (cutoffFreq(0) + cutoffFreq(1)), scale);
            (*m_pNumerCoef) /= scale;
        }
    }
    else if (!strcmp(type, "DC-1"))
    {
        int mid = order / 2;

        for (int i = 0; i < numBands; ++i)
        {
            if (i < numBands - 1)
            {
                lowCutoffFreq = cutoffFreq(2*i);
                highCutoffFreq = cutoffFreq(2*i+1);
            }
            else if (cutoffFreq.Size()%2 == 0)  // last band is band stop
            {
                lowCutoffFreq = cutoffFreq(2*i);
                highCutoffFreq = cutoffFreq(2*i+1);
            }
            else  // last band is low pass
            {
                lowCutoffFreq = cutoffFreq(cutoffFreq.Size()-1);
                highCutoffFreq = 0.0;
            }
#if 0
            /*
            // This is the Octave method. It adds the individual unnormalized band stop filters
            // together.
            hwDigitalFilter_FIR band(order, highCutoffFreq, lowCutoffFreq, window, false);
            */
#endif
            // In this method, the normalized band stop filters are added 
            // together. The normalization is done on the corresponding bandpass 
            // filter so that better rejection is acheived in the stop band.
            hwDigitalFilter_FIR band(order, highCutoffFreq, lowCutoffFreq, 
                                     window, normalize, true);

            m_status = band.Status();

            if (!m_status.IsOk())
            {
                return;
            }

            (*m_pNumerCoef) += (*band.GetNumerCoefs());
        }

        /* The multiband band stop filter is created by converting the individual filters
        to bandpass filters. The multiband bandpass filter is the converted back to a
        band stop filter. The arithmetic simplifies to summing the individual stop band
        filters and then modifying the center value. */
        (*m_pNumerCoef)(mid) -= numBands - 1;

        if (normalize)
        {
            double scale;
            RespMag(0.0, scale);
            (*m_pNumerCoef) /= scale;
        }
    }
    else
    {
        m_status(HW_MATH_ERR_FILTERTYPE, 4);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwDigitalFilter_FIR::~hwDigitalFilter_FIR()
{
}
