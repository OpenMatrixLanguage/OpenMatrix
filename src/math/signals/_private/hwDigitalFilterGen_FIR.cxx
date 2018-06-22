/**
* @file  hwDigitalFilterGen_FIR.cxx
* @date April 2009
* Copyright (C) 2009-2018 Altair Engineering, Inc.  
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
#include "hwMathStatus.h"
#include "hwDigitalFilterGen_FIR.h"
#include "hwDigitalFilter.h"
#include "hwFilterSpecs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwDigitalFilterGen_FIR::hwDigitalFilterGen_FIR(const hwFilterSpecs& filterSpecs,
                                               hwDigitalFilter&     digitalFilter,
                                               const hwMatrix*      window)
    : m_pWindow (nullptr)
{
    if (filterSpecs.Order() < 1)
    {
        m_status(HW_MATH_ERR_FILTERORDER, 1);
        return;
    }

    if (filterSpecs.BandType() == HighPass || filterSpecs.BandType() == BandStop)
    {
        if (filterSpecs.Order()%2 != 0)
        {
            m_status(HW_MATH_ERR_FILTERORDERODD, 1);
            return;
        }
    }

    m_status = filterSpecs.Status();
    if (!m_status.IsOk())   // checked after order for API reasons
    {
        return;
    }

    digitalFilter.SetSize(filterSpecs.Order() + 1, false);

    if (!digitalFilter.Status().IsOk())
    {
        m_status = digitalFilter.Status();
        return;
    }

    m_pDigitalFilter = &digitalFilter;
    m_pFilterSpecs = &filterSpecs;
    
    if (window)
    {
        if (!window->IsReal())
        {
            m_status(HW_MATH_ERR_COMPLEX, 3);
            return;
        }

        if (!window->IsEmptyOrVector())
        {
            m_status(HW_MATH_ERR_VECTOR, 3);
            return;
        }

        if (window->Size() != m_pFilterSpecs->Order() + 1)
        {
            m_status(HW_MATH_ERR_FIRWINDOW, 3);
            return;
        }

        m_pWindow = new hwMatrix(*window);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwDigitalFilterGen_FIR::~hwDigitalFilterGen_FIR()
{
    delete m_pWindow;
    m_pWindow = nullptr;
}
//------------------------------------------------------------------------------
// Compute transfer function polynomial coefficients for a low pass filter
//------------------------------------------------------------------------------
void hwDigitalFilterGen_FIR::CreateLowPassFilter(bool normalize)
{
    int numNumerCoefs = m_pDigitalFilter->GetNumerCoefs()->Size();
    double a = 0.5 * m_pFilterSpecs->Order();
    double omegaC = 2.0 * PI * m_pFilterSpecs->GetUpperCorner();        // normalized angular freq
    double* numerCoef = m_pDigitalFilter->GetNumerCoefs()->GetRealData();

    int i = 0;
    if (m_pWindow)
    {
        for (i = 0; i < numNumerCoefs/2; ++i)
        {
            numerCoef[i] = sin(omegaC*(i-a)) / (PI*(i-a)) * (*m_pWindow)(i);
        }

        if (m_pFilterSpecs->Order()%2 == 0)  // even order, odd number of coefficients
        {
            numerCoef[i] = omegaC / PI * (*m_pWindow)(i);
            ++i;
        }

        for ( ; i < numNumerCoefs; ++i)
        {
            numerCoef[i] = sin(omegaC*(i-a)) / (PI*(i-a)) * (*m_pWindow)(i);
        }
    }
    else
    {
        for (i = 0; i < numNumerCoefs/2; ++i)
        {
            numerCoef[i] = sin(omegaC*(i-a)) / (PI*(i-a));
            numerCoef[m_pFilterSpecs->Order()-i] = numerCoef[i];        // set symmetrical coefficient
        }

        if (m_pFilterSpecs->Order()%2 == 0)  // even order, odd number of coefficients
        {
            numerCoef[i] = omegaC / PI;
        }
    }

    // normalize
    if (normalize)
    {
        double scale;

        m_pDigitalFilter->RespMag(0.0, scale);

        for (i = 0; i < numNumerCoefs; ++i)
        {
            numerCoef[i] /= scale;
        }
    }
}
//------------------------------------------------------------------------------
// Compute transfer function polynomial coefficients for a high pass filter
//------------------------------------------------------------------------------
void hwDigitalFilterGen_FIR::CreateHighPassFilter(bool normalize)
{
    int i = 0;
    int numNumerCoefs = m_pDigitalFilter->GetNumerCoefs()->Size();
    double a = 0.5 * m_pFilterSpecs->Order();
    double omegaC = 2.0 * PI * m_pFilterSpecs->GetLowerCorner();        // normalized angular freq
    double* numerCoef = m_pDigitalFilter->GetNumerCoefs()->GetRealData();

    if (m_pWindow)
    {
        for (i = 0; i < numNumerCoefs/2; ++i)
        {
            numerCoef[i] = -sin(omegaC*(i-a)) / (PI*(i-a)) * (*m_pWindow)(i);
        }

        numerCoef[i] = (1.0 - omegaC / PI) * (*m_pWindow)(i);
        ++i;

        for ( ; i < numNumerCoefs; ++i)
        {
            numerCoef[i] = -sin(omegaC*(i-a)) / (PI*(i-a)) * (*m_pWindow)(i);
        }
    }
    else
    {
        for (i = 0; i < numNumerCoefs/2; ++i)
        {
            numerCoef[i] = -sin(omegaC*(i-a)) / (PI*(i-a));
            numerCoef[m_pFilterSpecs->Order()-i] = numerCoef[i];        // set symmetrical coefficient
        }

        numerCoef[i] = 1.0 - omegaC / PI;
    }

    // normalize
    if (normalize)
    {
        double scale;

        m_pDigitalFilter->RespMag(PI, scale);

        for (i = 0; i < numNumerCoefs; ++i)
        {
            numerCoef[i] /= scale;
        }
    }
}
//------------------------------------------------------------------------------
// Compute transfer function polynomial coefficients for a band pass filter
//------------------------------------------------------------------------------
void hwDigitalFilterGen_FIR::CreateBandPassFilter(bool normalize)
{
    int i = 0;
    int numNumerCoefs = m_pDigitalFilter->GetNumerCoefs()->Size();
    double a = 0.5 * m_pFilterSpecs->Order();
    double omegaC1 = 2.0 * PI * m_pFilterSpecs->GetLowerCorner();       // normalized angular freq
    double omegaC2 = 2.0 * PI * m_pFilterSpecs->GetUpperCorner();       // normalized angular freq
    double* numerCoef = m_pDigitalFilter->GetNumerCoefs()->GetRealData();

    if (m_pWindow)
    {
        for (i = 0; i < numNumerCoefs/2; ++i)
        {
            numerCoef[i] = (sin(omegaC2*(i-a)) - sin(omegaC1*(i-a))) / (PI*(i-a)) * (*m_pWindow)(i);
        }

        if (m_pFilterSpecs->Order()%2 == 0)                             // even order, odd number of coefficients
        {
            numerCoef[i] = (omegaC2 - omegaC1) / PI * (*m_pWindow)(i);
            ++i;
        }

        for ( ; i < numNumerCoefs; ++i)
        {
            numerCoef[i] = (sin(omegaC2*(i-a)) - sin(omegaC1*(i-a))) / (PI*(i-a)) * (*m_pWindow)(i);
        }
    }
    else
    {
        for (i = 0; i < numNumerCoefs/2; ++i)
        {
            numerCoef[i] = (sin(omegaC2*(i-a)) - sin(omegaC1*(i-a))) / (PI*(i-a));
            numerCoef[m_pFilterSpecs->Order()-i] = numerCoef[i];        // set symmetrical coefficient
        }

        if (m_pFilterSpecs->Order()%2 == 0)  // even order, odd number of coefficients
        {
            numerCoef[i] = (omegaC2 - omegaC1) / PI;
        }
    }

    // normalize
    if (normalize)
    {
        double scale;

        m_pDigitalFilter->RespMag(0.5 * (omegaC1 + omegaC2), scale);

        for (i = 0; i < numNumerCoefs; ++i)
        {
            numerCoef[i] /= scale;
        }
    }
}
//------------------------------------------------------------------------------
// Compute transfer function polynomial coefficients for a band stop filter
// -----------------------------------------------------------------------------
void hwDigitalFilterGen_FIR::CreateBandStopFilter(bool normalize, bool multiband)
{
#if 0 // Commented code
/*
    // previous code that creates the band stop filter directly without explicit
    // creation of the corresponding unnormalized bandpass filter
    int i;
    int numNumerCoefs = m_pDigitalFilter->GetNumerCoefs()->Size();
    double a = 0.5 * m_pFilterSpecs->Order();
    double omegaC1 = 2.0 * PI * m_pFilterSpecs->GetLowerCorner();       // normalized angular freq
    double omegaC2 = 2.0 * PI * m_pFilterSpecs->GetUpperCorner();       // normalized angular freq
    double* numerCoef = m_pDigitalFilter->GetNumerCoefs()->GetRealData();

    if (m_pWindow)
    {
        for (i = 0; i < numNumerCoefs/2; ++i)
            numerCoef[i] = -(sin(omegaC2*(i-a)) - sin(omegaC1*(i-a))) / (PI*(i-a)) * (*m_pWindow)(i);

        numerCoef[i] = 1.0 - ((omegaC2 - omegaC1) / PI) * (*m_pWindow)(i);
        ++i;

        for ( ; i < numNumerCoefs; ++i)
            numerCoef[i] = -(sin(omegaC2*(i-a)) - sin(omegaC1*(i-a))) / (PI*(i-a)) * (*m_pWindow)(i);
    }
    else
    {
        for (i = 0; i < numNumerCoefs/2; ++i)
        {
            numerCoef[i] = -(sin(omegaC2*(i-a)) - sin(omegaC1*(i-a))) / (PI*(i-a));
            numerCoef[m_pFilterSpecs->Order()-i] = numerCoef[i];        // set symmetrical coefficient
        }

        numerCoef[i] = 1.0 - (omegaC2 - omegaC1) / PI;
    }
*/
#endif

    // create corresponding passband filter, without normalization of the passband
    CreateBandPassFilter(normalize && multiband);

    // convert to stop band filter
    int numNumerCoefs = m_pDigitalFilter->GetNumerCoefs()->Size();
    double* numerCoef = m_pDigitalFilter->GetNumerCoefs()->GetRealData();
    int mid = numNumerCoefs/2;

    for (int i = 0; i < numNumerCoefs; ++i)
    {
        numerCoef[i] = -numerCoef[i];
    }

    if (m_pWindow)
    {
        numerCoef[mid] += (*m_pWindow)(mid);
    }
    else
    {
        numerCoef[mid] += 1.0;
    }
    // normalize stop band filter at DC
    if (normalize)
    {
        double scale;
        m_pDigitalFilter->RespMag(0.0, scale);

        for (int i = 0; i < numNumerCoefs; ++i)
        {
            numerCoef[i] /= scale;
        }
    }
}
