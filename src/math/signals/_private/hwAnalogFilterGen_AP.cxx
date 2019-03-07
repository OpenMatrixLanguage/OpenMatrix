/**
* @file hwAnalogFilterGen_AP.cxx
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
#include "hwAnalogFilterGen_AP.h"

#include "hwMatrix.h"
#include "hwAnalogFilter.h"
#include "hwFilterSpecs.h"
#include "GeneralFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwAnalogFilterGen_AP::hwAnalogFilterGen_AP(const hwFilterSpecs& filterSpecs,
                                           hwAnalogFilter&      analogFilter)
	: hwAnalogFilterGen_IIR(filterSpecs, analogFilter)
{
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwAnalogFilterGen_AP::~hwAnalogFilterGen_AP()
{
}
//------------------------------------------------------------------------------
// Compute transfer function polynomial coefficients for a low pass filter
//------------------------------------------------------------------------------
void hwAnalogFilterGen_AP::CreateLowPassFilter()
{
    // coefficients are stored in descending powers of s
    double omegaC = m_pFilterSpecs->GetUpperCorner();
    double sPlanePoleReal;
    double* numerCoef = m_pAnalogFilter->GetNumerCoefs()->GetRealData();
    double* denomCoef = m_pAnalogFilter->GetDenomCoefs()->GetRealData();

    // apply transformation to low pass prototype
    if (m_pFilterSpecs->Order() == 1)
    {
        // for one pole compute 1/(s-p)
        m_pAnalogFilter->GetSPlaneInfo(sPlanePoleReal);
        numerCoef[0] = 0.0;
        numerCoef[1] = -omegaC * sPlanePoleReal;
        denomCoef[0] = 1.0;
        denomCoef[1] = -omegaC * sPlanePoleReal;
    }
    else
    {
        // for each complex conjugate pair of poles compute [1/(s-p)][1/(s-p~)]
        double omegaCSq = omegaC * omegaC;
        double sPlanePoleMagSq;

        // first pole pair
        m_pAnalogFilter->GetSPlaneInfo(0, sPlanePoleReal, sPlanePoleMagSq);
        numerCoef[0] = 0.0;
        numerCoef[1] = 0.0;
        numerCoef[2] = omegaCSq * sPlanePoleMagSq;
        denomCoef[0] = 1.0;
        denomCoef[1] = -2.0 * omegaC * sPlanePoleReal;
        denomCoef[2] = omegaCSq * sPlanePoleMagSq;

        // this loop compiles the transfer function by sequentially multiplying the remaining
        // factors of the numerator and the denominator for each combined complex conjugate pair
        int k = 0;
        double nfc[3];          // numerator factor coefficients for a pole pair
        double dfc[3];          // denominator factor coefficients for a pole pair

        for (int i = 1; i < m_pFilterSpecs->Order()/2; ++i)
        {
            m_pAnalogFilter->GetSPlaneInfo(i, sPlanePoleReal, sPlanePoleMagSq);
            nfc[0] = 0.0;
            nfc[1] = 0.0;
            nfc[2] = omegaCSq * sPlanePoleMagSq;
            dfc[0] = 1.0;
            dfc[1] = -2.0 * omegaC * sPlanePoleReal;
            dfc[2] = omegaCSq * sPlanePoleMagSq;

            k += 2;
            polyNomMultAccum(numerCoef, k, nfc, 2);
            polyNomMultAccum(denomCoef, k, dfc, 2);
        }

        if (m_pFilterSpecs->Order()%2 == 1)
        {
            // last pole when there is an odd number
            m_pAnalogFilter->GetSPlaneInfo(sPlanePoleReal);
            nfc[0] = 0.0;
            nfc[1] = -omegaC * sPlanePoleReal;
            dfc[0] = 1.0;
            dfc[1] = -omegaC * sPlanePoleReal;

            polyNomMultAccum(numerCoef, m_pFilterSpecs->Order()-1, nfc, 1);
            polyNomMultAccum(denomCoef, m_pFilterSpecs->Order()-1, dfc, 1);
        }
    }

    ScaleCoefs();
}
//------------------------------------------------------------------------------
// Compute transfer function polynomial coefficients for a high pass filter
//------------------------------------------------------------------------------
void hwAnalogFilterGen_AP::CreateHighPassFilter()
{
    // coefficients are stored in descending powers of s
    double omegaC = m_pFilterSpecs->GetLowerCorner();
    double sPlanePoleReal;
    double* numerCoef = m_pAnalogFilter->GetNumerCoefs()->GetRealData();
    double* denomCoef = m_pAnalogFilter->GetDenomCoefs()->GetRealData();

    // apply transformation to low pass prototype
    if (m_pFilterSpecs->Order() == 1)
    {
        // for one pole compute 1/(s-p)
        m_pAnalogFilter->GetSPlaneInfo(sPlanePoleReal);
        numerCoef[0] = -sPlanePoleReal;
        numerCoef[1] = 0.0;
        denomCoef[0] = -sPlanePoleReal;
        denomCoef[1] = omegaC;
    }
    else
    {
        // for each complex conjugate pair of poles compute [1/(s-p)][1/(s-p~)]
        double omegaCSq = omegaC * omegaC;
        double sPlanePoleMagSq;

        // first pole pair
        m_pAnalogFilter->GetSPlaneInfo(0, sPlanePoleReal, sPlanePoleMagSq);
        numerCoef[0] = sPlanePoleMagSq;
        numerCoef[1] = 0.0;
        numerCoef[2] = 0.0;
        denomCoef[0] = sPlanePoleMagSq;
        denomCoef[1] = -2.0 * omegaC * sPlanePoleReal;
        denomCoef[2] = omegaCSq;

        // this loop compiles the transfer function by sequentially multiplying the remaining
        // factors of the numerator and the denominator for each combined complex conjugate pair
        int k = 0;
        double nfc[3];          // numerator factor coefficients for a pole pair
        double dfc[3];          // denominator factor coefficients for a pole pair

        for (int i = 1; i < m_pFilterSpecs->Order()/2; ++i)
        {
            m_pAnalogFilter->GetSPlaneInfo(i, sPlanePoleReal, sPlanePoleMagSq);
            nfc[0] = sPlanePoleMagSq;
            nfc[1] = 0.0;
            nfc[2] = 0.0;
            dfc[0] = sPlanePoleMagSq;
            dfc[1] = -2.0 * omegaC * sPlanePoleReal;
            dfc[2] = omegaCSq;

            k += 2;
            polyNomMultAccum(numerCoef, k, nfc, 2);
            polyNomMultAccum(denomCoef, k, dfc, 2);
        }

        if (m_pFilterSpecs->Order()%2 == 1)
        {
            // last pole when there is an odd number
            m_pAnalogFilter->GetSPlaneInfo(sPlanePoleReal);
            nfc[0] = -sPlanePoleReal;
            nfc[1] = 0.0;
            dfc[0] = -sPlanePoleReal;
            dfc[1] = omegaC;

            polyNomMultAccum(numerCoef, m_pFilterSpecs->Order()-1, nfc, 1);
            polyNomMultAccum(denomCoef, m_pFilterSpecs->Order()-1, dfc, 1);
        }
    }

    ScaleCoefs();
}
//------------------------------------------------------------------------------
// Compute transfer function polynomial coefficients for a band pass filter
//------------------------------------------------------------------------------
void hwAnalogFilterGen_AP::CreateBandPassFilter()
{
    // coefficients are stored in descending powers of s
    double bandWidth = m_pFilterSpecs->GetUpperCorner() - m_pFilterSpecs->GetLowerCorner();
    double bandCenterSq = m_pFilterSpecs->GetUpperCorner() * m_pFilterSpecs->GetLowerCorner();
    double sPlanePoleReal;
    double* numerCoef = m_pAnalogFilter->GetNumerCoefs()->GetRealData();
    double* denomCoef = m_pAnalogFilter->GetDenomCoefs()->GetRealData();

    // apply transformation to low pass prototype
    if (m_pFilterSpecs->Order() == 1)
    {
        // for one pole compute 1/(s-p)
        m_pAnalogFilter->GetSPlaneInfo(sPlanePoleReal);
        numerCoef[0] = 0.0;
        numerCoef[1] = -sPlanePoleReal * bandWidth;
        numerCoef[2] = 0.0;
        denomCoef[0] = 1.0;
        denomCoef[1] = -sPlanePoleReal * bandWidth;
        denomCoef[2] = bandCenterSq;
    }
    else
    {
        // for each complex conjugate pair of poles compute [1/(s-p)][1/(s-p~)]
        double sPlanePoleMagSq;

        // first pole pair
        m_pAnalogFilter->GetSPlaneInfo(0, sPlanePoleReal, sPlanePoleMagSq);
        numerCoef[0] = 0.0;
        numerCoef[1] = 0.0;
        numerCoef[2] = sPlanePoleMagSq * bandWidth * bandWidth;
        numerCoef[3] = 0.0;
        numerCoef[4] = 0.0;
        denomCoef[0] = 1.0;
        denomCoef[1] = -2.0 * sPlanePoleReal * bandWidth;
        denomCoef[2] = 2.0 * bandCenterSq + sPlanePoleMagSq * bandWidth * bandWidth;
        denomCoef[3] = -2.0 * sPlanePoleReal * bandCenterSq * bandWidth;
        denomCoef[4] = bandCenterSq * bandCenterSq;

        // this loop compiles the transfer function by sequentially multiplying the remaining
        // factors of the numerator and the denominator for each combined complex conjugate pair
        int k = 0;
        double nfc[5];          // numerator factor coefficients
        double dfc[5];          // denominator factor coefficients

        for (int i = 1; i < m_pFilterSpecs->Order()/2; ++i)
        {
            m_pAnalogFilter->GetSPlaneInfo(i, sPlanePoleReal, sPlanePoleMagSq);
            nfc[0] = 0.0;
            nfc[1] = 0.0;
            nfc[2] = sPlanePoleMagSq * bandWidth * bandWidth;
            nfc[3] = 0.0;
            nfc[4] = 0.0;
            dfc[0] = 1.0;
            dfc[1] = -2.0 * sPlanePoleReal * bandWidth;
            dfc[2] = 2.0 * bandCenterSq + sPlanePoleMagSq * bandWidth * bandWidth;
            dfc[3] = -2.0 * sPlanePoleReal * bandCenterSq * bandWidth;
            dfc[4] = bandCenterSq * bandCenterSq;

            k += 4;
            polyNomMultAccum(numerCoef, k, nfc, 4);
            polyNomMultAccum(denomCoef, k, dfc, 4);
        }

        if (m_pFilterSpecs->Order()%2 == 1)
        {
            // last pole when there is an odd number
            m_pAnalogFilter->GetSPlaneInfo(sPlanePoleReal);
            nfc[0] = 0.0;
            nfc[1] = -sPlanePoleReal * bandWidth;
            nfc[2] = 0.0;
            dfc[0] = 1.0;
            dfc[1] = -sPlanePoleReal * bandWidth;
            dfc[2] = bandCenterSq;

            polyNomMultAccum(numerCoef, 2*(m_pFilterSpecs->Order()-1), nfc, 2);
            polyNomMultAccum(denomCoef, 2*(m_pFilterSpecs->Order()-1), dfc, 2);
        }
    }

    ScaleCoefs();
}
//------------------------------------------------------------------------------
// Compute transfer function polynomial coefficients for a band stop filter
//------------------------------------------------------------------------------
void hwAnalogFilterGen_AP::CreateBandStopFilter()
{
    // coefficients are stored in descending powers of s
    double bandWidth = m_pFilterSpecs->GetUpperCorner() - m_pFilterSpecs->GetLowerCorner();
    double bandCenterSq = m_pFilterSpecs->GetUpperCorner() * m_pFilterSpecs->GetLowerCorner();
    double sPlanePoleReal;
    double* numerCoef = m_pAnalogFilter->GetNumerCoefs()->GetRealData();
    double* denomCoef = m_pAnalogFilter->GetDenomCoefs()->GetRealData();

    // apply transformation to low pass prototype
    if (m_pFilterSpecs->Order() == 1)
    {
        // for one pole compute 1/(s-p)
        m_pAnalogFilter->GetSPlaneInfo(sPlanePoleReal);
        numerCoef[0] = -sPlanePoleReal;
        numerCoef[1] = 0.0;
        numerCoef[2] = -sPlanePoleReal * bandCenterSq;
        denomCoef[0] = -sPlanePoleReal;
        denomCoef[1] = bandWidth;
        denomCoef[2] = -sPlanePoleReal * bandCenterSq;
    }
    else
    {
        // for each complex conjugate pair of poles compute [1/(s-p)][1/(s-p~)]
        double sPlanePoleMagSq;

        // first pole pair
        m_pAnalogFilter->GetSPlaneInfo(0, sPlanePoleReal, sPlanePoleMagSq);
        numerCoef[0] = sPlanePoleMagSq;
        numerCoef[1] = 0.0;
        numerCoef[2] = 2.0 * sPlanePoleMagSq * bandCenterSq;
        numerCoef[3] = 0.0;
        numerCoef[4] = sPlanePoleMagSq * bandCenterSq * bandCenterSq;
        denomCoef[0] = sPlanePoleMagSq;
        denomCoef[1] = -2.0 * sPlanePoleReal * bandWidth;
        denomCoef[2] = 2.0 * sPlanePoleMagSq * bandCenterSq + bandWidth * bandWidth;
        denomCoef[3] = -2.0 * sPlanePoleReal * bandCenterSq * bandWidth;
        denomCoef[4] = sPlanePoleMagSq * bandCenterSq * bandCenterSq;

        // this loop compiles the transfer function by sequentially multiplying the remaining
        // factors of the numerator and the denominator for each combined complex conjugate pair
        int k = 0;
        double nfc[5];          // numerator factor coefficients
        double dfc[5];          // denominator factor coefficients

        for (int i = 1; i < m_pFilterSpecs->Order()/2; ++i)
        {
            m_pAnalogFilter->GetSPlaneInfo(i, sPlanePoleReal, sPlanePoleMagSq);
            nfc[0] = sPlanePoleMagSq;
            nfc[1] = 0.0;
            nfc[2] = 2.0 * sPlanePoleMagSq * bandCenterSq;
            nfc[3] = 0.0;
            nfc[4] = sPlanePoleMagSq * bandCenterSq * bandCenterSq;
            dfc[0] = sPlanePoleMagSq;
            dfc[1] = -2.0 * sPlanePoleReal * bandWidth;
            dfc[2] = 2.0 * sPlanePoleMagSq * bandCenterSq + bandWidth * bandWidth;
            dfc[3] = -2.0 * sPlanePoleReal * bandCenterSq * bandWidth;
            dfc[4] = sPlanePoleMagSq * bandCenterSq * bandCenterSq;

            k += 4;
            polyNomMultAccum(numerCoef, k, nfc, 4);
            polyNomMultAccum(denomCoef, k, dfc, 4);
        }

        if (m_pFilterSpecs->Order()%2 == 1)
        {
            // last pole when there is an odd number
            m_pAnalogFilter->GetSPlaneInfo(sPlanePoleReal);
            nfc[0] = -sPlanePoleReal;
            nfc[1] = 0.0;
            nfc[2] = -sPlanePoleReal * bandCenterSq;
            dfc[0] = -sPlanePoleReal;
            dfc[1] = bandWidth;
            dfc[2] = -sPlanePoleReal * bandCenterSq;

            polyNomMultAccum(numerCoef, 2*(m_pFilterSpecs->Order()-1), nfc, 2);
            polyNomMultAccum(denomCoef, 2*(m_pFilterSpecs->Order()-1), dfc, 2);
        }
    }

    ScaleCoefs();
}
