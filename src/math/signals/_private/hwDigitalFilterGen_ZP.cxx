/**
* @file  hwDigitalFilterGen_ZP.cxx
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
#include "hwDigitalFilterGen_ZP.h"

#include "hwMatrix.h"
#include "hwDigitalFilter.h"
#include "hwFilterSpecs.h"
#include "GeneralFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwDigitalFilterGen_ZP::hwDigitalFilterGen_ZP(const hwFilterSpecs& filterSpecs,
                                             hwDigitalFilter&     digitalFilter)
	: hwDigitalFilterGen_IIR(filterSpecs, digitalFilter)
{
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwDigitalFilterGen_ZP::~hwDigitalFilterGen_ZP()
{
}
//------------------------------------------------------------------------------
// Compute transfer function polynomial coefficients for a low pass filter
//------------------------------------------------------------------------------
void hwDigitalFilterGen_ZP::CreateLowPassFilter()
{
    // coefficients are stored in descending powers of z
    double omegaC = tan(PI * m_pFilterSpecs->GetUpperCorner());  // prewarped angular cutoff freq
    double sPlanePoleReal;
    double* numerCoef = m_pDigitalFilter->GetNumerCoefs()->GetRealData();
    double* denomCoef = m_pDigitalFilter->GetDenomCoefs()->GetRealData();

    // apply bilinear transform s = (z'-1)/(z'+1) to low pass filter
    // see Sanjit K Mitra, "Digital Signal Processing"
    if (m_pFilterSpecs->Order() == 1)
    {
        // for one pole compute 1/(s-p)
        m_pDigitalFilter->GetSPlaneInfo(omegaC, sPlanePoleReal);
        numerCoef[0] = -sPlanePoleReal;
        numerCoef[1] = -sPlanePoleReal;
        denomCoef[0] = 1.0 - sPlanePoleReal;
        denomCoef[1] = -(1.0 + sPlanePoleReal);
    }
    else
    {
        // for each complex conjugate pair of poles compute [1/(s-p)][1/(s-p~)]
        double omegaCSq = omegaC * omegaC;

        double sPlaneZeroMagSq;
        double sPlanePoleMagSq;

        // first pole/zero pair
        m_pDigitalFilter->GetSPlaneInfo(0, omegaC, omegaCSq, sPlanePoleReal, 
                                        sPlanePoleMagSq, sPlaneZeroMagSq);
        
        double gainFactor = sPlanePoleMagSq / sPlaneZeroMagSq;
        numerCoef[0] = (sPlaneZeroMagSq + 1.0) * gainFactor;
        numerCoef[1] = 2.0 * (sPlaneZeroMagSq - 1.0) * gainFactor;
        numerCoef[2] = numerCoef[0];
        denomCoef[0] = sPlanePoleMagSq - 2.0 * sPlanePoleReal + 1.0;
        denomCoef[1] = 2.0 * (sPlanePoleMagSq - 1.0);
        denomCoef[2] = sPlanePoleMagSq + 2.0 * sPlanePoleReal + 1.0;

        // this loop compiles the transfer function by sequentially multiplying the remaining
        // factors of the numerator and the denominator for each combined complex conjugate pair
        int k = 0;
        double nfc[3];          // numerator factor coefficients for a pole/zero pair
        double dfc[3];          // denominator factor coefficients for a pole/zero pair

        for (int i = 1; i < m_pFilterSpecs->Order()/2; ++i)
        {
            m_pDigitalFilter->GetSPlaneInfo(i, omegaC, omegaCSq, sPlanePoleReal, sPlanePoleMagSq, sPlaneZeroMagSq);
            gainFactor = sPlanePoleMagSq / sPlaneZeroMagSq;
            nfc[0] = (sPlaneZeroMagSq + 1.0) * gainFactor;
            nfc[1] = 2.0 * (sPlaneZeroMagSq - 1.0) * gainFactor;
            nfc[2] = nfc[0];
            dfc[0] = sPlanePoleMagSq - 2.0 * sPlanePoleReal + 1.0;
            dfc[1] = 2.0 * (sPlanePoleMagSq - 1.0);
            dfc[2] = sPlanePoleMagSq + 2.0 * sPlanePoleReal + 1.0;

            k += 2;
            polyNomMultAccum(numerCoef, k, nfc, 2);
            polyNomMultAccum(denomCoef, k, dfc, 2);
        }

        if (m_pFilterSpecs->Order()%2 == 1)
        {
            // last pole when there is an odd number
            m_pDigitalFilter->GetSPlaneInfo(omegaC, sPlanePoleReal);
            nfc[0] = -sPlanePoleReal;
            nfc[1] = -sPlanePoleReal;
            dfc[0] = 1.0 - sPlanePoleReal;
            dfc[1] = -(1.0 + sPlanePoleReal);

            polyNomMultAccum(numerCoef, m_pFilterSpecs->Order()-1, nfc, 1);
            polyNomMultAccum(denomCoef, m_pFilterSpecs->Order()-1, dfc, 1);
        }
    }

    ScaleCoefs();
}
//------------------------------------------------------------------------------
// Compute transfer function polynomial coefficients for a high pass filter
//------------------------------------------------------------------------------
void hwDigitalFilterGen_ZP::CreateHighPassFilter()
{
    // coefficients are stored in descending powers of z
    double omegaC = tan(PI * m_pFilterSpecs->GetLowerCorner());        // prewarped angular cutoff freq
    double alpha  = -cos(2.0 * PI * m_pFilterSpecs->GetLowerCorner());
    double beta   = (1.0 + alpha) / (1.0 - alpha);
    double sPlanePoleReal;
    double* numerCoef = m_pDigitalFilter->GetNumerCoefs()->GetRealData();
    double* denomCoef = m_pDigitalFilter->GetDenomCoefs()->GetRealData();

    // apply bilinear transform s = (z'-1)/(z'+1) to low pass filter
    // set z' = -(z+alpha)/(alpha*z+1)
    // see Sanjit K Mitra, "Digital Signal Processing"
    if (m_pFilterSpecs->Order() == 1)
    {
        // for one pole compute 1/(s-p)
        m_pDigitalFilter->GetSPlaneInfo(omegaC, sPlanePoleReal);
        numerCoef[0] = -sPlanePoleReal;
        numerCoef[1] =  sPlanePoleReal;
        denomCoef[0] = beta - sPlanePoleReal;
        denomCoef[1] = beta + sPlanePoleReal;
    }
    else
    {
        // for each complex conjugate pair of poles compute [1/(s-p)][1/(s-p~)]
        double omegaCSq = omegaC * omegaC;
        double sPlanePoleMagSq;
        double sPlaneZeroMagSq;
        double betaSq = beta * beta;

        // first pole/zero pair
        m_pDigitalFilter->GetSPlaneInfo(0, omegaC, omegaCSq, sPlanePoleReal, sPlanePoleMagSq, sPlaneZeroMagSq);
        double gainFactor = sPlanePoleMagSq / sPlaneZeroMagSq;
        numerCoef[0] = (betaSq + sPlaneZeroMagSq) * gainFactor;
        numerCoef[1] = 2.0 * (betaSq - sPlaneZeroMagSq) * gainFactor;
        numerCoef[2] = numerCoef[0];
        denomCoef[0] = sPlanePoleMagSq - 2.0 * beta * sPlanePoleReal + betaSq;
        denomCoef[1] = 2.0 * (betaSq - sPlanePoleMagSq);
        denomCoef[2] = sPlanePoleMagSq + 2.0 * beta * sPlanePoleReal + betaSq;

        // this loop compiles the transfer function by sequentially multiplying the remaining
        // factors of the numerator and the denominator for each combined complex conjugate pair
        int k = 0;
        double nfc[3];          // numerator factor coefficients for a pole/zero pair
        double dfc[3];          // denominator factor coefficients for a pole/zero pair

        for (int i = 1; i < m_pFilterSpecs->Order()/2; ++i)
        {
            m_pDigitalFilter->GetSPlaneInfo(i, omegaC, omegaCSq, sPlanePoleReal, sPlanePoleMagSq, sPlaneZeroMagSq);
            gainFactor = sPlanePoleMagSq / sPlaneZeroMagSq;
            nfc[0] = (betaSq + sPlaneZeroMagSq) * gainFactor;
            nfc[1] = 2.0 * (betaSq - sPlaneZeroMagSq) * gainFactor;
            nfc[2] = nfc[0];
            dfc[0] = sPlanePoleMagSq - 2.0 * beta * sPlanePoleReal + betaSq;
            dfc[1] = 2.0 * (betaSq - sPlanePoleMagSq);
            dfc[2] = sPlanePoleMagSq + 2.0 * beta * sPlanePoleReal + betaSq;

            k += 2;
            polyNomMultAccum(numerCoef, k, nfc, 2);
            polyNomMultAccum(denomCoef, k, dfc, 2);
        }

        if (m_pFilterSpecs->Order()%2 == 1)
        {
            // last pole when there is an odd number
            m_pDigitalFilter->GetSPlaneInfo(omegaC, sPlanePoleReal);
            nfc[0] = -sPlanePoleReal;
            nfc[1] = sPlanePoleReal;
            dfc[0] = beta - sPlanePoleReal;
            dfc[1] = beta + sPlanePoleReal;

            polyNomMultAccum(numerCoef, m_pFilterSpecs->Order()-1, nfc, 1);
            polyNomMultAccum(denomCoef, m_pFilterSpecs->Order()-1, dfc, 1);
        }
    }

    ScaleCoefs();
}
//------------------------------------------------------------------------------
// Compute transfer function polynomial coefficients for a band pass filter
//------------------------------------------------------------------------------
void hwDigitalFilterGen_ZP::CreateBandPassFilter()
{
    // coefficients are stored in descending powers of z
    double bandWidth = m_pFilterSpecs->GetUpperCorner() - 
                       m_pFilterSpecs->GetLowerCorner();
    double bandCenter = 0.5 * (m_pFilterSpecs->GetLowerCorner() + 
                               m_pFilterSpecs->GetUpperCorner());

    double lowPassFC = bandWidth;            // cutoff of low pass filter to be converted
    double omegaC    = tan(PI * lowPassFC);  // prewarped angular cutoff freq
    double alpha     = cos(2.0 * PI * bandCenter) / cos(PI * bandWidth);
    double sPlanePoleReal;
    double* numerCoef = m_pDigitalFilter->GetNumerCoefs()->GetRealData();
    double* denomCoef = m_pDigitalFilter->GetDenomCoefs()->GetRealData();

    // apply bilinear transform s = (z'-1)/(z'+1) to low pass filter
    // set z' = z(z-alpha)/(alpha*z-1)
    // see Sanjit K Mitra, "Digital Signal Processing"
    if (m_pFilterSpecs->Order() == 1)
    {
        // for one pole compute 1/(s-p)
        m_pDigitalFilter->GetSPlaneInfo(omegaC, sPlanePoleReal);
        numerCoef[0] = -sPlanePoleReal;
        numerCoef[1] =  0.0;
        numerCoef[2] = sPlanePoleReal;
        denomCoef[0] = 1.0 - sPlanePoleReal;
        denomCoef[1] = -2.0 * alpha;
        denomCoef[2] = 1.0 + sPlanePoleReal;
    }
    else
    {
        // for each complex conjugate pair of poles compute [1/(s-p)][1/(s-p~)]
        double omegaCSq = omegaC * omegaC;
        double sPlanePoleMagSq;
        double sPlaneZeroMagSq;

        // first pole/zero pair
        m_pDigitalFilter->GetSPlaneInfo(0, omegaC, omegaCSq, sPlanePoleReal, sPlanePoleMagSq, sPlaneZeroMagSq);
        double gainFactor = sPlanePoleMagSq / sPlaneZeroMagSq;
        numerCoef[0] = (1.0 + sPlaneZeroMagSq) * gainFactor;
        numerCoef[1] = -4.0 * alpha * gainFactor;
        numerCoef[2] = 2.0 * (1.0 + 2.0 * alpha * alpha - sPlaneZeroMagSq) * gainFactor;
        numerCoef[3] = numerCoef[1];
        numerCoef[4] = numerCoef[0];
        denomCoef[0] = sPlanePoleMagSq - 2.0 * sPlanePoleReal + 1;
        denomCoef[1] = -4.0 * alpha * (1.0 - sPlanePoleReal);
        denomCoef[2] = -2.0 * (sPlanePoleMagSq - 2.0 * alpha * alpha - 1.0);
        denomCoef[3] = -4.0 * alpha * (1.0 + sPlanePoleReal);
        denomCoef[4] = sPlanePoleMagSq + 2.0 * sPlanePoleReal + 1;

        // this loop compiles the transfer function by sequentially multiplying the remaining
        // factors of the numerator and the denominator for each combined complex conjugate pair
        int k = 0;
        double nfc[5];          // numerator factor coefficients for a pole/zero pair
        double dfc[5];          // denominator factor coefficients for a pole/zero pair

        for (int i = 1; i < m_pFilterSpecs->Order()/2; ++i)
        {
            m_pDigitalFilter->GetSPlaneInfo(i, omegaC, omegaCSq, sPlanePoleReal, sPlanePoleMagSq, sPlaneZeroMagSq);
            gainFactor = sPlanePoleMagSq / sPlaneZeroMagSq;
            nfc[0] = (1.0 + sPlaneZeroMagSq) * gainFactor;
            nfc[1] = -4.0 * alpha * gainFactor;
            nfc[2] = 2.0 * (1.0 + 2.0 * alpha * alpha - sPlaneZeroMagSq) * gainFactor;
            nfc[3] = nfc[1];
            nfc[4] = nfc[0];
            dfc[0] = sPlanePoleMagSq - 2.0 * sPlanePoleReal + 1;
            dfc[1] = -4.0 * alpha * (1.0 - sPlanePoleReal);
            dfc[2] = -2.0 * (sPlanePoleMagSq - 2.0 * alpha * alpha - 1.0);
            dfc[3] = -4.0 * alpha * (1.0 + sPlanePoleReal);
            dfc[4] = sPlanePoleMagSq + 2.0 * sPlanePoleReal + 1;

            k += 4;
            polyNomMultAccum(numerCoef, k, nfc, 4);
            polyNomMultAccum(denomCoef, k, dfc, 4);
        }

        if (m_pFilterSpecs->Order()%2 == 1)
        {
            // last pole when there is an odd number
            m_pDigitalFilter->GetSPlaneInfo(omegaC, sPlanePoleReal);
            nfc[0] = -sPlanePoleReal;
            nfc[1] =  0.0;
            nfc[2] = sPlanePoleReal;
            dfc[0] = 1.0 - sPlanePoleReal;
            dfc[1] = -2.0 * alpha;
            dfc[2] = 1.0 + sPlanePoleReal;

            polyNomMultAccum(numerCoef, 2*(m_pFilterSpecs->Order()-1), nfc, 2);
            polyNomMultAccum(denomCoef, 2*(m_pFilterSpecs->Order()-1), dfc, 2);
        }
    }

    ScaleCoefs();
}
//------------------------------------------------------------------------------
// Compute transfer function polynomial coefficients for a band stop filter
//------------------------------------------------------------------------------
void hwDigitalFilterGen_ZP::CreateBandStopFilter()
{
    // coefficients are stored in descending powers of z
    double bandWidth = m_pFilterSpecs->GetUpperCorner() - m_pFilterSpecs->GetLowerCorner();
    double bandCenter = 0.5 * (m_pFilterSpecs->GetLowerCorner() + m_pFilterSpecs->GetUpperCorner());
    double lowPassFC = 0.5 - bandWidth;             // cutoff of low pass filter to be converted
    double omegaC = tan(PI * lowPassFC);            // prewarped angular cutoff freq
    double alpha = cos(2.0 * PI * bandCenter) / cos(PI * bandWidth);
    double sPlanePoleReal;
    double* numerCoef = m_pDigitalFilter->GetNumerCoefs()->GetRealData();
    double* denomCoef = m_pDigitalFilter->GetDenomCoefs()->GetRealData();

    // apply bilinear transform s = (z'-1)/(z'+1) to low pass filter
    // set z' = -z(z-alpha)/(alpha*z-1)
    // see Sanjit K Mitra, "Digital Signal Processing"
    if (m_pFilterSpecs->Order() == 1)
    {
        // for one pole compute 1/(s-p)
        m_pDigitalFilter->GetSPlaneInfo(omegaC, sPlanePoleReal);
        numerCoef[0] = -sPlanePoleReal;
        numerCoef[1] = 2.0 * alpha * sPlanePoleReal;
        numerCoef[2] = -sPlanePoleReal;
        denomCoef[0] = 1.0 - sPlanePoleReal;
        denomCoef[1] = 2.0 * alpha * sPlanePoleReal;
        denomCoef[2] = -(1.0 + sPlanePoleReal);
    }
    else
    {
        // for each complex conjugate pair of poles compute [1/(s-p)][1/(s-p~)]
        double omegaCSq = omegaC * omegaC;
        double sPlanePoleMagSq;
        double sPlaneZeroMagSq;

        // first pole/zero pair
        m_pDigitalFilter->GetSPlaneInfo(0, omegaC, omegaCSq, sPlanePoleReal, sPlanePoleMagSq, sPlaneZeroMagSq);
        double gainFactor = sPlanePoleMagSq / sPlaneZeroMagSq;
        numerCoef[0] = (1.0 + sPlaneZeroMagSq) * gainFactor;
        numerCoef[1] = -4.0 * alpha * sPlaneZeroMagSq * gainFactor;
        numerCoef[2] = 2.0 * (sPlaneZeroMagSq * (1.0 + 2.0 * alpha * alpha) - 1.0) * gainFactor;
        numerCoef[3] = numerCoef[1];
        numerCoef[4] = numerCoef[0];
        denomCoef[0] = sPlanePoleMagSq - 2.0 * sPlanePoleReal + 1.0;
        denomCoef[1] = -4.0 * alpha * (sPlanePoleMagSq - sPlanePoleReal);
        denomCoef[2] = 2.0 * ((1.0 + 2.0 * alpha * alpha) * sPlanePoleMagSq - 1.0);
        denomCoef[3] = -4.0 * alpha * (sPlanePoleMagSq + sPlanePoleReal);
        denomCoef[4] = sPlanePoleMagSq + 2.0 * sPlanePoleReal + 1.0;

        // this loop compiles the transfer function by sequentially multiplying the remaining
        // factors of the numerator and the denominator for each combined complex conjugate pair
        int k = 0;
        double nfc[5];          // numerator factor coefficients for a pole/zero pair
        double dfc[5];          // denominator factor coefficients for a pole/zero pair

        for (int i = 1; i < m_pFilterSpecs->Order()/2; ++i)
        {
            m_pDigitalFilter->GetSPlaneInfo(i, omegaC, omegaCSq, sPlanePoleReal, sPlanePoleMagSq, sPlaneZeroMagSq);
            gainFactor = sPlanePoleMagSq / sPlaneZeroMagSq;
            nfc[0] = (1.0 + sPlaneZeroMagSq) * gainFactor;
            nfc[1] = -4.0 * alpha * sPlaneZeroMagSq * gainFactor;
            nfc[2] = 2.0 * (sPlaneZeroMagSq * (1.0 + 2.0 * alpha * alpha) - 1.0) * gainFactor;
            nfc[3] = nfc[1];
            nfc[4] = nfc[0];
            dfc[0] = sPlanePoleMagSq - 2.0 * sPlanePoleReal + 1.0;
            dfc[1] = -4.0 * alpha * (sPlanePoleMagSq - sPlanePoleReal);
            dfc[2] = 2.0 * ((1.0 + 2.0 * alpha * alpha) * sPlanePoleMagSq - 1.0);
            dfc[3] = -4.0 * alpha * (sPlanePoleMagSq + sPlanePoleReal);
            dfc[4] = sPlanePoleMagSq + 2.0 * sPlanePoleReal + 1.0;

            k += 4;
            polyNomMultAccum(numerCoef, k, nfc, 4);
            polyNomMultAccum(denomCoef, k, dfc, 4);
        }

        if (m_pFilterSpecs->Order()%2 == 1)
        {
            // last pole when there is an odd number
            m_pDigitalFilter->GetSPlaneInfo(omegaC, sPlanePoleReal);
            nfc[0] = -sPlanePoleReal;
            nfc[1] = 2.0 * alpha * sPlanePoleReal;
            nfc[2] = -sPlanePoleReal;
            dfc[0] = 1.0 - sPlanePoleReal;
            dfc[1] = 2.0 * alpha * sPlanePoleReal;
            dfc[2] = -(1.0 + sPlanePoleReal);

            polyNomMultAccum(numerCoef, 2*(m_pFilterSpecs->Order()-1), nfc, 2);
            polyNomMultAccum(denomCoef, 2*(m_pFilterSpecs->Order()-1), dfc, 2);
        }
    }

    ScaleCoefs();
}
