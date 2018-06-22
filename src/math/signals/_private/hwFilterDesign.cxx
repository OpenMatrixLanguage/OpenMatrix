/**
* @file  hwFilterDesign.cxx
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
#include "hwFilterDesign.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwFilterDesign::hwFilterDesign()
{
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwFilterDesign::~hwFilterDesign()
{
}
//------------------------------------------------------------------------------
// Design an analog filter (low or high pass)
//------------------------------------------------------------------------------
hwMathStatus hwFilterDesign::Analog(double    passBandFreq, 
                                    double    stopBandFreq,
                                    double    passEdgeDb,
                                    double    stopEdgeDb,
                                    int&      order, 
                                    hwMatrix& freqC)
{
    // omegaC contains the permissible range for the corner freqency

    // check filter definition
    double nearZero = 1.0e-5;
    Type   filterType;

    if (passBandFreq < nearZero)
    {
        if (passBandFreq < 0.0)
        {
            return m_status(HW_MATH_ERR_FILTERFREQ_A, 1);
        }
        return m_status(HW_MATH_ERR_FILTERSPEC_E, 1);
    }

    if (stopBandFreq < nearZero)
    {
        if (stopBandFreq < 0.0)
        {
            return m_status(HW_MATH_ERR_FILTERFREQ_A, 2);
        }
        return m_status(HW_MATH_ERR_FILTERSPEC_E, 2);
    }

    if (passBandFreq < stopBandFreq)
    {
        filterType = LowPass;
    }
    else if (stopBandFreq < passBandFreq)
    {
        filterType = HighPass;
    }
    else
    {
        return m_status(HW_MATH_ERR_FILTERFREQS_EQ, 1, 2);
    }

    if (passEdgeDb < nearZero)
    {
        return m_status(HW_MATH_ERR_DB_SIGN, 3);
    }

    if (stopEdgeDb < nearZero)
    {
        return m_status(HW_MATH_ERR_DB_SIGN, 4);
    }
    if (stopEdgeDb < passEdgeDb + nearZero)
    {
        return m_status(HW_MATH_ERR_FILTERRIPPLE, 3, 4);
    }

    // size outputs
    m_status = freqC.Dimension(2, hwMatrix::REAL);

    if (!m_status.IsOk())
    {
        if (m_status.GetArg1() == 0)
        {
            m_status.SetArg1(6);
        }
        return m_status;
    }

    // compute lowpass prototype order and corner frequencies
    double omegaP = passBandFreq;
    double omegaS = stopBandFreq;

    double epsilon = sqrt(exp(passEdgeDb * log(10.0) / 10.0) - 1.0);
    double delta = sqrt(exp(stopEdgeDb * log(10.0) / 10.0) - 1.0);
    double omegaC1p;            // lowpass prototype normalized corner angular freq
    double omegaC2p;            // lowpass prototype normalized corner angular freq
    double omegaSP;             // lowpass prototype omegaS/omegaP ratio

    if (filterType == LowPass)
    {
        omegaSP = omegaS / omegaP;
    }
    else // (filterType == HighPass)
    {
        omegaSP = omegaP / omegaS;
    }

    DesignPrototype(epsilon, delta, omegaSP, order, omegaC1p, omegaC2p);

    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
        return m_status;
    }

	// transform range of prototype corner frequency
    if (filterType == LowPass)
    {
        freqC(0) = omegaP * omegaC1p;
        freqC(1) = omegaP * omegaC2p;
    }
    else // (filterType == HighPass)
    {
        freqC(0) = omegaP / omegaC1p;
        freqC(1) = omegaP / omegaC2p;
    }

    return m_status;
}
//------------------------------------------------------------------------------
// Design an analog filter (bandpass or bandstop)
//!
// -------------------------------------------------------------------
hwMathStatus hwFilterDesign::Analog(const hwMatrix& passBandFreq, 
                                    const hwMatrix& stopBandFreq,
                                    double          passEdgeDb, 
                                    double          stopEdgeDb,
                                    int&            order, 
                                    hwMatrix&       freqC1, 
                                    hwMatrix&       freqC2,
                                    bool            cheby2)
{
    // omegaC1 contains the permissible range for the first corner freqency
    // omegaC2 contains the permissible range for the second corner freqency

    // check filter definition
    double nearZero = 1.0e-5;
    Type   filterType;

    if (!passBandFreq.IsReal())
    {
        return m_status(HW_MATH_ERR_COMPLEX, 1);
    }

    if (passBandFreq.Size() != 2)
    {
        return m_status(HW_MATH_ERR_VECTOR2, 1);
    }

    if (passBandFreq(0) < nearZero)
    {
        if (passBandFreq(0) < 0.0)
        {
            return m_status(HW_MATH_ERR_FILTERFREQ_A, 1);
        }
        return m_status(HW_MATH_ERR_FILTERSPEC_E, 1);
    }

    if (passBandFreq(1) < nearZero)
    {
        if (passBandFreq(1) < 0.0)
        {
            return m_status(HW_MATH_ERR_FILTERFREQ_A, 1);
        }
        return m_status(HW_MATH_ERR_FILTERSPEC_E, 1);
    }

    if (passBandFreq(0) >= passBandFreq(1))
    {
        return m_status(HW_MATH_ERR_NONINCREASE, 1);
    }

    if (!stopBandFreq.IsReal())
    {
        return m_status(HW_MATH_ERR_COMPLEX, 2);
    }
    if (stopBandFreq.Size() != 2)
    {
        return m_status(HW_MATH_ERR_VECTOR2, 2);
    }
    if (stopBandFreq(0) < nearZero)
    {
        if (stopBandFreq(0) < 0.0)
        {
            return m_status(HW_MATH_ERR_FILTERFREQ_A, 2);
        }
        return m_status(HW_MATH_ERR_FILTERSPEC_E, 2);
    }

    if (stopBandFreq(0) >= stopBandFreq(1))
    {
        return m_status(HW_MATH_ERR_NONINCREASE, 2);
    }

    if (stopBandFreq(0) < passBandFreq(0) && passBandFreq(1) < stopBandFreq(1))
    {
        filterType = BandPass;
    }
    else if (passBandFreq(0) < stopBandFreq(0) && stopBandFreq(1) < passBandFreq(1))
    {
        filterType = BandStop;
    }
    else
    {
        return m_status(HW_MATH_ERR_FILTERBANDCONF, 1, 2);
    }

    if (passEdgeDb < nearZero)
    {
        return m_status(HW_MATH_ERR_DB_SIGN, 3);
    }

    if (stopEdgeDb < nearZero)
    {
        return m_status(HW_MATH_ERR_DB_SIGN, 4);
    }
    if (stopEdgeDb < passEdgeDb + nearZero)
    {
        return m_status(HW_MATH_ERR_FILTERTYPE, 3, 4);
    }

    // size outputs
    m_status = freqC1.Dimension(passBandFreq.M(), passBandFreq.N(), hwMatrix::REAL);

    if (!m_status.IsOk())
    {
        if (m_status.GetArg1() == 0)
        {
            m_status.SetArg1(6);
        }
        return m_status;
    }

    m_status = freqC2.Dimension(passBandFreq.M(), passBandFreq.N(), hwMatrix::REAL);

    if (!m_status.IsOk())
    {
        if (m_status.GetArg1() == 0)
        {
            m_status.SetArg1(7);
        }
        return m_status;
    }

    // compute lowpass prototype order and corner frequencies
    // work with steepest transition band and reassign the other stop
    // frequency for convenience
    double omegaP[2];
    double omegaS[2];

    omegaP[0] = passBandFreq(0);
    omegaP[1] = passBandFreq(1);
    omegaS[0] = stopBandFreq(0);
    omegaS[1] = stopBandFreq(1);

    double epsilon = sqrt(exp(passEdgeDb * log(10.0) / 10.0) - 1.0);
    double delta = sqrt(exp(stopEdgeDb * log(10.0) / 10.0) - 1.0);
    double omegaC1p;            // lowpass prototype normalized corner angular freq
    double omegaC2p;            // lowpass prototype normalized corner angular freq
    double omegaSP;             // lowpass prototype omegaS/omegaP ratio

    if (filterType == BandPass)
    {
        if (omegaS[0] * omegaS[1] < omegaP[0] * omegaP[1])
        {
            if (!cheby2)
            {
                omegaS[0] = omegaP[0] * omegaP[1] / omegaS[1];
            }
            else
            {
                omegaP[0] = omegaS[0] * omegaS[1] / omegaP[1];
            }
        }
        else
        {
            if (!cheby2)
            {
                omegaS[1] = omegaP[0] * omegaP[1] / omegaS[0];
            }
            else
            {
                omegaP[1] = omegaS[0] * omegaS[1] / omegaP[0];
            }
        }

        omegaSP = (omegaS[1] - omegaS[0]) / (omegaP[1] - omegaP[0]);
    }
    else // (filterType == BandStop)
    {
        if (omegaS[0] * omegaS[1] > omegaP[0] * omegaP[1])
        {
            if (!cheby2)
            {
                omegaS[0] = omegaP[0] * omegaP[1] / omegaS[1];
            }
            else
            {
                omegaP[0] = omegaS[0] * omegaS[1] / omegaP[1];
            }
        }
        else
        {
            if (!cheby2)
            {
                omegaS[1] = omegaP[0] * omegaP[1] / omegaS[0];
            }
            else
            {
                omegaP[1] = omegaS[0] * omegaS[1] / omegaP[0];
            }
        }

        omegaSP = (omegaP[1] - omegaP[0]) / (omegaS[1] - omegaS[0]);
    }

    DesignPrototype(epsilon, delta, omegaSP, order, omegaC1p, omegaC2p);

    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
        return m_status;
    }

    if (cheby2)
    {
        omegaC2p /= omegaSP;    // normalize to stopband
    }

    // transform range of first prototype corner frequency
    double a;
    double b;
    double c;

    if (filterType == BandPass)
    {
        if (!cheby2)
        {
            a = 1.0;
            b = -(omegaP[1] - omegaP[0]) * omegaC1p;
            c = -omegaP[0] * omegaP[1];
        }
        else
        {
            a = 1.0;
            b = -(omegaS[1] - omegaS[0]) * omegaC1p;
            c = -omegaS[0] * omegaS[1];
        }
    }
    else // (filterType == BandStop)
    {
        if (!cheby2)
        {
            a = omegaC1p;
            b = omegaP[1] - omegaP[0];
            c = -omegaP[0] * omegaP[1] * omegaC1p;
        }
        else
        {
            a = omegaC1p;
            b = omegaS[1] - omegaS[0];
            c = -omegaS[0] * omegaS[1] * omegaC1p;
        }
    }

    quadraticRoots(a, b, c, freqC1(0), freqC1(1));

    freqC1(0) = fabs(freqC1(0));
    freqC1(1) = fabs(freqC1(1));

    if (freqC1(0) > freqC1(1))
    {
        double value = freqC1(0);
        freqC1(0) = freqC1(1);
        freqC1(1) = value;
    }

    // transform range of second prototype corner frequency
    if (filterType == BandPass)
    {
        if (!cheby2)
        {
            a = 1.0;
            b = -(omegaP[1] - omegaP[0]) * omegaC2p;
            c = -omegaP[0] * omegaP[1];
        }
        else
        {
            a = 1.0;
            b = -(omegaS[1] - omegaS[0]) * omegaC2p;
            c = -omegaS[0] * omegaS[1];
        }
    }
    else // (filterType == BandStop)
    {
        if (!cheby2)
        {
            a = omegaC2p;
            b = omegaP[1] - omegaP[0];
            c = -omegaP[0] * omegaP[1] * omegaC2p;
        }
        else
        {
            a = omegaC2p;
            b = omegaS[1] - omegaS[0];
            c = -omegaS[0] * omegaS[1] * omegaC2p;
        }
    }

    quadraticRoots(a, b, c, freqC2(0), freqC2(1));

    freqC2(0) = fabs(freqC2(0));
    freqC2(1) = fabs(freqC2(1));

    if (freqC2(0) > freqC2(1))
    {
        double value = freqC2(0);
        freqC2(0) = freqC2(1);
        freqC2(1) = value;
    }

    return m_status;
}
//------------------------------------------------------------------------------
// Design a digital filter (low or high pass)
//------------------------------------------------------------------------------
hwMathStatus hwFilterDesign::Digital(double    passBandFreq, 
                                     double    stopBandFreq,
                                     double    passEdgeDb, 
                                     double    stopEdgeDb,
                                     int&      order,
                                     hwMatrix& freqC)
{
    // check filter definition
    double nearZero = 1.0e-5;
    Type   filterType;

    if (passBandFreq < nearZero || passBandFreq > 1.0 - nearZero)
    {
        if (passBandFreq < 0.0 || passBandFreq > 1.0)
        {
            return m_status(HW_MATH_ERR_FILTERFREQ_D, 1);
        }
        return m_status(HW_MATH_ERR_FILTERSPEC_E, 1);
    }

    if (stopBandFreq < nearZero || stopBandFreq > 1.0 - nearZero)
    {
        if (stopBandFreq < 0.0 || stopBandFreq > 1.0)
        {
            return m_status(HW_MATH_ERR_FILTERFREQ_D, 2);
        }
        return m_status(HW_MATH_ERR_FILTERSPEC_E, 2);
    }

    if (passBandFreq < stopBandFreq)
    {
        filterType = LowPass;
    }
    else if (stopBandFreq < passBandFreq)
    {
        filterType = HighPass;
    }
    else
    {
        return m_status(HW_MATH_ERR_FILTERFREQS_EQ, 1, 2);
    }

    if (passEdgeDb < nearZero)
    {
        return m_status(HW_MATH_ERR_DB_SIGN, 3);
    }
    if (stopEdgeDb < nearZero)
    {
        return m_status(HW_MATH_ERR_DB_SIGN, 4);
    }
    if (stopEdgeDb < passEdgeDb + nearZero)
    {
        return m_status(HW_MATH_ERR_FILTERRIPPLE, 3, 4);
    }

    // size outputs
    m_status = freqC.Dimension(2, hwMatrix::REAL);

    if (!m_status.IsOk())
    {
        if (m_status.GetArg1() == 0)
        {
            m_status.SetArg1(6);
        }
        return m_status;
    }

    // prewarp frequencies
    // The input frequencies are assumed to be normalized so that
    // F_normalized = F_cutoff / (F_sampling/2). However the algorithm
    // uses the F_normalized = F_cutoff / F_sampling convention. Therefore
    // the inputs must be multiplied by a factor of 1/2.
    double omegaP = tan(0.5 * PI * passBandFreq);
    double omegaS = tan(0.5 * PI * stopBandFreq);

    // compute lowpass prototype order and corner frequencies
    double epsilon = sqrt(exp(passEdgeDb * log(10.0) / 10.0) - 1.0);
    double delta = sqrt(exp(stopEdgeDb * log(10.0) / 10.0) - 1.0);
    double omegaC1p;            // lowpass prototype normalized corner angular freq
    double omegaC2p;            // lowpass prototype normalized corner angular freq
    double omegaSP;             // lowpass prototype omegaS/omegaP ratio

    if (filterType == LowPass)
    {
        omegaSP = omegaS / omegaP;
    }
    else // (filterType == HighPass)
    {
        omegaSP = omegaP / omegaS;
    }

    DesignPrototype(epsilon, delta, omegaSP, order, omegaC1p, omegaC2p);

    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
        return m_status;
    }

    // transform range of prototype corner frequency and reverse prewarp
    if (filterType == LowPass)
    {
        freqC(0) = 2.0 * atan(omegaP * omegaC1p) / PI;
        freqC(1) = 2.0 * atan(omegaP * omegaC2p) / PI;
    }
    else // (filterType == HighPass)
    {
        freqC(0) = 2.0 * atan(omegaP / omegaC1p) / PI;
        freqC(1) = 2.0 * atan(omegaP / omegaC2p) / PI;
    }

    return m_status;
}
//------------------------------------------------------------------------------
// Design a digital filter (bandpass or bandstop)
//------------------------------------------------------------------------------
hwMathStatus hwFilterDesign::Digital(const hwMatrix& passBandFreq, 
                                     const hwMatrix& stopBandFreq,
                                     double          passEdgeDb, 
                                     double          stopEdgeDb,
                                     int&            order, 
                                     hwMatrix&       freqC1, 
                                     hwMatrix&       freqC2, 
                                     bool            cheby2)
{
    // freqC1 contains the permissible range for the first corner freqency
    // freqC2 contains the permissible range for the second corner freqency

    // check filter definition
    double nearZero = 1.0e-5;
    Type filterType;

    if (!passBandFreq.IsReal())
    {
        return m_status(HW_MATH_ERR_COMPLEX, 1);
    }
    if (passBandFreq.Size() != 2)
    {
        return m_status(HW_MATH_ERR_VECTOR2, 1);
    }
    if (passBandFreq(0) < nearZero || passBandFreq(0) > 1.0 - nearZero)
    {
        if (passBandFreq(0) < 0.0 || passBandFreq(0) > 1.0)
        {
            return m_status(HW_MATH_ERR_FILTERFREQ_D, 1);
        }
        return m_status(HW_MATH_ERR_FILTERSPEC_E, 1);
    }

    if (passBandFreq(1) < nearZero || passBandFreq(1) > 1.0 - nearZero)
    {
        if (passBandFreq(1) < 0.0 || passBandFreq(1) > 1.0)
        {
            return m_status(HW_MATH_ERR_FILTERFREQ_D, 1);
        }
        return m_status(HW_MATH_ERR_FILTERSPEC_E, 1);
    }

    if (passBandFreq(0) >= passBandFreq(1))
    {
        return m_status(HW_MATH_ERR_NONINCREASE, 1);
    }
    if (!stopBandFreq.IsReal())
    {
        return m_status(HW_MATH_ERR_COMPLEX, 2);
    }
    if (stopBandFreq.Size() != 2)
    {
        return m_status(HW_MATH_ERR_VECTOR2, 2);
    }
    if (stopBandFreq(0) < nearZero || stopBandFreq(0) > 1.0 - nearZero)
    {
        if (stopBandFreq(0) < 0.0 || stopBandFreq(0) > 1.0)
        {
            return m_status(HW_MATH_ERR_FILTERFREQ_D, 2);
        }
        return m_status(HW_MATH_ERR_FILTERSPEC_E, 2);
    }

    if (stopBandFreq(1) < nearZero || stopBandFreq(1) > 1.0 - nearZero)
    {
        if (stopBandFreq(1) < 0.0 || stopBandFreq(1) > 1.0)
        {
            return m_status(HW_MATH_ERR_FILTERFREQ_D, 2);
        }
        return m_status(HW_MATH_ERR_FILTERSPEC_E, 2);
    }

    if (stopBandFreq(0) >= stopBandFreq(1))
    {
        return m_status(HW_MATH_ERR_NONINCREASE, 2);
    }
        
    if (stopBandFreq(0) < passBandFreq(0) && passBandFreq(1) < stopBandFreq(1))
    {
        filterType = BandPass;
    }
    else if (passBandFreq(0) < stopBandFreq(0) && stopBandFreq(1) < passBandFreq(1))
    {
        filterType = BandStop;
    }
    else
    {
        return m_status(HW_MATH_ERR_FILTERBANDCONF, 1, 2);
    }

    if (passEdgeDb < nearZero)
    {
        return m_status(HW_MATH_ERR_DB_SIGN, 3);
    }

    if (stopEdgeDb < nearZero)
    {
        return m_status(HW_MATH_ERR_DB_SIGN, 4);
    }

    if (stopEdgeDb < passEdgeDb + nearZero)
    {
        return m_status(HW_MATH_ERR_FILTERTYPE, 3, 4);
    }

    // size outputs
    m_status = freqC1.Dimension(passBandFreq.M(), passBandFreq.N(), hwMatrix::REAL);

    if (!m_status.IsOk())
    {
        if (m_status.GetArg1() == 0)
        {
            m_status.SetArg1(6);
        }
        return m_status;
    }

    m_status = freqC2.Dimension(passBandFreq.M(), passBandFreq.N(), hwMatrix::REAL);

    if (!m_status.IsOk())
    {
        if (m_status.GetArg1() == 0)
        {
            m_status.SetArg1(7);
        }
        return m_status;
    }

    // prewarp frequencies
    // The input frequencies are assumed to be normalized so that
    // F_normalized = F_cutoff / (F_sampling/2). However the algorithm
    // uses the F_normalized = F_cutoff / F_sampling convention. Therefore
    // the inputs must be multiplied by a factor of 1/2.
    double omegaP[2];
    double omegaS[2];

    omegaP[0] = tan(0.5 * PI * passBandFreq(0));
    omegaP[1] = tan(0.5 * PI * passBandFreq(1));
    omegaS[0] = tan(0.5 * PI * stopBandFreq(0));
    omegaS[1] = tan(0.5 * PI * stopBandFreq(1));

    // compute lowpass prototype order and corner frequencies
    // work with steepest transition band and reassign the other stop
    // frequency for convenience
    double epsilon = sqrt(exp(passEdgeDb * log(10.0) / 10.0) - 1.0);
    double delta = sqrt(exp(stopEdgeDb * log(10.0) / 10.0) - 1.0);
    double omegaC1p;            // lowpass prototype normalized corner angular freq
    double omegaC2p;            // lowpass prototype normalized corner angular freq
    double omegaSP;             // lowpass prototype omegaS/omegaP ratio

    if (filterType == BandPass)
    {
        if (omegaS[0] * omegaS[1] < omegaP[0] * omegaP[1])
        {
            if (!cheby2)
            {
                omegaS[0] = omegaP[0] * omegaP[1] / omegaS[1];
            }
            else
            {
                omegaP[0] = omegaS[0] * omegaS[1] / omegaP[1];
            }
        }
        else
        {
            if (!cheby2)
            {
                omegaS[1] = omegaP[0] * omegaP[1] / omegaS[0];
            }
            else
            {
                omegaP[1] = omegaS[0] * omegaS[1] / omegaP[0];
            }
        }

        omegaSP = (omegaS[1] - omegaS[0]) / (omegaP[1] - omegaP[0]);
    }
    else // (filterType == BandStop)
    {
        if (omegaS[0] * omegaS[1] > omegaP[0] * omegaP[1])
        {
            if (!cheby2)
            {
                omegaS[0] = omegaP[0] * omegaP[1] / omegaS[1];
            }
            else
            {
                omegaP[0] = omegaS[0] * omegaS[1] / omegaP[1];
            }
        }
        else
        {
            if (!cheby2)
            {
                omegaS[1] = omegaP[0] * omegaP[1] / omegaS[0];
            }
            else
            {
                omegaP[1] = omegaS[0] * omegaS[1] / omegaP[0];
            }
        }

        omegaSP = (omegaP[1] - omegaP[0]) / (omegaS[1] - omegaS[0]);
    }

    DesignPrototype(epsilon, delta, omegaSP, order, omegaC1p, omegaC2p);

    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
        return m_status;
    }

    if (cheby2)
    {
        omegaC2p /= omegaSP;    // normalize to stopband
    }

    // transform range of first prototype corner frequency and reverse prewarp
    double a;
    double b;
    double c;
    double omegaCa;             // prewarped corner freq
    double omegaCb;             // prewarped corner freq

    if (filterType == BandPass)
    {
        if (!cheby2)
        {
            a = 1.0;
            b = -(omegaP[1] - omegaP[0]) * omegaC1p;
            c = -omegaP[0] * omegaP[1];
        }
        else
        {
            a = 1.0;
            b = -(omegaS[1] - omegaS[0]) * omegaC1p;
            c = -omegaS[0] * omegaS[1];
        }
    }
    else // (filterType == BandStop)
    {
        if (!cheby2)
        {
            a = omegaC1p;
            b = omegaP[1] - omegaP[0];
            c = -omegaP[0] * omegaP[1] * omegaC1p;
        }
        else
        {
            a = omegaC1p;
            b = omegaS[1] - omegaS[0];
            c = -omegaS[0] * omegaS[1] * omegaC1p;
        }
    }

    quadraticRoots(a, b, c, omegaCa, omegaCb);

    freqC1(0) = fabs(2.0 * atan(omegaCa) / PI);
    freqC1(1) = fabs(2.0 * atan(omegaCb) / PI);

    if (freqC1(0) > freqC1(1))
    {
        double value = freqC1(0);
        freqC1(0) = freqC1(1);
        freqC1(1) = value;
    }

    // transform range of second prototype corner frequency and reverse prewarp
    if (filterType == BandPass)
    {
        if (!cheby2)
        {
            a = 1.0;
            b = -(omegaP[1] - omegaP[0]) * omegaC2p;
            c = -omegaP[0] * omegaP[1];
        }
        else
        {
            a = 1.0;
            b = -(omegaS[1] - omegaS[0]) * omegaC2p;
            c = -omegaS[0] * omegaS[1];
        }
    }
    else // (filterType == BandStop)
    {
        if (!cheby2)
        {
            a = omegaC2p;
            b = omegaP[1] - omegaP[0];
            c = -omegaP[0] * omegaP[1] * omegaC2p;
        }
        else
        {
            a = omegaC2p;
            b = omegaS[1] - omegaS[0];
            c = -omegaS[0] * omegaS[1] * omegaC2p;
        }
    }

    quadraticRoots(a, b, c, omegaCa, omegaCb);

    freqC2(0) = fabs(2.0 * atan(omegaCa) / PI);
    freqC2(1) = fabs(2.0 * atan(omegaCb) / PI);

    if (freqC2(0) > freqC2(1))
    {
        double value = freqC2(0);
        freqC2(0) = freqC2(1);
        freqC2(1) = value;
    }

    return m_status;
}
